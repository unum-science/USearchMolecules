"""Readers for the sharded Parquet collections and the USearch graphs published beside them.

Small molecules are keyed by an InChIKey UUIDv8 and searched through the ECFP4 and Chamfer graphs,
proteins and antibodies by the MD5 of their sequence and searched through the Needleman-Wunsch graph.

Example usage:
    dataset = FingerprintedDataset.open("data/pubchem")
    results = dataset.search("CC(=O)OC1=CC=CC=C1C(=O)O", count=10)  # Aspirin
"""

from __future__ import annotations

import bisect
import hashlib
import os
import random
from dataclasses import dataclass, field
from itertools import islice

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import stringzilla as sz
from tqdm import tqdm

from usearchmolecules.to_fingerprint import (
    FingerprintShape,
    shape_mixed,
    smiles_to_maccs_ecfp4_fcfp4,
)

SEED = 42  # For reproducibility
SHARD_SIZE = 1_000_000  # This would result in files between 150 and 300 MB
BATCH_SIZE = 100_000  # A good threshold to split insertions


@dataclass
class FingerprintedEntry:
    """SMILES string augmented with a potentially hybrid fingerprint of known `FingerprintShape` shape."""

    smiles: str
    fingerprint: np.ndarray
    key: int | None = None

    @staticmethod
    def from_table_row(table: pa.Table, row: int, shape: FingerprintShape) -> FingerprintedEntry:
        fingerprint = np.zeros(shape.nbytes, dtype=np.uint8)
        progress = 0
        if shape.include_maccs:
            fingerprint[progress : progress + 21] = table["maccs"][row].as_buffer()
            progress += 21
        if shape.nbytes_padding:
            progress += shape.nbytes_padding
        if shape.include_ecfp4:
            fingerprint[progress : progress + 256] = table["ecfp4"][row].as_buffer()
            progress += 256
        if shape.include_fcfp4:
            fingerprint[progress : progress + 256] = table["fcfp4"][row].as_buffer()
            progress += 256
        return FingerprintedEntry(smiles=table["smiles"][row], fingerprint=fingerprint)

    @staticmethod
    def from_parts(
        smiles: str | None,
        maccs: np.ndarray | None,
        ecfp4: np.ndarray | None,
        fcfp4: np.ndarray | None,
        shape: FingerprintShape,
    ) -> FingerprintedEntry:
        # Each modality is read only when the shape includes it, so an omitted
        # modality may be passed as None; likewise `smiles` is optional for index-only rows.
        fingerprint = np.zeros(shape.nbytes, dtype=np.uint8)
        progress = 0
        if shape.include_maccs:
            assert maccs is not None
            fingerprint[progress : progress + 21] = maccs
            progress += 21
        if shape.nbytes_padding:
            progress += shape.nbytes_padding
        if shape.include_ecfp4:
            assert ecfp4 is not None
            fingerprint[progress : progress + 256] = ecfp4
            progress += 256
        if shape.include_fcfp4:
            assert fcfp4 is not None
            fingerprint[progress : progress + 256] = fcfp4
            progress += 256
        return FingerprintedEntry(smiles=smiles, fingerprint=fingerprint)


def inchikey_of(key: bytes) -> tuple[str, int]:
    """InChIKey and enhanced-stereo flag packed into a small-molecule `key`.

    The key is an RFC 9562 UUIDv8 as 16 big-endian bytes. Its 122 custom bits hold a 117-bit payload
    `P << 2 | flag`, zero-extended and placed big-endian around the version and variant fields. `P` is the
    22 InChIKey block letters as one base-26 number, then 1 bit for the standard flag, 5 for the version
    letter and 5 for protonation. `flag` marks CXSMILES enhanced-stereo groups: 0 for none or only absolute,
    1 for AND groups, 2 for OR groups, 3 for both.
    """
    value = int.from_bytes(key, "big")
    payload = (value >> 80) << 74 | ((value >> 64) & 0xFFF) << 62 | value & ((1 << 62) - 1)
    flag, packed = payload & 3, payload >> 2
    protonation, packed = packed & 31, packed >> 5
    version, packed = packed & 31, packed >> 5
    standard, packed = packed & 1, packed >> 1
    letters = ""
    for _ in range(22):
        packed, letter = divmod(packed, 26)
        letters = chr(65 + letter) + letters
    return f"{letters[:14]}-{letters[14:]}{'S' if standard else 'N'}{chr(65 + version)}-{chr(65 + protonation)}", flag


def open_graph(kind, path: os.PathLike, residency: str) -> tuple:
    """A `.safegraphs` graph of USearch class `kind` and its executor, loaded into RAM or served from NVMe in place."""
    import usearch

    if residency == "ram":
        executor = usearch.Executor.for_ram()
        return kind.load_ram(path, executor=executor), executor
    executor = usearch.Executor.for_nvme()
    return kind.view_nvme(path, executor=executor), executor


def is_fingerprint_shard(filename: str) -> bool:
    """Whether a `parquet/` entry is a 2D fingerprint shard rather than a `*.3D.parquet` conformer shard.

    Both families share the directory and the stem, so a bare `*.parquet` glob picks up conformer shards
    too — and their stems parse as a duplicate `first_key` pointing at a `.3D.smi` that does not exist.
    """
    return filename.endswith(".parquet") and ".3D." not in filename


def shard_name(directory: str, from_index: int, to_index: int, kind: str):
    """Generate standardized shard filename with zero-padded indices."""
    return os.path.join(directory, kind, f"{from_index:0>10}-{to_index:0>10}.{kind}")


def write_table(table: pa.Table, path_out: os.PathLike):
    """Write PyArrow table to Parquet file with optimized settings for molecular data."""
    return pq.write_table(
        table,
        path_out,
        # Without compression the file size may be too large.
        # compression="NONE",
        write_statistics=False,
        store_schema=True,
        use_dictionary=False,
    )


@dataclass
class ShardConformers:
    """Lowest-energy stored conformer of every molecule in one shard, in row order.

    Row `i` occupies `coords[offsets[i] : offsets[i + 1]]` as tightly packed xyz triplets in the
    centered frame. An empty slice means the pipeline produced no geometry for that molecule.
    """

    offsets: np.ndarray
    coords: np.ndarray
    energies: np.ndarray

    def __getitem__(self, row: int) -> tuple[np.ndarray, float] | None:
        start, end = int(self.offsets[row]), int(self.offsets[row + 1])
        if start == end:
            return None
        return self.coords[start:end].reshape(-1, 3), float(self.energies[row])


@dataclass
class FingerprintedShard:
    """Potentially cached table and smiles path containing up to `SHARD_SIZE` entries."""

    first_key: int
    name: str

    table_path: os.PathLike
    smiles_path: os.PathLike
    conformers_path: os.PathLike
    table_cached: pa.Table | None = None
    smiles_caches: sz.Strs | None = None
    conformers_cached: ShardConformers | None = None

    @property
    def is_complete(self) -> bool:
        return os.path.exists(self.table_path) and os.path.exists(self.smiles_path)

    @property
    def table(self) -> pa.Table:
        return self.load_table()

    @property
    def smiles(self) -> sz.Strs:
        return self.load_smiles()

    def load_table(self, columns=None, view=False) -> pa.Table:
        """Load Parquet table into memory with optional memory-mapping and column filtering."""
        if not self.table_cached:
            self.table_cached = pq.read_table(
                self.table_path,
                memory_map=view,
                columns=columns,
            )
        return self.table_cached

    def load_smiles(self) -> sz.Strs:
        """Load SMILES strings using StringZilla for fast line-by-line access."""
        if not self.smiles_caches:
            self.smiles_caches = sz.Str(sz.File(self.smiles_path)).splitlines()
        return self.smiles_caches

    def load_conformers(self) -> ShardConformers:
        """Read the `*.3D.parquet` sibling, keeping each molecule's lowest-energy conformer.

        Read whole rather than sought into: the shards carry no Parquet page index, so reaching one
        molecule decodes every row before it anyway. Costs a few hundred megabytes per shard.
        """
        if self.conformers_cached is None:
            chunks: list[np.ndarray] = []
            lengths: list[np.ndarray] = []
            energies: list[np.ndarray] = []
            for batch in pq.ParquetFile(self.conformers_path).iter_batches(
                batch_size=BATCH_SIZE,
                columns=["conformer_index", "conformer_energy", "conformer_coords"],
            ):
                lowest = batch.filter(pc.equal(batch["conformer_index"], 0))
                coords = lowest["conformer_coords"]
                # `flatten` drops the null lists entirely, so lengths must be filled independently.
                chunks.append(np.asarray(coords.flatten(), dtype=np.float16))
                lengths.append(pc.fill_null(pc.list_value_length(coords), 0).to_numpy(zero_copy_only=False))
                energies.append(pc.fill_null(lowest["conformer_energy"], np.nan).to_numpy(zero_copy_only=False))

            offsets = np.zeros(sum(len(part) for part in lengths) + 1, dtype=np.int64)
            np.cumsum(np.concatenate(lengths), dtype=np.int64, out=offsets[1:])
            self.conformers_cached = ShardConformers(
                offsets=offsets,
                coords=np.concatenate(chunks),
                energies=np.concatenate(energies).astype(np.float32),
            )
        return self.conformers_cached


@dataclass
class FingerprintedDataset:
    """Sharded small-molecule collection with its fingerprints, conformers and graphs."""

    directory: os.PathLike
    shards: list[FingerprintedShard]
    shape: FingerprintShape | None = None
    graphs: dict[tuple[str, str], tuple] = field(default_factory=dict, repr=False)

    @staticmethod
    def open(
        directory: os.PathLike,
        shape: FingerprintShape | None = None,
        max_shards: int | None = None,
    ) -> FingerprintedDataset:
        """List the shards of a collection; graphs load on the first `graph` or `search` call.

        Args:
            directory: Collection root holding the `parquet/` subdirectory and the graphs
            shape: Default fingerprint layout for `head`
            max_shards: Optional limit on number of shards to load (useful for testing)
        """

        if directory is None:
            return FingerprintedDataset(directory=None, shards=[], shape=shape)

        shards = []
        filenames = sorted(f for f in os.listdir(os.path.join(directory, "parquet")) if is_fingerprint_shard(f))
        if max_shards:
            filenames = filenames[:max_shards]

        for filename in tqdm(filenames, unit="shard"):
            stem = filename.removesuffix(".parquet")
            shards.append(
                FingerprintedShard(
                    first_key=int(stem.partition("-")[0]),
                    name=stem,
                    table_path=os.path.join(directory, "parquet", stem + ".parquet"),
                    smiles_path=os.path.join(directory, "smiles", stem + ".smi"),
                    conformers_path=os.path.join(directory, "parquet", stem + ".3D.parquet"),
                )
            )

        print(f"Fetched {len(shards)} shards")
        return FingerprintedDataset(directory=directory, shards=shards, shape=shape)

    def shard_containing(self, key: int) -> FingerprintedShard:
        """Find the shard containing a given key/index, or raise if none does."""
        for shard in self.shards:
            if shard.first_key <= key < shard.first_key + SHARD_SIZE:
                return shard
        raise KeyError(f"No shard contains key {key}")

    def head(
        self,
        max_rows: int,
        shape: FingerprintShape | None = None,
        shuffle: bool = False,
    ) -> tuple[list[str], list[int], np.ndarray]:
        """Load the first N rows from the dataset for preview and testing.

        Args:
            max_rows: Maximum number of molecules to load
            shape: Fingerprint configuration (uses dataset default if None)
            shuffle: Whether to randomly shuffle the results

        Returns:
            Tuple of (SMILES strings, keys, fingerprint matrix)
        """

        if self.directory is None:
            return None

        shape = shape or self.shape
        entries = list(islice(self._iter_entries(shape), max_rows))

        smiles = np.array([entry.smiles for entry in entries], dtype=object)
        keys = np.arange(len(entries), dtype=np.uint64)
        fingers = np.vstack([entry.fingerprint for entry in entries])
        if shuffle:
            permutation = np.random.default_rng(SEED).permutation(len(entries))
            smiles, keys, fingers = smiles[permutation], keys[permutation], fingers[permutation]
        return smiles, keys, fingers

    def _iter_entries(self, shape: FingerprintShape):
        """Yield every entry across all shards in on-disk order."""
        for shard in self.shards:
            table = shard.load_table()
            for row in range(len(table)):
                yield FingerprintedEntry.from_table_row(table, row, shape)

    def graph(self, name: str = "ecfp4", residency: str = "ram") -> tuple:
        """The `index-<name>.safegraphs` graph, `ecfp4` or `chamfer`, and its executor, cached per residency."""
        import usearch

        if (name, residency) not in self.graphs:
            kind = {"ecfp4": usearch.DenseIndex, "chamfer": usearch.ShapesIndex}[name]
            path = os.path.join(self.directory, f"index-{name}.safegraphs")
            self.graphs[name, residency] = open_graph(kind, path, residency)
        return self.graphs[name, residency]

    def search(
        self,
        smiles: str,
        count: int = 10,
        expansion: int = 0,
        residency: str = "ram",
    ) -> list[tuple[bytes, float]]:
        """Molecules nearest to `smiles` in the ECFP4 graph as `(key, Hamming distance)`."""
        import usearch

        index, executor = self.graph("ecfp4", residency)
        _, ecfp4, _ = smiles_to_maccs_ecfp4_fcfp4(smiles)
        arena = usearch.Arena(executor)
        matches = index.on(executor).search(ecfp4[np.newaxis], count, arena, expansion_factor=expansion)
        return [(match[0].bytes, float(match[1])) for match in matches[0]]

    def search_shapes(
        self,
        points: np.ndarray,
        count: int = 10,
        expansion: int = 0,
        residency: str = "ram",
    ) -> list[tuple[bytes, float]]:
        """Molecules whose lowest-energy conformer is nearest to one `(n, 3)` centered cloud, under Chamfer distance."""
        import usearch

        index, executor = self.graph("chamfer", residency)
        cloud = np.ascontiguousarray(points, dtype=np.float16)
        arena = usearch.Arena(executor)
        matches = index.on(executor).search([cloud], count, arena, expansion_factor=expansion)
        return [(match[0].bytes, float(match[1])) for match in matches[0]]

    def conformer(self, key: int) -> tuple[np.ndarray, float] | None:
        """Lowest-energy stored conformer for `key` as `(xyz, energy)`, or None when none was produced.

        Coordinates follow the atom order of the molecule's largest fragment after `AddHs`, so they
        line up with an RDKit molecule rebuilt from the `smiles` column.
        """
        shard = self.shard_containing(key)
        return shard.load_conformers()[int(key - shard.first_key)]

    def __len__(self) -> int:
        return sum(pq.read_metadata(shard.table_path).num_rows for shard in self.shards)

    def random_smiles(self) -> str:
        """Get a random SMILES string from the dataset."""
        shard = random.choice(self.shards)
        row = random.randrange(pq.read_metadata(shard.table_path).num_rows)
        return read_row(shard.table_path, row)["smiles"][0].as_py()


def sequence_key(sequence: str) -> bytes:
    """MD5 of the uppercase ASCII sequence, the `key` column of the protein and antibody subsets."""
    return hashlib.md5(sequence.upper().encode("ascii")).digest()


def pair_key(heavy: str, light: str) -> bytes:
    """Key of a paired antibody cell, the MD5 of the heavy chain key followed by the light chain key."""
    return hashlib.md5(sequence_key(heavy) + sequence_key(light)).digest()


def read_row(path: os.PathLike, row: int) -> pa.Table:
    """Single-row table at `row` of a Parquet file, decoding only the row group that holds it."""
    file = pq.ParquetFile(path)
    for group in range(file.num_row_groups):
        group_rows = file.metadata.row_group(group).num_rows
        if row < group_rows:
            return file.read_row_group(group).slice(row, 1)
        row -= group_rows
    raise IndexError(f"Row past the end of {path}")


@dataclass
class SequenceDataset:
    """Row-addressed reader for the `alphafold`, `oas-unpaired` and `oas-paired` sequence subsets."""

    directory: os.PathLike
    stems: list[str]
    starts: list[int]
    graphs: dict[str, tuple] = field(default_factory=dict, repr=False)

    @staticmethod
    def open(directory: os.PathLike) -> SequenceDataset:
        """List the 2D stems under `directory/parquet` and their row counts from the Parquet footers."""
        filenames = sorted(f for f in os.listdir(os.path.join(directory, "parquet")) if is_fingerprint_shard(f))
        stems = [filename.removesuffix(".parquet") for filename in filenames]
        starts = [0]
        for filename in filenames:
            starts.append(starts[-1] + pq.read_metadata(os.path.join(directory, "parquet", filename)).num_rows)
        return SequenceDataset(directory=directory, stems=stems, starts=starts)

    def __len__(self) -> int:
        return self.starts[-1]

    def locate(self, index: int) -> tuple[str, int]:
        """Stem holding the dataset-wide row `index` and the row's position inside it."""
        if not 0 <= index < len(self):
            raise IndexError(f"Row {index} outside {len(self)} rows")
        position = bisect.bisect_right(self.starts, index) - 1
        return self.stems[position], index - self.starts[position]

    def row(self, index: int) -> dict:
        """All 2D columns of one row, with dictionary columns decoded to strings."""
        stem, row = self.locate(index)
        return read_row(os.path.join(self.directory, "parquet", stem + ".parquet"), row).to_pylist()[0]

    def structure(self, index: int) -> tuple[np.ndarray, np.ndarray | None]:
        """Cα coordinates of an `alphafold` row as `(n, 3)` float16 and its per-residue pLDDT, None when absent."""
        stem, row = self.locate(index)
        table = read_row(os.path.join(self.directory, "parquet", stem + ".3D.parquet"), row)
        points = np.asarray(pc.list_flatten(pc.list_flatten(table["points"])), dtype=np.float16).reshape(-1, 3)
        if "plddt" not in table.column_names:
            return points, None
        return points, np.asarray(pc.list_flatten(table["plddt"]), dtype=np.uint8)

    @property
    def graph_path(self) -> str:
        """Needleman-Wunsch graph under BLOSUM62 over the distinct keys, opened by `graph`.

        Large subsets keep a JSON manifest at this path beside its `index-nw-0000k-of-0000n.safegraphs` parts.
        """
        return os.path.join(self.directory, "index-nw.safegraphs")

    def graph(self, residency: str = "ram") -> tuple:
        """The Needleman-Wunsch graph and its executor, cached per residency."""
        import usearch

        if residency not in self.graphs:
            self.graphs[residency] = open_graph(usearch.SequenceIndex, self.graph_path, residency)
        return self.graphs[residency]

    def search(
        self,
        sequences: list[str],
        count: int = 10,
        expansion: int = 0,
        residency: str = "ram",
    ) -> list[list[tuple[bytes, float]]]:
        """Nearest distinct sequences to each query as `(key, Needleman-Wunsch distance)`, one list per query."""
        import usearch

        index, executor = self.graph(residency)
        arena = usearch.Arena(executor)
        matches = index.on(executor).search(sequences, count, arena, expansion_factor=expansion)
        return [[(match[0].bytes, float(match[1])) for match in matches[query]] for query in range(len(sequences))]


if __name__ == "__main__":
    dataset = FingerprintedDataset.open("data/example")
    print(dataset.search("C"))
