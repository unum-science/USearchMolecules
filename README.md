# USearchMolecules

USearchMolecules is an open, searchable atlas of chemical and biological entities: the small molecules people make, screen and buy, and the proteins and antibodies those molecules are designed against.
Every entity is one row under a permanent key derived from the entity itself, so a compound or a sequence lands on the same key in every collection that holds it.
Everything computed about an entity, from fingerprints and conformers to predicted structures and germline calls, is added beside the same rows, together with the graphs that link each entity to its nearest neighbors.
New collections, representations and annotations extend the atlas without renumbering what is already there.

![USearchMolecules atlas thumbnail](https://github.com/ashvardanian/ashvardanian/raw/master/repositories/USearchMolecules.jpg?raw=true)

It lets a chemist, a biologist, an automated pipeline or an agent ask, in milliseconds:

- _"Show me 100 molecules shaped like this one but without acrylates."_
- _"Show me everything in Enamine that swaps a methyl-pyridine for an ethyl-imidazole on this scaffold."_
- _"Cluster my 50K screening hits by their 3D shape similarity."_
- _"Find all molecules with the same Bemis-Murcko scaffold as paclitaxel that aren't PAINS-flagged."_
- _"Which human proteins align closest to this sequence, and which of them carry a kinase domain?"_
- _"Which observed antibodies sit nearest to this candidate, and how far has each matured from its germline?"_

The atlas includes three small-molecule collections, the predicted structures of UniProt, and the observed antibody repertoires:

- __Example__, `example`: 2'000'000 compounds drawn from the three below, the one to start with.
- __PubChem__, `pubchem`: 115'627'256 compounds deposited in NCBI PubChem.
- __GDB-13__, `gdb13`: 977'468'267 molecules enumerated by the University of Bern.
- __Enamine REAL__, `real`: 6'039'411'650 make-on-demand compounds.
- __AlphaFold DB__, `alphafold`: 241'070'489 predicted protein structures, synced to UniProt 2025_03.
- __OAS Unpaired__, `oas-unpaired`: 2'427'646'880 antibody chains from the Observed Antibody Space.
- __OAS Paired__, `oas-paired`: 3'087'576 antibody cells, one heavy and one light chain each.

Every entity carries what its kind is searched and filtered by:

- __Small molecules__: the SMILES graph; MACCS and PubChem substructure keys; ECFP4 and FCFP4 circular fingerprints; up to three MMFF94-relaxed conformers with their energies and USRCAT shape vectors.
- __Proteins__: the amino-acid sequence; the Cα trace with AlphaFold's per-residue confidence; organism and domain names; UniProt curation, gene names, enzyme classes, Gene Ontology terms and Pfam domains.
- __Antibodies__: the variable-domain sequence and its CDR3; framework and CDR lengths; V, D and J germline calls with V and J identity; the study, donor and run each chain was observed in.

[Unum's USearch](https://github.com/unum-cloud/USearch) graphs link every entity to its nearest neighbors under the same keys as the rows: compounds by ECFP4 fingerprints, PubChem compounds also by the Chamfer distance between their conformers, and proteins and antibodies by Needleman-Wunsch alignment.
USearch takes a user-defined metric, and [NumKong](https://github.com/ashvardanian/NumKong) supplies the geometric kernels one can call — Kabsch and Umeyama alignment over the conformer coordinates themselves.
The conformers come from [SmallField](https://github.com/unum-science/SmallField), which runs ETKDG embedding and MMFF94 relaxation entirely on the GPU.
It is published on [AWS Open Data](https://registry.opendata.aws/usearch-molecules/), downloadable without an AWS account, as covered under [Usage](#usage).

## Dataset Structure

Each collection holds `parquet/` with one row per entity, cut into 1M-row stems named by their first and end rows, and its graphs beside it as `index-<column>.safegraphs`.
A `<stem>.3D.parquet` file is row-aligned with the `<stem>.parquet` of the same stem, so row `i` of both describes one entity and the two need no join; across stems and collections, rows meet on `key`.
A graph split into parts keeps a JSON manifest under the bare name, beside its `-0000k-of-0000n` parts.
S3 verifies content on transfer, so a completed download needs no separate manifest.

```sh
.
└── data
    ├── example                                  # 3.5 GB, the same shape as pubchem without the Chamfer graph
    ├── pubchem
    │   ├── index-ecfp4.safegraphs               # 62.0 GB
    │   ├── index-chamfer.safegraphs             # ~52 GB, over the lowest-energy conformers
    │   └── parquet                              # 22.5 GB 2D + 129.8 GB 3D
    │       ├── 0000000000-0001000000.parquet    # 195 MB, SMILES and fingerprints
    │       ├── 0000000000-0001000000.3D.parquet # ~1.1 GB, conformers and shape vectors
    │       ├── ...
    │       └── 0115000000-0116000000.3D.parquet
    ├── gdb13                                    # 142.3 GB 2D + 681.1 GB 3D, the same shape as example
    ├── real
    │   ├── index-ecfp4.safegraphs               # manifest of the six parts below
    │   ├── index-ecfp4-00001-of-00006.safegraphs # ~540 GB, over 1B compounds
    │   ├── ...
    │   └── parquet                              # 1.1 TB 2D + 7.1 TB 3D
    ├── alphafold
    │   ├── index-nw.safegraphs                  # manifest of its index-nw-0000k-of-0000n parts
    │   └── parquet                              # 242 stems
    │       ├── 0000000000-0001000000.parquet    # ~225 MB, sequences and annotations
    │       ├── 0000000000-0001000000.3D.parquet # ~2.1 GB, Cα traces and pLDDT
    │       ├── ...
    │       └── 0241000000-0241070489.3D.parquet
    ├── oas-unpaired                             # 163.2 GB in 2'428 stems, the same shape as oas-paired, graph in parts
    └── oas-paired
        ├── index-nw.safegraphs                  # over heavy and light chains
        └── parquet                              # 374 MB, no 3D family
            ├── 0000000000-0001000000.parquet    # 121 MB, cells with both chains
            ├── ...
            └── 0003000000-0003087576.parquet
```

Every shard documents its own schema and provenance in Parquet metadata, readable without downloading the rows:

```py
import pyarrow.parquet as pq

pq.read_schema("data/example/parquet/0000000000-0001000000.parquet")      # 2D columns, key layout, schema version
pq.read_schema("data/example/parquet/0000000000-0001000000.3D.parquet")   # conformer columns, producer, max_atoms
```

One small-molecule row looks like:

| Column    | Value                                                      |
| :-------- | :--------------------------------------------------------- |
| `key`     | `0484b0ae-0a7a-8c49-bdc8-c6f845ff9034`                     |
| `smiles`  | `CNCC(C)NC(=O)C1(C(C)(C)OC)CC1`                            |
| `maccs`   | `0x00000200000002002021227C488B9C02100615FFCC`             |
| `pubchem` | `0x00733000000000000000000000001800…` truncated, 111 bytes |
| `ecfp4`   | `0x40000000000000000000800000002400…` truncated, 256 bytes |
| `fcfp4`   | `0xE0001400000000000000000000000000…` truncated, 256 bytes |

## Installation

Install straight from GitHub with [uv](https://github.com/astral-sh/uv), no clone required.

```sh
uv pip install "git+https://github.com/unum-science/USearchMolecules.git"                          # read and search the dataset
uv pip install "usearchmolecules[dev] @ git+https://github.com/unum-science/USearchMolecules.git"  # build fingerprints and indexes
uv pip install "usearchmolecules[viz] @ git+https://github.com/unum-science/USearchMolecules.git"  # StreamLit visualizer
uv pip install "usearchmolecules[all] @ git+https://github.com/unum-science/USearchMolecules.git"  # everything
```

To hack on the code, clone it and install editable instead.

```sh
git clone https://github.com/unum-science/USearchMolecules.git
cd USearchMolecules
uv venv --python 3.12
source .venv/bin/activate
uv pip install -e ".[all]"
```

The pip extras are CPU-only; for GPU-accelerated conformer generation, set up the environment via pixi as shown below.

### Via Pixi for GPU Acceleration

For GPU acceleration with nvMolKit, we recommend using [pixi](https://pixi.sh) which handles conda dependencies such as RDKit and nvMolKit seamlessly:

```sh
pixi install
pixi run python -m usearchmolecules.prep_conformers --datasets example --use-gpu --conformers 20 --batch-size 20
```

## Usage

### Downloading

Every collection sits under `s3://usearch-molecules/data/` on AWS Open Data and downloads without an account.
Start with `example`, then fetch whichever collections you need:

```sh
for collection in example pubchem gdb13 real alphafold oas-unpaired oas-paired; do
    aws s3 sync --no-sign-request s3://usearch-molecules/data/$collection/ data/$collection/
done
```

Add `--exclude "index-*"` to fetch only the shards and leave the graphs behind, or drop the path to browse: `aws s3 ls --no-sign-request s3://usearch-molecules/data/`.
The Nebius and Hugging Face copies still hold the previous schema until they are refreshed.

### Searching from Python

Small molecules are searched by their ECFP4 fingerprint, and every hit comes back as a key with its Hamming distance:

```py
from usearchmolecules.dataset import FingerprintedDataset, inchikey_of

data = FingerprintedDataset.open("data/example")

# No inspiration? Pick a random molecule with `data.random_smiles()`
for key, distance in data.search("CC(O)C(CN)=NNCC(C)(C)C", 100):
    inchikey, stereo_flag = inchikey_of(key)   # the InChIKey PubChem or ChEMBL hold for the same compound
```

The key decodes back to a standard InChIKey, so a hit joins any other InChIKey source, and the same key finds the compound in every collection that holds it.

PubChem compounds can also be searched by shape, with any conformer as the query, such as a row of the shards themselves:

```py
import numpy as np
import pyarrow.parquet as pq

row = pq.read_table("data/pubchem/parquet/0000000000-0001000000.3D.parquet").slice(0, 1).to_pylist()[0]
cloud = np.asarray(row["points"][: row["n_atoms"]], dtype=np.float16)   # the lowest-energy conformer

pubchem = FingerprintedDataset.open("data/pubchem")
shaped_alike = pubchem.search_shapes(cloud, 10)                         # keys and Chamfer distances
```

Proteins and antibodies open with one reader for all three collections, and search by global alignment:

```py
from usearchmolecules.dataset import SequenceDataset, sequence_key

proteins = SequenceDataset.open("data/alphafold")
entry = proteins.row(0)                   # sequence, accession, organism, gene, GO terms, Pfam domains
trace, plddt = proteins.structure(0)      # (n_residues, 3) float16 Cα coordinates, uint8 confidence
aligned = proteins.search([entry["sequence"]], 10)[0]

antibodies = SequenceDataset.open("data/oas-paired")
cell = antibodies.row(0)                  # both chains, their CDR3 loops, germline calls and identities
assert cell["key_heavy"] == sequence_key(cell["sequence_heavy"])
```

Every search takes `residency="nvme"` to read a graph in place from the drive rather than load it, so a multi-terabyte graph answers queries on a machine with a fraction of that memory.
A manifest such as `data/real/index-ecfp4.safegraphs` is searched like a single graph, its parts in turn.

### Searching from the Command Line

The same graphs are served by the `usearch` command line, which takes one sequence per line or one point cloud per query:

```sh
usearch sequence search data/alphafold/index-nw.safegraphs --query-texts my-proteins.txt --top-k 10
usearch sequence search data/oas-unpaired/index-nw.safegraphs --query-texts my-chains.txt --residency nvme
usearch shapes search data/pubchem/index-chamfer.safegraphs --query-points clouds.f16bin --query-counts clouds.ibin
```

Every match prints its `key`, so it joins straight back to the Parquet rows and, across collections, to every other row with the same compound or sequence.

## Exploring Dataset via Graphical Interface

The dataset also comes with a graphical sandbox for small molecules, implemented with StreamLit and 3DMol.js to help visualize similarities between them.

```sh
streamlit run streamlit_app.py
```

![USearchMolecules StreamLit demo preview](/assets/USearchMoleculesStreamLitPreview.gif)

## Methodology

### Dataset Sources

Original data came from:

- __PubChem__: [CID-SMILES](https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz)
- __GDB-13__: [Bern archive](https://zenodo.org/record/5172018/files/gdb13.tgz?download=1)
- __Enamine REAL__: [CXSMILES](https://enamine.net/compound-collections/real-compounds/real-database) by heavy-atom count
- __Example__: sampled from the three by `prep_sample`
- __AlphaFold DB v6__: [Foldseek](https://foldseek.steineggerlab.workers.dev/afdb.tar.gz) and [Foldcomp](https://foldcomp.steineggerlab.workers.dev/afdb_uniprot_v6) archives, CC BY 4.0
- __UniProt__: [2025_03](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2025_03/) and InterPro 106.0, CC BY 4.0 and CC0
- __OAS__: [unpaired and paired units](https://opig.stats.ox.ac.uk/webapps/oas/), CC BY 4.0

### Pipeline

The small-molecule shards come from the package's scripts and SmallField:

- `prep_parquet.py` converts each source into shuffled Parquet shards of SMILES strings.
- `prep_encode.py` adds the MACCS, PubChem, ECFP4 and FCFP4 fingerprints.
- [SmallField](https://github.com/unum-science/SmallField) embeds and relaxes the conformers and computes their USRCAT vectors into the `.3D.parquet` shards.
- `prep_sample.py` draws the `example` collection from the three sources.

Every script works with bigger-than-memory data in a "gliding-window" fashion, splitting the heavy parts between processes and threads, so processing 1 TB of molecules doesn't require 1 TB of RAM.

```sh
uv run python -m usearchmolecules.prep_parquet --datasets example
uv run python -m usearchmolecules.prep_encode --datasets example
```

Multiple passes over the finished shards then key every row and regroup the conformers:

- Each 2D row is keyed by RDKit's standard InChI of its `smiles`, parsed with sanitization and, when RDKit refuses that, without it.
  The InChIKey and the row's CXSMILES enhanced-stereo groups are packed into the UUIDv8 that `inchikey_of` decodes.
  Rows that InChI cannot name are dropped rather than given a substitute key.
- Each 3D shard is regrouped from one row per conformer into one row per molecule, in the row order of its 2D shard, keeping only embedded conformers.

The protein and antibody collections are converted from their archives in one pass each, natively rather than through these scripts:

- Every row is kept, including identical sequences and partial domains, and shuffled once with seed 42 before sharding.
- Each row is keyed by the MD5 of its uppercase amino-acid sequence, and a paired cell by the MD5 of its two chain keys.
- OAS run metadata is inlined into every row as dictionary columns rather than kept in a side table.
- UniProt and InterPro annotations join AlphaFold rows by accession, and only where the UniProt entry holds the same sequence.

The graphs are built with the `usearch` command line over the finished shards, inserting every row and then refining the whole graph until 99 % of queries find their exact nearest neighbor at a search expansion that keeps queries cheaper than inserts:

- ECFP4 graphs use connectivity 64 and Hamming distance over the 2048-bit fingerprint.
- The PubChem Chamfer graph is built in parts over each molecule's lowest-energy conformer with hydrogens, then absorbed into one graph and refined.
- Needleman-Wunsch graphs index each distinct sequence once, scoring global alignments under BLOSUM62 with linear gaps, and are built in parts across machines.

Once completed, every collection is uploaded to S3 with `cp` rather than `sync`, so a shard of the same size as its predecessor is never skipped:

```sh
for collection in example pubchem gdb13 real alphafold oas-unpaired oas-paired; do
    aws s3 cp --recursive --storage-class INTELLIGENT_TIERING data/$collection/ s3://usearch-molecules/data/$collection/
done
```

### Understanding Fingerprints

Not all per-molecule vectors are interchangeable, and using one where another is required produces silently wrong answers.
ECFP4 and FCFP4 hash local environments into bits with no chemical meaning, for whole-molecule similarity; MACCS and PubChem assign one bit per named pattern, for substructure screening; USRCAT describes each conformer's 3D shape.
A useful pipeline chains all three: a similarity prefilter, a substructure gate, and a shape ranking at the end.
The [`to_fingerprint`](https://github.com/unum-science/USearchMolecules/blob/main/usearchmolecules/to_fingerprint.py) module documents each family and its query semantics.

### What's Persisted

The guiding principle is: __cache what's expensive to recompute, skip what's cheap to reconstruct__.
Conformer generation with ETKDG + MMFF costs 60-600 ms per molecule, while parsing a SMILES string back into a full molecular graph takes under 0.2 ms, so the coordinates are persisted and the molecular graph is not.

__Small-molecule shard__ `parquet/<range>.parquet`, one row per compound:

- `key`, `binary(16)`: an RFC 9562 UUIDv8 packing the standard InChIKey and the CXSMILES enhanced-stereo flag, decoded by `inchikey_of`; sorting the bytes sorts by InChIKey.
- `smiles`, `utf8`, ~50 B: the canonical graph with atoms, bonds, charges and stereochemistry.
- `maccs`, `binary(21)`: MACCS structural keys, 166 bits.
- `ecfp4`, `binary(256)`: extended-connectivity fingerprint, radius 2, 2048 bits.
- `fcfp4`, `binary(256)`: functional-class fingerprint, radius 2, 2048 bits.
- `pubchem`, `binary(111)`: PubChem substructure fingerprint, 881 bits, computed with CDK.

__Conformer shard__ `parquet/<range>.3D.parquet`, up to three conformers per compound, lowest energy first:

- `key`, `binary(16)`: the 2D row's key.
- `n_atoms`, `n_heavy_atoms`, `n_bonds`, `uint16`: counts of the conformed fragment, hydrogens included in `n_atoms`.
- `molecular_weight`, `float32`: exact mass in Daltons.
- `energies`, `list<float32>`: MMFF94 energy per conformer in kcal/mol, ascending.
- `points`, `list<fixed_size_list<float16, 3>>`: every conformer's atoms end to end, centroid-centered.
- `usrcat`, `list<fixed_size_list<float16, 60>>`: one USRCAT shape vector per conformer, whose first 12 values are USR.

Conformer `c` of a row is `points[c * n_atoms : (c + 1) * n_atoms]`, and `len(energies)` counts the row's conformers:

```python
import numpy as np

energies = row["energies"]                                                   # MMFF94 energies, kcal/mol
xyz = np.asarray(row["points"], dtype=np.float16).reshape(len(energies), row["n_atoms"], 3)
usrcat = np.asarray(row["usrcat"], dtype=np.float16)                        # (len(energies), 60); USR is usrcat[:, :12]
```

Coordinates follow the atom order of the largest fragment after `AddHs`, so they line up with a molecule rebuilt from the matching `smiles`, and orientation is not canonicalized.
They are stored as `float16` rather than `bfloat16`, since Parquet encodes the former natively, and the rounding stays under 0.002 Å within 8 Å of the centroid and under 0.008 Å beyond it, well below thermal motion of about 0.1 Å.
Bond topology, atom types, formal charges and stereochemistry are not stored, because the `smiles` column encodes all four and RDKit's `MolFromSmiles` with `AddHs` recovers them losslessly.

__Protein shard__ `alphafold/parquet/<range>.parquet`, one row per AlphaFold DB entry:

- `key`, `binary(16)`: MD5 of the uppercase amino-acid sequence.
- `sequence`, `utf8`: one-letter amino-acid codes.
- `uniprot_accession` and `afdb_entry`, `utf8`: like `A0A0N7M9J1` and `AF-A0A0N7M9J1-F1-model_v6`.
- `fragment`, `uint16`: the fragment number of proteins split into several models.
- `taxonomy_id`, `uint32`: the NCBI taxon, null when unknown.
- `organism` and `domain`, `dictionary`: the taxon's scientific name, and Bacteria, Archaea, Eukaryota or a viral realm.
- `description`, `utf8`: the recommended protein name from UniProt.
- `reviewed`, `bool`: curated in Swiss-Prot rather than predicted in TrEMBL.
- `gene`, `dictionary`: the primary gene name.
- `ec_numbers`, `list<utf8>`: Enzyme Commission classes.
- `go_terms`, `list<uint32>`: Gene Ontology term numbers, `GO:0005524` stored as `5524`.
- `pfam_domains`, `list<utf8>`: Pfam families matched along the sequence, in order of start.

__Structure shard__ `alphafold/parquet/<range>.3D.parquet`, row-aligned with the protein shard:

- `key`, `binary(16)`: the protein row's key.
- `n_residues`, `uint16`: residues, one Cα each.
- `points`, `list<fixed_size_list<float16, 3>>`: Cα coordinates in Å, centered on their centroid.
- `plddt`, `list<uint8>`: AlphaFold's confidence per residue, from 0 to 100.

__Antibody chain shard__ `oas-unpaired/parquet/<range>.parquet`, one row per chain as OAS lists it:

- `key`, `binary(16)`: MD5 of the uppercase amino-acid sequence.
- `sequence` and `cdr3`, `utf8`: the variable domain and its third complementarity-determining region.
- `regions`, `binary(7)`: lengths of FWR1, CDR1, FWR2, CDR2, FWR3, CDR3 and FWR4.
- `v_call`, `d_call`, `j_call` and `locus`, `dictionary`: germline gene assignments and `IGH`, `IGK` or `IGL`.
- `redundancy`, `uint32`: times the sequence was observed in its run.
- `v_identity` and `j_identity`, `float32`: percent identity to the assigned V and J germline genes.
- `study`, `file`, `run`, `author`, `link`, `species`, `age`, `b_source`, `b_type`, `vaccine`, `disease`, `subject`, `longitudinal`, `isotype` and `chain`, `dictionary`: the run's metadata, inlined.

`oas-paired` stores the same per-chain columns twice, suffixed `_heavy` and `_light`, beside `key_heavy` and `key_light`, and keys the cell as `MD5(key_heavy ‖ key_light)`.

### Limits and Caveats

__Two atom counts matter for small molecules, and they are far apart.__
The `smiles` column implies a heavy-atom count, but any force field or 3D pipeline works on the molecule after explicit hydrogens are added, which is roughly twice as large.
Measured on a stratified sample — shards spread across each subset's full index range, a random row group per shard, uniform rows within it, largest fragment only:

| Subset    | Heavy p50 / p95 / p99 / max | With H p50 / p95 / p99 / max | Over 192 Atoms |
| :-------- | :-------------------------- | :--------------------------- | -------------: |
| `pubchem` | 24 / 49 / 79 / 202          | 46 / 90 / 168 / 412          |         0.68 % |
| `gdb13`   | 13 / 13 / 13 / 13           | 27 / 32 / 34 / 39            |         0.00 % |
| `real`    | 27 / 29 / 30 / 34           | 52 / 61 / 64 / 70            |         0.00 % |

PubChem is the only subset with a large-molecule tail, and its p99 is where the two counts diverge most — 79 heavy atoms become 168 with hydrogens.
GDB-13 is a spike rather than a distribution, since every molecule is enumerated at exactly 13 heavy atoms.
Enamine REAL is lead-like and narrow, which is why it never approaches a 192-atom ceiling.
Shards are not ordered by molecule size, so the first N shards of a subset are a fair sample of its size distribution.

__Conformer yield is not 100 %, and the shortfall is a size effect.__
Measured over a stratified sample of the three sources, and over all of `example`:

| Subset    | 3 Conformers |      2 |      1 |      0 |
| :-------- | -----------: | -----: | -----: | -----: |
| `example` |      99.71 % | 0.20 % | 0.09 % |      — |
| `pubchem` |      97.93 % | 0.21 % | 0.09 % | 1.77 % |
| `gdb13`   |     100.00 % |      — |      — |      — |
| `real`    |      99.99 % | 0.00 % | 0.00 % | 0.01 % |

`example` reads better than its sources because it is drawn from molecules whose first conformer succeeded, so it inherits none of their outright failures.
PubChem holds the only meaningful shortfall, and it is concentrated in the tail.
Yield stays near 99.95 % below 75 atoms and falls to 77.6 % in the 175-to-192-atom band.
A conformer that fails to embed or polish is left out, and a molecule refused before embedding, above the 192-atom ceiling, with an element MMFF94 cannot type, or with more rings than the pipeline handles, has empty lists and null counts.

__An embedded conformer's energy is finite, not necessarily low.__
Conformers far above their own molecule's minimum are structures the minimiser never relaxed, and they are still stored.
Filter on the gap to the molecule's lowest conformer rather than on absolute energy, which has no absolute zero in MMFF94:

```python
# energies are ascending, so energies[0] is the molecule's minimum
usable = [c for c, energy in enumerate(row["energies"]) if energy - row["energies"][0] < 50.0]
```

At a 100 kcal/mol gap this affects roughly 7 % of GDB13 molecules, whose small bridged cages ETKDG struggles to embed, against about 1 % of PubChem and Enamine REAL.
Because conformers are stored energy-sorted, the affected conformer is almost always the last of the three, and the first is the safest single geometry to consume.
The conformers of one molecule are diverse enough that deduplication would discard little: median inter-conformer Kabsch RMSD is 1.45 Å for PubChem, 1.08 Å for GDB13 and 1.50 Å for Enamine REAL.

__Multi-fragment SMILES keep a conformer for one fragment only.__
Where a SMILES names more than one disconnected fragment — a salt, a counter-ion, a mixture — the conformer covers only the fragment the pipeline kept, while the `smiles` column retains the whole string.
About 9 % of PubChem rows are affected; GDB13 and Enamine REAL are not.
Reproduce the pipeline's choice with RDKit's `LargestFragmentChooser` under `preferOrganic`, which keeps the organic parent even when a metal cluster carries more atoms:

```python
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
mol = chooser.choose(Chem.MolFromSmiles(row["smiles"]))
Chem.SanitizeMol(mol)   # ring perception is stale once atoms are removed
mol = Chem.AddHs(mol)   # now mol.GetNumAtoms() == row["n_atoms"]
```

__A key names a compound, not a row, so several rows can share one.__
Standard InChI merges different spellings of one structure, mobile-hydrogen tautomers and some charge and protonation forms, and every such row stays in the dataset:

| Subset    |          Rows | Shared-key groups | Rows in groups | Largest |
| :-------- | ------------: | ----------------: | -------------: | ------: |
| `pubchem` |   115'627'256 |           297'068 |        636'297 |      27 |
| `example` |     2'000'000 |                53 |            106 |       2 |
| `gdb13`   |   977'468'267 |        10'190'883 |     21'636'255 |      10 |
| `real`    | 6'039'411'650 |         9'553'952 |     19'470'151 |       … |

Most GDB13 groups are tautomers, which GDB13 enumerates as separate molecules.
PubChem groups mix salts beside their free bases, tautomers, and identical structures registered under several CIDs.
Enamine REAL lists about 93 million structures twice, once with an `&` stereo group and once without, and the enhanced-stereo flag in the key keeps those pairs apart.
Recomputing the full InChI of all 265 million grouped rows found no InChIKey hash collisions.
Nineteen PubChem groups share a key only because the protonation letter saturates at `A` from 13 protons on, so `/p-13` and `/p-15` read alike.
The same key also links subsets, so about 14.6 million PubChem compounds can be found in Enamine REAL by key alone.

__Eleven PubChem rows are not in the dataset.__
Seven are a bare `*` dummy atom, three are mixtures with a `*` fragment, and one is a polymer, and standard InChI cannot name any of them, so they have no key.
Each of the eleven shards that held them lost that one row.

__A sequence key names a sequence, not an entry, and the graphs index each key once.__
AlphaFold DB models 241'070'489 entries but only 220'322'986 distinct sequences, and the largest group of identical sequences spans 11'140 entries across strains and species.
OAS lists the same chain in many runs, donors and studies, and every listing stays a row.
A Needleman-Wunsch match is therefore a sequence, and joining its key back to the shards returns every entry, run or cell that holds it.

__UniProt annotations are present only where they still describe the modeled sequence.__
AlphaFold DB follows UniProt 2025_03, and an annotation is attached only when that release holds the accession with an identical sequence.
Rows whose entry was merged, retired or edited keep null annotation columns rather than borrowing another sequence's function.

__pLDDT qualifies the coordinates beside it.__
Residues below about 50 are commonly disordered, and their Cα positions are a placeholder rather than a structure, so mask them before aligning or measuring distances.

__Antibody rows keep partial and unusual chains.__
About 0.59 % of OAS chains do not tile cleanly into the seven framework and CDR regions, and their `regions` is null while the sequence and calls remain.
Germline identity is missing for 1 chain on V and 7'147 chains on J, where OAS left the field empty.
Neither collection carries antigen labels, so what an antibody binds is not part of the dataset.

## Related Work

PubChem3D is the closest analogue, offering public 3D shape similarity as a standing service rather than a download.
It covers about 85 million molecules, and being a service it answers only the question it was built to answer — you cannot point it at your own corpus, cluster a set you bring, or join shape to anything else.

ZINC-22 and VirtualFlow ship more 3D than that and no search at all.
Retrieving ZINC-22's 4.5 billion structures means moving about a petabyte first, and VirtualFlow's 68.7 billion sit behind a credentialed bucket.
Both are distribution rather than retrieval.

OpenEye's ROCS and FastROCS score Gaussian volume overlap over the whole geometry rather than a descriptor of it, but the work scales with the corpus on every query, the scoring is closed, and the alignment is a local optimum over whatever conformers were generated.
BioSolveIT's SpaceGrow never enumerates, searching billions on one CPU core by growing molecules from synthons, but works only inside a reaction space and cannot cluster an arbitrary set.

For proteins and antibodies, Foldseek and MMseqs2 search AlphaFold DB by structure and sequence through their own prefilters, and the OAS web interface searches its chains by sequence identity.
Neither ships the shards and the graph together, so a result cannot be joined to annotations, conformers or another corpus without re-indexing it.

USearchMolecules keeps the corpus indexed across modalities, so a substructure screen, a similarity sweep and a shape query all run against the same shards.

## Citation

If USearchMolecules helps your research or product, please cite it:

```bibtex
@dataset{Vardanian_USearchMolecules,
  author = {Vardanian, Ash},
  title = {{USearchMolecules: An Open, Searchable Atlas of Chemical and Biological Entities}},
  doi = {10.5281/zenodo.21613663},
  url = {https://github.com/unum-science/USearchMolecules},
  license = {Apache-2.0}
}
```

A machine-readable [`CITATION.cff`](CITATION.cff) is provided at the repository root.
