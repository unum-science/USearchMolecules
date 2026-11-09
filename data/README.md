---
license: apache-2.0
pretty_name: USearch Molecules
size_categories:
- 10B<n<100B
task_categories:
- feature-extraction
tags:
- chemistry
- biology
- molecules
- proteins
- antibodies
- drug-discovery
- cheminformatics
- smiles
- conformers
- molecular-fingerprints
- similarity-search
- vector-search
- 3d
- tabular
- pandas
- dask
- usearch
- rdkit
- numkong
configs:
- config_name: example
  data_files: data/example/parquet/??????????-??????????.parquet
- config_name: example-3d
  data_files: data/example/parquet/*.3D.parquet
- config_name: pubchem
  data_files: data/pubchem/parquet/??????????-??????????.parquet
- config_name: pubchem-3d
  data_files: data/pubchem/parquet/*.3D.parquet
- config_name: gdb13
  data_files: data/gdb13/parquet/??????????-??????????.parquet
- config_name: gdb13-3d
  data_files: data/gdb13/parquet/*.3D.parquet
- config_name: real-3d
  data_files: data/real/parquet/*.3D.parquet
---

# USearch Molecules

An open, searchable atlas of chemical and biological entities, from small molecules to the proteins they target and the antibodies that bind them.
Every entity is one row under a permanent key derived from the entity itself, with its representations and search graphs added beside the same rows over time.

Start with `example`: 2 million molecules drawn from all three sources, 4 GB, carrying everything the larger subsets do.

| Config    |     Molecules | Source                                   |
| :-------- | ------------: | :--------------------------------------- |
| `example` |     2'000'000 | drawn from the three below, not additive |
| `pubchem` |   115'627'256 | NCBI PubChem                             |
| `gdb13`   |   977'468'267 | University of Bern GDB13                 |
| `real`    | 6'039'411'650 | Enamine REAL                             |

Each subset has a `-3d` companion config holding up to three conformers per molecule, row for row with its fingerprints.

Enamine REAL is the exception on this mirror: its conformers are here and still growing, while its fingerprint shards live on the S3 mirrors rather than the Hub.

## Loading

The shards are plain Parquet, and `pyarrow`, `pandas`, `polars` and `dask` read both families directly:

```python
import pyarrow.parquet as pq

molecules = pq.read_table("data/example/parquet/0000000000-0001000000.parquet")
geometry = pq.read_table("data/example/parquet/0000000000-0001000000.3D.parquet")
```

The `datasets` library loads both families, with the key and the fingerprints as `fixed_size_binary` values:

```python
from datasets import load_dataset

geometry = load_dataset("unum-cloud/USearchMolecules", "example-3d", split="train")
```

Searching the fingerprints is what [USearch](https://github.com/unum-cloud/USearch) is for, rebuilding chemistry from a SMILES string is [RDKit](https://www.rdkit.org), and the Kabsch and Umeyama kernels for comparing conformers come from [NumKong](https://github.com/ashvardanian/NumKong).

## Columns

Fingerprint configs carry one row per molecule:

| Column    | Type          | Description                                             |
| :-------- | :------------ | :------------------------------------------------------ |
| `key`     | `binary(16)`  | Compound identity, an InChIKey-derived UUIDv8           |
| `smiles`  | `utf8`        | Canonical graph: atoms, bonds, charges, stereochemistry |
| `maccs`   | `binary(21)`  | MACCS structural keys, 166 bits                         |
| `ecfp4`   | `binary(256)` | Extended-connectivity fingerprint, radius 2, 2048 bits  |
| `fcfp4`   | `binary(256)` | Functional-class fingerprint, radius 2, 2048 bits       |
| `pubchem` | `binary(111)` | PubChem substructure fingerprint, 881 bits              |

The `-3d` configs carry one row per molecule, in the same order as its fingerprint shard, with up to three conformers, lowest energy first:

| Column                                | Type                                 | Description                                    |
| :------------------------------------ | :----------------------------------- | :--------------------------------------------- |
| `key`                                 | `binary(16)`                         | The fingerprint row's key                      |
| `n_atoms`, `n_heavy_atoms`, `n_bonds` | `uint16`                             | Counts, with and without hydrogens             |
| `molecular_weight`                    | `float32`                            | Exact mass in Daltons                          |
| `energies`                            | `list<float32>`                      | MMFF94 energy per conformer in kcal/mol        |
| `points`                              | `list<fixed_size_list<float16, 3>>`  | All conformers end to end, centroid-centered   |
| `usrcat`                              | `list<fixed_size_list<float16, 60>>` | Shape vector per conformer; USR is first 12    |

Conformer `c` is `points[c * n_atoms : (c + 1) * n_atoms]`, and only embedded conformers are stored.
Conformers come from ETKDG embedding over RDKit's experimental torsion preferences, followed by MMFF94 relaxation.
The key packs the standard InChIKey and the CXSMILES enhanced-stereo groups, and decodes back to the InChIKey that PubChem or ChEMBL holds; the [GitHub repository](https://github.com/unum-science/USearchMolecules) documents the bit layout.

## Caveats

Conformer yield is not complete: PubChem reaches three conformers for 97.9 % of molecules, and the shortfall is concentrated above 100 atoms.
An embedded conformer's energy is finite rather than necessarily low, so filter on the gap to a molecule's own lowest conformer rather than on absolute energy.

Where a SMILES names more than one fragment, the geometry covers only the largest, while the `smiles` column keeps the whole string.
Reproduce the choice with RDKit's `LargestFragmentChooser` under `preferOrganic`.

A key names a compound, not a row: tautomers, salts and repeated structures share one, so several rows can carry the same key.
Eleven PubChem rows that standard InChI cannot name, `*` dummy atoms and one polymer, are not part of the dataset.

## More

Pre-built USearch indexes, mirror choices, the SMARTS catalogs and the full methodology live in the [GitHub repository](https://github.com/unum-science/USearchMolecules).

## Citation

```bibtex
@dataset{Vardanian_USearchMolecules,
  author = {Vardanian, Ash},
  title = {{USearchMolecules: An Open, Searchable Atlas of Chemical and Biological Entities}},
  doi = {10.5281/zenodo.21613663},
  url = {https://github.com/unum-science/USearchMolecules},
  license = {Apache-2.0}
}
```
