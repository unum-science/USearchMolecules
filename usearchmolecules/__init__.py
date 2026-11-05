"""USearchMolecules - an open, searchable atlas of chemical and biological entities.

This package reads and builds the atlas: small molecules with their fingerprints and
3D conformers, AlphaFold proteins with their Cα traces, and OAS antibodies with their
germline calls, each keyed permanently and searched through USearch graphs.
"""

from importlib.metadata import PackageNotFoundError, version

from usearchmolecules.dataset import (
    BATCH_SIZE,
    SHARD_SIZE,
    FingerprintedDataset,
    FingerprintedEntry,
    FingerprintedShard,
    SequenceDataset,
    inchikey_of,
    pair_key,
    sequence_key,
)
from usearchmolecules.to_fingerprint import (
    FingerprintShape,
    shape_maccs,
    shape_mixed,
    smiles_to_maccs_ecfp4_fcfp4,
)

try:
    __version__ = version("usearchmolecules")
except PackageNotFoundError:  # running from a source checkout without an install
    __version__ = "0.0.0"

__all__ = [
    "BATCH_SIZE",
    "SHARD_SIZE",
    "FingerprintShape",
    "FingerprintedDataset",
    "FingerprintedEntry",
    "FingerprintedShard",
    "SequenceDataset",
    "inchikey_of",
    "pair_key",
    "sequence_key",
    "shape_maccs",
    "shape_mixed",
    "smiles_to_maccs_ecfp4_fcfp4",
]
