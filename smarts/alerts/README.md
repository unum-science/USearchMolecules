# Structural Alerts

Reactive-group, PAINS, toxicophore, and screening-deck filters from pharma and academic sources.
Used as a substructure-screen layer on similarity-search results — e.g. _"find ECFP4-similar molecules to X that are PAINS-clean and don't match any AbbVie ALARM-NMR alert."_

Every record is one self-contained JSON object on a single line, validated with RDKit's `Chem.MolFromSmarts`.
Per-source provenance, complexity numbers, and matching gotchas follow below.

## Sources

The _Patterns_ column and parse-rate are computed; everything else is provenance.

### Pharma-Attributed

| File                     | Patterns | Origin                                         | Citation                                               | License               |
| ------------------------ | -------: | ---------------------------------------------- | ------------------------------------------------------ | --------------------- |
| `lilly.jsonl`            |      104 | Eli Lilly inline-SMARTS subset                 | Bruns & Watson 2012, _J. Med. Chem._ 55, 9763          | Apache-2.0            |
| `gsk.jsonl`              |       55 | GSK reactive / undesirable                     | Hann et al. 1999, _J. Chem. Inf. Comput. Sci._ 39, 897 | BSD-3 (RDKit)         |
| `bms.jsonl`              |      180 | BMS HTS triage filters                         | Pearce et al. 2006, _J. Chem. Inf. Model._ 46, 1060    | BSD-3 (RDKit)         |
| `nibr.jsonl`             |      444 | Novartis Screening Deck 2019                   | Schuffenhauer et al. 2020, _J. Med. Chem._ 63, 14425   | BSD-3 (RDKit Contrib) |
| `pfizer_lint.jsonl`      |       57 | Pfizer LINT (community attribution via ChEMBL) | folkloric, no primary DOI                              | BSD-3 (RDKit)         |
| `pfizer_kalgutkar.jsonl` |       20 | Pfizer / reactive metabolites                  | Kalgutkar et al. 2005, _Curr. Drug Metab._ 6, 161      | hand-encoded          |
| `pfizer_stepan.jsonl`    |       21 | Pfizer / idiosyncratic-toxicity alerts         | Stepan et al. 2011, _Chem. Res. Toxicol._ 24, 1345     | hand-encoded          |
| `abbvie_alarm.jsonl`     |       11 | Abbott / AbbVie ALARM-NMR thiol-reactives      | Huth et al. 2005, _J. Am. Chem. Soc._ 127, 217         | hand-encoded          |

### Academic and Aggregated

| File                      | Patterns | Origin                                        | Citation                                           | License            |
| ------------------------- | -------: | --------------------------------------------- | -------------------------------------------------- | ------------------ |
| `kazius.jsonl`            |       16 | Kazius mutagenicity toxicophores              | Kazius et al. 2005, _J. Med. Chem._ 48, 312        | hand-encoded       |
| `brenk.jsonl`             |      105 | DNDi / Dundee unwanted-functional-groups      | Brenk et al. 2008, _ChemMedChem_ 3, 435            | BSD-3 (RDKit)      |
| `pains_a.jsonl`           |       16 | PAINS family A (high-confidence)              | Baell & Holloway 2010, _J. Med. Chem._ 53, 2719    | BSD-3 (RDKit)      |
| `pains_b.jsonl`           |       55 | PAINS family B                                | same                                               | BSD-3 (RDKit)      |
| `pains_c.jsonl`           |      409 | PAINS family C (largest, lowest-confidence)   | same                                               | BSD-3 (RDKit)      |
| `nih.jsonl`               |      180 | NIH MLSMR rejects (RDKit FilterCatalog "NIH") | Jadhav 2010 / Doveston 2014                        | BSD-3 (RDKit)      |
| `zinc.jsonl`              |       50 | UCSF ZINC reactivity screen                   | Irwin et al. 2012, _J. Chem. Inf. Model._ 52, 1757 | BSD-3 (RDKit)      |
| `dundee.jsonl`            |      105 | ChEMBL "Dundee" aggregate                     | ChEMBL `structural_alerts`                         | BSD-3 (RDKit)      |
| `inpharmatica.jsonl`      |       91 | ChEMBL "Inpharmatica"                         | same                                               | BSD-3 (RDKit)      |
| `mlsmr.jsonl`             |      116 | ChEMBL "MLSMR" alt                            | same                                               | BSD-3 (RDKit)      |
| `surechembl.jsonl`        |      166 | ChEMBL "SureChEMBL" non-MedChem-friendly      | same                                               | BSD-3 (RDKit)      |
| `benigni_bossa.jsonl`     |       77 | Benigni-Bossa carc / muta rulebase            | Benigni & Bossa 2008, _Mutat. Res._ 659, 248       | LGPL-2.1 (Toxtree) |
| `iss_ames.jsonl`          |       67 | ISS Ames in-vitro mutagenicity                | via Toxtree `toxtree-ames`                         | LGPL-2.1 (Toxtree) |
| `iss_micronucleus.jsonl`  |       26 | ISS in-vivo micronucleus                      | via Toxtree `toxtree-mic`                          | LGPL-2.1 (Toxtree) |
| `cramer3.jsonl`           |       22 | Cramer Class III high-concern decision tree   | Cramer et al. 1978; via Toxtree                    | LGPL-2.1 (Toxtree) |
| `michael_acceptors.jsonl` |        9 | Michael acceptor structural alerts            | via Toxtree `toxtree-michaelacceptors`             | LGPL-2.1 (Toxtree) |
| `eye_irritation.jsonl`    |       46 | Eye-irritation alert plugin                   | Gerner et al.; via Toxtree                         | LGPL-2.1 (Toxtree) |

## On the LGPL-2.1 Entries

Six files (`benigni_bossa`, `iss_ames`, `iss_micronucleus`, `cramer3`, `michael_acceptors`, `eye_irritation`) were derived from Toxtree (https://toxtree.sourceforge.net/) by regex-extracting SMARTS strings from the project's Java rule-class source files.
Toxtree is __LGPL-2.1__, and each derived record carries `"license": "LGPL-2.1"` so the constraint travels with the data even if files are reorganized.

What this license actually requires for a downstream consumer of these `.jsonl` files:

- __Attribution.__ Keep the `license`, `citation`, and `rule_id` fields in each record so the chain of provenance to Toxtree's authors is preserved.
- __Modifications stay LGPL.__ If you edit any LGPL-2.1 record (change a SMARTS, a name, etc.) and redistribute the modified file, the modified file remains LGPL-2.1.
- __Linking is not affected.__ LGPL is fundamentally library-linking-friendly.
  A program that loads these `.jsonl` files at runtime to do substructure screening does not become LGPL-encumbered itself, and does not need to be open-source.
  The license only constrains _redistribution of the data files themselves_, not the code that reads them.

What it does __not__ limit:

- Reading and matching SMARTS in commercial pipelines.
  USearchMolecules is free to call `Chem.SubstructMatch` on these patterns inside a closed-source product.
- Combining LGPL-tagged records with permissive-licensed ones in the same in-memory catalog at runtime.
- Reporting that a molecule matches an LGPL-tagged alert.
  Match outcomes are facts, not derivative works.

The earlier draft put these files under a `_lgpl/` subdirectory to flag them.
That was unnecessary because the per-record `license` field already carries the constraint and travels with the data; the directory split has been removed.

## Caveats

__BMS and NIH overlap heavily.__
RDKit's `chembl_bms.in` and `nih.in` both contain 180 named patterns, and __all 180 names are identical between them__.
Of the underlying SMARTS strings, __176 of the 180 (name, SMARTS) pairs are exact duplicates__, and __174 of the 178 unique SMARTS in each set appear in both__.
RDKit ships them as two named catalogs but the data is essentially the same set with different attributions.
A naive sum "BMS + NIH = 360" double-counts almost everything; the true union is ≈186 unique SMARTS.

__Lilly is a partial extract.__
Lilly's `.qry` files are written in a custom S-expression-style query language, not SMARTS.
We extract only the 131 queries that include an inline `(A C smarts "...")` field; of those, 27 use Lilly-only predicates (`G0`, `T1`, `/IWfss1c`, etc.) that RDKit cannot parse, leaving __104 records__.
The other 162 `.qry` files use Lilly's structured atom/bond constraint format and are not represented here; converting them to SMARTS faithfully would require implementing the Lilly query language.

__Hand-encoded entries are conservative.__
Kazius, Pfizer-Kalgutkar, Pfizer-Stepan, and AbbVie-ALARM had no machine-readable mirrors at the time of building, so each pattern was transcribed by hand from the primary paper's figures and tables, validated against RDKit, and spot-checked against a positive-control SMILES.
The counts (16, 20, 21, 11) are deliberately smaller than the comprehensive prose lists in those papers because we kept only patterns we were confident encoded the cited chemistry.

__Pharmas with no public SMARTS catalog__ — verified across rd_filters, ChEMBL `structural_alerts`, ToxAlerts, OpenEye FILTER, RDKit Contrib, Toxtree, useful_rdkit_utils, and primary supplementary info: Roche, AstraZeneca, Bayer, Boehringer Ingelheim, Merck & Co., Merck KGaA, J&J / Janssen, Sanofi, Takeda, Genentech, Amgen, Gilead, Biogen, Regeneron, Schering-Plough.
Their alert lists are internal.
Vertex/Schering-Plough's REOS (Walters & Murcko 2002) was never released as a canonical SMARTS file; the closest public artifact is the bundle in Pat Walters' `rd_filters`, which is a re-bundle of the public sets we already have.

## JSONL Record Schema

Every `.jsonl` line is one self-contained JSON object.
The required core is three fields; everything else is optional and only present when the source provides the data.

__Core (always present)__:

- `smarts`: the SMARTS string.
- `name`: human-readable label, sourced from the original catalog or the paper's figure / table caption.
- `source`: short identifier matching the file basename (e.g. `gsk`, `nibr`, `pfizer_kalgutkar`).

__Provenance__:

- `citation`: a short citation key (e.g. `Bruns2012`, `Schuffenhauer2020`).
- `license`: SPDX-style identifier, or `hand-encoded` for primary-source transcriptions.
- `family`: sub-grouping within the source (e.g. PAINS `A` / `B` / `C`, NIBR's `NIBR_Screeningdeck_2019`).

__Severity (when graded)__:

- `severity`: textual category (`REJECT`, `DEMERIT`, `EXCLUDE`, `FLAG`, `OTHER`).
- `severity_score`: integer score where the source provides one (NIBR uses 1 / 2; Lilly is binary in the inline-SMARTS subset).
- `covalent`: boolean — does the alert describe a covalent reactive group? (NIBR only.)
- `min_count`: minimum match count required to trigger (NIBR only).

__Examples (when bundled)__:

- `examples`: array of `{smiles, pubchem_cid?}` objects — example matching molecules.
  NIBR ships 5 PubChem CIDs per pattern; hand-encoded sources include 1–2 canonical positives.

__Toxtree-derived__:

- `rule_id`, `rule_title`, `description`: the originating rule's ID, full title, and prose explanation from the Java source.

__Lilly-only__:

- `qry_file`: the originating `.qry` filename in the Lilly repo.

__Source-specific extras__:

- NIBR: `original_set`, `num_matches_pubchem_2018`.
- Hand-encoded sources: `paper_section` referencing the figure or table number.

## Per-Source Complexity

Numbers are __median / p95 / max__ across all patterns in the source.

_Recur_ counts `$()` sub-patterns; _Depth max_ is the deepest `$()` nesting; _OR-arity_ is the widest `[A,B,C,...]` cardinality.
_Chirality_ is `@` inside `[...]` summed across the source; _Ring-bond_ is `@` outside `[...]` (e.g. `!@`, `-@`); _Multi-frag_ is `.` outside `[...]`; _Ring closures_ is digit / `%NN` ring labels.

| Source              | Patterns | Parse OK | Recur (med/p95/max) | Depth max | OR-arity (med/p95/max) | Chirality | Ring-bond | Multi-frag | Ring closures | Atoms (med/p95/max) |
| ------------------- | -------: | -------: | :-----------------: | :-------: | :--------------------: | --------: | --------: | ---------: | ------------: | :-----------------: |
| `lilly`             |      104 |  104/104 |      0 / 0 / 8      |     1     |       0 / 3 / 3        |         0 |         0 |          0 |           102 |     5 / 13 / 21     |
| `gsk`               |       55 |    55/55 |      0 / 0 / 0      |     0     |       0 / 2 / 4        |         0 |         0 |          3 |            38 |     5 / 12 / 21     |
| `bms`               |      180 |  180/180 |     0 / 13 / 88     |     2     |       1 / 3 / 11       |         0 |         9 |         37 |           159 |     5 / 14 / 24     |
| `nibr`              |      444 |  444/444 |     0 / 16 / 36     |     2     |       0 / 3 / 6        |       433 |       195 |         12 |          1013 |     6 / 17 / 51     |
| `pfizer_lint`       |       57 |    57/57 |      0 / 2 / 5      |     1     |       0 / 3 / 5        |         0 |        14 |          5 |            28 |     4 / 10 / 13     |
| `pfizer_kalgutkar`  |       20 |    20/20 |      0 / 2 / 3      |     1     |       0 / 2 / 3        |         0 |         0 |          0 |            24 |     4 / 9 / 10      |
| `pfizer_stepan`     |       21 |    21/21 |      0 / 1 / 2      |     1     |       0 / 1 / 2        |         0 |         0 |          0 |            26 |      5 / 8 / 9      |
| `abbvie_alarm`      |       11 |    11/11 |      0 / 1 / 1      |     1     |       1 / 3 / 3        |         0 |         0 |          0 |             8 |      4 / 8 / 8      |
| `kazius`            |       16 |    16/16 |      0 / 3 / 3      |     1     |       0 / 2 / 3        |         0 |         0 |          0 |             4 |      3 / 6 / 7      |
| `brenk`             |      105 |  105/105 |      0 / 0 / 6      |     1     |       0 / 3 / 9        |         0 |         7 |          2 |           104 |     4 / 14 / 19     |
| `pains_a`           |       16 |    16/16 |      0 / 3 / 5      |     1     |       0 / 0 / 0        |         0 |         1 |          0 |            32 |    10 / 15 / 16     |
| `pains_b`           |       55 |    55/55 |      0 / 5 / 8      |     1     |       0 / 0 / 3        |         0 |         9 |          0 |           142 |    12 / 20 / 24     |
| `pains_c`           |      409 |  409/409 |     0 / 2 / 12      |     2     |       0 / 0 / 3        |         0 |        11 |          0 |          1564 |    19 / 33 / 49     |
| `nih`               |      180 |  180/180 |     0 / 13 / 88     |     2     |       1 / 3 / 11       |         0 |         9 |         37 |           159 |     5 / 14 / 24     |
| `zinc`              |       50 |    50/50 |      0 / 0 / 0      |     0     |       0 / 3 / 4        |         0 |         0 |          0 |             6 |     4 / 10 / 14     |
| `dundee`            |      105 |  105/105 |      0 / 0 / 6      |     1     |       0 / 3 / 9        |         0 |         7 |          2 |           104 |     4 / 14 / 19     |
| `inpharmatica`      |       91 |    91/91 |      0 / 2 / 6      |     2     |       0 / 4 / 14       |        18 |        24 |          0 |            42 |     4 / 11 / 19     |
| `mlsmr`             |      116 |  116/116 |      0 / 0 / 6      |     1     |       0 / 3 / 4        |         0 |        25 |          6 |           110 |     5 / 16 / 21     |
| `surechembl`        |      166 |  166/166 |      0 / 1 / 7      |     2     |       0 / 3 / 72       |         2 |        29 |          2 |           170 |     5 / 16 / 19     |
| `benigni_bossa`     |       77 |    77/77 |      0 / 4 / 5      |     2     |       0 / 3 / 4        |         3 |         4 |          0 |           112 |     6 / 14 / 22     |
| `iss_ames`          |       67 |    67/67 |      0 / 4 / 8      |     2     |       0 / 3 / 3        |         3 |         2 |          1 |           156 |     7 / 21 / 23     |
| `iss_micronucleus`  |       26 |    26/26 |      0 / 3 / 4      |     1     |       0 / 3 / 3        |         0 |         0 |          0 |            16 |     4 / 11 / 13     |
| `cramer3`           |       22 |    22/22 |     0 / 4 / 14      |     1     |       0 / 3 / 6        |         0 |         5 |          0 |             2 |      4 / 6 / 6      |
| `eye_irritation`    |       46 |    46/46 |     0 / 4 / 19      |     1     |       1 / 3 / 5        |         0 |         0 |          1 |            70 |     8 / 14 / 20     |
| `michael_acceptors` |        9 |      9/9 |      0 / 0 / 0      |     0     |       0 / 5 / 5        |         0 |         0 |          0 |             6 |     5 / 12 / 12     |

## Extreme Patterns

The single largest pattern in each category, across all 25 sources.
SMARTS strings truncated to 100 characters; full strings are in the source `.jsonl` files.

| Metric                                |    Value | Source        | Name                          | SMARTS                                                                  |
| ------------------------------------- | -------: | ------------- | ----------------------------- | ----------------------------------------------------------------------- |
| Deepest `$()` nesting                 |    __2__ | `bms` / `nih` | `halo_acrylate`               | `[$([C;H2]),$([C&H1;$(C-F)]),$([C&H1;$(C-Cl)]),$([C&H1;$(C-Br)]),...`   |
| Most `$()` sub-patterns in one SMARTS |   __88__ | `bms` / `nih` | `tris_activated_aryl_ester`   | `[$(O=[C,S]Oc1a([$(S(=O)(=O)),F,$(C(F)(F)(F)),$(C#N),...])...`          |
| Widest `[A,B,...]` OR-list            |   __72__ | `surechembl`  | `Undesirable_Elements_Salts`  | `[Ac,Ag,Am,Ar,As,At,Au,Ba,Be,Bi,Bk,Cd,Ce,Cf,Cm,Cr,Cs,Dy,...]`           |
| Longest raw SMARTS string             | __1255__ | `nibr`        | `peptide_analogue_10`         | `N[C!r3!r4!r5!r6!r7!r8!r9](=O)CCN[...]...`                              |
| Largest parsed query (atoms)          |   __51__ | `nibr`        | `cyclosporine`                | `O=C(NCC(NC(C(C(CC=,-C)C)=,-[O,N])C(NCC(NCC(NCC(NCC(NCC(...)`           |
| Most ring-closure labels              |   __24__ | `nibr`        | `enolether_enolester`         | `[C;!$([#6]1=[#6]-[#6,#7&v3]=,:[#6,#7&v3]-[#6](=...)-1);!$(...);...]`   |
| Most ring-bond annotations            |   __17__ | `nibr`        | `polycyclic_systems_18_atoms` | `a@;:a@;:a@;:a@;:a@;:a@;:a@;:a@;:a@;:a@;:a@;:a@;:a@;:a@;:a@;:a@;:a@;:a` |
| Most `.` fragment delimiters          |    __6__ | `bms`         | `gte_7_aliphatic_OH`          | `C[O;D1].C[O;D1].C[O;D1].C[O;D1].C[O;D1].C[O;D1].C[O;D1]`               |

## On Depth and Directionality

Two complexity features that one might _expect_ to find in a SMARTS catalog turn out to be effectively absent across all 2'448 patterns.

__Maximum recursion depth across the whole catalog is 2.__
This holds even after surveying all 14 Toxtree plugins (only the 6 useful ones are included above; the other 8 — functional_groups, biodegradation, sicret, verhaar, kroes, dnabinding, proteinbinding, moa — were either out-of-scope or stored SMARTS in indirect form).
A handful of patterns have `$()` _width_ in the dozens (BMS / NIH `tris_activated_aryl_ester` packs 88 sub-patterns into a single SMARTS, eye_irritation has one with 19, cramer3 has one with 14), but the nesting never goes deeper than two levels.
This is intrinsic to the alert-cataloging task: deeper patterns are too brittle and over-specific to be generally useful.

__Directional bonds (`/`, `\`) appear zero times across the whole catalog.__
Structural alerts care about which atoms and bonds are present, not which geometric isomer is present.
A reactive functional group is reactive whether it is _E_ or _Z_; a Michael acceptor doesn't care.
Catalogs that DO use directional SMARTS — reaction templates, stereoisomer-specific binding queries, drug-substance specifications — are not structural-alert sets and are out of scope here.

If a downstream pipeline needs deeper-recursion or directional patterns, it should look at __reaction SMARTS__ (RDKit reaction templates use both heavily) or __pharmacophore feature definitions__ (`.fdef` files), not at alert catalogs.

## What to Look Out for When Matching

These are the cost drivers and correctness pitfalls a substructure matcher will hit, drawn from the metrics above.

__Recursive width, not depth, is the real cost driver.__
Maximum recursion depth across all 2'448 patterns is only 2, but BMS / NIH's `tris_activated_aryl_ester` packs 88 `$()` sub-patterns into a single SMARTS, mostly as alternatives in a wide OR-list of substituents.
Width is just as expensive as depth for matchers that don't share state across alternatives.
Plan for the matcher to memoize sub-pattern results per (atom, sub-pattern) pair.

__Wide OR-lists__ (`[A,B,C,...]`).
Median arity is small (0–1) but max arity reaches 11 in BMS / NIH and 72 in SureChEMBL's element-blacklist.
Naive linear comparison vs. bitset representation of allowed atomic numbers matters for hot loops.

__Ring queries and ring-bond annotations__ (`r`, `R`, `x`, `X` inside `[...]`; `@` outside `[...]`).
Both require SSSR (smallest-set-of-smallest-rings) on the target before matching.
NIBR uses ring-bond annotations heavily — 195 occurrences across 444 patterns, often as `!@` to require a non-ring bond at a specific position.
Skipping ring perception breaks these silently.

__Ring-closure labels__ (digit / `%NN` outside `[...]`).
PAINS-C alone has 1'564 ring-closure labels across 409 patterns — they are the dominant complexity feature there, not recursion.
Each label introduces a constraint that two atoms must be bonded, which interacts with the matcher's atom-pair enumeration.

__Chirality__ (`@` inside `[...]`).
Concentrated in NIBR (433 chirality specifiers) and used lightly by Inpharmatica, MLSMR, SureChEMBL, and the Toxtree-derived sets.
Requires chirality perception on the target SMILES; targets without explicit stereodescriptors will silently miss hits.

__Multi-fragment patterns__ (`.` outside `[...]`).
The query is disconnected; the matcher must enumerate fragment assignments.
BMS's `gte_7_aliphatic_OH` (7 disconnected hydroxyl fragments) is the extreme; BMS / NIH each have 37 multi-fragment patterns, NIBR has 12.
Bitset-fingerprint screens cannot reason across the `.` boundary.

__Lilly extensions__ (`G0`, `T1`, `/IWfss1c`, etc.).
Lilly's native query language uses predicates RDKit doesn't recognize.
27 of 131 inline-SMARTS rules were dropped during extraction.
A faithful Lilly matcher would need a separate query-language implementation.

## Implications for USearchMolecules

USearchMolecules indexes molecules by similarity fingerprints (MACCS, ECFP4, FCFP4, PubChem dictionary).
For substructure-screen workflows the relevant question is which patterns admit a fingerprint pre-filter and which do not.

__Dictionary fingerprints__ (MACCS, PubChem) are valid prefilters for the _non-recursive, single-fragment, non-stereo_ subset.
This covers roughly: GSK (55), Brenk (105), ZINC (50), Dundee (105), Lilly (104, ignoring the dropped extensions), and the manually-encoded Pfizer / AbbVie / Kazius alerts (~70 combined).
For these, a bitset-subset test (`fp(query) & fp(target) == fp(query)`) is a sound over-approximating screen.

__The recursive / ring-annotated majority__ — PAINS-C (409), NIH (180, ≈ BMS), most of NIBR (444), and the deeper ChEMBL aggregates — cannot be soundly screened by dictionary fingerprints.
The recursive sub-patterns describe _atom environments_ rather than connected substructures, and the ring-bond annotations encode topology that the precomputed fingerprint bits don't preserve.
These should go directly to RDKit's `Chem.SubstructMatch`, optionally with ECFP4 as a heuristic (lossy) prefilter — accepting that ECFP4 prefiltering can miss legitimate hits.

__Multi-fragment patterns__ require special handling regardless of catalog.
Fingerprint screening cannot reason across the `.` boundary; the matcher must enumerate fragment assignments on the target molecule.

__Stereo-bearing patterns__ (NIBR, Inpharmatica, MLSMR, SureChEMBL, the Toxtree-derived sets) require the target SMILES to carry chirality.
USearchMolecules' upstream pipelines should not strip stereodescriptors before substructure screening if these alerts matter.
