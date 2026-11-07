# Dataset catalog

Snapshot: graph-rag input commit `a204c52b63a5140415a80576971842396e224590`, checked 2026-11-07. A current pointer identifies a version; it does not prove that every file, external download or scientific claim has passed acceptance.

## Starter v1

The standalone starter selects activity records whose protein and compound endpoints both exist and whose DOI or PMID is present. Endpoint tables contain only referenced entities. Selection is deterministic from the existing ChEMBL 36 sample; it is not random or representative sampling.

| Table | Current candidate rows | Meaning | Source / license |
|---|---:|---|---|
| `data/proteins.tsv` | 87 | Distinct referenced protein IDs | UniProtKB/Swiss-Prot snapshot, fetched 2025-10-26; CC BY 4.0 |
| `data/molecules.tsv` | 102 | Distinct referenced InChIKeys | ChEMBL `chembl_36.db` identity records; CC BY-SA 3.0 |
| `data/activity_evidence.tsv` | 359 | Activity records, not unique entity pairs or publications | ChEMBL 36 sample and existing project-derived IDs/normalization; CC BY-SA 3.0 |

The actual built manifest is the authority for a package's counts, field lists and hashes. Input digests are:

| Input in this repository | SHA256 |
|---|---|
| `protein/release/external/protein-v6/tables/sequence/protein_master_v6_clean.tsv` | `e2e7b29d22216f6d0b9111cfcbafde8f9ab33ff8ca5960759acb772a1f324217` |
| `molecule/release/external/molecule-v1/tables/identifiers/molecule_structure_identifiers_v1.tsv` | `0533dfcb74d723dd13553084d6fc01ee503c950cb96a83e07f1fca18a6c22396` |
| `molecule/release/external/molecule-v1/tables/activity/molecule_activity_evidence_v2.sample.tsv` | `831d943dd65b0dcc1365bcc660d31c632d9578489638d28da95ae535222d4383` |

The package exports minimal identity, activity and provenance fields. It excludes structure strings from other sources, heuristic confidence scores and full-table aggregate statistics. See [the examples](../examples/README.md) and the package's `NOTICE.md` for the complete field/source map.

## Five-domain inventory

Counts below have different units and must not be added together. **Manifest-reported** means a stored manifest declares the count; **counted** means rows were read from the named TSV in the input snapshot. Neither designation proves scientific correctness or remote availability.

| Domain / version | Representative assets and scale | Distribution and current boundary |
|---|---|---|
| Protein / `protein-v6` | `protein_master_v6_clean.tsv`: 19,135 protein rows; `protein_edges.tsv`: 884,555 edge rows; `pdb_structures_v1.tsv`: 109,531 structure rows. These are manifest-reported full-table counts, separate from the starter. | Repository snapshot. Starter fields are source-checked; complete protein annotations and all legacy file downloads are outside starter acceptance. |
| RNA / `rna-v1` | Public manifest lists expression evidence and canonical ID mapping, plus explicit `.sample.tsv` assets and an RNA master sample. The public manifest does not report table row counts. | Repository package plus [legacy `rna-l1-v1` release](https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1). Release identity and five asset names were observed through GitHub API; full legacy bytes have not been re-downloaded or accepted here. |
| Molecule / `molecule-v1` | Identity TSV: 14,775 rows; activity evidence sample: 5,000 rows (both counted). Other structure, property, PK/toxicity and semantic assets are in the public manifest. | Activity files marked `.sample.tsv` are samples. The 359 starter records are a further connected subset; they do not measure full-domain coverage. Other fields need their own source/license checks. |
| Interaction / `interaction-l2-v1` | Public manifest lists PPI/PSI/RPI result assets, including `edges_ppi_v1.tsv.zst`; it does not give a common row-count field. | [Legacy release](https://github.com/hazelian0619/protian-entity/releases/tag/interaction-l2-v1) is recorded in the delivery report. Some large results use chunks. Complete chunk-set and whole-file integrity must be established before merging; historical reports alone do not prove present download success. |
| Disease / `disease-v1` | Disease nodes: 29,868 rows; phenotype nodes: 19,836 rows; disease–phenotype relations: 269,518 rows (counted). Other result tables include disease–gene/target, disease–drug and biomarkers. | `product.json` says `in-development`. The stored manifest contains machine-local absolute paths and `git_commit: null`; it is not a portable, accepted release manifest. No stable download promise is made here. |

## Canonical local entry points

| Domain | Version pointer | Public package | Manifest |
|---|---|---|---|
| Protein | [current.json](../protein/products/protein/current.json) | [protein-v6](../protein/release/external/protein-v6/) | [protein public manifest](../protein/release/external/protein-v6/manifests/protein_public_manifest_v1.json) |
| RNA | [current.json](../rna/products/rna/current.json) | [rna-v1](../rna/release/external/rna-v1/) | [RNA public manifest](../rna/release/external/rna-v1/manifests/rna_public_manifest_v1.json) |
| Molecule | [current.json](../molecule/products/molecule/current.json) | [molecule-v1](../molecule/release/external/molecule-v1/) | [molecule public manifest](../molecule/release/external/molecule-v1/manifests/molecule_public_manifest_v1.json) |
| Interaction | [current.json](../interaction/products/interaction/current.json) | [interaction-l2-v1](../interaction/release/external/interaction-l2-v1/) | [interaction public manifest](../interaction/release/external/interaction-l2-v1/manifests/interaction_public_manifest_v1.json) |
| Disease | [current.json](../disease/products/disease/current.json) | [disease-v1](../disease/release/external/disease-v1/) | [disease manifest](../disease/release/external/disease-v1/manifests/disease_v1_manifest.json) |

Protein uses a `latest` object. The other four pointers currently use `current/path`; the old index builder does not fully adapt them. This catalog explains the snapshot and does not create a second machine index or override those pointers.

## Interpreting evidence and licenses

The starter preserves ChEMBL record IDs, metrics, comparison symbols, original units/values, existing normalized values, and DOI/PMID. Two records with different values may refer to the same publication. Their experimental conditions have not been shown equivalent; examples explicitly mark comparability unknown.

The starter's field attribution is based on [UniProt's license statement](https://rest.uniprot.org/help/license) and [ChEMBL's license statement](https://www.ebi.ac.uk/chembldb/), checked 2026-11-07. UniProt fields remain CC BY 4.0; ChEMBL-derived tables remain CC BY-SA 3.0. Project code is MIT. This field-specific check does not clear all historical RNA, interaction, disease, structure, annotation or third-party data for redistribution.

## Reproducing the catalog check

Read each canonical pointer and its named public manifest at the input commit. Protein counts above come from `assets[].rows`; molecular/disease counts are TSV rows excluding the header. RNA and interaction are left uncounted where no compatible count was verified. Upstream data versions and local release labels are recorded separately.

For a current package, run the builder and examples described in the repository README, then inspect its `manifest.json`. Release tests cover their declared cases; they do not certify that every historical domain matrix passes or that every legacy external asset is available.
