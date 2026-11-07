# Architecture

## Repository Role

**graph-rag is the public-facing dataset entry point.** It contains a standalone starter candidate and a catalog of existing results. Starter publication and remote verification are pending. Some historical large assets remain in protian-entity Releases; their links and verification boundaries are explicit in the catalog.

Sister repository [protian-entity](https://github.com/hazelian0619/protian-entity) is the development engine where pipelines are built and data is produced. Only curated, QA-passed results are promoted to graph-rag for public delivery.

## Five Domains

| Domain | Scope | Latest version |
|---|---|---|
| `protein/` | Human protein L1 entity tables | protein-v6 |
| `rna/` | Human RNA L1/L2 entity tables | rna-v1 |
| `molecule/` | Small molecule L1/L2 entity tables | molecule-v1 |
| `interaction/` | Cross-entity L2 edges (PPI / PSI / RPI) | interaction-l2-v1 |
| `disease/` | Disease/phenotype relations; in development | disease-v1 |

## Per-Domain Layout

Each domain follows the same four-section structure:

```
<domain>/
  data/           # gitkeep + sample files only (raw/intermediate data is NOT stored here)
  pipelines/      # contract schemas, build scripts, QA reports
    <pipeline>/
      contracts/  # JSON data contracts (column schemas, validation rules)
      scripts/    # Python build and QA scripts
      reports/    # validation gates, QA metrics, manifests
  products/       # product definition (product.json) and version pointer (current.json)
  release/        # public and internal release packages
    external/     # curated public deliverables (tables, manifests, reports)
    internal/     # process evidence and intermediate QA
    schema/       # domain-specific schema (if any)
```

## Shared Layer (repository root)

```
schemas/           # shared JSON schemas (release index schema, etc.)
scripts/           # release tooling (build_release_index, validate_release_index, etc.)
tools/             # table validation utilities (kg_validate_table, kg_make_manifest)
tests/             # shared test suite
docs/              # repository standards and documentation
.github/workflows/ # CI: data-qa + release-metadata (coverage described below)
examples/          # standalone starter examples and field provenance
```

## Data Promotion Path

```
protian-entity (dev)                graph-rag (public)
  pipeline produces                   curated results committed to
  QA gates pass            ------>    <domain>/release/external/<version>/
  manifest generated                  product pointer updated (current.json)
                                      release index rebuilt
```

Rules:
- Large tables (L1/L2 result data) live only in `release/external/`.
- `data/` directories hold `.gitkeep` and sample files, not full datasets.
- Super-sized assets use GitHub Releases; the repo holds metadata only.

## Starter Delivery

`scripts/build_starter_dataset.py` selects connected protein/compound activity records from existing tables and writes a small package under ignored `build/`. Its manifest records input hashes, selection, fields, licenses and exported file hashes. `examples/run_examples.py` reads that manifest and the three tables offline; it has no dependency on repository modules or online services.

The planned public asset pair is a versioned `starter-v1.tar.gz` and `SHA256SUMS.txt`. The local candidate contains its own runner and attribution. Direct Release links will serve this first-use path after publication and remote verification; they do not require changing the five domain pointers or adapting every historical index.

## Historical Download Entry Points

Protein's `latest.artifacts[].download_url` fields use graph-rag raw URLs:

```
https://raw.githubusercontent.com/hazelian0619/graph-rag/main/<domain>/release/external/<version>/tables/...
```

The other four current pointers use `current/path`, without the same `latest` shape. The existing downloader supports protein/latest. Legacy RNA and interaction releases in protian-entity remain explicit catalog entries; this version does not migrate their assets or claim a uniform five-domain download CLI.

## CI Coverage

The existing workflows have narrower coverage than the five-domain catalog:

- **data-qa**: runs the protein contract validator, with selected path triggers. It is not five-domain scientific validation.
- **release-metadata**: retains the protein/RNA/molecule/interaction matrix and release tests. The non-protein pointer mismatch remains a known historical gap; consistency's `--allow-issues` is report mode, not acceptance.

Starter acceptance separately requires its release tests, repeated build hashes, clean-directory execution and actual remote download/checksum/extraction. Passing these checks does not imply that the historical domain matrix is green.
