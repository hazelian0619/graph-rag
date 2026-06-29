# Architecture

## Repository Role

**graph-rag is the sole public-facing dataset repository.** External users download all datasets from this repository. It does not serve as a development workspace or store raw pipeline inputs.

Sister repository [protian-entity](https://github.com/hazelian0619/protian-entity) is the development engine where pipelines are built and data is produced. Only curated, QA-passed results are promoted to graph-rag for public delivery.

## Four Domains

| Domain | Scope | Latest version |
|---|---|---|
| `protein/` | Human protein L1 entity tables | protein-v6 |
| `rna/` | Human RNA L1/L2 entity tables | rna-v1 |
| `molecule/` | Small molecule L1/L2 entity tables | molecule-v1 |
| `interaction/` | Cross-entity L2 edges (PPI / PSI / RPI) | interaction-l2-v1 |

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
.github/workflows/ # CI: data-qa + release-metadata (covers all four domains)
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

## Download Entry Point

All `download_url` fields in `products/*/current.json` point to graph-rag raw URLs:

```
https://raw.githubusercontent.com/hazelian0619/graph-rag/main/<domain>/release/external/<version>/tables/...
```

No download_url points to protian-entity or any other repository.

## CI Coverage

Two GitHub Actions workflows run on all four domains:

- **data-qa**: validates data contracts when any domain's `data/` or `pipelines/` changes
- **release-metadata**: builds release index, validates schema, checks consistency when any domain's `products/` or `release/` changes
