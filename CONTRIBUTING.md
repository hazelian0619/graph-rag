# Contributing

Thank you for contributing to Graph-RAG. This document explains development, validation, and release practices used by the project.

## Development workflow

- Branching: create a branch using `chore/` or `feat/` prefixes (e.g., `chore/standardize-docs-ci`, `feat/rna-v2`).
- Tests & validation: run validation before opening a PR:

```bash
# Protein example
python3 tools/kg_validate_table.py \
  --contract protein/pipelines/protein/contracts/protein_master_v6.json \
  --table protein/data/processed/protein_master_v6_clean.tsv \
  --out build/validate/protein_master_v6_report.json
```

- Commit messages: use clear, imperative messages. Follow conventional commits if possible.

## Pull requests

- Open PRs against `main`. Describe the change, test steps, and link to any data releases.
- Include validation reports for any changes to entity tables.
- CI must pass (data-qa + release-metadata workflows cover all four domains).

## Releases

- Curated result tables live in `<domain>/release/external/<version>/tables/`.
- Release index must validate after any product change.
- Large assets (L1 tables) use GitHub Releases; the repo holds metadata only.
- Attach `manifest.json` with checksums, row counts, git commit SHA, and QA reports.

## Contacts

- For maintenance and code ownership see `CODEOWNERS` or raise an issue.
- Development happens in [protian-entity](https://github.com/hazelian0619/protian-entity); this repo is the public deliverable.
