# RNA L1/L2 v2

This release publishes the RNA dataset artifacts that are intentionally excluded from git (`data/output/**`, `pipelines/**/reports/**`).

## Included
- 16 compressed RNA tables (`*.tsv.gz`)
- 15 non-sample validation reports (all PASS)
- key audit artifacts for tRNA v2 and ENST↔URS coverage
- `manifest.json` (row counts, checksums, commit)
- `SHA256SUMS.txt`

## Validation summary
- non-sample validations: **15/15 PASS**

## Notes
- Structural low-coverage items (ENST↔URS, Rfam full coverage, PDB sparsity) are data-source constraints.
- Active optimization line remains tRNA anticodon (v2 included in this release).
