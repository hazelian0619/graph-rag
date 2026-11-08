# Biomedical starter: three offline examples

The starter is a filtered, single-source ChEMBL activity subset linked to
UniProtKB/Swiss-Prot proteins. It contains three TSV tables and a standalone
Python script. No checkout, repository imports, network, API key, database or
third-party Python dependency is needed to run the extracted package.

From the repository root, build and run with an installed Python 3:

```bash
python3 scripts/build_starter_dataset.py --repo-root . --out build/starter-v1
python3 examples/run_examples.py --data-dir build/starter-v1/data --out build/example-results
```

For a downloaded package, first verify the archive against the separately
obtained checksum, then extract and run it:

```bash
shasum -a 256 -c SHA256SUMS.txt
tar -xzf starter-v1.tar.gz
python3 starter-v1/examples/run_examples.py --data-dir starter-v1/data --out example-results
```

The runner also defaults to the `data/` directory beside its packaged `README.md`.
The output directory is relative to the current working directory. A custom
builder `--out build/my-starter` creates `my-starter/`, `my-starter.tar.gz`, and
a sibling `SHA256SUMS.txt` containing only that archive's digest. It replaces
that checksum file on each successful build; keep different packages in
different parent directories if their checksum files must coexist.

## Outputs and interpretation

- `01_protein.json`: one recorded protein's exact UniProt ID, name, taxon,
  provenance and all activity evidence linked to it **within this subset**.
- `02_molecule.json`: one recorded molecule's InChIKey, ChEMBL ID and provenance,
  plus its linked targets, assay type, measurement type, comparator, raw values,
  normalized values, source record IDs and DOI/PMID from the subset.
- `03_record_comparison.json`: a deterministic pair with the same molecule,
  protein, assay type, measurement type, comparator and unit but numerically
  different `normalized_nM` values. All record values are read from the TSV,
  unchanged. It explicitly warns **实验可比性未知**. Numeric differences do not
  establish scientific conflict or clinical effects; the same PMID is not
  independent replication. Different references alone do not establish it either.

Defaults select a deterministic comparison pair and use its endpoints for the
first two examples. `--protein-id P56373` and `--molecule-id FIAFUQMPZJWCLV-UHFFFAOYSA-N`
override the first and second queries independently; the comparison remains the
deterministically selected package example. Unknown IDs produce explicit
`status: "unrecorded"`, a null entity and empty evidence, not a claim of no effect.

Before writing results, the runner checks all payload SHA256 values against
`manifest.json`, exact table schemas, row counts, unique IDs, source restrictions,
matching molecule ChEMBL IDs and endpoint closure. Missing files, malformed
tables, mismatched hashes or absence of a usable comparison pair exit nonzero.
The manifest is an integrity inventory, not a signature: archive verification
against a trusted external checksum is still required for authenticity.

## Data ownership, fields and modifications

`manifest.json` records exact exported fields, per-table rows and SHA256,
input repository-relative paths and actual hashes, sources, versions, license
links and filtering counts. Nothing infers an upstream release from a fetch date.
Only referenced endpoints are retained; source values and activity values are
not rewritten. No structures, confidence scores or aggregate fields are exported.

| Table | Origin and version | License and modifications |
| --- | --- | --- |
| `data/proteins.tsv` | UniProt Consortium, UniProtKB/Swiss-Prot; snapshot fetch dates retained, upstream release not recorded | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), [official source](https://rest.uniprot.org/help/license). Only approved ID/name/taxon/source/date fields and referenced proteins are selected. |
| `data/molecules.tsv` | ChEMBL team, EMBL-EBI; compound_structures, `ChEMBL:chembl_36.db` | [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/), [official source](https://www.ebi.ac.uk/chembldb/). Only approved identifiers/provenance and referenced molecules are selected; no structure text is exported. |
| `data/activity_evidence.tsv` | ChEMBL team, EMBL-EBI; activity sample, `chembl_36` | [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/), same official source. Existing project-generated evidence IDs and existing `normalized_nM` derivation are retained without recomputation; approved fields, exact endpoint links and DOI/PMID are selected. |

Keep attribution, source links, license links and modification notices when
redistributing; ChEMBL-derived adaptations remain under CC BY-SA 3.0. The source
organizations do not endorse this package. Code is separately MIT licensed;
the full MIT text and per-table attribution are in `NOTICE.md`. MIT does not
replace the data licenses. This package offers traceable records, not new
experimental results, complete biological coverage or validated clinical advice.
