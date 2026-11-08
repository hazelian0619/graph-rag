# Graph-RAG · Traceable bioactivity data

Explore protein–compound measurements with source identifiers, publication references, and runnable offline examples.

从一个小包开始：查蛋白、查小分子相关实验记录，再并排查看同一实体对的两条测量记录。每条记录保留来源和文献标识，方便继续核查。

## Start with the small package

The `starter-v1` subset connects **87 proteins**, **102 compounds**, and **359 activity records** selected from the existing 5,000-record ChEMBL 36 sample. It is an introduction to the data, not the size of the full project or a representative benchmark.

Download [dataset-starter-v1.0.0](https://github.com/hazelian0619/graph-rag/releases/tag/dataset-starter-v1.0.0). Anonymous download, checksum verification, extraction and all three examples passed on 2026-11-07.

```bash
curl --fail --location --output starter-v1.tar.gz https://github.com/hazelian0619/graph-rag/releases/download/dataset-starter-v1.0.0/starter-v1.tar.gz
curl --fail --location --output SHA256SUMS.txt https://github.com/hazelian0619/graph-rag/releases/download/dataset-starter-v1.0.0/SHA256SUMS.txt
shasum -a 256 -c SHA256SUMS.txt
tar -xzf starter-v1.tar.gz
python3 starter-v1/examples/run_examples.py --data-dir starter-v1/data --out example-results
```

The package runs with Python's standard library; no repository clone, API key, LLM, or database service is needed. See [example instructions and field provenance](examples/README.md). The archive is **32,619 bytes (about 32 KiB)**. Acceptance covers macOS with Python **3.9.21** and **3.12.6**, including isolated execution of the remotely downloaded package without installed site packages; Linux/Windows have not been independently accepted here.

| Example | What you get |
|---|---|
| Protein lookup | Protein identity and the activity records connected to it in this subset |
| Compound lookup | Targets, assay type, measurement, relation, unit, source record and DOI/PMID |
| Record comparison | Two measurements for the same entity pair, metric, relation and unit, retaining their original values and references |

The default comparison reads ChEMBL activity records **51795 / 51794** for target **P37231**, with EC50 values **8.6 / 54 nM**, both referring to PMID **14592510**. An excerpt from the actual output is:

```json
{
  "warning": "实验可比性未知",
  "match_fields": {
    "target_uniprot_accession": "P37231",
    "standard_type": "EC50",
    "standard_relation": "=",
    "standard_unit": "nM"
  },
  "same_pubmed_id": true
}
```

This is a record difference: experimental comparability is unknown, and a shared publication is not independent replication.

To rebuild from an existing checkout:

```bash
python3 scripts/build_starter_dataset.py --repo-root . --out build/starter-v1
python3 examples/run_examples.py --data-dir build/starter-v1/data --out build/example-results
```

The manifest records input hashes, exported fields, selection rules, table counts, file hashes and licenses. The builder retains both endpoints and a DOI or PMID for every exported activity record. It does not regenerate upstream measurements.

## Existing data beyond the starter

The repository also contains earlier protein, RNA, compound, interaction and disease results. Their maturity and distribution paths differ. [The dataset catalog](docs/DATASET_CATALOG.md) distinguishes checked files, manifest-reported counts, samples and legacy releases.

| Domain | Version pointer | Scope |
|---|---|---|
| Protein | [protein-v6](protein/products/protein/current.json) | Protein identity, sequence, structure and annotation tables |
| RNA | [rna-v1](rna/products/rna/current.json) | RNA tables and samples; older large assets have a separate release |
| Molecule | [molecule-v1](molecule/products/molecule/current.json) | Identifiers, structure and activity results; activity tables include samples |
| Interaction | [interaction-l2-v1](interaction/products/interaction/current.json) | PPI, protein–small molecule and RNA–protein results, including chunked legacy assets |
| Disease | [disease-v1](disease/products/disease/current.json) | Disease and phenotype relations; currently in development |

The existing downloader supports the protein/latest interface. The starter uses the fixed direct links above; a uniform five-domain CLI and full historical release verification are not delivered by this version.

## Sources, reuse and limits

- Starter protein fields originate from UniProtKB/Swiss-Prot under [CC BY 4.0](https://rest.uniprot.org/help/license). Starter compound and activity tables derive from ChEMBL 36 under [CC BY-SA 3.0](https://www.ebi.ac.uk/chembldb/). Attribution, modifications and table-specific terms travel with the package in `NOTICE.md` and `manifest.json`.
- Code uses [MIT](LICENSE). This does not relicense upstream data or establish a blanket license for every historical table. See the catalog before reusing other domains.
- This is a filtered, single-source bioactivity subset, joined to UniProt identity fields. It does not establish clinical efficacy, independent cross-source validation, calibrated confidence, or model performance.

When reusing the starter, cite Graph-RAG's package version and SHA256, retain source record IDs and references, and acknowledge UniProt and ChEMBL under their respective terms. Upstream citation details are included in the package instructions.

See [architecture](ARCHITECTURE.md), [repository standards](docs/REPO_STANDARDS.md), and [changes](CHANGELOG.md) for maintenance and release boundaries.
