# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]
- Add a source-attributed protein/compound starter with offline examples, input/file hashes and a versioned package build.
- Tighten release download/consistency integrity checks and repair moved schema/protein-table paths.
- Document all five domains, sample/full-table distinctions and unverified historical delivery boundaries.
- Formalize dual-repo architecture: graph-rag as sole public entry point
- Extend CI to all four domains (data-qa + release-metadata)
- Promote release schema to shared `schemas/` directory
- Rewrite shared docs for four-domain standard

## [protein-v6] - 2025-10-26
- Primary protein entity table `protein_master_v6_clean.tsv` (19,135 rows; see the table header for its current fields)
- Added gene ID fields and AlphaFold v6 updates

## [rna-v1] - 2026-04-16
- First public RNA package: sequence, structure, annotation, evidence, reference tables

## [molecule-v1] - 2026-04-16
- First public molecule package: structure, identifiers, activity, physicochemical, pk_toxicology, semantic tables

## [interaction-l2-v1] - 2026-04-16
- First public cross-entity interaction package: PPI, PSI, RPI edges and evidence
- Chunked distribution for large integrated tables
