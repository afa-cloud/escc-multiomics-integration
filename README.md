# ESCC layered public multi-omics integration

This repository contains the analysis code and processed, publication-facing result tables for a layered public multi-omics study of esophageal squamous cell carcinoma (ESCC).

The workflow integrates public genomic, transcriptomic, DNA-methylation, protein, metabolite and tissue-microbiome resources while preserving cohort identity, analysis-unit boundaries and negative validation results. The repository is a reproducibility archive, not a clinical software product.

## Repository contents

- `scripts/`: 29 ordered R, Python and shell scripts covering data inventory, ingestion, TCGA multi-assay construction, driver and heterogeneity analyses, cross-layer calibration, metabolomics, tissue 16S analysis, representation-overlap auditing, ECMS projection, external continuous-state calibration, visualization and manuscript-table assembly.
- `data/datasets.tsv`: versioned public-data inventory using `${RESEARCH_DATA_ROOT}` rather than machine-specific paths.
- `results/processed_tables/`: 71 cleaned, publication-facing TSV tables.
- `results/result_index.tsv`: machine-readable index linking processed tables to their scientific roles and interpretation boundaries.
- `results/escc_manuscript_table*.tsv`: machine-readable versions of the three main manuscript tables.
- `environment/`: software-package and external-tool inventory.
- `MANIFEST.tsv`: release file sizes and SHA-256 checksums.

## Scientific scope and boundaries

- Public datasets are used for hypothesis generation and calibration; associations are not causal proof.
- The nine retained genomic event–continuous state relationships remain within the TCGA source cohort after representation-overlap sensitivity analysis and are not independent replications.
- Factor1 and Factor3 RNA proxies showed strong held-out internal fidelity, but external survival and ECMS incremental associations were negative; the 28-patient pCR analysis is exploratory.
- Non-matched metabolomics, microbiome and host cohorts are retained as separate orthogonal modules and are not joined into a patient-level causal axis.
- Source data are not redistributed. Obtain them from the accessions and repositories listed in `data/datasets.tsv`.

## Reproduction

Set the shared data root before running the ordered scripts:

```bash
export RESEARCH_DATA_ROOT=/path/to/ResearchDataHub
```

Additional optional tool locations can be supplied through `AMP_ENV_PREFIX`, `ESCC_AMPLICON_ENV`, `BEDTOOLS` and `ARIA2C`. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for execution order, inputs and validation guidance.

## Citation

Citation metadata are provided in `CITATION.cff`. A version-specific Zenodo DOI is associated with the `v1.0.1` archival release.

## Licenses

- Analysis code: MIT License.
- Author-generated processed tables: CC BY 4.0 to the extent that rights are held by the authors.
- Underlying public datasets, reference resources and third-party software remain subject to their original repository terms and licenses.
