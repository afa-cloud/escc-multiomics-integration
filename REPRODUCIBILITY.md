# Reproducibility guide

## 1. Data location

The pipeline expects a reusable data warehouse at `RESEARCH_DATA_ROOT`. If the variable is not set, scripts use `~/ResearchDataHub`.

The warehouse is organized by source, stable accession and version. Raw source bytes are not included in this repository. Use `data/datasets.tsv` to obtain the exact accession, version, analysis role and limitations for each dataset.

## 2. Execution order

Scripts are numbered in dependency order:

1. `01`–`12`: inventory, source-data ingestion and TCGA multi-assay construction.
2. `13`–`19`: driver, continuous heterogeneity, cross-layer, metabolomics, microbiome and initial event–state analyses.
3. `20`–`23`: layered evidence integration, representation-overlap auditing and ECMS calibration.
4. `24`–`26`: main-figure generation and external continuous-state validation.
5. `27`–`29`: manuscript-table assembly, submission-readiness auditing and publication-facing supplement generation.

The scripts intentionally stop when required identities, versions, manifests or analysis-unit constraints are not satisfied. Do not bypass these checks by substituting unmatched samples or by treating overlapping cohorts as independent validation.

## 3. Software

The R packages detected in the released scripts and their versions in the reference execution environment are listed in `environment/r_packages.tsv`. External command-line dependencies are listed in `environment/system_tools.tsv`.

Python script `09_ingest_cao2020_proteomics_supplement.py` uses only Python standard-library modules.

## 4. Lightweight checks

These commands check syntax without downloading source data:

```bash
for f in scripts/*.R; do Rscript -e "parse(file='$f')" >/dev/null; done
for f in scripts/*.sh; do bash -n "$f"; done
python3 -m py_compile scripts/*.py
```

Processed-table integrity can be checked against `MANIFEST.tsv` with any SHA-256 implementation. Scientific regeneration requires the versioned source data and package/tool environment described above.

## 5. Interpretation limits

The repository preserves non-significant external results and cohort-dependence labels. Re-running or extending the workflow must not:

- label TCGA representation-overlap sensitivity results as external replication;
- reinterpret small calibration cohorts as causal mechanism validation;
- merge non-patient-matched metabolite, microbiome and host datasets into a causal chain;
- present the exploratory pCR association as a trained or validated clinical predictor.

