# CoupledMD

## Scientific Data descriptor

**CoupledMD: molecular dynamics trajectories of active GPCR–G-protein complexes**

The descriptor release contains 207 systems and 621 computational repeats (three 500-ns runs per system).

- [Descriptor code, methods and reuse guidance](descriptor_revision/README.md)
- [Final main Figures 1–6 and Supplementary Figures S1–S4](descriptor_revision/figures/README.md)
- [Supplementary Data S1–S8](descriptor_revision/supplementary_data/)
- [Stationarity results, including the interface assessment](descriptor_revision/revision/equilibration_out/)
- [Verify the packaged data and checksums](descriptor_revision/verify_tables.py)
- Dataset: [all versions](https://doi.org/10.5281/zenodo.21395291), [version used for this revision](https://doi.org/10.5281/zenodo.23106855)

The repository-root URL is the entry point for the manuscript's code citation:
https://github.com/HuangJianxiang-SJTU/coupledmd_database/

## Separate GPCR-atlas package

The legacy atlas scripts and S01–S11 bundle below concern a separate analysis. They are not the descriptor's Figures 1–6 or Supplementary Data S1–S8.

### Atlas figure package

This repository contains the figure-generation scripts and source data for the
GPCR atlas manuscript (5 main-text figures, 4 supplementary figures) and the
machine-readable supplementary data bundle.

## Repository structure

```
scripts/          Figure-generation scripts (one per figure)
  figure_common.py    Shared colour palette, portal definitions, helpers
  figure1.py          Fig 1 — Measurement problem
  figure2.py          Fig 2 — Replica filtering
  figure3.py          Fig 3 — Recurrent core
  figure4.py          Fig 4 — Lateral-opening hierarchy
  figure5.py          Fig 5 — Cavity association & evidence collapse
  figure_s1.py        Fig S1 — Cohort design
  figure_s2.py        Fig S2 — Gateway variance components
  figure_s3.py        Fig S3 — Within-receptor contrasts
  figure_s4.py        Fig S4 — Interface descriptors

data/             Source data consumed by the figure scripts
  source_data/        Input CSV tables (gateway metrics, interface descriptors,
                      cohort composition, known-site concordance, etc.)
  downstream/         Replica-analysis outputs (pocket fingerprints, cluster
                      recurrences, stability, permutation tests, etc.)
  audits/             G-protein α5-hook sequence audit
  manifests/          Receptor partition tables
  gateway_rigor_summary.json
  cavity_validation_rigor_summary.json
  gateway_static_geometry_summary.json

figures/          Generated SVG figures (output after running scripts)

supplementary_data/  Machine-readable data files backing the SI tables
  S01–S10/           Supplementary Data S1–S10 delivery bundle
  S11_statistical_inventory/  Statistical-test inventory artifacts (SI Table S11)
```

## Dependencies

- Python ≥ 3.10
- matplotlib
- numpy
- pandas
- scipy

## Usage

```bash
cd scripts
python3 figure1.py     # → ../figures/Fig1.svg
python3 figure2.py     # → ../figures/Fig2.svg
# ... etc for each figure
```

All scripts save their output to the adjacent `figures/` directory.
Each SVG is accompanied by a `_values.json` file with panel metadata and
numerical sources.

## Figure inventory

| Figure | Description |
|--------|-------------|
| Fig 1  | Effective-sample funnel, variance partition, three-tier distances, portal reliability |
| Fig 2  | Pooled → replica-filtered atlas, split-half recovery, effect sharpening |
| Fig 3  | Stable-core recurrence, null-calibrated site concordance, conservation dissociation |
| Fig 4  | Seven-portal hierarchy, rank stability, prepared-geometry control, stationarity |
| Fig 5  | 486-cell screen, R6–Gq robustness, construct evidence |
| Fig S1 | Cohort composition, receptor-by-family contrast design |
| Fig S2 | Per-portal mixed-model variance components |
| Fig S3 | Gateway context differences, pocket context distances, signal/noise |
| Fig S4 | Interface-descriptor replica reliability, descriptor–outcome correlations |

## Zenodo

A complete archive of the paper inputs, pipeline outputs, and software is
available at the companion Zenodo record (reserved DOI).

## Citation and licences

Please use the version-specific CITATION.cff in the manuscript package (see
Zenodo record, reserved DOI) when citing CoupledMD. Data are licensed under
CC BY 4.0, and code is licensed under the MIT License.
