# CoupledMD GPCR–G-protein Database — Paper 1 Figure Package

This repository contains the figure-generation scripts and source data for the
GPCR atlas manuscript (5 main-text figures, 4 supplementary figures).

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
  figure_s2.py        Fig S2 — Within-receptor contrasts
  figure_s3.py        Fig S3 — Gateway variance components
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
| Fig S2 | Gateway context differences, pocket context distances, signal/noise |
| Fig S3 | Per-portal mixed-model variance components |
| Fig S4 | Interface-descriptor replica reliability, descriptor–outcome correlations |

## Zenodo

A complete archive of the paper inputs, pipeline outputs, and software is
available at the companion Zenodo record (reserved DOI).

## Repository layout

- [`manuscript_figure_package_v14/`](manuscript_figure_package_v14/) — the
  complete manuscript release package:
  - [`figures/`](manuscript_figure_package_v14/figures/) contains the rendered
    manuscript figures.
  - [`source_data/`](manuscript_figure_package_v14/source_data/) contains the
    plotting inputs.
  - [`scripts/`](manuscript_figure_package_v14/scripts/) contains the figure
    scripts and analysis workflows.
  - [`supplementary_files/`](manuscript_figure_package_v14/supplementary_files/)
    contains Supplementary Data S1–S6.
- [`licences_and_citation/`](licences_and_citation/) — citation metadata and
  separate licences for data and code.

Start with the
[`v14 package README`](manuscript_figure_package_v14/README.md) for figure
generation instructions, or the
[`workflow README`](manuscript_figure_package_v14/scripts/workflows/README.md)
for reduction, quality-control, structural-validation, persistent-pocket, and
worked-example analyses.

The full molecular-dynamics trajectory archive is distributed separately as
described in the manuscript's Data Availability statement.

## Citation and licences

Please use [`CITATION.cff`](licences_and_citation/CITATION.cff) when citing
CoupledMD. Data are licensed under
[`CC BY 4.0`](licences_and_citation/LICENSE), and code is licensed under the
[`MIT License`](licences_and_citation/LICENSE-CODE).
