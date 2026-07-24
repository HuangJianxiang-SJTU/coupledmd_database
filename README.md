# CoupledMD database

This repository contains the data, figures, and reproducibility code accompanying
the CoupledMD manuscript and its frozen 207-system GPCR–G-protein complex cohort.

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
