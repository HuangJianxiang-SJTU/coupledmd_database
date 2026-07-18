# CoupledMD database

Publication figures, figure source data and Supplementary Data for **CoupledMD**, a curated cohort of **207 active-state GPCR–G-protein complexes** (95 Gi/o, 65 Gs, 41 Gq/11, 6 G12/13; 174 receptor names, 173 mapped UniProt accessions). The source simulation campaign comprises three 500-ns production replicas per system (621 replicas; 310.5 µs total). The distributed molecular dataset is a reduced derivative: one replica-1-derived retained-component PDB/XTC pair per system (414 files; 2,500 frames at approximately 200-ps intervals, 499.8 ns per trajectory).

This repository is the data/figure companion to the CoupledMD code repository and web portal:

- Portal and REST API: <https://www.coupledmd.cn> (`/api/docs`)
- Code repository (portal, API, analysis): <https://github.com/HuangJianxiang-SJTU/coupledmd>
- Molecular dataset: Zenodo DOI `10.5281/zenodo.21395292` (private draft; public upon publication)

## Layout

| Path | Contents |
|---|---|
| `figures/` | Publication figures 1-6 and Supplementary Figures S1-S2, PNG (600 dpi) and PDF |
| `source_data/` | Path-neutral, panel-level source tables (16 CSV) with a row/hash manifest and independent audit |
| `scripts/` | One plotting script per figure, shared helpers (`_common.py`), and the independent source-data validator |
| `assets/` | Representative structural rendering used by Figure 6 panel A |
| `supplementary_data/` | Supplementary Data S1-S10 (207-system cohort), package README and package audit |

## Reproducing the figures

Requirements: Python 3.10+ with `matplotlib`, `pandas`, `numpy`; Poppler's `pdftoppm` only for the Supplementary Figure S2 raster export.

Each main figure has its own script. From the repository root:

```bash
python scripts/figure1_composition.py
python scripts/figure2_organization.py
python scripts/figure3_pocket_validation.py
python scripts/figure4_gateways.py
python scripts/figure5_technical_validation.py
python scripts/figure6_reuse_atlas.py
python scripts/supplementary_figure_s1_receptor_profiles.py
python scripts/supplementary_figure_s2.py
```

Every script accepts `--data-dir` (default `source_data/`), `--assets-dir` (default `assets/`) and `--output-dir` (default `figures/`). All plots use the portal palette (Gi/o `#2F8F6B`, Gs `#2C6FB3`, Gq/11 `#C0741A`, G12/13 `#8A4AA0`) and prefer Helvetica when installed.

Re-run the independent source-data audit with:

```bash
python scripts/validate_figure_source_data.py
```

It recomputes row counts, key sets and checksums against `source_data/figure_source_data_manifest.csv` and fails on any internal-path or legacy-version token.

## Notes

- **Figure 6 panel A** was further post-processed by the authors for the manuscript (layout and annotation of the structural rendering), so the published panel may differ slightly from the programmatic rendering reproduced here.
- **Supplementary Data S10** inventories the code and build environment and is refreshed when the code repository release tag is frozen.
- Supplementary Figures S1 and S2 are unaffected by the Gq_7DWC exclusion: the receptor-dendrogram profiles (174 receptors) and the gateway method dictionary are cohort-count independent.
- Gq_7DWC is excluded from the release cohort as a duplicate source identity of Gq_8DWC; Gq_7E9W is excluded as a non-GPCR duplicate/mislabel. Supplementary Data S2 records this boundary.

## Licenses

- Code (everything under `scripts/`): MIT, see `LICENSE`.
- Data and figures (`figures/`, `source_data/`, `supplementary_data/`, `assets/`): CC BY 4.0, see `LICENSE-DATA`.

## Citation

Please cite the CoupledMD Scientific Data article (in preparation) and, once the dataset record is published, the exact dataset DOI and version. Until then, identify the portal access date when referencing derived records.
