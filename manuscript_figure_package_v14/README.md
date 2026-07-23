# CoupledMD manuscript figure package (v14)

This directory separates plotting inputs, executable plotting scripts, rendered
figures, and the six submission-facing supplementary data files.

## Layout

- `source_data/`: authoritative plotting inputs retrieved from the manuscript
  server, plus the current S1 system inventory used by `figure_s1.py`.
- `scripts/`: one runnable entry point per manuscript figure. Shared plotting
  functions are kept in `_plot_v12_core.py` to prevent code duplication.
- `figures/`: generated PDF and 600-dpi PNG files.
- `supplementary_files/`: current Supplementary Data S1–S6 and their README,
  refreshed from the `main` branch of `HuangJianxiang-SJTU/coupledmd_database`.

## Generate figures

Run from this package directory:

```bash
python scripts/figure1.py
python scripts/figure2.py
python scripts/figure3.py
python scripts/figure4.py
python scripts/figure5.py
python scripts/figure_s1.py
python scripts/figure_s3.py
```

Each command writes a PDF/PNG pair to `figures/`.

## Provenance

- Supplementary files: GitHub `main` commit
  `559d1f2824dfe2effe417723ab4f5527066633bb`, retrieved 2026-07-23.
- Main/S3 plotting implementation: server `plot_v12_figures.py`, SHA-256
  `f25738fc19456bd5504e017e15fdd68a1359a56fae0cb576f93db23b78b52c19`.
- S1 implementation: server `supplementary_s1_plot_207_timestamped.py`,
  SHA-256
  `b9c8b53562224ae602e9149a1f6b368a862e8fd896b64c4fab8c7327e6e6854e`,
  adapted only for the package-relative input and output paths.

The two `.npy` files in `source_data/` are required by Figure 5 and accompany
the CSV plotting inputs from the same authoritative server directory.
