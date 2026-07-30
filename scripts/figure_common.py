#!/usr/bin/env python3
"""Shared constants, paths, and helpers for the GPCR atlas figure scripts."""

from __future__ import annotations

import json
import math
import os
from decimal import Decimal
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/coupledmd_figure_mpl")

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# ---- Paths ----
HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
SRC = DATA / "source_data"
DOWN = DATA / "downstream"
CONSERVATION = DOWN / "classA_conservation"
OUT = HERE.parent / "figures"

# ---- Constants ----
MM_PER_IN = 25.4
FIG_WIDTH_IN = 180.0 / MM_PER_IN

COLORS = {
    "Gi": "#4477AA",
    "Gq": "#CC6677",
    "Gs": "#DDCC77",
    "G12-13": "#AAAAAA",
    "Class A": "#332288",
    "Class B": "#88CCEE",
    "Receptor identity": "#117733",
    "Residual replica": "#882255",
    "Administrative family": "#777777",
    "System / context": "#999933",
    "GPCR class": "#DDDDDD",
    "Significant": "#CC6677",
    "Non-significant": "#BBBBBB",
}

PORTALS = [
    "TM1-TM2", "TM2-TM3", "TM3-TM4", "TM4-TM5",
    "TM5-TM6", "TM6-TM7", "TM7-TM1",
]
PORTAL_LABELS = {portal: portal.replace("-", "–") for portal in PORTALS}
PORTAL_COLORS = {
    portal: mpl.colors.to_hex(plt.get_cmap("viridis")(x))
    for portal, x in zip(PORTALS, np.linspace(0.08, 0.92, 7))
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "lines.linewidth": 0.75,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "axes.unicode_minus": False,
})


# ---- Helpers ----
def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, **kwargs)


def read_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, pd.Series):
        return jsonable(value.tolist())
    if np.isscalar(value) and pd.isna(value):
        return None
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(jsonable(value), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def panel_label(ax, letter: str) -> None:
    ax.text(
        -0.13, 1.05, letter,
        transform=ax.transAxes,
        fontsize=8, fontweight="bold",
        va="bottom", ha="left", clip_on=False,
    )


def panel_title(ax, title: str) -> None:
    ax.set_title(title, loc="left", fontweight="bold", pad=3)


def placeholder(ax, note: str, aspect_width_over_height: float) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_box_aspect(1 / aspect_width_over_height)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linestyle((0, (4, 4)))
        spine.set_linewidth(0.5)
        spine.set_color("#999999")
    ax.text(0.5, 0.55, "AUTHOR-PREPARED PANEL",
            ha="center", va="center", fontsize=6, color="#999999",
            transform=ax.transAxes)
    ax.text(0.5, 0.42, note,
            ha="center", va="center", fontsize=5, color="#999999",
            style="italic", transform=ax.transAxes, wrap=True)


def report_placeholder(fig, ax, name: str) -> dict:
    fig.canvas.draw()
    bbox = ax.get_position()
    width_mm = bbox.width * fig.get_figwidth() * MM_PER_IN
    height_mm = bbox.height * fig.get_figheight() * MM_PER_IN
    info = {
        "panel": name,
        "width_mm": width_mm,
        "height_mm": height_mm,
        "aspect_width_over_height": width_mm / height_mm,
    }
    print(f"{name} placeholder: {width_mm:.2f} × {height_mm:.2f} mm; "
          f"aspect {width_mm / height_mm:.4f}:1")
    return info


def save_figure(fig, number: int, values: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    width_mm = fig.get_figwidth() * MM_PER_IN
    height_mm = fig.get_figheight() * MM_PER_IN
    if width_mm > 180.0001 or height_mm > 240.0001:
        raise AssertionError(
            f"Fig{number}: canvas {width_mm:.2f} × {height_mm:.2f} mm exceeds limit"
        )
    values["canvas"] = {
        "width_mm": width_mm,
        "height_mm": height_mm,
        "minimum_text_pt": 5,
    }
    fig.savefig(OUT / f"Fig{number}.svg", dpi=600)
    write_json(OUT / f"Fig{number}_values.json", values)
    plt.close(fig)


def save_supplementary(fig, stem: str, values: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.svg", dpi=600)
    values["canvas"] = {
        "width_mm": fig.get_figwidth() * MM_PER_IN,
        "height_mm": fig.get_figheight() * MM_PER_IN,
    }
    write_json(OUT / f"{stem}_values.json", values)
    plt.close(fig)
