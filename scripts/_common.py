#!/usr/bin/env python3
"""Shared constants, styling and helpers for the CoupledMD figure scripts."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgb  # noqa: F401  (re-exported for figure scripts)
from matplotlib.lines import Line2D  # noqa: F401
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: F401
import numpy as np  # noqa: F401
import pandas as pd

FAMILIES = ["Gi", "Gs", "Gq", "G12-13"]
FAMILY_LABELS = {
    "Gi": "Gi/o",
    "Gs": "Gs",
    "Gq": "Gq/11",
    "G12-13": "G12/13",
}
FAMILY_COLORS = {
    "Gi": "#2F8F6B",
    "Gs": "#2C6FB3",
    "Gq": "#C0741A",
    "G12-13": "#8A4AA0",
}
POCKET_COLORS = {
    "orthosteric": "#8172B2",
    "extracellular_vestibule": "#64B5CD",
    "tm_core_allosteric": "#2A9D8F",
    "intracellular_allosteric": "#DD8452",
}
INK = "#20242B"
MUTED = "#68707A"
GREY = "#B8B8B8"
PALE = "#E8EBEF"


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8.4,
        "axes.titlesize": 9.4,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.7,
        "xtick.labelsize": 7.7,
        "ytick.labelsize": 7.7,
        "legend.fontsize": 7.4,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def read(data_dir: Path, filename: str) -> pd.DataFrame:
    return pd.read_csv(data_dir / filename)


def panel_label(ax, value: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        value,
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=10.4,
        ha="left",
        va="bottom",
    )


def save(fig, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def family_sort(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["_family"] = pd.Categorical(
        result.g_protein_family, FAMILIES, ordered=True
    )
    return result.sort_values(["_family", "system_id"]).drop(columns="_family")




def cli(description: str) -> argparse.Namespace:
    """Standard three-directory command line for one-figure scripts."""
    repository_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=repository_root / "source_data",
    )
    parser.add_argument(
        "--assets-dir",
        type=Path,
        default=repository_root / "assets",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repository_root / "figures",
    )
    return parser.parse_args()
