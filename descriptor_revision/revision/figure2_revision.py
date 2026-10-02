#!/usr/bin/env python3
"""Shared plotting implementation for the organized CoupledMD figure scripts.

Publication-quality redesign: restrained Nature Communications / Scientific Data
palette, subtle dashed grids behind data, consistent bold panel labels, proper
typography (×, Å, µ, en-dash), and per-panel layout polish described in the
design brief.
"""
from __future__ import annotations

import re
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, Rectangle
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parent
SOURCE = PACKAGE_ROOT / "supplementary_data"
OUT = PACKAGE_ROOT / "figures"

FAMILIES = ["Gi", "Gs", "Gq", "G12-13"]
FAMILY_LABELS = {
    "Gi": "Gi/o",
    "Gs": "Gs",
    "Gq": "Gq/11",
    "G12-13": "G12/13",
}
COLORS = {
    "Gi": "#2F8F6B",
    "Gs": "#2C6FB3",
    "Gq": "#C0741A",
    "G12-13": "#8A4AA0",
}
INK = "#1E252B"
MUTED = "#626B73"
PALE = "#E8ECEF"
TEAL = "#287C7E"
BLUE = "#3F6FA0"
GREEN = "#3D8466"
AMBER = "#C4862C"
RED = "#B4514A"

# Sequential contact colormap: white -> saturated teal (for Fig 5A)
CONTACT_CMAP = LinearSegmentedColormap.from_list(
    "contact", ["#FFFFFF", "#D7ECEC", "#7FBDBE", "#2F8F8F", "#1E5E5E"]
)


def style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial",
                "Liberation Sans",
                "Helvetica",
                "DejaVu Sans",
            ],
            "font.size": 7.8,
            "axes.titlesize": 9.0,
            "axes.titleweight": "bold",
            "axes.labelsize": 8.2,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 6.8,
            "axes.linewidth": 0.65,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
        }
    )


def panel(ax: plt.Axes, letter: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=INK,
    )

def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ["pdf", "png"]:
        fig.savefig(
            OUT / f"{stem}.{extension}",
            dpi=600 if extension == "png" else None,
            bbox_inches="tight",
            pad_inches=0.03,
        )
    plt.close(fig)

def box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    subtitle: str = "",
    face: str = "#F7F8F9",
    edge: str = "#69737B",
    title_size: float = 7.6,
) -> None:
    x, y = xy
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.012,rounding_size=0.016",
            facecolor=face,
            edgecolor=edge,
            linewidth=0.8,
        )
    )
    ax.text(
        x + width / 2,
        y + height * 0.64,
        title,
        ha="center",
        va="center",
        fontsize=title_size,
        fontweight="bold",
        color=INK,
    )
    if subtitle:
        ax.text(
            x + width / 2,
            y + height * 0.29,
            subtitle,
            ha="center",
            va="center",
            fontsize=6.3,
            color=MUTED,
            linespacing=1.15,
        )

def arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = MUTED,
    dashed: bool = False,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=0.9,
            linestyle="--" if dashed else "-",
            color=color,
            shrinkA=1,
            shrinkB=1,
        )
    )

def figure2() -> None:
    ligand = pd.read_csv(SOURCE / "Supplementary_Data_S1_system_inventory.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.6), constrained_layout=True)
    ax = axes.ravel()
    for item in ax:
        item.axis("off")
        item.set_xlim(0, 1)
        item.set_ylim(0, 1)

    # ── Panel A: molecular-record boundary ──
    box(
        ax[0],
        (0.08, 0.73),
        0.84,
        0.16,
        "FULL-SYSTEM SOURCE SIMULATIONS",
        "protein complex + bound ligands + membrane + ions\n(solvent coordinates not archived)",
        face="#F1F2F3",
        edge=MUTED,
    )
    box(
        ax[0],
        (0.18, 0.43),
        0.64,
        0.15,
        "REDUCTION AND ALIGNMENT",
        "retain selected components · 200-ps frame spacing",
        face="#EDF4F4",
        edge=TEAL,
    )
    box(
        ax[0],
        (0.08, 0.10),
        0.84,
        0.17,
        "MATCHED REDUCED PDB/XTC RECORDS",
        "protein complex + bound ligands retained;\nnot a full-system archive",
        face="#EAF0F7",
        edge=BLUE,
    )
    # arrows between stages; callouts sit in the gaps so they label the
    # transition rather than overlapping the box content
    arrow(ax[0], (0.5, 0.72), (0.5, 0.59))
    arrow(ax[0], (0.5, 0.42), (0.5, 0.28))
    ax[0].text(
        0.62,
        0.655,
        "drop membrane\nand ions",
        ha="left",
        va="center",
        fontsize=6.4,
        color=RED,
    )
    ax[0].text(
        0.62,
        0.355,
        "align · 200-ps\nframe spacing",
        ha="left",
        va="center",
        fontsize=6.4,
        color=TEAL,
    )
    # ligand retention summary bar
    n_small = int(ligand.retained_small_molecule.sum())
    n_none = int((~ligand.retained_small_molecule).sum())
    total_width = 0.45
    bar_x = 0.08
    bar_y = 0.045
    bar_h = 0.035
    sm_frac = n_small / 207
    ax[0].add_patch(Rectangle(
        (bar_x, bar_y), total_width * sm_frac, bar_h,
        facecolor=AMBER, edgecolor="white", linewidth=0.5
    ))
    ax[0].add_patch(Rectangle(
        (bar_x + total_width * sm_frac, bar_y),
        total_width * (1 - sm_frac), bar_h,
        facecolor="#D8DDE1", edgecolor="white", linewidth=0.5
    ))
    ax[0].text(
        bar_x + total_width / 2, bar_y - 0.007,
        f"Retained small molecule: {n_small} yes / {n_none} no",
        ha="center", va="top", fontsize=5.8, color=MUTED,
    )
    ax[0].set_title("Molecular-record boundary")
    panel(ax[0], "A", x=-0.04)

    # ── Panel B: support matrix ──
    rows = [
        ("Reduced trajectory", ["yes", "yes", "no", "no", "no"]),
        ("Gateway summary", ["no", "no", "summary", "no", "no"]),
        ("Full-system source", ["yes", "yes", "yes", "yes", "no"]),
    ]
    columns = [
        "protein\ncomplex",
        "ligand\ncontact",
        "lipid\ngateways",
        "ion\ncoords",
        "solvent\ncoords",
    ]
    x_positions = [0.365, 0.4875, 0.61, 0.7325, 0.855]
    for x, label in zip(x_positions, columns):
        ax[1].text(x, 0.86, label, ha="center", va="bottom", fontsize=5.8,
                   fontweight="bold")
    for row_index, (label, values) in enumerate(rows):
        y = 0.66 - row_index * 0.22
        ax[1].text(0.02, y, label, ha="left", va="center", fontsize=7)
        for x, value in zip(x_positions, values):
            if value in {"yes", "summary"}:
                face = GREEN if value == "yes" else BLUE
                txt = "Y" if value == "yes" else "S"
            else:
                face = "#D8DDE1"
                txt = "—"
            ax[1].add_patch(
                FancyBboxPatch(
                    (x - 0.055, y - 0.052),
                    0.11,
                    0.104,
                    boxstyle="round,pad=.004,rounding_size=.010",
                    facecolor=face,
                    edgecolor="white",
                    linewidth=0.8,
                )
            )
            ax[1].text(
                x,
                y,
                txt,
                ha="center",
                va="center",
                color="white" if value in {"yes", "summary"} else MUTED,
                fontsize=7.4,
                fontweight="bold",
            )
        # row separator
        ax[1].plot([0.02, 0.92], [y - 0.072, y - 0.072], color=PALE,
                   linewidth=0.5)
    legend_items = [
        (0.07, 0.13, "Y", "supported", GREEN),
        (0.39, 0.45, "S", "summary only", BLUE),
        (0.72, 0.78, "—", "not supported", "#D8DDE1"),
    ]
    for symbol_x, label_x, symbol, label, color in legend_items:
        ax[1].add_patch(FancyBboxPatch(
            (symbol_x - 0.038, 0.035), 0.076, 0.082,
            boxstyle="round,pad=.003,rounding_size=.009",
            facecolor=color, edgecolor="white", linewidth=0.7,
        ))
        ax[1].text(symbol_x, 0.076, symbol, ha="center", va="center",
                   color="white" if symbol != "—" else MUTED,
                   fontsize=8.2, fontweight="bold")
        ax[1].text(label_x, 0.076, label, ha="left", va="center",
                   fontsize=6.6, color=INK, fontweight="bold")
    ax[1].set_title("What each record can support")
    panel(ax[1], "B", x=-0.04)

    # ── Panel C: system→replica→file joins ──
    box(
        ax[2],
        (0.18, 0.72),
        0.64,
        0.17,
        "SYSTEM",
        "207 system IDs\nmetadata join key",
        face="#F5F6F7",
    )
    box(
        ax[2],
        (0.18, 0.43),
        0.64,
        0.17,
        "REPLICA",
        "three computational\nrepeats per system",
        face="#EDF4F4",
        edge=TEAL,
    )
    box(
        ax[2],
        (0.18, 0.14),
        0.64,
        0.17,
        "FILES",
        "PDB + XTC\nmanifest record",
        face="#EAF0F7",
        edge=BLUE,
    )
    arrow(ax[2], (0.50, 0.71), (0.50, 0.61))
    arrow(ax[2], (0.50, 0.42), (0.50, 0.32))
    # circled ratio annotations
    for y, ratio in [(0.655, "1:3"), (0.365, "1:2")]:
        ax[2].add_patch(
            FancyBboxPatch(
                (0.55, y - 0.023),
                0.10,
                0.045,
                boxstyle="round,pad=.004,rounding_size=.020",
                facecolor="white",
                edgecolor=TEAL,
                linewidth=0.8,
            )
        )
        ax[2].text(0.60, y, ratio, ha="center", va="center", fontsize=6.4,
                   color=TEAL, fontweight="bold")
    ax[2].text(
        0.5,
        0.04,
        "Manifest joins sizes, atom counts, frame counts\nand component scope.",
        ha="center",
        va="center",
        fontsize=6.6,
        color=MUTED,
    )
    ax[2].set_title("Replica-aware organization")
    panel(ax[2], "C", x=-0.04)

    # ── Panel D: access and reuse layers ──
    access = [
        (
            "Supplementary Data",
            "inventories · annotations\nQC and source tables",
            GREEN,
        ),
        (
            "REST API",
            "coupledmd.cn · /api/v1\nmetadata · pockets · gateway records",
            TEAL,
        ),
        (
            "GitHub · CoupledMD",
            "descriptor workflows + portal code\nsee Code Availability",
            BLUE,
        ),
        (
            "Zenodo · replicas 1–3",
            "10.5281/zenodo.21395291\n621 reduced PDB/XTC record pairs",
            AMBER,
        ),
    ]
    for index, (title, subtitle, color) in enumerate(access):
        y = 0.72 - index * 0.22
        box(ax[3], (0.10, y), 0.80, 0.17, title, subtitle,
            face="#F7F8F9", edge=color, title_size=7.1)
    ax[3].set_title("Access and reuse layers")
    panel(ax[3], "D", x=-0.04)
    save(fig, "Figure_2")

if __name__ == "__main__":
    style(); figure2()
