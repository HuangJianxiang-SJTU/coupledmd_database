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
SOURCE = PACKAGE_ROOT / "source_data"
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


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(SOURCE / name)


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


def soft_grid(ax: plt.Axes, axis: str = "y") -> None:
    ax.grid(
        axis=axis,
        color=PALE,
        linestyle="--",
        linewidth=0.5,
        alpha=0.85,
        zorder=0,
    )
    ax.set_axisbelow(True)


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


# ---------------------------------------------------------------------------
# Schematic helpers (Figure 2, S1)
# ---------------------------------------------------------------------------

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


def card(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    number: str,
    label: str,
    accent: str = TEAL,
    number_size: float = 16,
) -> None:
    """Publication data card: subtle fill, thin coloured top border, bold number."""
    x, y = xy
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.0,rounding_size=0.012",
            facecolor="#F2F5F6",
            edgecolor="#D3D9DC",
            linewidth=0.6,
        )
    )
    # accent top border
    ax.add_patch(
        Rectangle(
            (x, y + height - 0.018),
            width,
            0.018,
            facecolor=accent,
            edgecolor="none",
        )
    )
    ax.text(
        x + width / 2,
        y + height * 0.54,
        number,
        ha="center",
        va="center",
        fontsize=number_size,
        fontweight="bold",
        color=INK,
    )
    ax.text(
        x + width / 2,
        y + height * 0.20,
        label,
        ha="center",
        va="center",
        fontsize=6.6,
        color=MUTED,
        linespacing=1.1,
    )


# ---------------------------------------------------------------------------
# Syntax-highlighted code rendering (Figure 4F, 5B)
# ---------------------------------------------------------------------------

_CODE_KEYWORDS = {
    "import", "from", "as", "assert", "def", "return", "for", "in", "if",
    "elif", "else", "with", "class", "None", "True", "False", "not", "and",
    "or", "while", "is", "lambda",
}

_TOKEN_COLOURS = {
    "comment": "#5C8C5C",
    "keyword": COLORS["Gs"],
    "string": COLORS["Gq"],
    "number": COLORS["G12-13"],
    "default": INK,
    "punct": MUTED,
    "ws": INK,
}


def _tokenize_line(line: str) -> list[tuple[str, str]]:
    tokens: list[tuple[str, str]] = []
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c == "#":
            tokens.append((line[i:], "comment"))
            break
        if c in '"\'':
            j = i + 1
            while j < n and line[j] != c:
                j += 1
            tokens.append((line[i : j + 1], "string"))
            i = j + 1
            continue
        if c.isdigit():
            j = i
            while j < n and (line[j].isdigit() or line[j] == "."):
                j += 1
            tokens.append((line[i:j], "number"))
            i = j
            continue
        if c.isalpha() or c == "_":
            j = i
            while j < n and (line[j].isalnum() or line[j] == "_"):
                j += 1
            word = line[i:j]
            kind = "keyword" if word in _CODE_KEYWORDS else "default"
            tokens.append((word, kind))
            i = j
            continue
        if c.isspace():
            j = i
            while j < n and line[j].isspace():
                j += 1
            tokens.append((line[i:j], "ws"))
            i = j
            continue
        tokens.append((c, "punct"))
        i += 1
    return tokens


def render_code(
    ax: plt.Axes,
    x: float,
    y_top: float,
    code: str,
    fontsize: float = 5.8,
    line_height: float = 0.030,
    char_width: float | None = None,
) -> float:
    """Render code one complete line at a time, returning y after the last line.

    Drawing every token as a separate text artist made glyph placement depend on
    backend-specific font metrics and caused the PDF/PNG text to overlap.  A
    single monospace artist per line is stable in both outputs.
    """
    y = y_top
    for line in code.split("\n"):
        stripped = line.lstrip()
        color = _TOKEN_COLOURS["comment"] if stripped.startswith("#") else INK
        ax.text(
            x,
            y,
            line,
            family="DejaVu Sans Mono",
            fontsize=fontsize,
            color=color,
            va="top",
            ha="left",
        )
        y -= line_height
    return y


# ---------------------------------------------------------------------------
# Figure 1 — cohort scope
# ---------------------------------------------------------------------------

def figure1() -> None:
    composition = read("Figure_1A_cohort_composition.csv")
    boundary = read("Figure_1B_release_boundary.csv")
    identifiers = read("Figure_1D_identifier_coverage.csv")
    fig = plt.figure(figsize=(7.15, 5.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.92])
    ax = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
    ]

    # ── Panel A: reserved for a representative MD snapshot ──
    ax[0].axis("off")
    ax[0].set_xlim(0, 1)
    ax[0].set_ylim(0, 1)
    panel(ax[0], "A", x=-0.04)

    # ── Panel B: cohort composition heatmap ──
    matrix = (
        composition.pivot(
            index="gpcr_class",
            columns="g_protein_family",
            values="systems",
        )
        .reindex(index=["A", "B"], columns=FAMILIES)
        .fillna(0)
    )
    maximum = matrix.to_numpy().max()
    ax[1].set_xlim(-0.5, 3.5)
    ax[1].set_ylim(-0.5, 1.5)
    for row in range(2):
        for column, family in enumerate(FAMILIES):
            count = int(matrix.iloc[row, column])
            strength = 0.30 + 0.60 * count / maximum
            base = np.array(to_rgb(COLORS[family]))
            fill = 1 - strength * (1 - base)
            ax[1].add_patch(
                Rectangle(
                    (column - 0.48, row - 0.48),
                    0.96,
                    0.96,
                    facecolor=fill,
                    edgecolor="white",
                    linewidth=1.6,
                    zorder=2,
                )
            )
            text_color = "white" if strength > 0.55 else INK
            ax[1].text(
                column,
                row,
                str(count),
                ha="center",
                va="center",
                fontweight="bold",
                fontsize=7.4,
                color=text_color,
                zorder=3,
            )
    ax[1].set_xticks(range(4), [FAMILY_LABELS[x] for x in FAMILIES])
    ax[1].set_yticks([0, 1], ["Class A", "Class B"])
    ax[1].tick_params(length=0)
    for spine in ax[1].spines.values():
        spine.set_visible(False)
    # family colour legend below
    handles = [
        Patch(facecolor=COLORS[x], edgecolor="none", label=FAMILY_LABELS[x])
        for x in FAMILIES
    ]
    ax[1].legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=4,
        frameon=False,
        handlelength=1.0,
        columnspacing=1.0,
    )
    ax[1].set_title("Frozen cohort composition")
    panel(ax[1], "B", x=-0.04)

    # ── Panel C: release boundary ──
    order = ["included", "unresolved", "excluded"]
    boundary = boundary.set_index("release_status").loc[order].reset_index()
    bar_colors = [TEAL, "#AAB1B7", RED]
    hatches = [None, "..", None]
    y_positions = np.arange(3)
    bars = ax[2].barh(y_positions, boundary.systems, color=bar_colors, height=0.60)
    for bar, hatch in zip(bars, hatches):
        if hatch:
            bar.set_hatch(hatch)
            bar.set_edgecolor("white")
            bar.set_linewidth(0.4)
            bar.set_alpha(0.85)
    ax[2].set_yticks(y_positions, ["Included", "Unresolved", "Excluded"])
    ax[2].invert_yaxis()
    ax[2].set_xlabel("working-inventory records")
    ax[2].set_xlim(0, 240)
    soft_grid(ax[2], "x")
    for value, ypos in zip(boundary.systems, y_positions):
        ax[2].text(
            value + 4,
            ypos,
            str(int(value)),
            va="center",
            fontweight="bold",
            fontsize=7.4,
        )
    ax[2].set_title("Release boundary")
    panel(ax[2], "C", x=-0.04)

    # ── Panel D: identifier coverage cards ──
    ax[3].axis("off")
    ax[3].set_xlim(0, 1)
    ax[3].set_ylim(0, 1)
    lookup = identifiers.set_index("metric").value.to_dict()
    cards = [
        ("systems", "systems", TEAL),
        ("production_replicas", "computational\nrepeats", BLUE),
        ("receptor_names", "receptor\nnames", GREEN),
        ("mapped_uniprot_accessions", "mapped UniProt\naccessions", "#7C6C9C"),
    ]
    accents = [TEAL, BLUE, GREEN, "#7C6C9C"]
    for index, (key, label, accent) in enumerate(zip(
        ["systems", "production_replicas", "receptor_names", "mapped_uniprot_accessions"],
        ["systems", "computational\nrepeats", "receptor\nnames", "mapped UniProt\naccessions"],
        accents,
    )):
        col = index % 2
        rowi = index // 2
        x = 0.04 + col * 0.49
        y = 0.54 - rowi * 0.45
        card(
            ax[3],
            (x, y),
            0.43,
            0.34,
            f"{int(lookup[key]):,}",
            label,
            accent=accent,
            number_size=10.5,
        )
    ax[3].set_title("Identifier and repeat coverage")
    panel(ax[3], "D", x=-0.04)
    save(fig, "figure1_cohort_scope")


# ---------------------------------------------------------------------------
# Figure 2 — record boundary and reuse model
# ---------------------------------------------------------------------------

def figure2() -> None:
    roles = read("Figure_2_record_roles.csv").set_index("record")
    ligand = read("Figure_2_retained_ligand_classification.csv")
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
        "protein complex + bound ligands + membrane + solvent + ions",
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
        "drop membrane,\nsolvent, ions",
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
    n_small = int((ligand.category == "small_molecule").sum())
    n_none = int((ligand.category == "none").sum())
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
        f"{n_small} with non-protein ligand  |  {n_none} protein-only",
        ha="center", va="top", fontsize=5.8, color=MUTED,
    )
    ax[0].set_title("Molecular-record boundary")
    panel(ax[0], "A", x=-0.04)

    # ── Panel B: support matrix ──
    rows = [
        ("Reduced trajectory", ["yes", "yes", "no", "no"]),
        ("Gateway summary", ["summary", "summary", "no", "yes"]),
        ("Full-system source", ["yes", "yes", "yes", "yes"]),
    ]
    columns = [
        "protein\ncomplex",
        "ligand\ncontact",
        "lipid / solvent\n/ ions",
        "recompute\nfrom source",
    ]
    x_positions = [0.40, 0.54, 0.69, 0.86]
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
        ax[1].plot([0.02, 0.89], [y - 0.072, y - 0.072], color=PALE,
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
            "GitHub · coupledmd",
            "HuangJianxiang-SJTU/coupledmd\ncode · schemas · workflows",
            BLUE,
        ),
        (
            "Zenodo · replica 1",
            "10.5281/zenodo.21395292\nreduced PDB/XTC record",
            AMBER,
        ),
    ]
    for index, (title, subtitle, color) in enumerate(access):
        y = 0.72 - index * 0.22
        box(ax[3], (0.10, y), 0.80, 0.17, title, subtitle,
            face="#F7F8F9", edge=color, title_size=7.1)
    ax[3].set_title("Access and reuse layers")
    panel(ax[3], "D", x=-0.04)
    save(fig, "figure2_records_reuse_boundary")


# ---------------------------------------------------------------------------
# Figure 3 — annotation validation
# ---------------------------------------------------------------------------

def _legacy_orthosteric_figure() -> None:
    controls = read("Figure_3A_positive_control.csv")
    completeness = read("Figure_3B_pocket_completeness.csv")
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.15, 3.25),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.35, 0.85]},
    )
    eligible = controls.sort_values(["g_family", "system_id"]).reset_index(drop=True)
    recovered = eligible.ortho_recovered.astype(bool)
    n_systems = len(eligible)

    # x-axis = best orthosteric-pocket frequency, y-axis = system index (58
    # peptide-ligand positive controls). The vertical dashed line at 0.85 is
    # the recovery threshold; recovered systems sit to its right.
    rng = np.random.default_rng(7)
    y_jitter = rng.normal(0, 0.32, n_systems)
    for family in FAMILIES:
        selected = eligible.g_family.eq(family)
        good = selected & recovered
        failed = selected & ~recovered
        axes[0].scatter(
            eligible.loc[good, "best_ortho_freq"],
            eligible.index[good] + y_jitter[good.to_numpy()],
            s=20,
            color=COLORS[family],
            label=FAMILY_LABELS[family],
            edgecolor="white",
            linewidth=0.3,
            alpha=0.85,
            zorder=3,
        )
        # not-recovered systems have no frequency; mark them at the left edge
        axes[0].scatter(
            np.full(failed.sum(), 0.805),
            eligible.index[failed] + y_jitter[failed.to_numpy()],
            s=26,
            marker="x",
            color="#9AA1A6",
            zorder=3,
        )
    axes[0].axvline(0.85, linestyle="--", linewidth=1.0, color=INK, zorder=2)
    axes[0].text(
        0.852,
        -3.2,
        "threshold 0.85",
        fontsize=6.4,
        color=INK,
        ha="left",
        va="top",
    )
    axes[0].set_xlim(0.78, 1.01)
    axes[0].set_ylim(n_systems, -3)  # invert so first system is at the top
    axes[0].set_yticks([])
    axes[0].set_xlabel("best orthosteric-pocket frequency")
    axes[0].set_ylabel("58 peptide-ligand positive controls")
    # recovery callout, placed in the lower-left away from the data cloud
    axes[0].text(
        0.812,
        n_systems - 4,
        "49/58 recovered\n(84.5%)",
        ha="left",
        va="top",
        fontweight="bold",
        fontsize=7.2,
        color=INK,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#F2F5F6",
                 edgecolor="#C8CFD3", linewidth=0.6),
    )
    legend_handles = [
        mpl.lines.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color=COLORS[family],
            label=FAMILY_LABELS[family],
            markersize=4.5,
        )
        for family in ["Gi", "Gs", "Gq"]
    ]
    legend_handles.append(
        mpl.lines.Line2D(
            [0],
            [0],
            marker="x",
            linestyle="",
            color="#9AA1A6",
            label="not recovered",
            markersize=5,
        )
    )
    axes[0].legend(
        handles=legend_handles,
        frameon=False,
        ncol=4,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.28),
    )
    axes[0].set_title("Orthosteric recovery as a technical positive control")
    panel(axes[0], "A")

    labels = ["Detected\npockets", "Valid zero-\npocket result"]
    values = [
        int(
            completeness.loc[
                completeness.record_status.eq("available_detected_pocket"),
                "systems",
            ].iloc[0]
        ),
        int(
            completeness.loc[
                completeness.record_status.eq("available_zero_pockets"),
                "systems",
            ].iloc[0]
        ),
    ]
    bars = axes[1].bar(
        [0, 1],
        values,
        color=[TEAL, "#B9C0C5"],
        width=0.62,
        zorder=3,
    )
    axes[1].set_xticks([0, 1], labels)
    axes[1].set_ylabel("systems")
    axes[1].set_ylim(0, 230)
    soft_grid(axes[1], "y")
    for bar, value in zip(bars, values):
        axes[1].text(
            bar.get_x() + bar.get_width() / 2,
            value + 5,
            str(value),
            ha="center",
            fontweight="bold",
        )
    axes[1].text(
        0.56,
        0.78,
        "2,149 detected-pocket rows",
        transform=axes[1].transAxes,
        ha="left",
        fontsize=7,
        color=MUTED,
    )
    axes[1].text(
        0.56,
        0.68,
        "205 + 2 = 207 complete records",
        transform=axes[1].transAxes,
        ha="left",
        fontsize=7,
        fontweight="bold",
    )
    # donut inset showing 205/2 split
    inset = axes[1].inset_axes([0.60, 0.08, 0.30, 0.30])
    inset.set_aspect("equal")
    inset.axis("off")
    inset.pie(
        values,
        colors=[TEAL, "#B9C0C5"],
        startangle=90,
        wedgeprops=dict(width=0.42, edgecolor="white", linewidth=1.0),
    )
    inset.text(0, 0, "207", ha="center", va="center", fontsize=7,
               fontweight="bold", color=INK)
    axes[1].set_title("Pocket-record completeness")
    panel(axes[1], "B")
    save(fig, "figure3_annotation_validation")


# ---------------------------------------------------------------------------
# Figure 4 — technical validation + QC
# ---------------------------------------------------------------------------

def violin_by_replica(
    ax: plt.Axes,
    data: pd.DataFrame,
    field: str,
    ylabel: str,
    title: str,
) -> None:
    replica_field = "replica" if "replica" in data.columns else "release_replica"
    values = [
        data.loc[data[replica_field].eq(replica), field].dropna().to_numpy()
        for replica in [1, 2, 3]
    ]
    plot = ax.violinplot(
        values,
        positions=[1, 2, 3],
        widths=0.72,
        showextrema=False,
        showmedians=False,
    )
    shades = ["#9CC0CE", "#5B91A7", "#2E5F76"]
    for body, color in zip(plot["bodies"], shades):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_linewidth(0.6)
        body.set_alpha(0.55)
    # strip-plot of individual points behind violins (light, jittered)
    rng = np.random.default_rng(0)
    for pos, vals in zip([1, 2, 3], values):
        jitter = rng.normal(0, 0.06, len(vals))
        ax.scatter(
            np.full(len(vals), pos) + jitter,
            vals,
            s=3.5,
            color=INK,
            alpha=0.25,
            zorder=2,
            edgecolor="none",
        )
    ax.boxplot(
        values,
        positions=[1, 2, 3],
        widths=0.18,
        showfliers=False,
        patch_artist=True,
        medianprops={"color": INK, "linewidth": 1.2},
        boxprops={"facecolor": "white", "edgecolor": TEAL, "linewidth": 0.8},
        whiskerprops={"color": INK, "linewidth": 0.6},
        capprops={"color": INK, "linewidth": 0.6},
    )
    ax.set_xticks([1, 2, 3], ["Rep 1", "Rep 2", "Rep 3"])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    soft_grid(ax, "y")


def _legacy_qc_figure() -> None:
    validation = read("Figure_4A_replica_validation.csv")
    qc_summary = read("Figure_4D_release_qc_summary.csv")
    sizes = read("Figure_4E_replica_file_sizes.csv")
    example = read("Figure_4F_worked_example.csv")
    available = validation[validation.validation_status.eq("available")].copy()
    n_avail = len(available)
    total_val = len(validation)
    fig, axes = plt.subplots(3, 2, figsize=(7.15, 8.2), constrained_layout=True)
    ax = axes.ravel()

    # ── Panels A–C: violin plots ──
    violin_by_replica(
        ax[0], available,
        "tm_core_rmsd_A_p95",
        "TM-core Cα RMSD P95 (Å)",
        f"Receptor TM-core displacement\n({n_avail}/{total_val} available; Gs_8HTI excluded)",
    )
    panel(ax[0], "A")

    violin_by_replica(
        ax[1], available,
        "galpha_interface_rmsd_A_p95",
        "Gα-interface Cα RMSD P95 (Å)",
        f"Gα interface-region displacement\n({n_avail}/{total_val} available; Gs_8HTI excluded)",
    )
    panel(ax[1], "B")

    violin_by_replica(
        ax[2], available,
        "contact_retention_p05",
        "initial contacts retained P05",
        f"Initial interface-contact retention\n({n_avail}/{total_val} available; Gs_8HTI excluded)",
    )
    ax[2].set_ylim(-0.03, 1.03)
    panel(ax[2], "C")

    # ── Panel D: reduced-release QC summary (clean heatmap table) ──
    ax[3].axis("off")
    ax[3].set_xlim(0, 1)
    ax[3].set_ylim(0, 1)

    col_xs = [0.42, 0.53, 0.64]
    col_labels = ["Rep 1", "Rep 2", "Rep 3"]
    row_labels = [
        ("2,500-frame target", "uniform_2500_frames"),
        ("PDB/XTC atom match", "pdb_xtc_atom_counts_match"),
        ("Finite coordinates", "finite_coordinates"),
        ("Monotonic time", "finite_monotonic_time"),
        ("Valid periodic box", "valid_periodic_box"),
        ("Chain break absent", "catastrophic_chain_break_absent"),
        ("Separation absent", "complex_separation_absent"),
        ("Scatter absent", "coordinate_scatter_absent"),
        ("Full-frame QC pass", "full_frame_qc_passed"),
    ]
    for x, label in zip(col_xs, col_labels):
        ax[3].text(x, 0.94, label, ha="center", fontsize=7.0, fontweight="bold")
    for row_idx, (label, field) in enumerate(row_labels):
        y = 0.85 - row_idx * 0.084
        ax[3].text(0.40, y, label, ha="right", va="center", fontsize=6.4)
        for col_idx in range(3):
            rep = col_idx + 1
            x = col_xs[col_idx]
            val = qc_summary.loc[qc_summary.release_replica == rep, field].iloc[0]
            status = "PASS" if val == 207 else f"{val}/207"
            ax[3].add_patch(
                FancyBboxPatch(
                    (x - 0.052, y - 0.030),
                    0.104,
                    0.054,
                    boxstyle="round,pad=.003,rounding_size=.008",
                    facecolor=GREEN,
                    edgecolor="white",
                    linewidth=0.6,
                )
            )
            ax[3].text(
                x - 0.030, y, "P", ha="center", va="center",
                fontsize=6.0, color="white", fontweight="bold",
            )
            ax[3].text(
                x + 0.012, y, status, ha="center", va="center",
                fontsize=5.3, color="white", fontweight="bold",
            )
    ax[3].text(
        0.5, 0.01,
        "All 207 systems × 3 replicas = 621 records pass every check.",
        ha="center", fontsize=6.2, color=MUTED,
    )
    ax[3].set_title("Reduced-release QC (full-frame)")
    panel(ax[3], "D", x=-0.04)

    # ── Panel E: per-replica XTC size and atom counts ──
    ax[4].axis("off")
    ax[4].set_xlim(0, 1)
    ax[4].set_ylim(0, 1)
    left = ax[4].inset_axes([0.02, 0.08, 0.46, 0.82])
    values_s = [
        sizes.loc[sizes.release_replica.eq(rep), "size_gb"].dropna().to_numpy()
        for rep in [1, 2, 3]
    ]
    left_v = left.violinplot(values_s, positions=[1, 2, 3], widths=0.72,
                              showextrema=False, showmedians=False)
    for body, c in zip(left_v["bodies"], ["#9CC0CE", "#5B91A7", "#2E5F76"]):
        body.set_facecolor(c); body.set_edgecolor(c); body.set_linewidth(0.6)
        body.set_alpha(0.55)
    left.boxplot(values_s, positions=[1, 2, 3], widths=0.18, showfliers=False,
                 patch_artist=True,
                 medianprops={"color": INK, "linewidth": 1.2},
                 boxprops={"facecolor": "white", "edgecolor": BLUE, "linewidth": 0.8},
                 whiskerprops={"color": INK, "linewidth": 0.6},
                 capprops={"color": INK, "linewidth": 0.6})
    rng = np.random.default_rng(1)
    for pos, vals in zip([1, 2, 3], values_s):
        jitter = rng.normal(0, 0.05, len(vals))
        left.scatter(np.full(len(vals), pos) + jitter, vals, s=3.5,
                     color=INK, alpha=0.25, edgecolor="none", zorder=2)
    left.set_xticks([1, 2, 3], ["Rep 1", "Rep 2", "Rep 3"], fontsize=6.5)
    left.set_ylabel("XTC size (GB)", fontsize=7.5)
    soft_grid(left, "y")
    left.tick_params(labelsize=6.5)
    right = ax[4].inset_axes([0.52, 0.08, 0.46, 0.82])
    values_a = [
        sizes.loc[sizes.release_replica.eq(rep), "n_atoms"].dropna().to_numpy()
        for rep in [1, 2, 3]
    ]
    right_v = right.violinplot(values_a, positions=[1, 2, 3], widths=0.72,
                                showextrema=False, showmedians=False)
    for body, c in zip(right_v["bodies"], ["#9CC0CE", "#5B91A7", "#2E5F76"]):
        body.set_facecolor(c); body.set_edgecolor(c); body.set_linewidth(0.6)
        body.set_alpha(0.55)
    right.boxplot(values_a, positions=[1, 2, 3], widths=0.18, showfliers=False,
                  patch_artist=True,
                  medianprops={"color": INK, "linewidth": 1.2},
                  boxprops={"facecolor": "white", "edgecolor": TEAL, "linewidth": 0.8},
                  whiskerprops={"color": INK, "linewidth": 0.6},
                  capprops={"color": INK, "linewidth": 0.6})
    for pos, vals in zip([1, 2, 3], values_a):
        jitter = rng.normal(0, 0.05, len(vals))
        right.scatter(np.full(len(vals), pos) + jitter, vals, s=3.5,
                      color=INK, alpha=0.25, edgecolor="none", zorder=2)
    right.set_xticks([1, 2, 3], ["Rep 1", "Rep 2", "Rep 3"], fontsize=6.5)
    right.set_ylabel("retained atoms", fontsize=7.5)
    soft_grid(right, "y")
    right.tick_params(labelsize=6.5)
    ax[4].set_title("Reduced XTC size and atom counts")
    panel(ax[4], "E", x=-0.04)

    # ── Panel F: worked checksum-and-load example (G12_7SF7) ──
    ax[5].axis("off")
    ax[5].set_xlim(0, 1)
    ax[5].set_ylim(0, 1)
    pdb = example[example.file_role.eq("structure_pdb")].iloc[0]
    xtc_ex = example[example.file_role.eq("trajectory_xtc")].iloc[0]
    code = (
        f"# G12_7SF7, release replica 1\n"
        f"sha256sum G12_7SF7/structure.pdb\n"
        f"# expected: {str(pdb.sha256)[:18]}…\n"
        f"sha256sum G12_7SF7/traj.xtc\n"
        f"# expected: {str(xtc_ex.sha256)[:18]}…\n\n"
        "import MDAnalysis as mda\n"
        "u = mda.Universe(\"G12_7SF7/structure.pdb\", \"G12_7SF7/traj.xtc\")\n"
        f"assert u.atoms.n_atoms == {int(pdb.n_atoms)}\n"
        "assert len(u.trajectory) == 2500"
    )
    ax[5].add_patch(
        FancyBboxPatch(
            (0.04, 0.10), 0.92, 0.75,
            boxstyle="round,pad=.018",
            facecolor="#F7F8F9", edgecolor="#C8CFD3", linewidth=0.8,
        )
    )
    render_code(ax[5], 0.08, 0.78, code, fontsize=5.8, line_height=0.058,
               char_width=0.0050)
    ax[5].text(
        0.5, 0.03,
        "Verify both digests before analysis; load the PDB as the XTC topology.",
        ha="center", fontsize=6.2, color=MUTED,
    )
    ax[5].set_title("Checksum-first reuse: G12_7SF7 (G12/13)")
    panel(ax[5], "F", x=-0.04)

    save(fig, "legacy_figure4_technical_validation_qc")


# ---------------------------------------------------------------------------
# Revised Figure 3 — QC and stability (main-text panels only)
# ---------------------------------------------------------------------------

def figure3() -> None:
    validation = read("Figure_4A_replica_validation.csv")
    available = validation[validation.validation_status.eq("available")].copy()
    n_avail = len(available)
    total_val = len(validation)
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 3.15), constrained_layout=True)

    specs = [
        (
            "tm_core_rmsd_A_p95",
            "TM-core Cα RMSD P95 (Å)",
            "Receptor TM-core\ndisplacement",
        ),
        (
            "galpha_interface_rmsd_A_p95",
            "Gα-interface Cα RMSD P95 (Å)",
            "Gα interface-region\ndisplacement",
        ),
        (
            "contact_retention_p05",
            "Initial contacts retained P05",
            "Interface-contact\nretention",
        ),
    ]
    for letter, ax_i, (field, ylabel, title) in zip("ABC", axes, specs):
        violin_by_replica(ax_i, available, field, ylabel, title)
        values = available[field].dropna().to_numpy()
        if field == "contact_retention_p05":
            ax_i.set_ylim(0.0, 1.04)
        else:
            lower = max(0.0, float(values.min()) - 0.06 * np.ptp(values))
            upper = float(values.max()) + 0.08 * np.ptp(values)
            ax_i.set_ylim(lower, upper)
        ax_i.tick_params(axis="both", labelsize=7.2)
        panel(ax_i, letter, x=-0.18)
    save(fig, "figure3_qc_stability")


# ---------------------------------------------------------------------------
# Revised Figure 4 — orthosteric recovery
# ---------------------------------------------------------------------------

def figure4() -> None:
    controls = read("Figure_3A_positive_control.csv")
    ranked = controls.copy()
    ranked["plot_frequency"] = ranked.best_ortho_freq.fillna(0.0)
    ranked = ranked.sort_values(
        ["plot_frequency", "g_family", "system_id"],
        ascending=[True, True, True],
    ).reset_index(drop=True)
    present_families = [
        family for family in FAMILIES
        if ranked.g_family.eq(family).any()
    ]

    fig = plt.figure(figsize=(7.15, 3.85))
    outer = fig.add_gridspec(
        1, 2, width_ratios=[1.05, 0.95],
        left=0.09, right=0.98, bottom=0.16, top=0.89, wspace=0.30,
    )

    # Frequencies occupy two separated ranges (zero and >=0.85).  Use the
    # same explicit 0.1--0.8 axis break in both panels so the nine misses stay
    # visible without devoting most of the figure to an empty interval.
    overall_grid = outer[0].subgridspec(
        1, 2, width_ratios=[0.38, 0.62], wspace=0.08,
    )
    ax_overall_left = fig.add_subplot(overall_grid[0])
    ax_overall_right = fig.add_subplot(
        overall_grid[1], sharey=ax_overall_left,
    )
    overall_axes = (ax_overall_left, ax_overall_right)

    family_grid = outer[1].subgridspec(
        len(present_families), 1, hspace=0.12,
    )
    family_axes = []
    for index in range(len(present_families)):
        row = family_grid[index].subgridspec(
            1, 2, width_ratios=[0.38, 0.62], wspace=0.08,
        )
        ax_left = fig.add_subplot(row[0])
        ax_right = fig.add_subplot(row[1], sharey=ax_left)
        family_axes.append((ax_left, ax_right))

    def format_broken_pair(ax_left, ax_right, *, bottom_marks=True):
        ax_left.set_xlim(-0.005, 0.105)
        ax_right.set_xlim(0.795, 1.005)
        ax_left.set_xticks([0.0, 0.1])
        ax_right.set_xticks([0.8, 0.9, 1.0])
        ax_left.spines["right"].set_visible(False)
        ax_right.spines["left"].set_visible(False)
        ax_right.tick_params(axis="y", left=False, labelleft=False)

        # Small diagonal marks make the omitted 0.1--0.8 interval explicit.
        mark = dict(
            color=INK, clip_on=False, linewidth=0.8,
            transform=ax_left.transAxes,
        )
        if bottom_marks:
            ax_left.plot((0.975, 1.025), (-0.018, 0.018), **mark)
            ax_right.plot(
                (-0.025, 0.025), (-0.018, 0.018),
                color=INK, clip_on=False, linewidth=0.8,
                transform=ax_right.transAxes,
            )

    # Panel A: overall ranked distribution. Not-recovered systems are retained
    # as zero-frequency points on the same 0–1 axis.
    for ax_overall in overall_axes:
        for family in present_families:
            selected = ranked.g_family.eq(family)
            ax_overall.scatter(
                ranked.loc[selected, "plot_frequency"],
                ranked.index[selected] + 1,
                s=22,
                color=COLORS[family],
                edgecolor="white",
                linewidth=0.35,
                zorder=3,
                clip_on=False,
            )
        ax_overall.set_ylim(0, len(ranked) + 1)
        soft_grid(ax_overall, "x")
    ax_overall_right.axvline(
        0.85, linestyle="--", linewidth=1.0, color=INK, zorder=2,
    )
    format_broken_pair(ax_overall_left, ax_overall_right)
    ax_overall_left.set_yticks([1, 10, 20, 30, 40, 50, 58])
    ax_overall_left.set_ylabel("Ranked peptide controls")
    legend_handles = [
        mpl.lines.Line2D(
            [0], [0], marker="o", linestyle="", markersize=4.8,
            markerfacecolor=COLORS[family], markeredgecolor="white",
            label=FAMILY_LABELS[family],
        )
        for family in present_families
    ]
    ax_overall_left.legend(
        handles=legend_handles, frameon=False, ncol=2,
        loc="upper left", columnspacing=0.9, handletextpad=0.35,
        bbox_to_anchor=(0.0, 1.01), fontsize=6.6,
    )
    ax_overall_left.set_title(
        "Overall recovery (n=58)", loc="left", fontsize=8.4, pad=7,
    )
    panel(ax_overall_left, "A", x=-0.25, y=1.10)

    # Panel B: dot-only family distributions on the same frequency axis.
    for index, (family, axes_family) in enumerate(
        zip(present_families, family_axes)
    ):
        ax_left, ax_right = axes_family
        values = np.sort(
            ranked.loc[ranked.g_family.eq(family), "plot_frequency"].to_numpy()
        )
        n_family = len(values)
        ranked_y = np.arange(1, n_family + 1)
        ax_right.axvline(
            0.85, linestyle="--", linewidth=0.9, color=INK, zorder=2,
        )
        for ax_family in axes_family:
            ax_family.scatter(
                values, ranked_y, s=20, color=COLORS[family],
                edgecolor="white", linewidth=0.35, zorder=4,
                clip_on=False,
            )
            ax_family.set_ylim(-0.5, n_family + 1.5)
            soft_grid(ax_family, "x")
        ax_left.set_yticks([1, n_family])
        ax_left.set_yticklabels(["1", str(n_family)])
        ax_left.tick_params(axis="y", labelsize=6.3, pad=2)
        format_broken_pair(ax_left, ax_right)
        ax_left.text(
            0.04, 0.88, f"{FAMILY_LABELS[family]}  n={n_family}",
            transform=ax_left.transAxes, ha="left", va="center",
            fontsize=6.8, fontweight="bold", color=COLORS[family],
        )
        if index < len(family_axes) - 1:
            ax_left.tick_params(axis="x", labelbottom=False)
            ax_right.tick_params(axis="x", labelbottom=False)
    family_axes[0][0].set_title(
        "Family distributions", loc="left", fontsize=8.4, pad=7,
    )
    panel(family_axes[0][0], "B", x=-0.25, y=1.13)

    fig.text(0.287, 0.055, "Best orthosteric-pocket frequency", ha="center")
    fig.text(0.760, 0.055, "Best orthosteric-pocket frequency", ha="center")
    fig.text(
        0.565, 0.52, "Within-family rank", rotation=90,
        ha="center", va="center",
    )
    save(fig, "figure4_orthosteric_recovery")


# ---------------------------------------------------------------------------
# Figure 5 — example data reuse
# ---------------------------------------------------------------------------

def figure5() -> None:
    """Example data reuse: contact fingerprint + worked reuse."""
    contact = np.load(SOURCE / "Figure_5A_contact_matrix_full.npy")
    res_info = read("Figure_5A_contact_params.csv")
    rmsd = read("Figure_5B_rmsd_trajectory.csv")
    release_manifest = read("Figure_5B_release_file_manifest.csv")

    fig = plt.figure(figsize=(7.15, 8.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.18, 0.82], width_ratios=[1.25, 0.75],
                          left=0.10, right=0.97, bottom=0.06, top=0.95,
                          hspace=0.30, wspace=0.30)

    # ── Panel A: top-30 contact-fingerprint heatmap + frequency bar ──
    gs_a = gs[0, :].subgridspec(1, 2, width_ratios=[5.0, 1.15], wspace=0.05)
    ax_a = fig.add_subplot(gs_a[0, 0])
    ax_a_freq = fig.add_subplot(gs_a[0, 1], sharey=ax_a)

    res_info = res_info.sort_values("order").reset_index(drop=True)
    top_res = res_info.sort_values(
        ["persistence", "order"], ascending=[False, True]
    ).head(30).copy()
    top_indices = top_res["order"].astype(int).to_numpy()
    contact = contact[:, top_indices]
    top_res = top_res.reset_index(drop=True)
    n_res = len(top_res)
    n_frames_plot = contact.shape[0]

    step = 10
    contact_sub = contact[::step, :].T
    extent = [0, n_frames_plot * 0.2, n_res, 0]

    ax_a.imshow(contact_sub, aspect="auto", cmap=CONTACT_CMAP,
                extent=extent, vmin=0, vmax=1, interpolation="nearest")

    ax_a.set_xlim(0, n_frames_plot * 0.2)
    ax_a.set_xlabel("Time (ns)", fontsize=8)
    ax_a.set_ylabel("Contacting protein residue", fontsize=8)
    ax_a.set_yticks(np.arange(n_res) + 0.5)
    ax_a.set_yticklabels(
        [f"{r.resname}{int(r.resid)}" for _, r in top_res.iterrows()],
        fontsize=6.4,
    )
    ax_a.tick_params(axis="x", labelsize=6.8)

    persistence = top_res.persistence.to_numpy()
    y_positions = np.arange(n_res) + 0.5
    ax_a_freq.barh(
        y_positions, persistence, height=0.72,
        color=[TEAL if value >= 0.5 else "#B9C0C5" for value in persistence],
        edgecolor="none", zorder=3,
    )
    ax_a_freq.set_xlim(0, 1.04)
    ax_a_freq.set_xlabel("Persistence", fontsize=7.2)
    ax_a_freq.set_xticks([0.5, 1.0])
    ax_a_freq.tick_params(axis="x", labelsize=6.3)
    ax_a_freq.tick_params(axis="y", left=False, labelleft=False)
    soft_grid(ax_a_freq, "x")

    ax_a.set_title("Top 30 ligand–protein contacts: S1P in G12_7T6B (2,500 frames, 500 ns)")

    # ── Panel B: worked reuse example ──
    ax_b_code = fig.add_subplot(gs[1, 0])
    ax_b_code.axis("off")
    ax_b_code.set_xlim(0, 1)
    ax_b_code.set_ylim(0, 1)

    manifest_rows = release_manifest[
        release_manifest.system_id.eq("G12_7T6B")
        & release_manifest.release_replica.eq(1)
    ]
    pdb_row = manifest_rows[manifest_rows.file_role.eq("structure_pdb")].iloc[0]
    xtc_row = manifest_rows[manifest_rows.file_role.eq("trajectory_xtc")].iloc[0]
    pdb_sha = str(pdb_row.sha256)
    xtc_sha = str(xtc_row.sha256)
    n_atoms = int(pdb_row.n_atoms)
    all_resids = res_info.sort_values("resid").resid.astype(int).tolist()
    resid_lines = [
        " ".join(str(value) for value in all_resids[i : i + 10])
        for i in range(0, len(all_resids), 10)
    ]
    resid_block = "\n".join(
        f'    "{line}{" " if index < len(resid_lines) - 1 else ")"}"'
        for index, line in enumerate(resid_lines)
    )

    code = (
        "# G12_7T6B, release replica 1\n"
        "# Verify the manifest SHA-256 values\n"
        "$ sha256sum structure.pdb\n"
        f"{pdb_sha}\n"
        "$ sha256sum traj.xtc\n"
        f"{xtc_sha}\n"
        "import MDAnalysis as mda\n"
        'u = mda.Universe("structure.pdb", "traj.xtc")\n'
        f"assert u.atoms.n_atoms == {n_atoms}\n"
        "assert len(u.trajectory) == 2500\n\n"
        "# Select binding-pocket Cα atoms\n"
        "# (all 48 residues that contact the S1P ligand)\n"
        'pocket = u.select_atoms(\n'
        '    "name CA and (resid "\n'
        f"{resid_block}\n"
        ")\n"
        "from MDAnalysis.analysis import rms\n"
        "rmsd_analysis = rms.RMSD(\n"
        "    u, pocket, select=pocket, ref_frame=0\n"
        ").run()\n"
    )
    ax_b_code.add_patch(
        FancyBboxPatch(
            (0.035, 0.03), 1.010, 0.93,
            boxstyle="round,pad=.015",
            facecolor="#F7F8F9", edgecolor="#C8CFD3", linewidth=0.7,
            clip_on=False,
        )
    )
    ax_b_code.text(
        0.060, 0.90, "Reuse code",
        fontsize=7.5, fontweight="bold", va="top", color=INK,
    )
    render_code(ax_b_code, 0.050, 0.84, code, fontsize=6.0,
               line_height=0.032)

    ax_b_plot = fig.add_subplot(gs[1, 1])
    step = max(1, len(rmsd) // 250)
    rmsd_plot = rmsd.iloc[::step]
    ax_b_plot.fill_between(
        rmsd_plot.time_ns, 0, rmsd_plot.ca_rmsd_fit_A,
        color=BLUE, alpha=0.12, zorder=1,
    )
    ax_b_plot.plot(
        rmsd_plot.time_ns, rmsd_plot.ca_rmsd_fit_A,
        color=BLUE, linewidth=0.4, alpha=0.6, zorder=2,
    )
    smooth = rmsd.ca_rmsd_fit_A.rolling(window=51, center=True).mean()
    ax_b_plot.plot(rmsd.time_ns, smooth, color=INK, linewidth=1.3, zorder=3)
    mean_rmsd = float(rmsd.ca_rmsd_fit_A.mean())
    ax_b_plot.axhline(mean_rmsd, linestyle="--", linewidth=0.8, color=AMBER,
                     zorder=2)
    ax_b_plot.set_xlabel("Time (ns)", fontsize=8)
    ax_b_plot.set_ylabel("Cα RMSD after fit (Å)", fontsize=8)
    ax_b_plot.set_ylim(0, 2.5)
    soft_grid(ax_b_plot, "y")
    ax_b_plot.tick_params(labelsize=6.5)
    from matplotlib.lines import Line2D
    rmsd_legend = [
        Line2D([0], [0], color=BLUE, linewidth=0.4, alpha=0.6, label="per-frame"),
        Line2D([0], [0], color=INK, linewidth=1.3, label="rolling mean (51-frame)"),
        Line2D([0], [0], color=AMBER, linewidth=0.8, linestyle="--",
               label=f"mean = {mean_rmsd:.2f} Å"),
    ]
    ax_b_plot.legend(handles=rmsd_legend, fontsize=5.8, frameon=False,
                     loc="lower right")
    ax_b_plot.set_title(
        "Binding-pocket Cα RMSD: G12_7T6B\n"
        "(fit to 48 pocket residues, mean=1.35 Å, P95=1.66 Å)",
        fontsize=8.5,
    )
    # Figure-level labels keep A/B vertically aligned and B/C horizontally
    # aligned despite the unequal panel widths and the two-line C title.
    for letter, x, y in [
        ("A", 0.075, 0.962),
        ("B", 0.075, 0.430),
        ("C", 0.615, 0.430),
    ]:
        fig.text(
            x, y, letter,
            fontsize=10.5, fontweight="bold",
            ha="left", va="bottom", color=INK,
        )

    save(fig, "figure5_example_data_reuse")


# ---------------------------------------------------------------------------
# Supplementary figures
# ---------------------------------------------------------------------------

def supplementary_figure1() -> None:
    data = read("Supplementary_Figure_1_gateway_provenance.csv").sort_values(
        "stage_order"
    )
    fig, ax = plt.subplots(figsize=(7.15, 3.25), constrained_layout=True)
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    stages = data.iloc[:3]
    positions = [0.03, 0.36, 0.69]
    for x, row in zip(positions, stages.itertuples(index=False)):
        box(
            ax,
            (x, 0.55),
            0.27,
            0.23,
            "\n".join(textwrap.wrap(str(row.record).upper(), 23)),
            "\n".join(textwrap.wrap(str(row.role), 34)),
            face="#F2F5F5",
            edge=TEAL,
            title_size=7.0,
        )
    arrow(ax, (0.30, 0.665), (0.355, 0.665))
    arrow(ax, (0.63, 0.665), (0.685, 0.665))
    reduced = data.iloc[3]
    box(
        ax,
        (0.24, 0.13),
        0.52,
        0.20,
        str(reduced.record).upper(),
        str(reduced.role),
        face="#EAF0F7",
        edge=BLUE,
    )
    arrow(ax, (0.50, 0.34), (0.50, 0.51), color=RED, dashed=True)
    ax.text(
        0.54,
        0.43,
        "cannot regenerate\nlipid-gateway values",
        color=RED,
        fontsize=6.7,
        va="center",
    )
    ax.set_title(
        "Gateway-summary provenance and the reduced-record reproducibility boundary",
        pad=10,
    )
    save(fig, "supplementary_figure_s1_gateway_provenance")


def supplementary_figure2() -> None:
    data = read("Supplementary_Figure_2_sampling_provenance.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.2), constrained_layout=True)
    campaign = data[data.scope.eq("full final-207 campaign")]
    axes[0].bar(
        np.arange(len(campaign)),
        campaign.sampling_us,
        color=[TEAL, AMBER],
        width=0.62,
        zorder=3,
    )
    axes[0].set_xticks(
        np.arange(len(campaign)),
        ["Nominal\nprotocol", "Original file-\nobserved"],
    )
    axes[0].set_ylim(310.35, 310.55)
    axes[0].set_ylabel("aggregate sampling (µs)")
    soft_grid(axes[0], "y")
    for index, value in enumerate(campaign.sampling_us):
        axes[0].text(
            index,
            value + 0.005,
            f"{value:.4f}",
            ha="center",
            fontweight="bold",
        )
    axes[0].text(
        0.5,
        0.05,
        "difference = 83.2 ns",
        transform=axes[0].transAxes,
        ha="center",
        color=MUTED,
    )
    axes[0].set_title("Final-207 campaign")
    panel(axes[0], "A")

    record = data[data.scope.eq("Gq_8ZPT replica 2")]
    order = [
        "nominal protocol",
        "original file-observed coordinates",
        "locally repaired reduced source",
    ]
    record = record.set_index("sampling_definition").loc[order].reset_index()
    bars = axes[1].bar(
        np.arange(3),
        record.sampling_us,
        color=[TEAL, AMBER, BLUE],
        width=0.62,
        zorder=3,
    )
    axes[1].set_xticks(
        np.arange(3),
        ["Nominal", "Original\nobserved", "Repaired\nreduced source"],
    )
    axes[1].set_ylabel("coordinate span (µs)")
    axes[1].set_ylim(0, 0.56)
    soft_grid(axes[1], "y")
    for bar, value in zip(bars, record.sampling_us):
        axes[1].text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.012,
            f"{value:.5f}".rstrip("0"),
            ha="center",
            fontsize=6.7,
            fontweight="bold",
        )
    axes[1].text(
        0.5,
        0.05,
        "Local repair does not establish distribution of a full-system archive.",
        transform=axes[1].transAxes,
        ha="center",
        fontsize=6.2,
        color=MUTED,
        wrap=True,
    )
    axes[1].set_title("Gq_8ZPT replica 2 provenance")
    panel(axes[1], "B")
    save(fig, "supplementary_figure_s2_sampling_provenance")


def supplementary_figure3() -> None:
    """QC matrix and reduced-file scale moved out of main Figure 3."""
    qc_summary = read("Figure_4D_release_qc_summary.csv")
    sizes = read("Figure_4E_replica_file_sizes.csv")
    fig = plt.figure(figsize=(7.15, 4.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1.0])
    ax_table = fig.add_subplot(gs[0, :])
    ax_size = fig.add_subplot(gs[1, 0])
    ax_atoms = fig.add_subplot(gs[1, 1])

    ax_table.axis("off")
    ax_table.set_xlim(0, 1)
    ax_table.set_ylim(0, 1)
    col_xs = [0.47, 0.62, 0.77]
    row_labels = [
        ("2,500-frame target", "uniform_2500_frames"),
        ("PDB/XTC atom match", "pdb_xtc_atom_counts_match"),
        ("Finite coordinates", "finite_coordinates"),
        ("Monotonic time", "finite_monotonic_time"),
        ("Valid periodic box", "valid_periodic_box"),
        ("Chain break absent", "catastrophic_chain_break_absent"),
        ("Separation absent", "complex_separation_absent"),
        ("Coordinate scatter absent", "coordinate_scatter_absent"),
        ("Full-frame QC pass", "full_frame_qc_passed"),
    ]
    for x, label in zip(col_xs, ["Rep 1", "Rep 2", "Rep 3"]):
        ax_table.text(x, 0.94, label, ha="center", va="center",
                      fontsize=7.3, fontweight="bold")
    for row_idx, (label, field) in enumerate(row_labels):
        y = 0.84 - row_idx * 0.086
        ax_table.text(0.39, y, label, ha="right", va="center", fontsize=6.5)
        for x, rep in zip(col_xs, [1, 2, 3]):
            value = int(qc_summary.loc[qc_summary.release_replica.eq(rep), field].iloc[0])
            status = "PASS" if value == 207 else f"{value}/207"
            ax_table.add_patch(FancyBboxPatch(
                (x - 0.055, y - 0.029), 0.11, 0.058,
                boxstyle="round,pad=.003,rounding_size=.007",
                facecolor=GREEN if value == 207 else AMBER,
                edgecolor="white", linewidth=0.6,
            ))
            ax_table.text(x, y, status, ha="center", va="center",
                          fontsize=5.8, color="white", fontweight="bold")
    ax_table.set_title("Full-frame reduced-release QC: 207 systems × 3 replicas")
    panel(ax_table, "A", x=-0.01, y=0.98)

    violin_by_replica(ax_size, sizes, "size_gb", "XTC size (GB)",
                      "Reduced trajectory size")
    panel(ax_size, "B")
    violin_by_replica(ax_atoms, sizes, "n_atoms", "Retained atoms",
                      "Protein-complex atom count")
    panel(ax_atoms, "C")
    save(fig, "supplementary_figure_s3_release_qc_scale")

