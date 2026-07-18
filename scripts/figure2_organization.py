#!/usr/bin/env python3
"""Generate CoupledMD Figure 2 (released molecular-data organization)."""

from __future__ import annotations

from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import FancyBboxPatch
import matplotlib.pyplot as plt

from _common import (
    cli,
    INK,
    MUTED,
    panel_label,
    read,
    save,
)


def diagram_node(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    subtitle: str = "",
    face: str = "#EEF4FA",
    edge: str = "#4C72B0",
    title_size: float = 8.2,
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=.014",
            fc=face,
            ec=edge,
            lw=0.85,
        )
    )
    ax.text(
        x + width / 2,
        y + height * 0.63,
        title,
        ha="center",
        va="center",
        fontweight="bold",
        fontsize=title_size,
    )
    ax.text(
        x + width / 2,
        y + height * 0.27,
        subtitle,
        ha="center",
        va="center",
        fontsize=6.4,
        color=MUTED,
        linespacing=1.15,
    )


def diagram_arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = MUTED,
    style: str = "-|>",
    width: float = 0.9,
    size: float = 8,
    connection: str = "arc3",
    linestyle: str = "-",
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle=style,
            mutation_scale=size,
            lw=width,
            color=color,
            connectionstyle=connection,
            linestyle=linestyle,
            shrinkA=0,
            shrinkB=0,
            clip_on=False,
        )
    )


def figure2(data_dir: Path, output_dir: Path, _assets_dir: Path) -> None:
    facts = read(data_dir, "Figure_2_release_architecture_source.csv").set_index("quantity")
    assert int(float(facts.loc["included systems", "value"])) == 207
    assert int(float(facts.loc["released molecular files", "value"])) == 414
    assert int(float(facts.loc["frames per released trajectory", "value"])) == 2500

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.8), constrained_layout=True)
    ax = axes.ravel()
    for item in ax:
        item.axis("off")
        item.set(xlim=(0, 1), ylim=(0, 1))

    stages = [
        (0.80, "release metadata + checksums", "cohort boundary - file manifest", "#E8EFF8"),
        (0.57, "four family ZIP archives", "Gi/o - Gs - Gq/11 - G12/13", "#EAF3EC"),
        (0.34, "207 system directories", "stable system_id join key", "#F3EEF8"),
        (
            0.11,
            "structure.pdb + traj.xtc",
            "matched retained-component pair - 2,500 frames",
            "#FCF0E8",
        ),
    ]
    for y, title, subtitle, face in stages:
        diagram_node(ax[0], 0.12, y, 0.76, 0.13, title, subtitle, face=face)
    for upper, lower in zip(stages[:-1], stages[1:]):
        diagram_arrow(ax[0], (0.5, upper[0] - 0.012), (0.5, lower[0] + 0.155), width=1, size=9)
    ax[0].text(
        0.5,
        0.035,
        "derived from replica 1; full-system source trajectories are not distributed",
        ha="center",
        va="center",
        fontsize=6.4,
        color=MUTED,
    )
    panel_label(ax[0], "A", -0.04, 1.02)
    ax[0].set_title("Released molecular-data hierarchy")

    diagram_node(
        ax[1], 0.02, 0.38, 0.27, 0.27, "SYSTEM", "system_id\nreceptor - class - family", face="#F7F8FA"
    )
    diagram_node(
        ax[1],
        0.46,
        0.62,
        0.45,
        0.22,
        "SOURCE CAMPAIGN",
        "replicas 1-3 - 500 ns each\nprovenance and derived analyses",
        face="#F3EEF8",
        edge="#8172B2",
    )
    diagram_node(
        ax[1],
        0.46,
        0.17,
        0.45,
        0.22,
        "RELEASED FILE PAIR",
        "replica-1-derived PDB + XTC\nfile role - size - SHA-256",
        face="#EAF3EC",
        edge="#55A868",
    )
    diagram_arrow(ax[1], (0.30, 0.56), (0.45, 0.72), color=INK, width=1, size=9)
    diagram_arrow(ax[1], (0.30, 0.47), (0.45, 0.28), color=INK, width=1, size=9)
    diagram_arrow(
        ax[1],
        (0.685, 0.605),
        (0.685, 0.405),
        color=MUTED,
        width=0.8,
        size=8,
        linestyle="--",
    )
    ax[1].text(0.71, 0.505, "replica 1 selected", fontsize=6.3, color=MUTED, va="center")
    panel_label(ax[1], "B", -0.04, 1.02)
    ax[1].set_title("Stable identities and source provenance")

    diagram_node(
        ax[2],
        0.05,
        0.65,
        0.40,
        0.20,
        "REDUCED DATASET",
        "PDB/XTC + manifest\nfull-frame QC",
        face="#EAF3EC",
        edge="#55A868",
    )
    diagram_node(
        ax[2],
        0.55,
        0.65,
        0.40,
        0.20,
        "PORTAL + REST API",
        "search + NGL\nmetadata + annotations",
        face="#E8EFF8",
    )
    diagram_node(
        ax[2],
        0.30,
        0.10,
        0.40,
        0.20,
        "CODE REPOSITORY",
        "processing + validation\nfigure + source-data scripts",
        face="#F3EEF8",
        edge="#8172B2",
    )
    diagram_node(
        ax[2],
        0.32,
        0.39,
        0.36,
        0.12,
        "STABLE LINKS",
        "system_id + checksums\ncitation",
        face="#F7F8FA",
        edge="#737A82",
        title_size=7.4,
    )
    diagram_arrow(ax[2], (0.25, 0.64), (0.40, 0.52), width=0.8, size=8)
    diagram_arrow(ax[2], (0.75, 0.64), (0.60, 0.52), width=0.8, size=8)
    diagram_arrow(ax[2], (0.50, 0.38), (0.50, 0.315), width=0.8, size=8)
    panel_label(ax[2], "C", -0.04, 1.02)
    ax[2].set_title("Linked release and access layers")

    diagram_node(ax[3], 0.08, 0.82, 0.84, 0.11, "SEARCH", "system ID - receptor - metadata", face="#F7F8FA")
    diagram_node(ax[3], 0.08, 0.61, 0.84, 0.11, "FILTER", "family - class - ligand context", face="#F7F8FA")
    diagram_node(ax[3], 0.08, 0.40, 0.84, 0.11, "SYSTEM RECORD", "metadata + derived annotations", face="#E8EFF8")
    diagram_node(ax[3], 0.02, 0.09, 0.29, 0.15, "NGL VIEW", "reduced PDB/XTC", face="#EAF3EC", edge="#55A868", title_size=7.5)
    diagram_node(ax[3], 0.355, 0.09, 0.29, 0.15, "API ENTRY", "machine-readable\nrecords", face="#EAF3EC", edge="#55A868", title_size=7.5)
    diagram_node(ax[3], 0.69, 0.09, 0.29, 0.15, "DATASET", "checksummed files", face="#EAF3EC", edge="#55A868", title_size=7.5)
    diagram_arrow(ax[3], (0.5, 0.805), (0.5, 0.735), width=1, size=9)
    diagram_arrow(ax[3], (0.5, 0.595), (0.5, 0.525), width=1, size=9)
    for endpoint in (0.165, 0.50, 0.835):
        diagram_arrow(ax[3], (0.5, 0.385), (endpoint, 0.255), width=0.8, size=8)
    panel_label(ax[3], "D", -0.04, 1.02)
    ax[3].set_title("Portal-to-reuse workflow")
    save(fig, output_dir, "figure2_organization")


def main() -> None:
    args = cli("Generate CoupledMD Figure 2 (released molecular-data organization).")
    figure2(args.data_dir, args.output_dir, args.assets_dir)
    print(f"generated {args.output_dir / 'figure2_organization.png'}")


if __name__ == "__main__":
    main()
