#!/usr/bin/env python3
"""Generate CoupledMD Supplementary Figure S1 (receptor representation)."""

from __future__ import annotations

from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

from _common import (
    cli,
    FAMILIES,
    FAMILY_LABELS,
    FAMILY_COLORS,
    INK,
    MUTED,
    read,
    save,
)


def supplementary_figure_s1(
    data_dir: Path, output_dir: Path, _assets_dir: Path
) -> None:
    from scipy.cluster.hierarchy import leaves_list, linkage, to_tree

    source = read(data_dir, "Supplementary_Figure_S1_receptor_profile_source.csv")
    source = source.sort_values("receptor_name").reset_index(drop=True)
    assert len(source) == 174
    profiles = {
        row.receptor_name: set(str(row.represented_families).split(";"))
        for row in source.itertuples(index=False)
    }
    names = sorted(profiles)
    values = np.asarray(
        [[family in profiles[name] for family in FAMILIES] for name in names], dtype=float
    )
    linkage_matrix = linkage(values, method="average", metric="jaccard")
    order = leaves_list(linkage_matrix)
    root = to_tree(linkage_matrix)
    angles = np.zeros(len(names))
    angles[order] = np.pi / 2 - 2 * np.pi * np.arange(len(names)) / len(names)
    node_info: dict[int, tuple[float, float, float]] = {}

    def walk(node):
        if node.is_leaf():
            node_info[id(node)] = (angles[node.id], angles[node.id], 1.0)
            return node_info[id(node)]
        left, right = walk(node.left), walk(node.right)
        value = (
            left[0],
            right[1],
            0.19 + 0.77 * (1 - node.dist / (root.dist or 1)),
        )
        node_info[id(node)] = value
        return value

    walk(root)
    fig, ax = plt.subplots(figsize=(13.5, 13.5))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set(xlim=(-1.33, 1.33), ylim=(-1.33, 1.33))

    def arc(radius, start, end):
        theta = np.linspace(start, end, 50)
        ax.plot(
            radius * np.cos(theta),
            radius * np.sin(theta),
            c="#A9AFB7",
            lw=0.55,
            zorder=1,
        )

    def draw(node):
        start, end, radius = node_info[id(node)]
        if node.is_leaf():
            return
        for child in (node.left, node.right):
            child_start, child_end, child_radius = node_info[id(child)]
            middle = (child_start + child_end) / 2
            ax.plot(
                [radius * np.cos(middle), child_radius * np.cos(middle)],
                [radius * np.sin(middle), child_radius * np.sin(middle)],
                c="#A9AFB7",
                lw=0.55,
                zorder=1,
            )
            draw(child)
        arc(radius, start, end)

    draw(root)
    indexed = source.set_index("receptor_name")
    for index, name in enumerate(names):
        angle = angles[index]
        represented = profiles[name]
        x, y = np.cos(angle), np.sin(angle)
        color = (
            "#4A4F55"
            if len(represented) > 1
            else FAMILY_COLORS[next(iter(represented))]
        )
        ax.plot(
            [0.985 * x, 1.035 * x],
            [0.985 * y, 1.035 * y],
            c=color,
            lw=2.0,
            solid_capstyle="butt",
            zorder=3,
        )
        row = indexed.loc[name]
        suffix = " *" if row.uniprot_mapping_status == "unmapped" else ""
        horizontal = "left" if x >= 0 else "right"
        rotation = np.degrees(angle) if x >= 0 else np.degrees(angle) + 180
        ax.text(
            1.058 * x,
            1.058 * y,
            str(row.representative_pdb) + suffix,
            fontsize=8.2,
            rotation=rotation,
            rotation_mode="anchor",
            ha=horizontal,
            va="center",
            color=color,
            fontweight="bold" if len(represented) > 1 or suffix else "normal",
        )
    ax.add_patch(plt.Circle((0, 0), 0.32, fc="white", ec="#D9DDE2", lw=0.8, zorder=2))
    ax.text(0, 0.055, "174", ha="center", va="center", fontsize=24, fontweight="bold", color=INK, zorder=4)
    ax.text(0, -0.025, "distinct receptor names", ha="center", va="center", fontsize=10.2, color=INK, zorder=4)
    ax.text(0, -0.09, "presence/absence profiles across\nfour G-protein families", ha="center", va="center", fontsize=8.4, color=MUTED, zorder=4)
    handles = [
        Line2D([0], [0], color=FAMILY_COLORS[family], lw=3, label=FAMILY_LABELS[family])
        for family in FAMILIES
    ] + [Line2D([0], [0], color="#4A4F55", lw=3, label=">1 family")]
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.01),
        ncol=5,
        frameon=False,
        fontsize=8.5,
        handlelength=1.7,
        columnspacing=1.5,
    )
    ax.text(-1.28, 1.27, "S1", fontweight="bold", fontsize=15.5, ha="left", va="top")
    ax.text(0, -1.22, "Leaf label = representative PDB identifier. * Gs_8HTI has no mapped UniProt accession.", ha="center", fontsize=8.2, color=INK)
    ax.text(0, -1.27, "Branch proximity reflects release-coverage profiles only; it does not represent sequence, structural, or mechanistic similarity.", ha="center", fontsize=8.2, color=MUTED)
    save(fig, output_dir, "supplementary_figure_s1_receptor_profiles")


def main() -> None:
    args = cli("Generate CoupledMD Supplementary Figure S1 (receptor representation).")
    supplementary_figure_s1(args.data_dir, args.output_dir, args.assets_dir)
    print(f"generated {args.output_dir / 'supplementary_figure_s1_receptor_profiles.png'}")


if __name__ == "__main__":
    main()
