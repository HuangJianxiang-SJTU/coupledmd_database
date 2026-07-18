#!/usr/bin/env python3
"""Generate CoupledMD Figure 6 (GPCR-centered pocket reuse atlas)."""

from __future__ import annotations

from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

from _common import (
    cli,
    FAMILIES,
    FAMILY_LABELS,
    FAMILY_COLORS,
    POCKET_COLORS,
    PALE,
    panel_label,
    read,
    save,
)


def figure6(data_dir: Path, output_dir: Path, assets_dir: Path) -> None:
    panel_b = read(data_dir, "Figure_6_panel_B_class_counts.csv")
    panel_c = read(data_dir, "Figure_6_panel_C_system_summary.csv")
    panel_d = read(data_dir, "Figure_6_panel_D_position_counts.csv")
    assert panel_b.pocket_records.sum() == 669 and len(panel_c) == 201

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(7.1, 5.75),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1, 1.12]},
    )
    ax = axes.ravel()
    snapshot_path = assets_dir / "figure6_panel_A_structure.png"
    if not snapshot_path.exists():
        snapshot_path = assets_dir / "figure6_1_panelA_structure_clean.png"
    snapshot = plt.imread(snapshot_path)
    rgb = snapshot[..., :3]
    mask = np.any(rgb < 0.985, axis=2)
    rows, columns = np.where(mask)
    snapshot = snapshot[
        max(0, rows.min() - 20) : min(snapshot.shape[0], rows.max() + 21),
        max(0, columns.min() - 20) : min(snapshot.shape[1], columns.max() + 21),
    ]
    ax[0].axis("off")
    panel_label(ax[0], "A", -0.22, 1.02)
    ax[0].set_title("Representative GPCR pocket classes on Gi_6K42", x=0.58)
    snapshot_ax = ax[0].inset_axes([-0.08, -0.02, 0.72, 1.05])
    snapshot_ax.imshow(snapshot)
    snapshot_ax.axis("off")
    legend_entries = [
        ("orthosteric", "C1 - orthosteric"),
        ("extracellular_vestibule", "C26 - extracellular vestibule"),
        ("tm_core_allosteric", "C22, C68 - TM-core allosteric"),
        ("intracellular_allosteric", "C11, C2 - intracellular allosteric"),
    ]
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            ls="none",
            ms=5,
            mfc=POCKET_COLORS[key],
            mec="none",
        )
        for key, _ in legend_entries
    ]
    ax[0].legend(
        handles,
        [label for _, label in legend_entries],
        frameon=False,
        ncol=1,
        loc="center left",
        bbox_to_anchor=(0.58, 0.46),
        labelspacing=0.85,
        handletextpad=0.35,
    )

    classes = panel_b.anatomical_class.tolist()
    y = np.arange(len(panel_b))
    bars = ax[1].barh(
        y,
        panel_b.pocket_records,
        color=[POCKET_COLORS[key] for key in classes],
        height=0.62,
    )
    ax[1].set(
        yticks=y,
        yticklabels=panel_b.panel_label,
        xlabel="GPCR-centered pocket records",
    )
    ax[1].invert_yaxis()
    ax[1].grid(axis="x", color=PALE, lw=0.5)
    ax[1].set_xlim(0, panel_b.pocket_records.max() * 1.22)
    for bar, count in zip(bars, panel_b.pocket_records):
        ax[1].text(
            count + 5,
            bar.get_y() + bar.get_height() / 2,
            f"{count:,}",
            va="center",
            fontsize=7,
            fontweight="bold",
        )
    panel_label(ax[1], "B", -0.14, 1.02)
    ax[1].set_title("GPCR-centered records by anatomical class")

    for family in FAMILIES:
        subset = panel_c[panel_c.g_protein_family.eq(family)]
        ax[2].scatter(
            subset.gpcr_centered_pockets,
            subset.mean_gpcr_pocket_occupancy,
            c=FAMILY_COLORS[family],
            s=20,
            label=FAMILY_LABELS[family],
            alpha=0.72,
            edgecolor="white",
            linewidth=0.25,
        )
    ax[2].set(
        xlabel="GPCR-centered pockets per system",
        ylabel="mean GPCR-pocket occupancy",
    )
    ax[2].grid(color=PALE, lw=0.5)
    ax[2].legend(frameon=False, ncol=2, loc="upper right")
    panel_label(ax[2], "C")
    ax[2].set_title("Per-system GPCR-pocket landscape")

    plotted = panel_d.sort_values("rank", ascending=False)
    bars = ax[3].barh(
        plotted.gpcrdb_generic_position,
        plotted.orthosteric_pocket_records,
        color="#4C72B0",
        height=0.68,
    )
    ax[3].set(xlabel="orthosteric GPCR-pocket records")
    ax[3].grid(axis="x", color=PALE, lw=0.5)
    for bar, count in zip(bars, plotted.orthosteric_pocket_records):
        ax[3].text(
            count + 0.5,
            bar.get_y() + bar.get_height() / 2,
            str(count),
            va="center",
            fontsize=6.6,
        )
    panel_label(ax[3], "D")
    ax[3].set_title("Recurrent orthosteric GPCRdb positions")
    save(fig, output_dir, "figure6_reuse_atlas")


def main() -> None:
    args = cli("Generate CoupledMD Figure 6 (GPCR-centered pocket reuse atlas).")
    figure6(args.data_dir, args.output_dir, args.assets_dir)
    print(f"generated {args.output_dir / 'figure6_reuse_atlas.png'}")


if __name__ == "__main__":
    main()
