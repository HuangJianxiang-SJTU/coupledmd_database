#!/usr/bin/env python3
"""Generate CoupledMD Figure 3 (pocket validation and availability)."""

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
    GREY,
    PALE,
    family_sort,
    panel_label,
    read,
    save,
)


def figure3(data_dir: Path, output_dir: Path, _assets_dir: Path) -> None:
    positive = family_sort(read(data_dir, "Figure_3_positive_control_source.csv"))
    availability = family_sort(read(data_dir, "Figure_3_pocket_availability_source.csv"))
    positive["ortho_recovered"] = positive.ortho_recovered.astype(str).str.lower().eq("true")
    assert len(positive) == 58 and int(positive.ortho_recovered.sum()) == 49

    fig, ax = plt.subplots(1, 2, figsize=(7.1, 3.35), constrained_layout=True)
    y_start = 0
    for family in FAMILIES:
        subset = positive[positive.g_protein_family.eq(family)]
        if subset.empty:
            continue
        y = np.arange(y_start, y_start + len(subset))
        recovered = subset.ortho_recovered.to_numpy()
        ax[0].scatter(
            subset.loc[recovered, "best_ortho_freq"],
            y[recovered],
            c=FAMILY_COLORS[family],
            s=18,
            label=FAMILY_LABELS[family],
        )
        ax[0].scatter(
            np.full((~recovered).sum(), 0.805),
            y[~recovered],
            c=GREY,
            s=20,
            marker="x",
            label="undefined / not recovered" if family == "Gi" else None,
        )
        y_start += len(subset)
    ax[0].axvline(0.85, ls="--", c=INK, lw=0.8)
    ax[0].set(
        xlabel="best orthosteric-pocket frequency",
        ylabel="eligible systems (n = 58)",
        xlim=(0.80, 1.01),
        ylim=(-1, 68),
        yticks=[],
    )
    ax[0].axhspan(59, 68, color="white", zorder=2)
    ax[0].axhline(59, color=PALE, lw=0.6, zorder=3)
    ax[0].text(0.805, 66.5, "49/58 (84.5%) recovered", fontsize=7.8, fontweight="bold", va="top", zorder=4)
    ax[0].text(0.805, 62.6, "Gi/o 16/18 - Gs 15/18 - Gq/11 18/22 - G12/13 N/A", fontsize=6.6, va="top", color=MUTED, zorder=4)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            ls="",
            color=FAMILY_COLORS[family],
            label=FAMILY_LABELS[family],
            markersize=4.5,
        )
        for family in FAMILIES
        if (positive.g_protein_family == family).any()
    ]
    handles.append(
        Line2D(
            [0],
            [0],
            marker="x",
            ls="",
            color=GREY,
            label="undefined / not recovered",
            markersize=5,
        )
    )
    ax[0].legend(
        handles=handles,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.155),
        ncol=4,
        columnspacing=0.8,
        handletextpad=0.35,
        borderaxespad=0,
    )
    panel_label(ax[0], "A")
    ax[0].set_title("Orthosteric positive control")

    values = [
        availability.loc[availability.g_protein_family.eq(family), "detected_pockets"]
        for family in FAMILIES
    ]
    for index, (family, family_values) in enumerate(zip(FAMILIES, values), start=1):
        rng = np.random.default_rng(index)
        ax[1].scatter(
            rng.normal(index, 0.065, len(family_values)),
            family_values,
            s=6,
            c=FAMILY_COLORS[family],
            alpha=0.42,
            edgecolors="none",
            zorder=1,
        )
    box = ax[1].boxplot(
        values,
        tick_labels=[FAMILY_LABELS[family] for family in FAMILIES],
        patch_artist=True,
        showfliers=False,
        widths=0.42,
        medianprops={"color": INK, "lw": 1.1},
        boxprops={"facecolor": "none", "lw": 1.05},
        whiskerprops={"color": INK, "lw": 0.75},
        capprops={"color": INK, "lw": 0.75},
    )
    for patch, family in zip(box["boxes"], FAMILIES):
        patch.set_facecolor("none")
        patch.set_edgecolor(FAMILY_COLORS[family])
    ax[1].set(ylabel="pockets per system")
    panel_label(ax[1], "B")
    ax[1].set_title("Pocket-record availability")
    save(fig, output_dir, "figure3_pocket_validation")


def main() -> None:
    args = cli("Generate CoupledMD Figure 3 (pocket validation and availability).")
    figure3(args.data_dir, args.output_dir, args.assets_dir)
    print(f"generated {args.output_dir / 'figure3_pocket_validation.png'}")


if __name__ == "__main__":
    main()
