#!/usr/bin/env python3
"""Generate CoupledMD Figure 4 (transmembrane-gateway opening)."""

from __future__ import annotations

from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import numpy as np

from _common import (
    cli,
    FAMILIES,
    FAMILY_LABELS,
    FAMILY_COLORS,
    INK,
    GREY,
    PALE,
    family_sort,
    panel_label,
    read,
    save,
)


def figure4(data_dir: Path, output_dir: Path, _assets_dir: Path) -> None:
    source = read(data_dir, "Figure_4_gateway_open_fraction_source.csv")
    interfaces = [
        "TM1-TM2",
        "TM2-TM3",
        "TM3-TM4",
        "TM4-TM5",
        "TM5-TM6",
        "TM6-TM7",
        "TM7-TM1",
    ]
    metadata = source[["system_id", "g_protein_family"]].drop_duplicates()
    metadata = family_sort(metadata)
    wide = source.pivot(index="system_id", columns="interface", values="mean").loc[
        metadata.system_id, interfaces
    ]
    matrix = wide.to_numpy(dtype=float)
    assert matrix.shape == (207, 7)

    fig, ax = plt.subplots(
        1,
        2,
        figsize=(7.1, 4.25),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [0.88, 1.22]},
    )
    cmap = plt.cm.viridis.copy()
    cmap.set_bad(GREY)
    image = ax[0].imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap=cmap, interpolation="nearest")
    ax[0].set(xticks=range(7), xticklabels=interfaces, yticks=[])
    ax[0].tick_params(axis="x", rotation=38)
    counts = [(metadata.g_protein_family == family).sum() for family in FAMILIES]
    starts = np.cumsum([0] + counts[:-1])
    ends = np.cumsum(counts)
    side = ax[0].inset_axes([-0.055, 0, 0.025, 1], transform=ax[0].transAxes)
    codes = np.concatenate([np.full(count, index) for index, count in enumerate(counts)])[:, None]
    side.imshow(
        codes,
        aspect="auto",
        cmap=ListedColormap([FAMILY_COLORS[family] for family in FAMILIES]),
        vmin=0,
        vmax=3,
        interpolation="nearest",
    )
    side.set(xticks=[], yticks=[])
    for start, end, family in zip(starts, ends, FAMILIES):
        if start:
            ax[0].axhline(start - 0.5, c="white", lw=0.8)
        ax[0].text(
            -0.085,
            1 - (start + end) / (2 * len(metadata)),
            f"{FAMILY_LABELS[family]}  {end - start}",
            transform=ax[0].transAxes,
            ha="right",
            va="center",
            fontsize=7.1,
            color=FAMILY_COLORS[family],
            fontweight="bold",
        )
    colorbar = fig.colorbar(image, ax=ax[0], orientation="vertical", fraction=0.052, pad=0.035)
    colorbar.set_ticks([0, 0.5, 1])
    colorbar.ax.tick_params(labelsize=6.6, pad=2)
    colorbar.ax.set_title("open\nfraction", fontsize=7, pad=5, linespacing=0.9, x=0.66)
    panel_label(ax[0], "A", -0.22, 1.03)
    ax[0].set_title("Per-system gateway records", loc="center", pad=8)

    values = [source.loc[source.interface.eq(interface), "mean"].to_numpy() for interface in interfaces]
    violin = ax[1].violinplot(values, showmeans=False, showmedians=False, showextrema=False)
    for body in violin["bodies"]:
        body.set_facecolor("#8DB6C8")
        body.set_edgecolor("#496A78")
        body.set_alpha(0.65)
        body.set_linewidth(0.5)
    ax[1].boxplot(
        values,
        widths=0.18,
        showfliers=False,
        patch_artist=True,
        medianprops={"color": INK, "lw": 1},
        boxprops={"facecolor": "white", "edgecolor": INK, "lw": 0.6},
        whiskerprops={"color": INK, "lw": 0.6},
        capprops={"color": INK, "lw": 0.6},
    )
    ax[1].set(xticks=range(1, 8), xticklabels=interfaces, ylabel="open fraction", ylim=(0, 1))
    ax[1].tick_params(axis="x", rotation=38)
    ax[1].grid(axis="y", color=PALE, lw=0.5)
    panel_label(ax[1], "B")
    ax[1].set_title("Interface distributions across 207 systems")
    save(fig, output_dir, "figure4_gateways")


def main() -> None:
    args = cli("Generate CoupledMD Figure 4 (transmembrane-gateway opening).")
    figure4(args.data_dir, args.output_dir, args.assets_dir)
    print(f"generated {args.output_dir / 'figure4_gateways.png'}")


if __name__ == "__main__":
    main()
