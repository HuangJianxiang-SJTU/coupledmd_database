#!/usr/bin/env python3
"""Generate CoupledMD Figure 5 (structural stability and interface retention)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from _common import (
    cli,
    FAMILIES,
    FAMILY_LABELS,
    FAMILY_COLORS,
    INK,
    PALE,
    panel_label,
    read,
    save,
)


def figure5(data_dir: Path, output_dir: Path, _assets_dir: Path) -> None:
    source = read(data_dir, "Figure_5_panel_source_data.csv")
    metrics = [
        "tm_core_rmsd_A_p95",
        "galpha_interface_rmsd_A_p95",
        "contact_retention_p05",
    ]
    assert len(source) == 618 and source.system_id.nunique() == 206
    fig, ax = plt.subplots(1, 3, figsize=(7.1, 3.2), constrained_layout=True)
    titles = [
        "Receptor TM-core stability",
        r"G$\alpha$ interface-region displacement",
        "Initial interface-contact retention",
    ]
    labels = [
        r"TM-core C$\alpha$ RMSD, P95 ($\AA$)",
        r"G$\alpha$ interface C$\alpha$ RMSD, P95 ($\AA$)",
        "initial contacts retained, P05",
    ]
    upper = [
        max(4, np.ceil(source[metrics[0]].max() + 0.5)),
        max(6, np.ceil(source[metrics[1]].max() + 1)),
        1.02,
    ]
    limits = [(0, upper[0]), (0, upper[1]), (-0.02, upper[2])]
    for index, (axis, metric, title, label, limit) in enumerate(
        zip(ax, metrics, titles, labels, limits)
    ):
        values = [
            source.loc[source.g_protein_family.eq(family), metric].to_numpy()
            for family in FAMILIES
        ]
        violin = axis.violinplot(
            values,
            positions=np.arange(4),
            widths=0.72,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for body, family in zip(violin["bodies"], FAMILIES):
            body.set_facecolor(FAMILY_COLORS[family])
            body.set_edgecolor(FAMILY_COLORS[family])
            body.set_alpha(0.18)
            body.set_linewidth(0.7)
        box = axis.boxplot(
            values,
            positions=np.arange(4),
            widths=0.30,
            showfliers=False,
            patch_artist=True,
            medianprops={"color": INK, "lw": 1.05},
            boxprops={"facecolor": "none", "lw": 0.85},
            whiskerprops={"color": INK, "lw": 0.65},
            capprops={"color": INK, "lw": 0.65},
        )
        for patch, family in zip(box["boxes"], FAMILIES):
            patch.set_facecolor("none")
            patch.set_edgecolor(FAMILY_COLORS[family])
        axis.set(
            xticks=np.arange(4),
            xticklabels=[FAMILY_LABELS[family] for family in FAMILIES],
            ylabel=label,
            ylim=limit,
        )
        axis.grid(axis="y", color=PALE, lw=0.5)
        axis.set_title(title)
        panel_label(axis, chr(65 + index), -0.17, 1.04)
    ax[2].set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax[2].set_yticklabels(["0", ".25", ".50", ".75", "1.00"])
    save(fig, output_dir, "figure5_technical_validation")


def main() -> None:
    args = cli("Generate CoupledMD Figure 5 (structural stability and interface retention).")
    figure5(args.data_dir, args.output_dir, args.assets_dir)
    print(f"generated {args.output_dir / 'figure5_technical_validation.png'}")


if __name__ == "__main__":
    main()
