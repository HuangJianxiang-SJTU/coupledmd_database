#!/usr/bin/env python3
"""Generate CoupledMD Figure 1 (cohort composition)."""

from __future__ import annotations

from matplotlib.colors import to_rgb
from matplotlib.patches import FancyBboxPatch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import (
    cli,
    FAMILIES,
    FAMILY_LABELS,
    FAMILY_COLORS,
    INK,
    MUTED,
    GREY,
    panel_label,
    read,
    save,
)


def figure1(data_dir: Path, output_dir: Path, _assets_dir: Path) -> None:
    systems = read(data_dir, "Figure_1_system_composition_source.csv")
    ligand = read(data_dir, "Figure_1_ligand_context_source.csv")
    assert len(systems) == 207 and len(ligand) == 207

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 6.2), constrained_layout=True)
    ax = axes.ravel()

    matrix = pd.crosstab(systems.gpcr_class, systems.g_protein_family).reindex(
        index=["A", "B"], columns=FAMILIES, fill_value=0
    )
    rgba = np.ones((2, 4, 4))
    peak = matrix.to_numpy().max()
    for column, family in enumerate(FAMILIES):
        for row in range(2):
            strength = 0.22 + 0.70 * matrix.iloc[row, column] / peak
            rgba[row, column, :3] = 1 - strength * (
                1 - np.asarray(to_rgb(FAMILY_COLORS[family]))
            )
    ax[0].imshow(rgba)
    ax[0].set(
        xticks=range(4),
        xticklabels=[FAMILY_LABELS[family] for family in FAMILIES],
        yticks=range(2),
        yticklabels=["Class A", "Class B"],
    )
    for row in range(2):
        for column in range(4):
            ax[0].text(
                column,
                row,
                str(matrix.iloc[row, column]),
                ha="center",
                va="center",
                fontweight="bold",
                color=INK,
            )
    panel_label(ax[0], "A")
    ax[0].set_title("System counts")

    profiles = systems.groupby("receptor_name").g_protein_family.agg(
        lambda values: sorted(set(values))
    )
    multiple = profiles[profiles.map(len) > 1].sort_values(
        key=lambda values: values.map(len), ascending=False
    )
    y = np.arange(len(multiple))
    left = np.zeros(len(multiple))
    for family in FAMILIES:
        values = multiple.map(lambda represented: family in represented).astype(int)
        ax[1].barh(
            y,
            values,
            left=left,
            color=FAMILY_COLORS[family],
            label=FAMILY_LABELS[family],
        )
        left += values
    ax[1].set(
        yticks=y,
        yticklabels=multiple.index,
        xlabel="represented G-protein families",
        xlim=(0, 3.65),
    )
    ax[1].invert_yaxis()
    ax[1].legend(frameon=False, ncol=1, loc="lower right", borderaxespad=0.4)
    panel_label(ax[1], "B", -0.17)
    ax[1].set_title("Receptors represented with >1 family")

    context = pd.crosstab(ligand.g_protein_family, ligand.ligand_context).reindex(
        FAMILIES, fill_value=0
    )
    left = np.zeros(4)
    for key, color, label in [
        ("peptide", "#64B5CD", "peptide eligible"),
        ("none", GREY, "no orthosteric-defining ligand"),
    ]:
        values = context[key].to_numpy()
        ax[2].bar(
            range(4), values, bottom=left, label=label, color=color, edgecolor="white", lw=0.4
        )
        for index, (bottom, count) in enumerate(zip(left, values)):
            if count:
                ax[2].text(
                    index,
                    bottom + count / 2,
                    str(int(count)),
                    ha="center",
                    va="center",
                    fontsize=7.2,
                    color="white" if key == "peptide" else INK,
                )
        left += values
    ax[2].set(
        xticks=range(4),
        xticklabels=[FAMILY_LABELS[family] for family in FAMILIES],
        ylabel="systems",
        ylim=(0, 104),
        xlim=(-0.55, 3.55),
    )
    ax[2].legend(frameon=False, loc="upper right", borderaxespad=0)
    panel_label(ax[2], "C")
    ax[2].set_title("Ligand context used for pocket validation", pad=8)

    ax[3].axis("off")
    for x, value, text in [
        (0.07, 174, "receptor names"),
        (0.55, 173, "mapped\nUniProt"),
    ]:
        ax[3].add_patch(
            FancyBboxPatch(
                (x, 0.47),
                0.38,
                0.39,
                boxstyle="round,pad=.02",
                fc="#F1F3F5",
                ec="#737A82",
                lw=0.8,
            )
        )
        ax[3].text(
            x + 0.19,
            0.69,
            str(value),
            ha="center",
            va="center",
            fontsize=21,
            fontweight="bold",
            color=INK,
        )
        ax[3].text(
            x + 0.19, 0.55, text, ha="center", va="center", fontsize=8.3, color=INK
        )
    ax[3].text(
        0.5,
        0.37,
        "one accession gap",
        ha="center",
        fontsize=7.8,
        fontweight="bold",
        color=MUTED,
    )
    ax[3].text(
        0.5,
        0.25,
        "Gs_8HTI - consensus OR52c model\nno canonical UniProt accession",
        ha="center",
        fontsize=7.2,
        color=MUTED,
    )
    panel_label(ax[3], "D", -0.04, 1.02)
    ax[3].set_title("Receptor identifier coverage")
    fig.canvas.draw()
    fig.set_layout_engine(None)
    right_shift = max(0, ax[1].get_position().y1 - ax[0].get_position().y1)
    for axis in (ax[1], ax[3]):
        position = axis.get_position()
        axis.set_position(
            [position.x0, position.y0 - right_shift, position.width, position.height]
        )
    c_position = ax[2].get_position()
    d_position = ax[3].get_position()
    new_x = c_position.x0 + 0.012
    new_width = min(c_position.width + 0.115, d_position.x0 - new_x - 0.045)
    ax[2].set_position(
        [new_x, c_position.y0, new_width, c_position.height]
    )
    d_position = ax[3].get_position()
    trim = 0.055
    ax[3].set_position(
        [
            d_position.x0,
            d_position.y0 + trim,
            d_position.width,
            d_position.height - trim,
        ]
    )
    save(fig, output_dir, "figure1_composition")


def main() -> None:
    args = cli("Generate CoupledMD Figure 1 (cohort composition).")
    figure1(args.data_dir, args.output_dir, args.assets_dir)
    print(f"generated {args.output_dir / 'figure1_composition.png'}")


if __name__ == "__main__":
    main()
