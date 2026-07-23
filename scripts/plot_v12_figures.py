#!/usr/bin/env python3
"""Generate the v12 Scientific Data figures from path-neutral source tables."""
from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "v12_figure_source_data"
OUT = HERE / "publication_figures_v12"

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


def figure1() -> None:
    composition = read("Figure_1A_cohort_composition.csv")
    boundary = read("Figure_1B_release_boundary.csv")
    protocols = read("Figure_1C_protocol_groups.csv")
    identifiers = read("Figure_1D_identifier_coverage.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.5), constrained_layout=True)
    ax = axes.ravel()

    matrix = (
        composition.pivot(
            index="gpcr_class",
            columns="g_protein_family",
            values="systems",
        )
        .reindex(index=["A", "B"], columns=FAMILIES)
        .fillna(0)
    )
    rgba = np.ones((2, 4, 4))
    maximum = matrix.to_numpy().max()
    for row in range(2):
        for column, family in enumerate(FAMILIES):
            strength = 0.18 + 0.72 * matrix.iloc[row, column] / maximum
            rgba[row, column, :3] = 1 - strength * (
                1 - np.array(to_rgb(COLORS[family]))
            )
    ax[0].imshow(rgba, aspect="auto")
    ax[0].set_xticks(range(4), [FAMILY_LABELS[x] for x in FAMILIES])
    ax[0].set_yticks(range(2), ["Class A", "Class B"])
    for row in range(2):
        for column in range(4):
            ax[0].text(
                column,
                row,
                str(int(matrix.iloc[row, column])),
                ha="center",
                va="center",
                fontweight="bold",
                fontsize=9,
            )
    ax[0].set_title("Frozen cohort composition")
    panel(ax[0], "A")

    order = ["included", "unresolved", "excluded"]
    boundary = boundary.set_index("release_status").loc[order].reset_index()
    colors = [TEAL, "#AAB1B7", RED]
    bars = ax[1].barh(
        np.arange(3),
        boundary.systems,
        color=colors,
        height=0.62,
    )
    ax[1].set_yticks(np.arange(3), ["Included", "Unresolved", "Excluded"])
    ax[1].invert_yaxis()
    ax[1].set_xlabel("working-inventory records")
    ax[1].set_xlim(0, 230)
    ax[1].grid(axis="x", color=PALE, linewidth=0.6)
    for bar, value in zip(bars, boundary.systems):
        ax[1].text(
            value + 4,
            bar.get_y() + bar.get_height() / 2,
            str(int(value)),
            va="center",
            fontweight="bold",
        )
    ax[1].text(
        0.99,
        0.02,
        "207 + 13 + 2 = 222",
        transform=ax[1].transAxes,
        ha="right",
        color=MUTED,
        fontsize=7,
    )
    ax[1].set_title("Release boundary")
    panel(ax[1], "B")

    protocols = protocols.sort_values("protocol_group")
    bars = ax[2].bar(
        protocols.protocol_group,
        protocols.systems,
        color=[BLUE, TEAL, "#7C6C9C"],
        width=0.62,
    )
    ax[2].set_ylabel("systems")
    ax[2].set_ylim(0, 198)
    ax[2].grid(axis="y", color=PALE, linewidth=0.6)
    for bar, systems, replicas in zip(
        bars, protocols.systems, protocols.replicas
    ):
        ax[2].text(
            bar.get_x() + bar.get_width() / 2,
            systems + 5,
            f"{int(systems)} systems\n{int(replicas)} repeats",
            ha="center",
            va="bottom",
            fontsize=6.7,
        )
    ax[2].set_title("Source-simulation protocol groups")
    panel(ax[2], "C")

    ax[3].axis("off")
    lookup = identifiers.set_index("metric").value.to_dict()
    cards = [
        ("systems", "systems"),
        ("production_replicas", "computational\nrepeats"),
        ("receptor_names", "receptor\nnames"),
        ("mapped_uniprot_accessions", "mapped UniProt\naccessions"),
    ]
    for index, (key, label) in enumerate(cards):
        x = 0.04 + (index % 2) * 0.49
        y = 0.54 - (index // 2) * 0.45
        box(
            ax[3],
            (x, y),
            0.43,
            0.34,
            str(int(lookup[key])),
            label,
            face="#F2F5F6",
            edge=TEAL,
            title_size=16,
        )
    ax[3].text(
        0.5,
        0.03,
        "Gs_8HTI retains an explicit null accession",
        transform=ax[3].transAxes,
        ha="center",
        color=MUTED,
        fontsize=6.8,
    )
    ax[3].set_title("Identifier and repeat coverage")
    panel(ax[3], "D", x=-0.04)
    save(fig, "figure1_cohort_scope")


def figure2() -> None:
    roles = read("Figure_2_record_roles.csv").set_index("record")
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.6), constrained_layout=True)
    ax = axes.ravel()
    for item in ax:
        item.axis("off")
        item.set_xlim(0, 1)
        item.set_ylim(0, 1)

    box(
        ax[0],
        (0.08, 0.73),
        0.84,
        0.16,
        "FULL-SYSTEM SOURCE SIMULATIONS",
        "protein complex + membrane + solvent + ions",
        face="#F1F2F3",
        edge=MUTED,
    )
    box(
        ax[0],
        (0.25, 0.43),
        0.50,
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
        "protein-complex reuse; not a full-system archive",
        face="#EAF0F7",
        edge=BLUE,
    )
    arrow(ax[0], (0.5, 0.72), (0.5, 0.59))
    arrow(ax[0], (0.5, 0.42), (0.5, 0.28))
    ax[0].text(
        0.96,
        0.49,
        "membrane · solvent · ions\nexcluded",
        ha="right",
        va="center",
        fontsize=6.7,
        color=RED,
    )
    ax[0].set_title("Molecular-record boundary")
    panel(ax[0], "A", x=-0.04)

    rows = [
        ("Reduced trajectory", ["yes", "no", "no", "no"]),
        ("Gateway summary", ["summary", "no", "no", "yes"]),
        ("Full-system source", ["yes", "yes", "yes", "yes"]),
    ]
    columns = ["protein", "lipid", "solvent\nions", "gateway\nrecalc."]
    x_positions = [0.42, 0.58, 0.74, 0.90]
    for x, label in zip(x_positions, columns):
        ax[1].text(x, 0.84, label, ha="center", va="bottom", fontsize=6.4)
    for row_index, (label, values) in enumerate(rows):
        y = 0.66 - row_index * 0.22
        ax[1].text(0.02, y, label, ha="left", va="center", fontsize=7)
        for x, value in zip(x_positions, values):
            color = GREEN if value in {"yes", "summary"} else "#D8DDE1"
            ax[1].add_patch(
                FancyBboxPatch(
                    (x - 0.055, y - 0.055),
                    0.11,
                    0.11,
                    boxstyle="round,pad=.008",
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.6,
                )
            )
            ax[1].text(
                x,
                y,
                "Y" if value == "yes" else ("S" if value == "summary" else "—"),
                ha="center",
                va="center",
                color="white" if value in {"yes", "summary"} else MUTED,
                fontsize=7.4,
                fontweight="bold",
            )
    ax[1].text(
        0.02,
        0.07,
        "S = processed summary supplied; molecular recalculation still requires the source trajectories",
        fontsize=6.2,
        color=MUTED,
        wrap=True,
    )
    ax[1].set_title("What each record can support")
    panel(ax[1], "B", x=-0.04)

    box(
        ax[2],
        (0.05, 0.38),
        0.24,
        0.25,
        "SYSTEM",
        "207 system IDs\nmetadata join key",
        face="#F5F6F7",
    )
    box(
        ax[2],
        (0.38, 0.38),
        0.24,
        0.25,
        "REPLICA",
        "three computational\nrepeats per system",
        face="#EDF4F4",
        edge=TEAL,
    )
    box(
        ax[2],
        (0.71, 0.38),
        0.24,
        0.25,
        "FILES",
        "PDB + XTC\nSHA-256",
        face="#EAF0F7",
        edge=BLUE,
    )
    arrow(ax[2], (0.30, 0.505), (0.37, 0.505))
    arrow(ax[2], (0.63, 0.505), (0.70, 0.505))
    ax[2].text(0.335, 0.57, "1:3", ha="center", fontsize=6.4, color=MUTED)
    ax[2].text(
        0.665, 0.57, "1:2", ha="center", fontsize=6.4, color=MUTED
    )
    ax[2].text(
        0.5,
        0.20,
        "File-level manifest joins sizes, checksums, atom counts,\nframe counts and component scope.",
        ha="center",
        va="center",
        fontsize=6.6,
        color=MUTED,
    )
    ax[2].set_title("Replica-aware organization")
    panel(ax[2], "C", x=-0.04)

    access = [
        ("Supplementary Data", "inventories · annotations · QC", GREEN),
        ("REST API", "metadata and processed records", TEAL),
        ("GitHub repositories", "code · paper source · schemas", BLUE),
        ("Zenodo records", "reduced PDB/XTC + checksums", AMBER),
    ]
    for index, (title, subtitle, color) in enumerate(access):
        y = 0.76 - index * 0.20
        box(
            ax[3],
            (0.10, y),
            0.80,
            0.13,
            title,
            subtitle,
            face="#F7F8F9",
            edge=color,
        )
    ax[3].text(
        0.5,
        0.06,
        "Public-release claims remain conditional on unauthenticated DOI access.",
        ha="center",
        fontsize=6.4,
        color=MUTED,
    )
    ax[3].set_title("Access and reuse layers")
    panel(ax[3], "D", x=-0.04)
    save(fig, "figure2_records_reuse_boundary")


def figure3() -> None:
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
    for family in FAMILIES:
        selected = eligible.g_family.eq(family)
        good = selected & recovered
        failed = selected & ~recovered
        axes[0].scatter(
            eligible.loc[good, "best_ortho_freq"],
            eligible.index[good],
            s=20,
            color=COLORS[family],
            label=FAMILY_LABELS[family],
            edgecolor="white",
            linewidth=0.25,
        )
        axes[0].scatter(
            np.full(failed.sum(), 0.805),
            eligible.index[failed],
            s=24,
            marker="x",
            color="#9AA1A6",
        )
    axes[0].axvline(0.85, linestyle="--", linewidth=0.8, color=INK)
    axes[0].set_xlim(0.80, 1.005)
    axes[0].set_ylim(-2, 60)
    axes[0].set_yticks([])
    axes[0].set_xlabel("best orthosteric-pocket frequency")
    axes[0].set_ylabel("58 peptide-ligand positive controls")
    axes[0].grid(axis="x", color=PALE, linewidth=0.5)
    axes[0].text(
        0.803,
        58.5,
        "49/58 recovered (84.5%)",
        va="top",
        fontweight="bold",
        fontsize=7.5,
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
    )
    axes[1].set_xticks([0, 1], labels)
    axes[1].set_ylabel("systems")
    axes[1].set_ylim(0, 230)
    axes[1].grid(axis="y", color=PALE, linewidth=0.5)
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
        0.73,
        "2,149 detected-pocket rows",
        transform=axes[1].transAxes,
        ha="left",
        fontsize=7,
        color=MUTED,
    )
    axes[1].text(
        0.56,
        0.64,
        "205 + 2 = 207 complete system records",
        transform=axes[1].transAxes,
        ha="left",
        fontsize=7,
        fontweight="bold",
    )
    axes[1].set_title("Pocket-record completeness")
    panel(axes[1], "B")
    save(fig, "figure3_annotation_validation")


def violin_by_replica(
    ax: plt.Axes,
    data: pd.DataFrame,
    field: str,
    ylabel: str,
    title: str,
) -> None:
    values = [
        data.loc[data.replica.eq(replica), field].dropna().to_numpy()
        for replica in [1, 2, 3]
    ]
    plot = ax.violinplot(
        values,
        positions=[1, 2, 3],
        widths=0.72,
        showextrema=False,
        showmedians=False,
    )
    for body, color in zip(plot["bodies"], ["#82AFC0", "#5B91A7", "#326F89"]):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.34)
    boxes = ax.boxplot(
        values,
        positions=[1, 2, 3],
        widths=0.25,
        showfliers=False,
        patch_artist=True,
        medianprops={"color": INK, "linewidth": 1.0},
        boxprops={"facecolor": "white", "edgecolor": TEAL, "linewidth": 0.8},
        whiskerprops={"color": INK, "linewidth": 0.6},
        capprops={"color": INK, "linewidth": 0.6},
    )
    for item in boxes["boxes"]:
        item.set_alpha(0.9)
    ax.set_xticks([1, 2, 3], ["Rep 1", "Rep 2", "Rep 3"])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", color=PALE, linewidth=0.5)


def figure4() -> None:
    validation = read("Figure_4A_replica_validation.csv")
    reduced = read("Figure_4B_reduced_record_qc.csv")
    available = validation[validation.validation_status.eq("available")].copy()
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.4), constrained_layout=True)
    ax = axes.ravel()
    violin_by_replica(
        ax[0],
        available,
        "tm_core_rmsd_A_p95",
        "TM-core Cα RMSD P95 (Å)",
        "Receptor TM-core displacement",
    )
    panel(ax[0], "A")
    violin_by_replica(
        ax[1],
        available,
        "galpha_interface_rmsd_A_p95",
        "Gα-interface Cα RMSD P95 (Å)",
        "Gα interface-region displacement",
    )
    panel(ax[1], "B")
    violin_by_replica(
        ax[2],
        available,
        "contact_retention_p05",
        "initial contacts retained P05",
        "Initial interface-contact retention",
    )
    ax[2].set_ylim(-0.03, 1.03)
    panel(ax[2], "C")

    ax[3].axis("off")
    ax[3].set_xlim(0, 1)
    ax[3].set_ylim(0, 1)
    checks = [
        ("Structural metric records", 618, 621, "206 systems; Gs_8HTI unavailable"),
        ("Reduced-record audit", int(reduced.status.eq("OK").sum()), 207, "exact final-207 cohort"),
        (
            "PDB/XTC atom-count match",
            int((~reduced.atom_count_mismatch.astype(bool)).sum()),
            207,
            "40 sampled frames per record",
        ),
    ]
    for index, (label, passed, total, note) in enumerate(checks):
        y = 0.82 - index * 0.25
        ax[3].text(0.02, y + 0.08, label, fontsize=7.3, fontweight="bold")
        ax[3].add_patch(
            Rectangle((0.02, y), 0.78, 0.055, facecolor=PALE, edgecolor="none")
        )
        ax[3].add_patch(
            Rectangle(
                (0.02, y),
                0.78 * passed / total,
                0.055,
                facecolor=GREEN,
                edgecolor="none",
            )
        )
        ax[3].text(
            0.83,
            y + 0.0275,
            f"{passed}/{total}",
            va="center",
            fontweight="bold",
            fontsize=7.2,
        )
        ax[3].text(0.02, y - 0.04, note, fontsize=6.2, color=MUTED)
    distribution = reduced.n_frames.value_counts().sort_index().to_dict()
    text = " · ".join(f"{count}×{frames}" for frames, count in distribution.items())
    ax[3].text(
        0.02,
        0.06,
        f"Frame-count distribution: {text}",
        fontsize=6.5,
        color=MUTED,
    )
    ax[3].set_title("Completeness and reduced-record QC")
    panel(ax[3], "D", x=-0.04)
    save(fig, "figure4_replica_technical_validation")


def figure5() -> None:
    status = read("Figure_5A_release_status.csv").sort_values("replica")
    manifest = read("Figure_5B_release_file_manifest.csv")
    example = read("Figure_5D_worked_example.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.55), constrained_layout=True)
    ax = axes.ravel()

    ax[0].axis("off")
    ax[0].set_xlim(0, 1)
    ax[0].set_ylim(0, 1)
    column_labels = ["Manifest", "Full-frame QC", "Private draft", "Public DOI"]
    xs = [0.38, 0.57, 0.76, 0.92]
    for x, label in zip(xs, column_labels):
        ax[0].text(x, 0.89, label, ha="center", fontsize=6.2)
    for index, row in enumerate(status.itertuples(index=False)):
        y = 0.69 - index * 0.23
        ax[0].text(0.02, y, f"Replica {int(row.replica)}", va="center", fontweight="bold")
        cells = [
            (
                row.local_manifest_rows == 207,
                f"{int(row.local_manifest_rows)}/207",
            ),
            (
                row.local_full_frame_qc_status == "passed",
                str(row.local_full_frame_qc_status),
            ),
            (
                row.remote_state == "done",
                f"{int(row.remote_uploaded_files)} files\n{row.remote_state}",
            ),
            (
                bool(row.publicly_resolvable),
                "resolves" if row.publicly_resolvable else "not public",
            ),
        ]
        for x, (passed, label) in zip(xs, cells):
            face = GREEN if passed else AMBER
            if "not public" in label:
                face = RED
            ax[0].add_patch(
                FancyBboxPatch(
                    (x - 0.068, y - 0.066),
                    0.136,
                    0.132,
                    boxstyle="round,pad=.006",
                    facecolor=face,
                    edgecolor="white",
                    linewidth=0.5,
                )
            )
            ax[0].text(
                x,
                y,
                label,
                ha="center",
                va="center",
                fontsize=5.7,
                color="white",
                linespacing=1.05,
            )
    ax[0].text(
        0.02,
        0.05,
        "Reserved identifiers are not treated as publication.",
        fontsize=6.4,
        color=MUTED,
    )
    ax[0].set_title("Machine-verified release state")
    panel(ax[0], "A", x=-0.04)

    xtc = manifest[manifest.file_role.eq("trajectory_xtc")].copy()
    replicas = sorted(xtc.release_replica.unique())
    if replicas:
        values = [
            xtc.loc[xtc.release_replica.eq(replica), "size_bytes"].dropna()
            / 1e9
            for replica in replicas
        ]
        plot = ax[1].violinplot(
            values,
            positions=replicas,
            widths=0.7,
            showextrema=False,
        )
        for body in plot["bodies"]:
            body.set_facecolor(BLUE)
            body.set_edgecolor(BLUE)
            body.set_alpha(0.3)
        ax[1].boxplot(
            values,
            positions=replicas,
            widths=0.25,
            showfliers=False,
            patch_artist=True,
            medianprops={"color": INK},
            boxprops={"facecolor": "white", "edgecolor": BLUE},
        )
    ax[1].set_xticks([1, 2, 3], ["Rep 1", "Rep 2", "Rep 3"])
    ax[1].set_ylabel("XTC size per system (GB)")
    ax[1].grid(axis="y", color=PALE, linewidth=0.5)
    ax[1].set_title("Reduced-trajectory record sizes")
    panel(ax[1], "B")

    if replicas:
        atoms = [
            xtc.loc[xtc.release_replica.eq(replica), "n_atoms"].dropna()
            for replica in replicas
        ]
        plot = ax[2].violinplot(
            atoms,
            positions=replicas,
            widths=0.7,
            showextrema=False,
        )
        for body in plot["bodies"]:
            body.set_facecolor(TEAL)
            body.set_edgecolor(TEAL)
            body.set_alpha(0.3)
        ax[2].boxplot(
            atoms,
            positions=replicas,
            widths=0.25,
            showfliers=False,
            patch_artist=True,
            medianprops={"color": INK},
            boxprops={"facecolor": "white", "edgecolor": TEAL},
        )
    ax[2].set_xticks([1, 2, 3], ["Rep 1", "Rep 2", "Rep 3"])
    ax[2].set_ylabel("retained atoms")
    ax[2].grid(axis="y", color=PALE, linewidth=0.5)
    ax[2].set_title("Matched PDB/XTC atom counts")
    panel(ax[2], "C")

    ax[3].axis("off")
    ax[3].set_xlim(0, 1)
    ax[3].set_ylim(0, 1)
    pdb = example[example.file_role.eq("structure_pdb")].iloc[0]
    xtc_example = example[example.file_role.eq("trajectory_xtc")].iloc[0]
    code = (
        f"# G12_7SF7, release replica 1\n"
        f"sha256sum {Path(str(pdb.archive_member_path)).name}\n"
        f"# expected: {str(pdb.sha256)[:18]}…\n"
        f"sha256sum {Path(str(xtc_example.archive_member_path)).name}\n"
        f"# expected: {str(xtc_example.sha256)[:18]}…\n\n"
        "import MDAnalysis as mda\n"
        "u = mda.Universe(pdb_path, xtc_path)\n"
        "assert u.atoms.n_atoms == expected_atoms\n"
        "assert len(u.trajectory) == 2500"
    )
    ax[3].add_patch(
        FancyBboxPatch(
            (0.02, 0.12),
            0.96,
            0.72,
            boxstyle="round,pad=.018",
            facecolor="#F4F6F7",
            edgecolor="#9BA3A8",
            linewidth=0.8,
        )
    )
    ax[3].text(
        0.05,
        0.78,
        code,
        family="monospace",
        fontsize=6.25,
        va="top",
        color=INK,
        linespacing=1.28,
    )
    ax[3].text(
        0.5,
        0.04,
        "Verify both digests before analysis; load the PDB as the XTC topology.",
        ha="center",
        fontsize=6.3,
        color=MUTED,
    )
    ax[3].set_title("Checksum-first worked reuse example")
    panel(ax[3], "D", x=-0.04)
    save(fig, "figure5_release_integrity_reuse")


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
    )
    axes[0].set_xticks(
        np.arange(len(campaign)),
        ["Nominal\nprotocol", "Original file-\nobserved"],
    )
    axes[0].set_ylim(310.35, 310.55)
    axes[0].set_ylabel("aggregate sampling (µs)")
    axes[0].grid(axis="y", color=PALE, linewidth=0.5)
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
    )
    axes[1].set_xticks(
        np.arange(3),
        ["Nominal", "Original\nobserved", "Repaired\nreduced source"],
    )
    axes[1].set_ylabel("coordinate span (µs)")
    axes[1].set_ylim(0, 0.56)
    axes[1].grid(axis="y", color=PALE, linewidth=0.5)
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


def main() -> int:
    style()
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    supplementary_figure1()
    supplementary_figure2()
    print(f"Wrote 7 PDF/PNG figure pairs to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
