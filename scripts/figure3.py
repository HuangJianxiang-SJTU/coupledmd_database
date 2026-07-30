#!/usr/bin/env python3
"""Figure 3 — The recurrent core is null-calibrated against prespecified sites but not sequence-conserved."""

from figure_common import *
from scipy.stats import gaussian_kde

# Must be defined before use in the figure3 plotting code
def ensure_union_null_draws() -> Path:
    """Export the deterministic union-null draws omitted by the producer."""
    derived_dir = OUT / "derived_source_data"
    path = derived_dir / "classA_stable_core_union_matched_null.csv"
    if path.exists():
        return path
    position = read_csv(CONSERVATION / "source_data" / "classA_generic_position_conservation.csv")
    clusters2 = read_csv(DOWN / "replica_classA_cluster_effective_recurrence.csv")
    summary2 = read_json(DOWN / "replica_pocket_downstream_summary.json")
    stable_ids = summary2["classA_stable_recurrent_clusters"]
    stable2 = clusters2.loc[clusters2["cluster"].isin(stable_ids)].copy()
    conservation_summary2 = read_json(CONSERVATION / "cavity_conservation_summary.json")
    n_null = int(conservation_summary2["matched_null_draws"])
    min_coverage = int(conservation_summary2["minimum_h40_coverage"])
    n_candidates = int(conservation_summary2["matched_candidates_per_position"])
    lookup = position.set_index("generic_number")

    def pool_for(generic: str) -> np.ndarray:
        target_row = lookup.loc[generic]
        pool = position.loc[
            position["topology_category"].eq(target_row["topology_category"])
            & position["n_h40_groups"].ge(min_coverage)
        ].copy()
        pool["coverage_difference"] = (pool["n_h40_groups"] - target_row["n_h40_groups"]).abs()
        pool = pool.sort_values(["coverage_difference", "generic_number"]).head(min(n_candidates, len(pool)))
        return pool["h40_equal_weight_conservation"].to_numpy(dtype=float)

    rng = np.random.default_rng(20260726 + 811)
    union_set: set[str] = set()
    for cluster in stable2.itertuples(index=False):
        core = list(dict.fromkeys(str(cluster.core_generic_numbers).split(";")))
        union_set.update(core)
        for generic in core:
            pool = pool_for(generic)
            rng.integers(0, len(pool), size=n_null)
    union_sorted = sorted(union_set)
    sampled = []
    for generic in union_sorted:
        pool = pool_for(generic)
        sampled.append(pool[rng.integers(0, len(pool), size=n_null)])
    null_draws = np.vstack(sampled).mean(axis=0)
    derived_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"draw": np.arange(1, n_null + 1), "mean_h40_conservation": null_draws, "n_union_positions": len(union_sorted)}).to_csv(path, index=False)
    return path


values: dict[str, Any] = {"figure": 3, "panels": {}}
summary = read_json(DOWN / "replica_pocket_downstream_summary.json")
recurrence = read_csv(DOWN / "replica_cluster_effective_recurrence.csv")
stable = recurrence.loc[recurrence["cluster"].isin(summary["stable_recurrent_clusters"])].copy()
stable = stable.sort_values("n_receptors", ascending=True)

fig = plt.figure(figsize=(FIG_WIDTH_IN, 8.05))
gs = fig.add_gridspec(3, 2, height_ratios=[0.88, 1.15, 1.0], width_ratios=[0.86, 1.14], hspace=0.47, wspace=0.62)

ax_a = fig.add_subplot(gs[0, 0])
panel_label(ax_a, "A")
panel_title(ax_a, "Replica-stable residue cores")
placeholder(ax_a, "Gi_7E33 — 5-HT1E–Gi,\nR0/R1/R2 residue cores", 1.0)
values["panels"]["A"] = {"placeholder_note": "Gi_7E33 — 5-HT1E–Gi, R0/R1/R2 residue cores", "required_aspect_width_over_height": 1.0}

ax = fig.add_subplot(gs[0, 1])
panel_label(ax, "B")
panel_title(ax, "Stable-core recurrence")
y = np.arange(len(stable))
ax.barh(y - 0.17, stable["n_receptors"], height=0.31, color=COLORS["Receptor identity"], label="Receptors")
ax.barh(y + 0.17, stable["n_H40_homology_groups"], height=0.31, color="#8DBF9A", label="H40 groups")
import re
rule = summary["stable_recurrent_cluster_definition"]
rec_match = re.search(r"n_receptors >= ([0-9]+)", rule)
rec_match2 = re.search(r"probability >= ([0-9.]+)", rule)
j_match = re.search(r"Jaccard >= ([0-9.]+)", rule)
min_receptors = int(rec_match.group(1)) if rec_match else 5
min_recovery = float(rec_match2.group(1)) if rec_match2 else 0.80
min_j = float(j_match.group(1)) if j_match else 0.50
ax.axvline(min_receptors, color=COLORS["Significant"], linestyle="--", linewidth=0.7)
ax.text(min_receptors + 0.8, len(stable) - 0.3, f"≥{min_receptors}-receptor threshold", fontsize=5, color=COLORS["Significant"])
for yi, row in enumerate(stable.itertuples(index=False)):
    if row.cluster in ("R0", "R1", "R2", "R6"):
        ax.text(row.n_receptors + 1, yi - 0.17, f"{row.n_systems} sys / {row.n_receptors} rec / {row.n_H40_homology_groups} H40", fontsize=5, va="center")
ax.set_yticks(y, stable["cluster"])
ax.set_xlabel("Effective recurrence")
ax.legend(frameon=False, loc="lower right")
ax.text(0.99, 0.01, f"stable: recovery ≥{min_recovery:.2f} at core J≥{min_j:.2f}", transform=ax.transAxes, ha="right", fontsize=5)
values["panels"]["B"] = {
    "clusters": stable.to_dict(orient="records"),
    "minimum_receptors": min_receptors, "minimum_split_half_recovery": min_recovery, "minimum_core_jaccard": min_j,
}

ax = fig.add_subplot(gs[1, 0])
panel_label(ax, "C")
panel_title(ax, "Cluster × generic segment")
segment = read_csv(SRC / "figure2_stable_core_segment_map_replica_final207.csv")
row_order = stable.sort_values("n_receptors", ascending=False)["cluster"].tolist()
col_order = ["TM1", "TM2", "TM3", "TM4", "TM5", "TM6", "TM7", "Loops", "H8"]
matrix = segment.pivot(index="cluster", columns="segment", values="core_positions").reindex(index=row_order, columns=col_order).fillna(0)
im = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="Blues", vmin=0)
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        value = int(matrix.iloc[i, j])
        if value:
            ax.text(j, i, str(value), ha="center", va="center", fontsize=5, color="white" if value >= matrix.to_numpy().max() * 0.55 else "#222222")
ax.set_xticks(np.arange(len(col_order)), col_order, rotation=45, ha="right")
ax.set_yticks(np.arange(len(row_order)), row_order)
ax.set_xlabel("GPCRdb generic segment")
ax.set_ylabel("Replica-stable cluster")
cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cbar.set_label("Core positions", fontsize=6)
cbar.ax.tick_params(labelsize=5)
values["panels"]["C"] = {"segments": col_order, "clusters": row_order, "core_position_counts": matrix.to_dict(orient="index")}

ax = fig.add_subplot(gs[1, 1])
panel_label(ax, "D")
panel_title(ax, "Null-calibrated site concordance")
site = read_csv(DOWN / "replica_known_site_null_calibrated_concordance.csv")
class_a = site.loc[site["atlas_scope"].eq("Class A only")].copy()
plot_aliases = ["NaP", "IC56", "TM35mid", "IC345", "UP124", "OBP", "IC23H8"]
plot_rows = class_a.set_index("alias").loc[plot_aliases].reset_index()
display_labels = {
    "NaP": "NaP · A44 †", "IC56": "IC56 · A1", "TM35mid": "TM345/TM35 · A3",
    "IC345": "IC345 · A2", "UP124": "UP124 · A8", "OBP": "OBP druggable · A0",
    "IC23H8": "IC23H8 · A42 ‡",
}
y = np.arange(len(plot_rows))[::-1]
for yi, row in zip(y, plot_rows.itertuples(index=False)):
    ax.plot([row.null_95th_percentile_max_jaccard, row.jaccard], [yi, yi], color="#A0A0A0", linewidth=1.0)
    ax.scatter(row.null_95th_percentile_max_jaccard, yi, s=20, facecolor="none", edgecolor="#555555", linewidth=0.8)
    significant = row.q_bh_prespecified_nonredundant_families < 0.05
    size = 20 + 34 * float(row.reference_recall)
    ax.scatter(row.jaccard, yi, s=size, facecolor=COLORS["Significant"] if significant else COLORS["Non-significant"],
               edgecolor="white", linewidth=0.5, zorder=3)
    result_label = (f"recall {row.reference_recall:.2f}; q={row.q_bh_prespecified_nonredundant_families:.4f}"
                    if significant else f"p={row.empirical_p_value:.3f}; n.s.; {int(row.cluster_n_receptors)} receptors")
    ax.text(max(row.jaccard, row.null_95th_percentile_max_jaccard) + 0.025, yi, result_label, va="center", fontsize=5,
            color="#555555" if not significant else "#222222")
ax.set_yticks(y, [display_labels[a] for a in plot_aliases])
ax.set_xlim(0, 0.96)
ax.set_xlabel("Jaccard concordance")
ax.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=COLORS["Significant"],markeredgecolor="white",label="Observed, q<0.05"),
                   Line2D([0],[0],marker="o",color="none",markerfacecolor=COLORS["Non-significant"],markeredgecolor="white",label="Observed, n.s."),
                   Line2D([0],[0],marker="o",color="none",markerfacecolor="none",markeredgecolor="#555555",label="Matched-null 95th")],
          frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.02), handletextpad=0.35)
ic23 = class_a.loc[class_a["alias"].eq("IC23H8")].iloc[0]
orth = class_a.loc[class_a["alias"].eq("OBP-orth-track")].iloc[0]
ax.text(0.01, -0.19, "†/‡ two-receptor clusters; both are flagged as underpowered.\nOBP orthosteric-track concordance is not estimable after the combined repair.",
        transform=ax.transAxes, fontsize=5, va="top")
values["panels"]["D"] = {
    "primary_scope": "Class A only", "null_iterations_per_definition": int(plot_rows["null_iterations"].iloc[0]),
    "plotted_definitions": plot_rows.to_dict(orient="records"),
    "underpowered_ic23h8": ic23.to_dict(), "orthosteric_track_obp": orth.to_dict(),
}

sub = gs[2, :].subgridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.34)
ax_left = fig.add_subplot(sub[0, 0])
panel_label(ax_left, "E")
panel_title(ax_left, "Conservation dissociation — stable-core union")
conservation_summary = read_json(CONSERVATION / "cavity_conservation_summary.json")
union = conservation_summary["union"]
null_path = ensure_union_null_draws()
null = read_csv(null_path)["mean_h40_conservation"].to_numpy()
density = gaussian_kde(null)
xgrid = np.linspace(null.min(), null.max(), 250)
ygrid = density(xgrid)
ax_left.fill_between(xgrid, 0, ygrid, color="#D9D9D9", edgecolor="#777777", linewidth=0.6)
ax_left.axvspan(union["matched_null_ci95_lo"], union["matched_null_ci95_hi"], color="#BBBBBB", alpha=0.28)
ax_left.axvline(union["observed_mean_h40_conservation"], color=COLORS["Receptor identity"], linewidth=1.6)
ax_left.text(union["observed_mean_h40_conservation"], ygrid.max() * 0.92,
             f" observed {union['observed_mean_h40_conservation']:.3f}", rotation=90, va="top", ha="right",
             fontsize=5.3, color=COLORS["Receptor identity"])
ax_left.text(0.02, 0.96, f"null mean {union['matched_null_mean']:.3f}; 95% {union['matched_null_ci95_lo']:.3f}–"
            f"{union['matched_null_ci95_hi']:.3f}\ntwo-sided p={union['empirical_two_sided_p_value']:.3f}",
            transform=ax_left.transAxes, va="top", fontsize=5.5)
ax_left.set_xlabel("Mean H40-weighted conservation")
ax_left.set_ylabel("Matched-null density")
ax_left.set_yticks([])

ax_right = fig.add_subplot(sub[0, 1])
panel_title(ax_right, "Recurrence versus mean conservation")
cons_table = read_csv(CONSERVATION / "source_data" / "cavity_classA_conservation_matched_null.csv")
for row in cons_table.itertuples(index=False):
    significant = row.q_bh_classA_stable_clusters_two_sided < 0.05
    color = COLORS["Significant"] if significant else COLORS["Non-significant"]
    ax_right.scatter(row.n_receptors_with_cluster, row.mean_h40_conservation, s=28 if significant else 18, color=color, edgecolor="white", linewidth=0.4)
    if row.cluster in ("A0", "A6", "A11"):
        ax_right.annotate(row.cluster, (row.n_receptors_with_cluster, row.mean_h40_conservation), xytext=(3,3), textcoords="offset points", fontsize=5.3, fontweight="bold")
rho = conservation_summary["recurrence_vs_mean_conservation_spearman_rho"]
p = conservation_summary["recurrence_vs_mean_conservation_naive_p"]
ax_right.text(0.98, 0.96, f"Spearman ρ={rho:.3f}, p={p:.3f}", transform=ax_right.transAxes, ha="right", va="top", fontsize=5.5)
ax_right.set_xlabel("Receptors with cluster")
ax_right.set_ylabel("Mean H40-weighted conservation")
footer = (f"{union['unique_classA_stable_core_positions']} unique positions; {conservation_summary['classA_receptor_accessions']} accessions; "
          f"{conservation_summary['classA_h40_groups']} H40 groups; positions with <{conservation_summary['minimum_h40_coverage']} H40 groups excluded")
fig.text(0.5, 0.004, footer, ha="center", va="bottom", fontsize=5)
values["panels"]["E"] = {
    "union": union, "union_null_draws": null.tolist(), "display_density_x": xgrid.tolist(), "display_density_y": ygrid.tolist(),
    "cluster_conservation": cons_table.to_dict(orient="records"),
    "recurrence_spearman_rho": rho, "recurrence_spearman_p": p,
    "classA_receptor_accessions": conservation_summary["classA_receptor_accessions"],
    "classA_h40_groups": conservation_summary["classA_h40_groups"],
    "minimum_h40_coverage": conservation_summary["minimum_h40_coverage"],
}
values["panels"]["A"]["reserved_box"] = report_placeholder(fig, ax_a, "3A")
save_figure(fig, 3, values)
print("Figure 3 done.")
