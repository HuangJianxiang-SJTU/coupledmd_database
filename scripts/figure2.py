#!/usr/bin/env python3
"""Figure 2 — Replica filtering removes pooled-only cavity cores and amplifies the surviving association estimate."""

from figure_common import *

values: dict[str, Any] = {"figure": 2, "panels": {}}
summary = read_json(DOWN / "replica_pocket_downstream_summary.json")
pooled_summary = read_json(DATA / "cavity_validation_rigor_summary.json")
recurrence = read_csv(DOWN / "replica_cluster_effective_recurrence.csv")
pooled = pooled_summary["stable_recurrent_core_ids"]
replica = summary["stable_recurrent_clusters"]
stable = recurrence.loc[recurrence["cluster"].isin(replica)].copy()

fig = plt.figure(figsize=(FIG_WIDTH_IN, 5.9))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.08, 1.0], hspace=0.42, wspace=0.40)

ax = fig.add_subplot(gs[0, 0])
panel_label(ax, "A")
panel_title(ax, "Pooled diagnostic → replica-filtered atlas")
stable_predecessors = set(stable["best_pooled_baseline_cluster"])
no_successor = [cluster for cluster in pooled if cluster not in stable_predecessors]
left_y = np.linspace(0.93, 0.12, len(pooled))
right_y = np.linspace(0.90, 0.15, len(replica))
pooled_y = dict(zip(pooled, left_y))
replica_y = dict(zip(replica, right_y))
predecessor = dict(zip(stable["cluster"], stable["best_pooled_baseline_cluster"]))
for r_cluster in replica:
    p_cluster = predecessor[r_cluster]
    if p_cluster in pooled_y:
        ax.plot([0.37, 0.63], [pooled_y[p_cluster], replica_y[r_cluster]], color="#BBBBBB", linewidth=0.7, zorder=0)
for cluster in pooled:
    color = COLORS["Significant"] if cluster in no_successor else "#444444"
    ax.text(0.34, pooled_y[cluster], cluster, ha="right", va="center", fontsize=5.2, color=color)
    if cluster in no_successor:
        ax.plot([0.28, 0.35], [pooled_y[cluster], pooled_y[cluster]], color=color, linewidth=0.8)
for cluster in replica:
    ax.text(0.66, replica_y[cluster], cluster, ha="left", va="center", fontsize=5.2, color="#222222", fontweight="bold")
ax.text(0.19, 1.01, f"Pooled diagnostic\n{pooled_summary['druggable_clusters']} clusters; {len(pooled)}-member core",
        ha="center", fontsize=5.5, fontweight="bold")
ax.text(0.81, 1.01, f"Replica-filtered\n{summary['druggable_consensus_clusters']} clusters; {len(replica)}-member core",
        ha="center", fontsize=5.5, fontweight="bold")
ax.text(0.5, 0.01,
    f"{len(no_successor)} of {len(pooled)} pooled core clusters do not survive\n"
    "receptor-disjoint replica reconstruction",
    ha="center", va="bottom", fontsize=5.4, color=COLORS["Significant"])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.08)
ax.axis("off")
values["panels"]["A"] = {
    "pooled_consensus_clusters": pooled_summary["druggable_clusters"],
    "pooled_core": pooled,
    "replica_consensus_clusters": summary["druggable_consensus_clusters"],
    "replica_core": replica,
    "pooled_without_replica_stable_successor": no_successor,
}

ax = fig.add_subplot(gs[0, 1])
panel_label(ax, "B")
panel_title(ax, "Stable-cluster correspondence")
predecessors = stable["best_pooled_baseline_cluster"].tolist()
unique_pred = list(dict.fromkeys(predecessors))
y_left = np.linspace(0.90, 0.10, len(unique_pred))
y_right = np.linspace(0.90, 0.10, len(replica))
pred_y = dict(zip(unique_pred, y_left))
rep_y = dict(zip(replica, y_right))
pooled_set = set(pooled)
highlight = {"R0", "R1", "R6"}
links = []
for row in stable.sort_values("cluster", key=lambda x: x.str.extract(r"(\d+)")[0].astype(int)).itertuples(index=False):
    promoted = row.best_pooled_baseline_cluster not in pooled_set
    color = COLORS["Significant"] if row.cluster in highlight else ("#332288" if promoted else "#999999")
    lw = 1.8 if row.cluster in highlight else 0.8
    ax.plot([0.32, 0.68], [pred_y[row.best_pooled_baseline_cluster], rep_y[row.cluster]], color=color, linewidth=lw, alpha=0.9)
    links.append({
        "replica_cluster": row.cluster,
        "pooled_predecessor": row.best_pooled_baseline_cluster,
        "core_jaccard": row.best_pooled_baseline_core_jaccard,
        "predecessor_in_pooled_core": not promoted,
    })
for p in unique_pred:
    promoted = p not in pooled_set
    ax.text(0.29, pred_y[p], p + ("  promoted" if promoted else ""), ha="right", va="center", fontsize=5.1,
            color="#332288" if promoted else "#333333", fontweight="bold" if promoted else "normal")
for r in replica:
    ax.text(0.71, rep_y[r], r, ha="left", va="center", fontsize=5.1, fontweight="bold" if r in highlight else "normal")
ax.text(0.20, 0.99, "Pooled predecessor", ha="center", fontsize=5.5, fontweight="bold")
ax.text(0.80, 0.99, "Replica stable cluster", ha="center", fontsize=5.5, fontweight="bold")
promoted_pairs = [(row.best_pooled_baseline_cluster, row.cluster) for row in stable.itertuples(index=False) if row.best_pooled_baseline_cluster not in pooled_set]
ax.text(0.5, 0.015,
    f"{' and '.join(pair[0] for pair in promoted_pairs)} were outside the provisional pooled core;\n"
    f"replica filtering promotes their successors {' and '.join(pair[1] for pair in promoted_pairs)}.",
    ha="center", va="bottom", fontsize=5.1, color="#332288")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.05)
ax.axis("off")
values["panels"]["B"] = {"links": links}

ax = fig.add_subplot(gs[1, 0])
panel_label(ax, "C")
panel_title(ax, "Split-half recovery")
stability = read_csv(DOWN / "replica_cluster_split_half_stability.csv")
stable_stability = stability.loc[stability["cluster"].isin(replica)].copy()
stable_stability["cluster_number"] = stable_stability["cluster"].str.extract(r"(\d+)").astype(int)
stable_stability = stable_stability.sort_values("cluster_number")
split_runs = read_csv(DOWN / "replica_cluster_split_half_runs.csv")
halves = split_runs.loc[split_runs["half"].isin([1, 2])]
# Reconstruct stable_rule from summary text
rule = summary["stable_recurrent_cluster_definition"]
import re
rec_match = re.search(r"n_receptors >= ([0-9]+)", rule)
rec_match2 = re.search(r"probability >= ([0-9.]+)", rule)
j_match = re.search(r"Jaccard >= ([0-9.]+)", rule)
min_receptors = int(rec_match.group(1)) if rec_match else 5
recovery_threshold = float(rec_match2.group(1)) if rec_match2 else 0.80
j_threshold = float(j_match.group(1)) if j_match else 0.50
x = np.arange(len(stable_stability))
ax.bar(x, stable_stability["split_half_recovery_probability_j_ge_0p5"], color="#6F6F6F", width=0.72)
ax.axhline(recovery_threshold, color=COLORS["Significant"], linestyle="--", linewidth=0.8)
ax.text(0.05, recovery_threshold + 0.008, f"prespecified {recovery_threshold:.2f}",
        color=COLORS["Significant"], ha="left", fontsize=5,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4, "alpha": 0.9})
ax.set_xticks(x, stable_stability["cluster"])
ax.set_ylim(0.75, 1.035)
ax.set_ylabel(f"Recovery probability (J≥{j_threshold:.2f})")
ax.text(0.01, 0.02, f"{split_runs['iteration'].nunique()} receptor-disjoint iterations\n"
        f"{len(halves)} half-atlas reconstructions", transform=ax.transAxes, fontsize=5, va="bottom")
inset = ax.inset_axes([0.51, 0.27, 0.47, 0.29])
inset.set_facecolor("white")
sensitivity = read_csv(DOWN / "replica_cluster_parameter_sensitivity.csv")
sensitivity_rows = []
for threshold in (0.66, 0.78):
    row = sensitivity.loc[sensitivity["linkage"].eq("average") & sensitivity["distance_threshold"].eq(threshold)].iloc[0]
    sensitivity_rows.append((f"Full\navg {threshold:.2f}", row["adjusted_rand_index_vs_baseline"], row["fraction_baseline_clusters_best_j_ge_0p5"]))
sensitivity_rows.insert(1, ("Split halves\navg 0.72", halves["adjusted_rand_index_vs_full"].mean(), halves["fraction_full_clusters_best_j_ge_0p5"].mean()))
single = sensitivity.loc[sensitivity["linkage"].eq("single")].iloc[0]
sensitivity_rows.append(("Full\nsingle 0.72", single["adjusted_rand_index_vs_baseline"], np.nan))
sx = np.arange(len(sensitivity_rows))
inset.scatter(sx, [r[1] for r in sensitivity_rows], marker="o", color=COLORS["Gi"], s=14, label="ARI")
for xi, row in zip(sx, sensitivity_rows):
    if np.isfinite(row[2]):
        inset.scatter(xi, row[2], marker="s", color=COLORS["Gq"], s=14)
inset.set_xticks(sx, ["Avg 0.66", "Split 0.72", "Avg 0.78", "Single 0.72"], fontsize=4.6, rotation=28, ha="right")
inset.tick_params(axis="x", pad=1, length=1.5)
inset.set_xlim(-0.45, len(sensitivity_rows) - 0.55)
inset.set_ylim(0, 1.03)
inset.tick_params(axis="y", labelsize=5)
inset.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=COLORS["Gi"],label="ARI"),
                       Line2D([0],[0],marker="s",color="none",markerfacecolor=COLORS["Gq"],label="Recovered")],
             frameon=False, fontsize=5, loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2, columnspacing=0.8, handletextpad=0.3)
inset.set_title("Parameter / receptor-split sensitivity", fontsize=5.5)
values["panels"]["C"] = {
    "stable_cluster_recovery": stable_stability.drop(columns="cluster_number").to_dict(orient="records"),
    "recovery_threshold": recovery_threshold, "core_jaccard_threshold": j_threshold,
    "split_iterations": summary["parameters"]["split_iterations"],
    "half_atlas_reconstructions": int(len(halves)),
    "sensitivity": [{"setting": label, "ari": ari, "fraction_recovered": recovered} for label, ari, recovered in sensitivity_rows],
}

ax = fig.add_subplot(gs[1, 1])
panel_label(ax, "D")
panel_title(ax, "Effect sharpening")
replica_sens = read_csv(DOWN / "legacy_C7_replica_sensitivity.csv")
pooled_sens = read_csv(SRC / "C7_Gq_sensitivity_variants_corrected.csv")
maxstat = read_csv(DOWN / "replica_cluster_family_maxstat_permutation.csv")
pooled_maxstat = read_csv(SRC / "C7_Gq_maxstat_all_classes_corrected.csv")
rows = []
pooled_row = pooled_sens.loc[pooled_sens["variant"].eq("all systems")].iloc[0]
rows.append({
    "label": "Pooled C7 × Gq\nsuperseded", "rr": pooled_row["risk_ratio"],
    "lo": pooled_row["risk_ratio_bootstrap_ci95_lo"], "hi": pooled_row["risk_ratio_bootstrap_ci95_hi"],
    "p": pooled_maxstat.loc[pooled_maxstat["analysis"].eq("all_classes") & pooled_maxstat["cluster"].eq("C7") & pooled_maxstat["family"].eq("Gq"), "max_stat_fwer_p_value"].iloc[0],
    "status": "superseded", "test_label": "FWER",
})
variants = [
    ("all systems", "Replica R6 × Gq\nprimary", "primary"),
    ("Class A systems", "R6 × Gq, Class A", "standard"),
    ("alpha5-hook-supported systems", "R6 × Gq, α5-hook", "standard"),
    ("receptor-family contexts", "R6 × Gq, receptor contexts", "standard"),
    ("single-family receptors", "R6 × Gq, single-family", "standard"),
    ("Direct G-alpha accession systems", "Direct-Gα accession\nunderpowered", "underpowered"),
]
for variant, label, status in variants:
    row = replica_sens.loc[replica_sens["variant"].eq(variant)].iloc[0]
    p = row["p_value_one_sided_fisher"]
    if variant in ("all systems", "Class A systems", "alpha5-hook-supported systems"):
        analysis_label = {"all systems": "replica atlas: all classes", "Class A systems": "replica atlas: Class A only",
                          "alpha5-hook-supported systems": "replica atlas: alpha5-hook-supported systems"}[variant]
        p = maxstat.loc[maxstat["analysis"].eq(analysis_label) & maxstat["cluster"].eq("R6") & maxstat["family"].eq("Gq"), "max_stat_fwer_p_value"].iloc[0]
    rows.append({
        "label": label, "rr": row["risk_ratio"],
        "lo": row["risk_ratio_bootstrap_ci95_lo"] if variant in ("all systems", "Class A systems", "alpha5-hook-supported systems") else np.nan,
        "hi": row["risk_ratio_bootstrap_ci95_hi"] if variant in ("all systems", "Class A systems", "alpha5-hook-supported systems") else np.nan,
        "katz_lo": row["risk_ratio_ci95_lo_katz"] if variant == "all systems" else np.nan,
        "katz_hi": row["risk_ratio_ci95_hi_katz"] if variant == "all systems" else np.nan,
        "p": p, "status": status,
        "test_label": "FWER" if variant in ("all systems", "Class A systems", "alpha5-hook-supported systems") else "Fisher",
    })
y = np.arange(len(rows))[::-1]
ax.axvline(1, color="#777777", linestyle="--", linewidth=0.7)
for yi, row in zip(y, rows):
    color = "#888888" if row["status"] == "superseded" else COLORS["Gq"]
    face = "none" if row["status"] == "underpowered" else color
    if np.isfinite(row.get("katz_lo", np.nan)):
        ax.hlines(yi, row["katz_lo"], row["katz_hi"], color="#D8A7B0", linewidth=3.0, alpha=0.8, zorder=1)
    if np.isfinite(row["lo"]):
        ax.hlines(yi, row["lo"], row["hi"], color=color, linewidth=1.2, zorder=2)
    ax.scatter(row["rr"], yi, s=25, facecolor=face, edgecolor=color, linewidth=0.8, zorder=3)
    p_label = f"FWER {row['p']:.1e}" if row["test_label"] == "FWER" else f"Fisher {row['p']:.3g}"
    ax.text(row["rr"] * 1.08, yi, f"{row['rr']:.2f}; {p_label}", va="center", fontsize=5)
fold = rows[1]["rr"] / rows[0]["rr"]
ax.annotate(f"{fold:.2f}× larger estimate", xy=(rows[1]["rr"], y[1]),
            xytext=(rows[0]["rr"] * 1.05, y[0] - 0.42), fontsize=5, color=COLORS["Gq"],
            arrowprops={"arrowstyle": "->", "color": COLORS["Gq"], "lw": 0.6})
ax.text(0.02, -0.14, "Direct-accession subset: underpowered — not evidence of absence", transform=ax.transAxes, fontsize=5, color="#666666")
ax.set_xscale("log")
ax.set_xlim(0.8, 120)
ax.set_yticks(y, [row["label"] for row in rows])
ax.set_xlabel("Risk ratio (log scale)")
values["panels"]["D"] = {
    "analyses": rows, "replica_to_pooled_fold_change": fold, "reference_risk_ratio": 1,
}
save_figure(fig, 2, values)
print("Figure 2 done.")
