#!/usr/bin/env python3
"""Figure 5 — One construct-confounded cavity association and no global family code."""

from figure_common import *
import re

values: dict[str, Any] = {"figure": 5, "panels": {}}
summary = read_json(DOWN / "replica_pocket_downstream_summary.json")
maxstat = read_csv(DOWN / "replica_cluster_family_maxstat_permutation.csv")
presence = read_csv(DOWN / "replica_cluster_presence_systems.csv")

fig = plt.figure(figsize=(FIG_WIDTH_IN, 7.25))
gs = fig.add_gridspec(2, 2, height_ratios=[1.18, 1.0], width_ratios=[0.82, 1.18], hspace=0.40, wspace=0.34)

ax_a = fig.add_subplot(gs[0, 0])
panel_label(ax_a, "A")
panel_title(ax_a, "R6–Gq construct context")
placeholder(ax_a, "Gq_7EZM — CCKAR, R6 core + Gα α5 helix", 1.0)
values["panels"]["A"] = {"placeholder_note": "Gq_7EZM — CCKAR, R6 core + Gα α5 helix", "required_aspect_width_over_height": 1.0}

ax_container = fig.add_subplot(gs[0, 1])
panel_label(ax_container, "B")
panel_title(ax_container, f"The {len(maxstat)}-cell screen")
ax_container.axis("off")
analyses = ["replica atlas: all classes", "replica atlas: Class A only", "replica atlas: alpha5-hook-supported systems"]
facet_titles = ["All classes", "Class A", "α5-hook supported"]
clusters = sorted(maxstat["cluster"].unique(), key=lambda item: int(re.search(r"\d+", item).group()))
families = ["Gi", "Gq", "Gs", "G12-13"]
facets = []
for facet_index, (analysis, title) in enumerate(zip(analyses, facet_titles)):
    x0 = 0.01 + facet_index * 0.333
    ax = ax_container.inset_axes([x0, 0.08, 0.31, 0.76])
    screen = maxstat.loc[maxstat["analysis"].eq(analysis)].copy()
    matrix = np.ones((len(clusters), len(families)))
    observed = np.zeros_like(matrix)
    for i, cluster in enumerate(clusters):
        for j, family in enumerate(families[:3]):
            row = screen.loc[screen["cluster"].eq(cluster) & screen["family"].eq(family)].iloc[0]
            matrix[i, j] = row["max_stat_fwer_p_value"]
            observed[i, j] = row["observed_family_cluster_systems"]
        g12_mask = presence["g_protein_family"].eq("G12-13")
        if analysis == "replica atlas: Class A only":
            g12_mask &= presence["analysis_gpcr_class"].eq("A")
        elif analysis == "replica atlas: alpha5-hook-supported systems":
            g12_mask &= presence["interface_family_sequence_supported"].eq(True)
        observed[i, 3] = presence.loc[g12_mask, cluster].sum()
    for i in range(len(clusters)):
        for j, family in enumerate(families):
            if family == "G12-13":
                rect = patches.Rectangle((j-0.5,i-0.5), 1, 1, facecolor=COLORS["G12-13"], edgecolor="white", linewidth=0.15, hatch="////")
            else:
                sig = matrix[i, j] < 0.05
                rect = patches.Rectangle((j-0.5,i-0.5), 1, 1, facecolor=COLORS["Significant"] if sig else COLORS["Non-significant"], edgecolor="white", linewidth=0.15)
            ax.add_patch(rect)
            if observed[i, j] == 0:
                ax.plot(j, i, marker="x", ms=1.6, mew=0.35, color="white" if family != "Gs" else "#555555")
    r6_i = clusters.index("R6")
    ax.add_patch(patches.Rectangle((-0.5, r6_i-0.5), 3, 1, fill=False, edgecolor="#222222", linewidth=0.8))
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylim(len(clusters) - 0.5, -0.5)
    ax.set_xticks(range(4), ["Gi", "Gq", "Gs", "G12/13"], rotation=45, ha="right", fontsize=5)
    tick_indices = sorted(set(list(range(0, len(clusters), 6)) + [r6_i, len(clusters) - 1]))
    ax.set_yticks(tick_indices, [clusters[i] for i in tick_indices], fontsize=5)
    ax.set_title(title, fontsize=5.5)
    ax.tick_params(length=1.5)
    facets.append({
        "analysis": analysis, "clusters": clusters, "families": families,
        "max_stat_fwer_p": matrix.tolist(), "observed_family_cluster_systems": observed.tolist(),
    })
ax_container.legend(handles=[patches.Patch(facecolor=COLORS["Significant"], edgecolor="white", label="max-stat FWER p<0.05"),
                              patches.Patch(facecolor=COLORS["Non-significant"], edgecolor="white", label="non-significant"),
                              Line2D([0],[0],marker="x",color="#555555",linestyle="none",label="zero observed systems"),
                              patches.Patch(facecolor=COLORS["G12-13"],edgecolor="#777777",hatch="////",label="G12/13 descriptive only")],
                    frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=4.8, columnspacing=0.9, handletextpad=0.35)
survivors = maxstat.loc[maxstat["max_stat_fwer_p_value"].lt(0.05)]
survivor_pairs = survivors[["cluster", "family"]].drop_duplicates()
survivor_description = ", ".join(f"{row.cluster} × {row.family}" for row in survivor_pairs.itertuples(index=False))
ax_container.text(0.5, 0.93, f"Exactly {len(survivors)} survivors — {survivor_description} in all {len(analyses)} prespecified screens",
                  transform=ax_container.transAxes, ha="center", fontsize=5.5, fontweight="bold", color=COLORS["Significant"])
values["panels"]["B"] = {
    "facets": facets, "cells_per_inferential_screen": int(len(maxstat) / 3), "total_inferential_cells": int(len(maxstat)),
    "fwer_survivors": survivors.to_dict(orient="records"), "g12_13_status": "descriptive only; excluded from 486-cell inferential universe",
}

ax = fig.add_subplot(gs[1, 0])
panel_label(ax, "C")
panel_title(ax, "Association robustness")
mapping = read_csv(DOWN / "legacy_C7_replica_cluster_mapping.csv")
selected = mapping.loc[mapping["selected_outcome_blind_candidate"]].iloc[0]
next_best = mapping.loc[mapping["anatomical_replication_eligible"] & ~mapping["selected_outcome_blind_candidate"]].sort_values("core_jaccard_to_legacy_C7", ascending=False).iloc[0]
sensitivity = read_csv(DOWN / "legacy_C7_replica_sensitivity.csv")
powered_variants = ["all systems", "Class A systems", "alpha5-hook-supported systems", "receptor-family contexts", "single-family receptors"]
powered = sensitivity.loc[sensitivity["variant"].isin(powered_variants)]
resamples = read_csv(DOWN / "legacy_C7_replica_one_system_per_receptor.csv")
influence = read_csv(DOWN / "legacy_C7_replica_influence.csv")
mapping_callout = ("Outcome-blind anatomical mapping (not a p value):\n"
    f"R6↔C7 J={selected['core_jaccard_to_legacy_C7']:.2f}, {int(selected['shared_core_positions'])}/{int(selected['legacy_core_size'])} positions; "
    f"next eligible {next_best['replica_cluster']}={next_best['core_jaccard_to_legacy_C7']:.3f}")
robustness_rows = [
    {"label": "Powered restrictions", "value": -math.log10(powered["p_value_one_sided_fisher"].max()),
     "display": f"RR {powered['risk_ratio'].min():.2f}–{powered['risk_ratio'].max():.2f}; max p={powered['p_value_one_sided_fisher'].max():.2e}"},
    {"label": "1 system / receptor", "value": -math.log10(resamples["p_value_one_sided_fisher"].max()),
     "display": f"{(resamples['p_value_one_sided_fisher'] < 0.05).sum():,}/{len(resamples):,} significant; max p={resamples['p_value_one_sided_fisher'].max():.5f}"},
    {"label": "Leave one receptor", "value": -math.log10(influence.loc[influence["leave_one_level"].eq("receptor"), "p_value_one_sided_fisher"].max()),
     "display": f"max p={influence.loc[influence['leave_one_level'].eq('receptor'), 'p_value_one_sided_fisher'].max():.2e}"},
    {"label": "Leave one H40", "value": -math.log10(influence.loc[influence["leave_one_level"].eq("H40 homology group"), "p_value_one_sided_fisher"].max()),
     "display": f"max p={influence.loc[influence['leave_one_level'].eq('H40 homology group'), 'p_value_one_sided_fisher'].max():.2e}"},
]
y = np.arange(len(robustness_rows))[::-1]
max_bar = max(row["value"] for row in robustness_rows)
for yi, row in zip(y, robustness_rows):
    ax.barh(yi, row["value"], color=COLORS["Significant"], height=0.55)
    ax.text(row["value"] + 0.12, yi, row["display"], va="center", fontsize=5)
ax.set_yticks(y, [row["label"] for row in robustness_rows])
ax.set_xlabel("−log₁₀ maximum one-sided Fisher p")
ax.set_xlim(0, max_bar * 1.65)
ax.set_ylim(-0.65, len(robustness_rows) + 0.35)
ax.text(0.01, 0.98, mapping_callout, transform=ax.transAxes, va="top", fontsize=5, color="#555555",
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "#F2F2F2", "edgecolor": "#999999", "linewidth": 0.4}, wrap=True)
ax.text(0.01, 0.01, "All powered and resampling checks retain the association.", transform=ax.transAxes, fontsize=5, fontweight="bold")
values["panels"]["C"] = {
    "successor_mapping": mapping.to_dict(orient="records"), "powered_sensitivity_variants": powered.to_dict(orient="records"),
    "one_system_per_receptor_draws": resamples.to_dict(orient="records"), "leave_one_out": influence.to_dict(orient="records"),
    "anatomical_mapping_callout": mapping_callout, "display_rows": robustness_rows,
}

ax = fig.add_subplot(gs[1, 1])
panel_label(ax, "D")
panel_title(ax, "The evidence collapse")
evidence = read_csv(SRC / "figure5_R6_legacy_C7_Gq_construct_evidence_replica_final207.csv")
order = ["Simulated hook supports Gq", "Deposited = simulated full hook", "Direct accession supports label"]
evidence = evidence.set_index("evidence_tier").loc[order].reset_index()
cohort_r6_gq_total = int(evidence["systems"].max())
direct_support = int(evidence.loc[evidence["evidence_tier"].eq("Direct accession supports label"), "systems"].iloc[0])
x = np.arange(len(evidence))
bar_colors = [COLORS["Gq"], "#C8A1A9", "#333333"]
bars = ax.bar(x, evidence["systems"], color=bar_colors, width=0.62, edgecolor="white", linewidth=0.5)
for bar, count in zip(bars, evidence["systems"]):
    ax.text(bar.get_x() + bar.get_width() / 2, count + 0.25, f"{int(count)} / {cohort_r6_gq_total}", ha="center", fontsize=7, fontweight="bold")
for i in range(len(evidence) - 1):
    ax.annotate("", xy=(i+1-0.31, evidence.loc[i+1, "systems"]+0.05), xytext=(i+0.31, evidence.loc[i, "systems"]-0.05),
                arrowprops={"arrowstyle": "->", "lw": 1.0, "color": "#555555"})
ax.set_xticks(x, ["Gq-like simulated\nα5 hook", "Deposited = simulated\nfull hook", "Direct Gα accession\nsupport"])
ax.set_ylim(0, cohort_r6_gq_total + 1.4)
ax.set_ylabel("R6–administrative-Gq systems")
presence2 = read_csv(DOWN / "replica_cluster_presence_systems.csv")
r6_gq = presence2.loc[presence2["R6"].eq(1) & presence2["g_protein_family"].eq("Gq")]
hook = read_csv(SRC / "gprotein_alpha5_hook_sequence_audit.csv")
hook_r6 = hook.loc[hook["system_id"].isin(r6_gq["system_id"])]
near = hook_r6.loc[hook_r6["simulation_deposition_hook_status"].eq("near_match_sequence_difference")]
hook_summary = read_json(DATA / "audits" / "gprotein_alpha5_hook_audit_summary.json")
ax.text(0.5, 0.57, f"Direct Gα-accession support\nin {direct_support} {'system' if direct_support == 1 else 'systems'}",
        transform=ax.transAxes, ha="center", va="center", fontsize=8, fontweight="bold", color="#222222",
        bbox={"boxstyle": "round,pad=0.20", "facecolor": "white", "edgecolor": "none", "alpha": 0.88})
ax.text(0.02, 0.86, "Two-position deposited/simulated differences, still Gq-like:\n" + ", ".join(near["system_id"].tolist()),
        transform=ax.transAxes, va="top", fontsize=5)
ax.text(0.98, 0.02, f"Cohort-wide: hook supports {hook_summary['interface_family_sequence_supported']}/{hook_summary['systems']} labels; {hook_summary['assigned_family_hook_conflicts']} conflicts retained",
        transform=ax.transAxes, ha="right", fontsize=5)
values["panels"]["D"] = {
    "evidence_tiers": evidence.to_dict(orient="records"),
    "r6_administrative_gq_systems": r6_gq["system_id"].tolist(),
    "r6_gq_hook_audit_rows": hook_r6.to_dict(orient="records"),
    "near_match_two_position_systems": near["system_id"].tolist(),
    "cohort_hook_supported": hook_summary["interface_family_sequence_supported"],
    "cohort_total": hook_summary["systems"],
    "cohort_conflicts_retained": hook_summary["assigned_family_hook_conflicts"],
}
values["panels"]["A"]["reserved_box"] = report_placeholder(fig, ax_a, "5A")
save_figure(fig, 5, values)
print("Figure 5 done.")
