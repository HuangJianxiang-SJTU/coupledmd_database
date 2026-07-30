#!/usr/bin/env python3
"""Figure 1 — The measurement problem: small inferential base, large replica noise."""

from figure_common import *

values: dict[str, Any] = {"figure": 1, "panels": {}}

s1 = read_csv(SRC / "provenance" / "Supplementary_Data_S1_system_inventory_corrected.csv")
homology = read_csv(DATA / "manifests" / "receptor_partitions.csv")
gateway_contrasts = read_csv(SRC / "gateway_context_noise_controlled_contrasts.csv")
pocket_context = read_csv(DOWN / "replica_pocket_context_receptor_inference.csv")
interface = read_csv(SRC / "interface_descriptor_outcome_correlations.csv")
summary = read_json(DOWN / "replica_pocket_downstream_summary.json")

receptor_keys = s1["receptor_uniprot"].fillna(s1["receptor_name"])
mapped = set(s1["receptor_uniprot"].dropna())
h40 = homology.loc[
    homology["homology_threshold"].eq(0.4)
    & homology["requested_accession"].isin(mapped)
]
stage_rows = [
    {"label": "Systems", "count": s1["system_id"].nunique(), "detail": ""},
    {"label": "Receptor entities", "count": receptor_keys.nunique(),
     "detail": f"{s1['receptor_uniprot'].nunique()} mapped accessions"},
    {"label": "H40 homology groups", "count": h40["homology_cluster_id"].nunique(), "detail": ""},
    {"label": "Gateway contrasts", "count": gateway_contrasts.shape[0],
     "detail": f"across {gateway_contrasts['receptor_uniprot'].nunique()} receptors"},
    {"label": "Pocket contrasts", "count": int(summary["distance_tier_inference"]["contrast_rows"]),
     "detail": f"across {pocket_context.shape[0]} receptors"},
    {"label": "Complete interface contrasts", "count": int(interface["n_receptors"].max()), "detail": ""},
]
values["panels"]["A"] = {"stages": stage_rows}

fig = plt.figure(figsize=(FIG_WIDTH_IN, 5.55))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.36)

ax = fig.add_subplot(gs[0, 0])
panel_label(ax, "A")
panel_title(ax, "Effective-sample funnel")
main = stage_rows[:3]
nested = stage_rows[3:]
y_main = np.array([5.0, 4.0, 3.0])
max_count = max(row["count"] for row in stage_rows)
for y, row in zip(y_main, main):
    ax.barh(y, row["count"], height=0.63, color="#D9D9D9", edgecolor="#666666", linewidth=0.5)
    text = f"{row['count']}  {row['label']}"
    if row["detail"]:
        text += f"\n({row['detail']})"
    ax.text(3, y, text, va="center", ha="left", fontsize=6)
nested_y = np.array([1.9, 1.0, 0.1])
nested_max = max(row["count"] for row in nested)
for y, row, color in zip(nested_y, nested, [COLORS["Gi"], COLORS["Gq"], COLORS["Gs"]]):
    ax.barh(y, row["count"], height=0.55, color=color, alpha=0.75, edgecolor="#555555", linewidth=0.4)
    text = f"{row['count']}  {row['label']}"
    if row["detail"]:
        text += f" ({row['detail']})"
    ax.text(row["count"] + 3, y, text, va="center", ha="left", fontsize=5.5)
ax.add_patch(patches.FancyBboxPatch(
    (-3, -0.32), nested_max + 118, 2.75,
    boxstyle="round,pad=0.02", facecolor="none", edgecolor="#777777",
    linewidth=0.5, linestyle=(0, (3, 2)),
))
ax.text(nested_max + 109, 2.12, "parallel contrast sets", fontsize=5, color="#666666", ha="right")
ax.set_xlim(0, max_count * 1.05)
ax.set_ylim(-0.5, 5.6)
ax.set_yticks([])
ax.set_xlabel("Effective units")
ax.spines["left"].set_visible(False)

ax = fig.add_subplot(gs[0, 1])
panel_label(ax, "B")
panel_title(ax, "Fingerprint variance partition")
variance = read_csv(DOWN / "replica_pocket_fingerprint_variance_partition.csv")
primary = variance.loc[variance["transform"].eq("Hellinger-normalized fingerprint")]
order = [
    ("GPCR class", "GPCR class"),
    ("receptor identity within class", "Receptor identity"),
    ("administrative G-protein family after receptor", "Administrative family"),
    ("system/context after receptor and family", "System / context"),
    ("residual replica variation within system", "Residual replica"),
]
left = 0.0
component_rows = []
for source_name, label in order:
    row = primary.loc[primary["component"].eq(source_name)].iloc[0]
    width = float(row["fraction_total_variation"])
    ax.barh(0.68, width, left=left, height=0.32, color=COLORS[label],
            edgecolor=("#666666" if label in ("GPCR class", "Administrative family") else "white"),
            linewidth=0.5, hatch="////" if label == "Administrative family" else None)
    if width > 0.04:
        ax.text(left + width / 2, 0.68, f"{width:.3f}", ha="center", va="center", fontsize=5,
                color="#333333" if label in ("GPCR class", "Administrative family") else "white")
    component_rows.append({
        "component": label, "value": width,
        "ci95": [row["receptor_bootstrap_ci95_lo"], row["receptor_bootstrap_ci95_hi"]],
    })
    left += width
for i, item in enumerate(component_rows):
    y = 0.40 - i * 0.060
    x = item["value"]
    lo, hi = item["ci95"]
    ax.errorbar(x, y, xerr=[[x - lo], [hi - x]], fmt="o", ms=2.8, color=COLORS[item["component"]],
                ecolor=COLORS[item["component"]], capsize=1.5, linewidth=0.7)
    label = item["component"]
    if label == "Administrative family":
        label += " (permutation n.s.)"
    ax.text(0.60, y, label, va="center", fontsize=5.0, color=COLORS[item["component"]])
family = primary.loc[primary["component"].eq("administrative G-protein family after receptor")].iloc[0]
family_p = summary["fingerprint_family_permutation"]["empirical_p_value"]
ax.text(0.02, 0.015,
    "Administrative-family increment beyond receptor identity:\n"
    f"{family['fraction_total_variation']:.4f} "
    f"({family['receptor_bootstrap_ci95_lo']:.4f}–{family['receptor_bootstrap_ci95_hi']:.4f})\n"
    f"permutation p={family_p:.3f}, not significant",
    fontsize=5, color="#777777", transform=ax.transAxes, va="bottom")
n_obs = summary["fingerprint_family_permutation"]["replica_observations"]
n_blocks = summary["fingerprint_family_permutation"]["receptor_blocks"]
ax.text(0.02, 0.96, f"n={n_obs} replicas; {n_blocks} receptor blocks", transform=ax.transAxes, fontsize=5.5, va="top")
ax.set_xlim(0, 1.0)
ax.set_ylim(-0.04, 1.02)
ax.set_yticks([])
ax.set_xlabel("Fraction of total Hellinger-normalised variation")
ax.spines["left"].set_visible(False)
values["panels"]["B"] = {
    "components": component_rows,
    "administrative_family_increment": {
        "value": family["fraction_total_variation"],
        "ci95": [family["receptor_bootstrap_ci95_lo"], family["receptor_bootstrap_ci95_hi"]],
        "permutation_p": family_p,
    },
    "n_replica_observations": n_obs,
    "n_receptor_blocks": n_blocks,
}

ax = fig.add_subplot(gs[1, 0])
panel_label(ax, "C")
panel_title(ax, "Three-tier fingerprint distances")
tiers = read_csv(DOWN / "replica_pocket_distance_tiers.csv")
tier_order = ["within-system replica noise", "same-receptor different-context", "different-receptor ceiling"]
tier_labels = ["Within-system\nreplica floor", "Same receptor,\ndifferent context", "Different receptor\n(matched)"]
tier_colors = [COLORS["Residual replica"], COLORS["System / context"], COLORS["Receptor identity"]]
rng = np.random.default_rng(20260729)
raw_values = {}
for i, (name, label, color) in enumerate(zip(tier_order, tier_labels, tier_colors)):
    data = tiers.loc[tiers["tier"].eq(name), "jaccard_distance"].to_numpy()
    raw_values[name] = data.tolist()
    violin_data = data
    if len(data) > 2500:
        violin_data = data[np.linspace(0, len(data) - 1, 2500, dtype=int)]
    parts = ax.violinplot(violin_data, positions=[i], widths=0.72, showmeans=False, showmedians=False, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor(color); body.set_edgecolor(color); body.set_alpha(0.28)
    strip = data if len(data) <= 300 else data[np.linspace(0, len(data) - 1, 300, dtype=int)]
    jitter = rng.uniform(-0.12, 0.12, len(strip))
    ax.scatter(np.full(len(strip), i) + jitter, strip, s=3, alpha=0.35, color=color, linewidths=0)
    mean = float(np.mean(data))
    median = float(np.median(data))
    ax.scatter([i], [mean], marker="D", s=17, color=color, edgecolor="white", linewidth=0.4, zorder=5)
    ax.hlines(median, i - 0.18, i + 0.18, color="#222222", linewidth=1.0, zorder=5)
    ax.text(i, 1.01, f"mean {mean:.3f}\nmedian {median:.3f}\nn={len(data):,}", ha="center", va="bottom", fontsize=5)
ax.set_xticks(range(3), tier_labels)
ax.set_ylabel("Jaccard distance")
ax.set_ylim(0, 1.12)
inference = summary["distance_tier_inference"]
within_system_median = float(np.median(
    tiers.loc[tiers["tier"].eq("within-system replica noise"), "jaccard_distance"]
))
ax.text(0.02, 0.03,
    f"{inference['receptors_context_exceeds_replica_noise']}/"
    f"{inference['receptors_with_complete_context_contrast']} receptors "
    "exceed replica noise; CRHR2 is the exception\n"
    f"median S/N={inference['median_receptor_signal_to_noise']:.3f} "
    f"({inference['median_receptor_signal_to_noise_bootstrap_ci95_lo']:.3f}–"
    f"{inference['median_receptor_signal_to_noise_bootstrap_ci95_hi']:.3f}); "
    f"one-sided sign p={inference['receptor_level_one_sided_sign_p']:.4f}",
    transform=ax.transAxes, fontsize=5, va="bottom")
pooled_median = summary["upstream_reducer_summary"]["median_within_system_replica_jaccard_distance"]
ax.text(0.02, -0.34,
    f"Pairwise pooled within-system median={pooled_median:.3f}; "
    f"the plotted {within_system_median:.3f} median is across per-system means.",
    transform=ax.transAxes, fontsize=5, color="#555555")
values["panels"]["C"] = {
    "raw_jaccard_distances": raw_values,
    "distance_tier_inference": inference,
    "pooled_pairwise_within_system_median": pooled_median,
}

ax = fig.add_subplot(gs[1, 1])
panel_label(ax, "D")
panel_title(ax, "Portal replica reliability")
icc = read_csv(SRC / "gateway_replica_icc.csv").sort_values("icc_1_1")
y = np.arange(len(icc))
ax.axvspan(0, 0.5, color="#F2F2F2", zorder=0)
ax.text(0.02, len(icc) - 0.35, "poor reliability", color="#777777", fontsize=5)
for i, row in enumerate(icc.itertuples(index=False)):
    color = PORTAL_COLORS[row.portal]
    ax.hlines(i, row.ci95_lo_system_bootstrap, row.ci95_hi_system_bootstrap, color=color, linewidth=1.2)
    ax.plot(row.icc_1_1, i, "o", color=color, ms=4)
    ax.text(row.icc_1_1 + 0.025, i, f"{row.icc_1_1:.3f}", va="center", fontsize=5)
ax.set_yticks(y, [PORTAL_LABELS[p] for p in icc["portal"]])
ax.set_xlim(0, 1.0)
ax.set_xlabel("ICC(1,1), system-bootstrap 95% CI")
ax.axvline(0.5, color="#999999", linestyle="--", linewidth=0.6)
replicas_per_system = int(summary["upstream_reducer_summary"]["replica_calls"] / summary["systems_complete"])
ax.text(0.98, 0.03, f"n={int(icc['n_systems'].iloc[0])} systems × {replicas_per_system} replicas",
        transform=ax.transAxes, ha="right", fontsize=5)
values["panels"]["D"] = {
    "n_systems": int(icc["n_systems"].iloc[0]),
    "portals": icc.to_dict(orient="records"),
    "poor_reliability_threshold": 0.5,
    "portal_colors": PORTAL_COLORS,
}
save_figure(fig, 1, values)
print("Figure 1 done.")
