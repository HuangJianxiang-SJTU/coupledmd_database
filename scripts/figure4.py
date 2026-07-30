#!/usr/bin/env python3
"""Figure 4 — TM3–TM4 leads a reproducible inter-helical-opening hierarchy that prepared geometry does not explain."""

from figure_common import *

values: dict[str, Any] = {"figure": 4, "panels": {}}
gateway_summary = read_json(DATA / "gateway_rigor_summary.json")
threshold = read_csv(SRC / "gateway_threshold_sensitivity.csv")
outcomes = read_csv(SRC / "gateway_static_geometry_outcomes.csv")
outcome_means = outcomes.groupby("portal")["mean_open_fraction"].mean()
depth_discrepancy = {
    float(depth): float(
        (rows.set_index("portal")["mean"].reindex(outcome_means.index) - outcome_means).abs().max()
    )
    for depth, rows in threshold.groupby("depth_threshold_angstrom")
}
primary_depth = min(depth_discrepancy, key=depth_discrepancy.get)
primary = threshold.loc[threshold["depth_threshold_angstrom"].eq(primary_depth)].copy()

fig = plt.figure(figsize=(FIG_WIDTH_IN, 7.65))
gs = fig.add_gridspec(3, 2, height_ratios=[0.88, 1.05, 1.05], width_ratios=[0.88, 1.12], hspace=0.48, wspace=0.36)

ax_a = fig.add_subplot(gs[0, 0])
panel_label(ax_a, "A")
panel_title(ax_a, "Seven lateral portals")
placeholder(ax_a, "Gi_6K42 — ADRA2B–Gi, seven portal definitions", 4 / 3)
values["panels"]["A"] = {"placeholder_note": "Gi_6K42 — ADRA2B–Gi, seven portal definitions", "required_aspect_width_over_height": 4 / 3}

ax = fig.add_subplot(gs[0, 1])
panel_label(ax, "B")
panel_title(ax, "Seven-portal schematic with cohort means")
theta = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, 7, endpoint=False)
coords = np.c_[np.cos(theta), np.sin(theta)]
mean_map = dict(zip(primary["portal"], primary["mean"]))
norm = Normalize(vmin=primary["mean"].min(), vmax=primary["mean"].max())
for i, portal in enumerate(PORTALS):
    j = (i + 1) % 7
    xy1, xy2 = coords[i], coords[j]
    value = mean_map[portal]
    ax.plot([xy1[0], xy2[0]], [xy1[1], xy2[1]], color=PORTAL_COLORS[portal], linewidth=0.8+5.0*value, solid_capstyle="round", alpha=0.86)
    edge_mid = (xy1 + xy2) / 2
    mid = edge_mid / np.linalg.norm(edge_mid) * 1.62
    ax.text(mid[0], mid[1], f"{PORTAL_LABELS[portal]}\n{value:.3f}", ha="center", va="center", fontsize=5,
            color=PORTAL_COLORS[portal], fontweight="bold" if portal in ("TM3-TM4", "TM7-TM1") else "normal")
for i, (x, y) in enumerate(coords, 1):
    ax.add_patch(patches.Circle((x, y), 0.18, facecolor="white", edgecolor="#444444", linewidth=0.6, zorder=3))
    ax.text(x, y, f"TM{i}", ha="center", va="center", fontsize=5.5, zorder=4)
ax.text(0, 0, f"cohort mean\nopen fraction\n{primary_depth:.1f} Å threshold", ha="center", va="center", fontsize=5.2)
ax.set_aspect("equal")
ax.set_xlim(-2.05, 2.05)
ax.set_ylim(-1.95, 1.95)
ax.axis("off")
values["panels"]["B"] = {
    "primary_depth_threshold_angstrom": primary_depth, "primary_threshold_source_discrepancy": depth_discrepancy,
    "portal_cohort_means": mean_map, "portal_colors": PORTAL_COLORS, "n_systems": gateway_summary["systems"],
}

ax = fig.add_subplot(gs[1, :])
panel_label(ax, "C")
panel_title(ax, "Rank stability across opening thresholds")
depths = sorted(threshold["depth_threshold_angstrom"].unique())
for portal in PORTALS:
    rows = threshold.loc[threshold["portal"].eq(portal)].sort_values("depth_threshold_angstrom")
    lead = portal in ("TM3-TM4", "TM7-TM1")
    color = PORTAL_COLORS[portal] if lead else "#C7C7C7"
    ax.plot(rows["depth_threshold_angstrom"], rows["rank"], marker="o", ms=4 if lead else 2.5,
            linewidth=2.0 if lead else 0.7, color=color, alpha=1.0 if lead else 0.85,
            label=PORTAL_LABELS[portal], zorder=3 if lead else 1)
    if lead:
        ax.text(depths[-1] + 0.025, rows.iloc[-1]["rank"], PORTAL_LABELS[portal], color=color, fontsize=6, va="center", fontweight="bold")
ax.set_xticks(depths, [f"{d:.1f}" for d in depths])
ax.set_xlim(min(depths) - 0.03, max(depths) + 0.16)
ax.set_ylim(7.35, 0.65)
ax.set_yticks(range(1, 8))
ax.set_xlabel("Depth threshold (Å)")
ax.set_ylabel("Portal rank")
ax.text(0.01, 0.04, "TM3–TM4 ranks first and the canonical TM7–TM1 portal second at every threshold.",
        transform=ax.transAxes, fontsize=5.5, fontweight="bold")
values["panels"]["C"] = {
    "threshold_rank_rows": threshold.to_dict(orient="records"),
    "leading_portals_at_every_threshold": gateway_summary["threshold_top_two"], "portal_colors": PORTAL_COLORS,
}

sub_d = gs[2, 0].subgridspec(1, 2, width_ratios=[1.1, 0.9], wspace=0.46)
ax_d1 = fig.add_subplot(sub_d[0, 0])
panel_label(ax_d1, "D")
panel_title(ax_d1, "Prepared-geometry control")
associations = read_csv(SRC / "gateway_static_geometry_associations.csv")
assoc = associations.loc[associations["analysis_status"].eq("primary inferential")].copy()
assoc["portal_order"] = assoc["portal"].map({p: i for i, p in enumerate(PORTALS)})
assoc = assoc.sort_values("portal_order", ascending=False)
y = np.arange(len(assoc))
for yi, row in zip(y, assoc.itertuples(index=False)):
    color = PORTAL_COLORS[row.portal]
    ax_d1.hlines(yi, row.receptor_bootstrap_ci95_lo, row.receptor_bootstrap_ci95_hi, color=color, linewidth=1)
    ax_d1.plot(row.spearman_rho, yi, "o", color=color, ms=3.5)
ax_d1.axvline(0, color="#777777", linestyle="--", linewidth=0.6)
ax_d1.set_yticks(y, [PORTAL_LABELS[p] for p in assoc["portal"]])
ax_d1.set_xlabel("Within-portal Spearman ρ")
n_significant = int(assoc["q_bh_7_primary_portals"].lt(0.05).sum())
ax_d1.text(0.02, 0.02, f"{n_significant}/{len(assoc)} significant after BH correction\nn={int(assoc['n_systems'].iloc[0])} systems; {int(assoc['n_receptors'].iloc[0])} receptor blocks",
           transform=ax_d1.transAxes, fontsize=5, va="bottom")

ax_d2 = fig.add_subplot(sub_d[0, 1])
portal_means = outcomes.groupby("portal", as_index=False).agg(mean_prepared_width=("ca_symmetric_nearest_mean_angstrom", "mean"), mean_open_fraction=("mean_open_fraction", "mean"))
for row in portal_means.itertuples(index=False):
    ax_d2.scatter(row.mean_prepared_width, row.mean_open_fraction, s=22, color=PORTAL_COLORS[row.portal])
    if row.portal in ("TM3-TM4", "TM7-TM1"):
        ax_d2.annotate(PORTAL_LABELS[row.portal], (row.mean_prepared_width, row.mean_open_fraction), xytext=(2,2), textcoords="offset points", fontsize=5, color=PORTAL_COLORS[row.portal])
static_summary = read_json(DATA / "gateway_static_geometry_summary.json")
ax_d2.text(0.98, 0.98, f"Between portals\nρ={static_summary['across_7_portal_mean_rank_spearman']:.3f}, p={static_summary['across_7_portal_mean_rank_naive_p']:.3f}",
           transform=ax_d2.transAxes, ha="right", va="top", fontsize=5.2)
ax_d2.set_xlabel("Prepared width (Å)")
ax_d2.set_ylabel("Mean open fraction")
values["panels"]["D"] = {
    "within_portal_associations": assoc.drop(columns="portal_order").to_dict(orient="records"),
    "between_portal_means": portal_means.to_dict(orient="records"),
    "between_portal_rank_spearman": static_summary["across_7_portal_mean_rank_spearman"],
    "between_portal_rank_p": static_summary["across_7_portal_mean_rank_naive_p"],
    "n_systems": static_summary["systems"], "n_receptors": static_summary["receptors"], "portal_colors": PORTAL_COLORS,
}

ax = fig.add_subplot(gs[2, 1])
panel_label(ax, "E")
panel_title(ax, "Stationarity: last minus first 100 ns")
stationarity = read_csv(SRC / "gateway_stationarity_tests.csv")
stationarity["portal_order"] = stationarity["portal"].map({p: i for i, p in enumerate(PORTALS)})
stationarity = stationarity.sort_values("portal_order", ascending=False)
y = np.arange(len(stationarity))
for yi, row in zip(y, stationarity.itertuples(index=False)):
    color = PORTAL_COLORS[row.portal]
    ax.hlines(yi, row.ci95_lo_receptor_bootstrap, row.ci95_hi_receptor_bootstrap, color=color, linewidth=1.1)
    face = color if row.q_bh_7_portals < 0.05 else "white"
    ax.plot(row.mean_last_minus_first, yi, marker="o", ms=4, markerfacecolor=face, markeredgecolor=color, markeredgewidth=0.8)
ax.axvline(0, color="#777777", linestyle="--", linewidth=0.6)
ax.set_yticks(y, [PORTAL_LABELS[p] for p in stationarity["portal"]])
ax.set_xlabel("Open-fraction change (last − first)")
rank_rho = gateway_summary["stationarity"]["first_last_portal_rank_spearman_rho"]
tm34 = stationarity.loc[stationarity["portal"].eq("TM3-TM4")].iloc[0]
ax.text(0.02, 0.98, f"Rank correlation ρ={rank_rho:.3f}\nTM3–TM4: Δ={tm34['mean_last_minus_first']:.3f}, 95% CI {tm34['ci95_lo_receptor_bootstrap']:.3f} to {tm34['ci95_hi_receptor_bootstrap']:.3f}; q={tm34['q_bh_7_portals']:.3f}",
        transform=ax.transAxes, va="top", fontsize=5.2)
ax.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor="#666666",markeredgecolor="#666666",label="Filled: BH q<0.05"),
                   Line2D([0],[0],marker="o",color="none",markerfacecolor="white",markeredgecolor="#666666",label="Open: BH q≥0.05")],
          frameon=False, loc="lower right", fontsize=5, handletextpad=0.3)
ax.text(0.98, -0.12, f"n={int(stationarity['n_systems'].iloc[0])} systems; {int(stationarity['n_receptors'].iloc[0])} receptors",
        transform=ax.transAxes, ha="right", fontsize=5)
values["panels"]["E"] = {
    "portal_changes": stationarity.drop(columns="portal_order").to_dict(orient="records"),
    "first_last_rank_spearman": rank_rho, "n_systems": int(stationarity['n_systems'].iloc[0]),
    "n_receptors": int(stationarity['n_receptors'].iloc[0]), "portal_colors": PORTAL_COLORS,
}
values["panels"]["A"]["reserved_box"] = report_placeholder(fig, ax_a, "4A")
save_figure(fig, 4, values)
print("Figure 4 done.")
