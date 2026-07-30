#!/usr/bin/env python3
"""Supplementary Figure S4 — Interface descriptor reliability and descriptor–outcome power bounds."""

from figure_common import *

values: dict[str, Any] = {"figure": "S4", "panels": {}}
icc = read_csv(SRC / "interface_descriptor_replica_icc.csv")
corr = read_csv(SRC / "interface_descriptor_outcome_correlations.csv")
systems = read_csv(SRC / "interface_descriptor_systems.csv")

fig, axes = plt.subplots(1, 2, figsize=(FIG_WIDTH_IN, 3.8))
fig.subplots_adjust(wspace=0.46)

ax = axes[0]
panel_label(ax, "A")
panel_title(ax, "Interface-descriptor replica reliability")
ordered = icc.sort_values("icc_1_1")
y = np.arange(len(ordered))
ax.barh(y, ordered["icc_1_1"], color=COLORS["Residual replica"])
ax.set_yticks(y, ordered["feature_label"], fontsize=5)
ax.set_xlim(0, 1)
ax.set_xlabel("ICC(1,1)")
n_systems_max = int(ordered["n_systems_complete_triplicates"].max())
n_systems_min = int(ordered["n_systems_complete_triplicates"].min())
replica_counts = systems["replica_means"].dropna().str.split("|").str.len()
n_replicas = int(replica_counts.mode().iloc[0])
ax.text(0.02, 0.98, f"range {ordered['icc_1_1'].min():.3f}–{ordered['icc_1_1'].max():.3f}\nn={n_systems_min}–{n_systems_max} systems × {n_replicas} replicas (up to {n_systems_max})",
        transform=ax.transAxes, va="top", fontsize=5)
values["panels"]["A"] = {
    "features": ordered.to_dict(orient="records"), "n_systems_min": n_systems_min,
    "n_systems_max": n_systems_max, "n_replicas": n_replicas,
}

ax = axes[1]
panel_label(ax, "B")
panel_title(ax, "Descriptor–outcome patterns (descriptive only)")
outcome_order = ["gateway_mean_absolute_difference", "pocket_set_jaccard_distance"]
matrix = corr.pivot(index="feature_label", columns="outcome", values="spearman_rho")
feature_order = icc["feature_label"].tolist()
matrix = matrix.reindex(index=feature_order, columns=outcome_order)
im = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        ax.text(j, i, f"{matrix.iloc[i, j]:.2f}", ha="center", va="center", fontsize=5)
ax.set_yticks(range(len(feature_order)), feature_order, fontsize=5)
ax.set_xticks(range(2), ["Gateway difference", "Pocket distance"], rotation=25, ha="right")
cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
cbar.set_label("Spearman ρ", fontsize=6)
cbar.ax.tick_params(labelsize=5)
min_detectable = corr["rho_for_80pct_power_nominal_alpha_0p05"].dropna().iloc[0]
ax.text(0.02, -0.32, f"Only 4 complete receptor contrasts; inferentially ineligible. Minimum detectable |ρ|≈{min_detectable:.3f} at 80% power.",
        transform=ax.transAxes, fontsize=5)
values["panels"]["B"] = {
    "correlations": corr.to_dict(orient="records"), "n_complete_contrasts": int(corr["n_contrasts"].max()),
    "minimum_detectable_abs_rho_nominal": min_detectable,
}
save_supplementary(fig, "FigS4_interface_descriptors", values)
print("Figure S4 done.")
