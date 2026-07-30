#!/usr/bin/env python3
"""Supplementary Figure S3 — Opening mixed-model variance components."""

from figure_common import *

values: dict[str, Any] = {"figure": "S3", "panels": {}}
variance = read_csv(SRC / "gateway_variance_components.csv")

fig, ax = plt.subplots(figsize=(FIG_WIDTH_IN, 3.5))
panel_label(ax, "A")
panel_title(ax, "Per-portal mixed-model variance components")
ordered = variance.set_index("portal").loc[PORTALS].reset_index()
components = [
    ("fixed_fraction_of_total", "Fixed family", COLORS["GPCR class"]),
    ("receptor_fraction_of_total", "Receptor", COLORS["Receptor identity"]),
    ("system_fraction_of_total", "System/context", COLORS["System / context"]),
    ("replica_fraction_of_total", "Replica", COLORS["Residual replica"]),
]
left = np.zeros(len(ordered))
y = np.arange(len(ordered))
for column, label, color in components:
    ax.barh(y, ordered[column], left=left, color=color, label=label)
    left += ordered[column].to_numpy()
ax.set_yticks(y, [PORTAL_LABELS[p] for p in ordered["portal"]])
ax.set_xlim(0, 1)
ax.set_xlabel("Fraction of modeled variance")
ax.legend(frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.27))
n_obs = int(ordered["n_observations"].iloc[0])
n_sys = int(ordered["n_systems"].iloc[0])
n_rec = int(ordered["n_receptors"].iloc[0])
ax.text(0.99, 0.98, f"n={n_obs} observations; {n_sys} systems; {n_rec} receptors", transform=ax.transAxes, ha="right", va="top", fontsize=5)
values["panels"]["A"] = {"portals": ordered.to_dict(orient="records"), "n_observations": n_obs, "n_systems": n_sys, "n_receptors": n_rec}
save_supplementary(fig, "FigS3_gateway_variance_components", values)
print("Figure S3 done.")
