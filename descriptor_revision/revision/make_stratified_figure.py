"""Stratified structural-validation metrics figure (Reviewer 1, comment 3).
Merges Supplementary Data S5 (per-replica metrics) with S1 (class/family) and
plots violin + jitter distributions by G-protein family and GPCR class.
"""
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/MDdata/data02/jxhuang/gpcr_g/coupledmd_database/manuscript_figure_package_v14/supplementary_files/"
OUT = "Figure_stratified_metrics"

s1 = pd.read_csv(BASE + "Supplementary_Data_S1_system_inventory.csv")
s5 = pd.read_csv(BASE + "Supplementary_Data_S5_structural_validation.csv")
df = s5.merge(s1[["system_id", "gpcr_class"]], on="system_id", how="left")
df = df.dropna(subset=["tm_core_rmsd_A_p95"])  # drops 3 null-metric Gs_8HTI replicas

fam_map = {"Gi": "Gi/o", "Gq": "Gq/11", "Gs": "Gs", "G12-13": "G12/13"}
df["family"] = df.g_protein_family.map(fam_map)
fam_order = ["Gs", "Gi/o", "Gq/11", "G12/13"]
metrics = [("tm_core_rmsd_A_p95", "TM-core Ca RMSD P95 (A)"),
           ("galpha_interface_rmsd_A_p95", "Ga-interface Ca RMSD P95 (A)"),
           ("contact_retention_p05", "Interface contact retention P05")]
colors = {"Gs": "#4C9BD6", "Gi/o": "#E07B39", "Gq/11": "#54A24B", "G12/13": "#B279A2"}

fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.5))
rng = np.random.default_rng(7)
for j, (col, label) in enumerate(metrics):
    for i, (split, title) in enumerate([("family", "by G-protein family"), ("gpcr_class", "by GPCR class")]):
        ax = axes[i, j]
        groups = fam_order if split == "family" else ["A", "B"]
        data = [df.loc[df[split] == g, col].values for g in groups]
        vp = ax.violinplot(data, showmedians=True, widths=0.8)
        for k, (body, g) in enumerate(zip(vp["bodies"], groups)):
            body.set_facecolor(colors[g] if split == "family" else "#7f7f7f")
            body.set_alpha(0.55)
        for k, (g, d) in enumerate(zip(groups, data)):
            ax.scatter(rng.normal(k + 1, 0.05, len(d)), d, s=4, alpha=0.25, color="k", zorder=1)
            ax.text(k + 1, ax.get_ylim()[0], f"n={len(d)}", ha="center", va="bottom", fontsize=7)
        ax.set_xticks(range(1, len(groups) + 1), groups)
        if i == 0:
            ax.set_title(label, fontsize=10)
        if j == 0:
            ax.set_ylabel(title, fontsize=10)
        ax.tick_params(labelsize=9)
axes[0, 0].set_ylim(bottom=0)
axes[1, 0].set_ylim(bottom=0)
fig.suptitle("Structural-validation metrics stratified by G-protein family and GPCR class "
             "(618 replica records; descriptive)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=300)
