#!/usr/bin/env python3
"""Supplementary Figure S1 — Cohort composition and receptor-by-family analytical design."""

from figure_common import *

values: dict[str, Any] = {"figure": "S1", "panels": {}}
cohort = read_csv(SRC / "figure1_cohort_by_family_class.csv")
design = read_csv(SRC / "figure1_pocket_contrast_design.csv")
design.loc[design["receptor_uniprot"].eq("P48960"), "receptor"] = "ADGRE5"

fig, axes = plt.subplots(1, 2, figsize=(FIG_WIDTH_IN, 3.3), gridspec_kw={"width_ratios": [0.85, 1.15]})
fig.subplots_adjust(wspace=0.36)

ax = axes[0]
panel_label(ax, "A")
panel_title(ax, "Cohort composition by family and GPCR class")
families = ["Gi", "Gq", "Gs", "G12-13"]
y = np.arange(len(families))
left = np.zeros(len(families))
for gpcr_class, label in (("A", "Class A"), ("B", "Class B")):
    counts = np.array([cohort.loc[cohort["g_protein_family"].eq(family) & cohort["analysis_gpcr_class"].eq(gpcr_class), "systems"].sum() for family in families])
    ax.barh(y, counts, left=left, color=COLORS[label], label=label)
    for yi, lft, count in zip(y, left, counts):
        if count:
            ax.text(lft + count / 2, yi, str(int(count)), ha="center", va="center", fontsize=5, color="white" if gpcr_class == "A" else "#222222")
    left += counts
ax.set_yticks(y, ["Gi", "Gq", "Gs", "G12/13"])
ax.set_xlabel("Systems")
ax.legend(frameon=False)
values["panels"]["A"] = cohort.to_dict(orient="records")

ax = axes[1]
panel_label(ax, "B")
panel_title(ax, "Receptor-by-family contrast design")
matrix = design.pivot_table(index=["receptor", "receptor_uniprot"], columns="family", values="n_systems", fill_value=0).reindex(columns=families, fill_value=0)
matrix["_span"] = (matrix[families] > 0).sum(axis=1)
matrix["_systems"] = matrix[families].sum(axis=1)
matrix = matrix.sort_values(["_span", "_systems"], ascending=[False, False])
display = matrix[families]
cmap = mpl.colors.ListedColormap(["#F4F4F4", COLORS["Gi"]])
ax.imshow((display.to_numpy() > 0).astype(int), aspect="auto", cmap=cmap, vmin=0, vmax=1)
for i in range(display.shape[0]):
    for j in range(display.shape[1]):
        count = int(display.iloc[i, j])
        if count:
            ax.text(j, i, str(count), ha="center", va="center", fontsize=5, color="white")
ax.set_xticks(range(4), ["Gi", "Gq", "Gs", "G12/13"])
ax.set_yticks(range(len(display)), [idx[0] for idx in display.index], fontsize=5)
ax.set_xlabel("Administrative family")
values["panels"]["B"] = {
    "receptors": [idx[0] for idx in display.index],
    "receptor_uniprots": [idx[1] for idx in display.index],
    "family_counts": display.to_dict(orient="split"),
}
save_supplementary(fig, "FigS1_cohort_design", values)
print("Figure S1 done.")
