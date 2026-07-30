#!/usr/bin/env python3
"""Supplementary Figure S2 — Within-receptor context contrasts and recorded confounds."""

from figure_common import *

values: dict[str, Any] = {"figure": "S2", "panels": {}}
contrasts = read_csv(SRC / "figure4_contrast_replica_context_replica_final207.csv")
receptor = read_csv(DOWN / "replica_pocket_context_receptor_inference.csv")
contrasts["receptor"] = contrasts["receptor"].replace({"P48960": "ADGRE5"})
receptor["receptor"] = receptor["receptor"].replace({"P48960": "ADGRE5"})

fig, axes = plt.subplots(2, 2, figsize=(FIG_WIDTH_IN, 5.7))
fig.subplots_adjust(hspace=0.46, wspace=0.38)

ax = axes[0, 0]
panel_label(ax, "A")
panel_title(ax, "Gateway context difference versus sampling floor")
significant = contrasts["gateway_sampling_null_q"].lt(0.05)
ax.scatter(contrasts["gateway_expected_triplicate_sampling_mad"], contrasts["gateway_mean_absolute_difference"],
           c=np.where(significant, COLORS["Significant"], COLORS["Non-significant"]), s=20)
limit = max(contrasts["gateway_expected_triplicate_sampling_mad"].max(), contrasts["gateway_mean_absolute_difference"].max())
ax.plot([0, limit], [0, limit], "--", color="#777777", linewidth=0.7)
for row in contrasts.loc[significant].itertuples(index=False):
    ax.annotate(row.receptor, (row.gateway_expected_triplicate_sampling_mad, row.gateway_mean_absolute_difference), xytext=(2,2), textcoords="offset points", fontsize=5)
ax.set_xlabel("Expected triplicate sampling MAD")
ax.set_ylabel("Observed context MAD")
ax.text(0.98, 0.02, f"{significant.sum()}/{len(contrasts)} exceed sampling null; BH q<0.05", transform=ax.transAxes, ha="right", fontsize=5)
values["panels"]["A"] = contrasts[["receptor","contrast","gateway_expected_triplicate_sampling_mad","gateway_mean_absolute_difference","gateway_sampling_null_p","gateway_sampling_null_q"]].to_dict(orient="records")

ax = axes[0, 1]
panel_label(ax, "B")
panel_title(ax, "Pocket context distance by contrast")
ordered = contrasts.sort_values("between_context_mean_replica_jaccard")
y = np.arange(len(ordered))
ax.barh(y, ordered["between_context_mean_replica_jaccard"], color=COLORS["System / context"])
ax.set_yticks(y, ordered["receptor"] + " " + ordered["contrast"], fontsize=5)
ax.set_xlabel("Between-context mean replica Jaccard distance")
values["panels"]["B"] = ordered[["receptor","contrast","between_context_mean_replica_jaccard"]].to_dict(orient="records")

ax = axes[1, 0]
panel_label(ax, "C")
panel_title(ax, "Receptor-level context signal / replica noise")
ordered_r = receptor.sort_values("mean_signal_to_noise")
y = np.arange(len(ordered_r))
colors = np.where(ordered_r["context_exceeds_replica_noise"], COLORS["Receptor identity"], COLORS["Non-significant"])
ax.barh(y, ordered_r["mean_signal_to_noise"], color=colors)
ax.axvline(1, color="#777777", linestyle="--", linewidth=0.7)
ax.set_yticks(y, ordered_r["receptor"], fontsize=5)
ax.set_xlabel("Mean signal / noise")
values["panels"]["C"] = ordered_r.to_dict(orient="records")

ax = axes[1, 1]
panel_label(ax, "D")
panel_title(ax, "Context matching audit")
match_columns = [("ligand_matched", "Ligand"), ("binder_matched", "Binder"),
                  ("receptor_construct_description_matched", "Receptor construct"),
                  ("all_systems_label_supported", "All labels supported"),
                  ("all_context_covariates_matched", "All covariates")]
matched = np.array([contrasts[col].fillna(False).astype(bool).sum() for col, _ in match_columns])
total = len(contrasts)
x = np.arange(len(match_columns))
ax.bar(x, matched, color=COLORS["Gi"], label="Matched")
ax.bar(x, total - matched, bottom=matched, color="#E5E5E5", label="Not matched")
ax.set_xticks(x, [label for _, label in match_columns], rotation=35, ha="right")
ax.set_ylabel("Contrasts")
ax.legend(frameon=False)
ax.text(0.98, 0.98, "Context-associated, not causal partner swaps", transform=ax.transAxes, ha="right", va="top", fontsize=5, fontweight="bold")
values["panels"]["D"] = {label: {"matched": int(count), "total": total} for (_, label), count in zip(match_columns, matched)}
save_supplementary(fig, "FigS2_within_receptor_contrasts", values)
print("Figure S2 done.")
