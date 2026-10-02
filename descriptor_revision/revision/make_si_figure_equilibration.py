#!/usr/bin/env python3
"""Merged SI equilibration figure (P1 of k3max_new_review.md).

Three panels, all denominators explicit (fixes E6):
  A  Per-replica t_eq histogram stacked by G-protein family; the 29
     nonstationary replicas shown as a separate grey category after an
     axis break (not plotted as a fake t_eq ~ 490 ns bin).
  B  Cumulative fraction: solid = of the 589 stationary replicas;
     dashed = of all 618 evaluated records (nonstationary counted as
     never plateauing; asymptote 589/618 = 95.3%).
  C  Per-system recommended discard (max t_eq over stationary replicas),
     n = 206 systems, median/P90 marked.

Reads equilibration_out/equilibration_per_replica.csv and
equilibration_per_system.csv. Writes Figure_SI_equilibration_merged.{png,pdf}
and prints the caption numbers.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
OUT = BASE / "equilibration_out"
FIG_BASE = "Figure_SI_equilibration_merged"

FAM_ORDER = ["Gs", "Gi/o", "Gq/11", "G12/13"]
FAM_MAP = {"Gs": "Gs", "Gi": "Gi/o", "Gq": "Gq/11", "G12-13": "G12/13"}
COLORS = {"Gs": "#4C9BD6", "Gi/o": "#E07B39", "Gq/11": "#54A24B", "G12/13": "#B279A2"}


def main():
    pr = pd.read_csv(OUT / "equilibration_per_replica.csv")
    ps = pd.read_csv(OUT / "equilibration_per_system.csv")

    ok = pr[pr.status.eq("ok")].copy()
    ok["t_eq_ns"] = pd.to_numeric(ok["t_eq_ns"], errors="coerce")
    ok["family"] = ok["g_protein_family"].map(FAM_MAP)
    stat = ok.dropna(subset=["t_eq_ns"])
    nonstat = ok[ok["t_eq_ns"].isna()]
    n_eval, n_stat, n_nonstat = len(ok), len(stat), len(nonstat)
    assert n_eval == 618 and n_nonstat == 29, (n_eval, n_nonstat)

    teq = stat["t_eq_ns"].values
    med, p90 = np.median(teq), np.quantile(teq, 0.9)
    f50s = (teq <= 50).mean()                       # of stationary
    f50a = (teq <= 50).sum() / n_eval               # of all 618
    f250s = (teq <= 250).mean()
    f250a = (teq <= 250).sum() / n_eval

    disc = ps["recommended_discard_ns"].dropna().values
    dmed, dp90 = np.median(disc), np.quantile(disc, 0.9)
    n_sys_flag = int((ps["n_nonstationary"] > 0).sum())

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13.8, 4.0))

    # ---- Panel A: stacked histogram + separated nonstationary category ----
    bins = np.arange(0, 510, 10)
    bottom = np.zeros(len(bins) - 1)
    for fam in FAM_ORDER:
        vals = stat.loc[stat.family.eq(fam), "t_eq_ns"]
        counts, _ = np.histogram(vals, bins=bins)
        ax1.bar(bins[:-1] + 5, counts, width=9, bottom=bottom,
                color=COLORS[fam], alpha=0.85, edgecolor="white", linewidth=0.4,
                label=fam)
        bottom += counts
    # nonstationary bar after an axis break at x = 545
    ns_x = 545
    ax1.bar(ns_x, n_nonstat, width=20, color="#7f7f7f", alpha=0.85,
            edgecolor="white", linewidth=0.4, label="nonstationary")
    ax1.set_xlim(-12, 580)
    ax1.set_xticks(list(range(0, 501, 100)) + [ns_x])
    ax1.set_xticklabels([str(v) for v in range(0, 501, 100)] + ["non-\nstationary"],
                        fontsize=8)
    # break marks between 500 and the nonstationary category
    for xb in (512, 522):
        ax1.plot([xb - 4, xb + 4], [-9, 7], color="k", lw=1.0, clip_on=False)
    ax1.axvline(med, color="k", ls="-", lw=1.2)
    ax1.text(med + 8, ax1.get_ylim()[1] * 0.80, f"median = {med:.0f} ns",
             fontsize=8, rotation=0)
    ax1.set_xlabel(r"Per-replica equilibration time $t_{eq}$ (ns)", fontsize=10)
    ax1.set_ylabel("Number of replica records", fontsize=10)
    ax1.set_title(f"A  Per-replica $t_{{eq}}$ ({n_stat} stationary of {n_eval} evaluated)",
                  fontsize=10, loc="left")
    ax1.legend(fontsize=8, frameon=False)

    # ---- Panel B: cumulative fractions, both denominators ----
    xs = np.sort(teq)
    ax2.plot(xs, np.arange(1, len(xs) + 1) / len(xs), color="k", lw=1.6,
             label=f"of {n_stat} stationary replicas")
    ax2.plot(xs, np.arange(1, len(xs) + 1) / n_eval, color="k", lw=1.6, ls="--",
             label=f"of all {n_eval} evaluated records")
    ax2.axhline(n_stat / n_eval, color="grey", lw=0.8, ls=":")
    ax2.text(12, 0.90,
             f"dashed curve asymptotes at {n_stat}/{n_eval} = {100*n_stat/n_eval:.1f}%\n"
             f"({n_nonstat} nonstationary never plateau)",
             fontsize=7.5, ha="left", va="top", color="#444444")
    for v, lab, ls in [(50, f"{100*f50s:.1f}% ≤ 50 ns", ":"),
                       (250, f"{100*f250s:.1f}% ≤ 250 ns", ":")]:
        ax2.axvline(v, color="grey", ls=ls, lw=1.0)
        ax2.text(v + 6, 0.06 if v == 50 else 0.30, lab, fontsize=8)
    ax2.set_xlim(0, 500)
    ax2.set_ylim(0, 1.02)
    ax2.set_xlabel(r"$t_{eq}$ (ns)", fontsize=10)
    ax2.set_ylabel("Cumulative fraction", fontsize=10)
    ax2.set_title("B  Cohort cumulative distribution", fontsize=10, loc="left")
    ax2.legend(fontsize=8, frameon=False, loc="lower right")

    # ---- Panel C: per-system recommended discard ----
    ax3.hist(disc, bins=bins, color="#E07B39", alpha=0.75, edgecolor="white",
             linewidth=0.4)
    for v, lab, ls, ha, dx in [(dmed, f"median = {dmed:.0f} ns", "-", "right", -8),
                               (dp90, f"P90 = {dp90:.0f} ns", "--", "left", 8)]:
        ax3.axvline(v, color="k", ls=ls, lw=1.2)
        ax3.text(v + dx, ax3.get_ylim()[1] * 0.97, lab, fontsize=8, va="top", ha=ha)
    ax3.set_xlim(0, 500)
    ax3.set_xlabel("Recommended discard per system (ns)", fontsize=10)
    ax3.set_ylabel("Number of systems", fontsize=10)
    ax3.set_title(f"C  Recommended exclusion window ({len(disc)} systems)",
                  fontsize=10, loc="left")
    ax3.text(0.97, 0.80,
             f"{n_sys_flag} systems contain ≥1 nonstationary replica;\n"
             "discard computed over stationary replicas only (S7)",
             transform=ax3.transAxes, fontsize=7.5, ha="right", va="top",
             color="#7f7f7f")

    for ax in (ax1, ax2, ax3):
        ax.tick_params(labelsize=9)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{FIG_BASE}.{ext}", dpi=300)
    print(f"Wrote {FIG_BASE}.png / .pdf")

    print("\n=== Caption numbers ===")
    print(f"evaluated {n_eval}, stationary {n_stat}, nonstationary {n_nonstat}")
    print(f"t_eq median {med:.0f} ns, P90 {p90:.0f} ns")
    print(f"<=50 ns : {100*f50s:.1f}% of stationary, {100*f50a:.1f}% of all evaluated")
    print(f"<=250 ns: {100*f250s:.1f}% of stationary, {100*f250a:.1f}% of all evaluated")
    print(f"discard median {dmed:.0f} ns, P90 {dp90:.0f} ns; systems with flag: {n_sys_flag}")


if __name__ == "__main__":
    main()
