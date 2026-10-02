#!/usr/bin/env python3
"""Cohort-wide equilibration assessment for CoupledMD (Reviewer 2, Fig. 5C).

Reuses the manuscript's structural-validation machinery
(build_core_interface_validation.build_tasks in the original scidata_figures
working directory) so the TM-core selection, trajectory routing (including
corrected/repaired replicas) and fitting protocol match the published metrics
exactly. For each of the 621 replica tasks it computes the full TM-core Calpha
RMSD time series (1,001 evenly spaced frames, least-squares fit on the replica
starting frame) and applies an objective windowed-plateau stationarity test.

Stationarity criterion (operational):
  series split into non-overlapping 10-ns windows; t_eq = start of the
  earliest window W such that (a) every window mean from W onward lies within
  +/-0.5 A of the mean of the final 100 ns and (b) the least-squares slope of
  window means from W onward is <= 0.1 A per 100 ns. Replicas with no such
  window are flagged nonstationary. The Galpha-interface RMSD series is
  assessed in parallel at a relaxed 1.0 A tolerance as a secondary indicator.
"""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

SCIDATA = Path("/MDdata/data02/jxhuang/gpcr_g/gpcr_g_server/manuscript_figures/scidata_figures")
sys.path.insert(0, str(SCIDATA))

import MDAnalysis as mda

from build_core_interface_validation import (  # noqa: E402
    N_SAMPLES, build_tasks, ca_atom_indices, kabsch, parse_ordinals,
    sample_indices,
)

OUTDIR = Path(__file__).resolve().parent / "equilibration_out"
PLATEAU_TOL_A = 0.5
SLOPE_TOL_A_PER_100NS = 0.1
IFACE_PLATEAU_TOL_A = 1.0
WINDOW_NS = 10.0
TOTAL_NS = 500.0


def equilibration_time(t, values, plateau_tol, slope_tol, window_ns=WINDOW_NS):
    edges = np.arange(0.0, t[-1] + 1e-9, window_ns)
    wmean, wstart = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (t >= lo) & (t < hi)
        if m.any():
            wmean.append(values[m].mean())
            wstart.append(lo)
    wmean = np.asarray(wmean)
    wstart = np.asarray(wstart)
    final_ref = wmean[wstart >= t[-1] - 100.0].mean()
    for k in range(len(wmean)):
        seg = wmean[k:]
        if np.all(np.abs(seg - final_ref) <= plateau_tol):
            slope = np.polyfit(wstart[k:], seg, 1)[0] * 100.0
            if abs(slope) <= slope_tol:
                return float(wstart[k]), float(final_ref), False
    return float("nan"), float(final_ref), True


def process_replica_series(task):
    sid, rep = task["system_id"], int(task["replica"])
    if task["selection_status"] != "available":
        return dict(system_id=sid, replica=rep, status="no_selection",
                    reason=task.get("selection_reason", ""))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(task["topology_path"], task["trajectory_path"])
        tm_idx = ca_atom_indices(u, parse_ordinals(task["tm_ordinals"]))
        ga_idx = ca_atom_indices(u, parse_ordinals(task["galpha_ordinals"]))
        rec_idx = ca_atom_indices(u, parse_ordinals(task["receptor_ordinals"]))
        indices = sample_indices(len(u.trajectory))

        u.trajectory[indices[0]]
        ref_tm = u.atoms[tm_idx].positions.astype(float).copy()
        rec0 = u.atoms[rec_idx].positions.astype(float)
        ga0 = u.atoms[ga_idx].positions.astype(float)
        from scipy.spatial.distance import cdist
        d0 = cdist(rec0, ga0)
        iface_local = np.where((d0 <= 12.0).any(axis=0))[0]
        iface_idx = ga_idx[iface_local]
        ref_iface = u.atoms[iface_idx].positions.astype(float).copy()

        tm = np.empty(N_SAMPLES)
        iface = np.empty(N_SAMPLES)
        for j, frame in enumerate(indices):
            u.trajectory[int(frame)]
            mob = u.atoms[tm_idx].positions.astype(float)
            rot, mc, rc = kabsch(mob, ref_tm)
            fitted = (mob - mc) @ rot + rc
            tm[j] = np.sqrt(np.mean(np.sum((fitted - ref_tm) ** 2, axis=1)))
            fi = (u.atoms[iface_idx].positions.astype(float) - mc) @ rot + rc
            iface[j] = np.sqrt(np.mean(np.sum((fi - ref_iface) ** 2, axis=1)))

        t = np.linspace(0.0, TOTAL_NS, N_SAMPLES)
        teq, plateau, nonstat = equilibration_time(t, tm, PLATEAU_TOL_A, SLOPE_TOL_A_PER_100NS)
        teq_i, plateau_i, nonstat_i = equilibration_time(
            t, iface, IFACE_PLATEAU_TOL_A, SLOPE_TOL_A_PER_100NS)
        np.save(OUTDIR / "timeseries" / f"{sid}_rep{rep}_tm_core_rmsd.npy",
                np.column_stack([t, tm, iface]))
        return dict(system_id=sid, replica=rep, status="ok",
                    g_protein_family=task["g_protein_family"],
                    t_eq_ns=teq, plateau_rmsd_A=plateau, nonstationary=nonstat,
                    iface_t_eq_ns=teq_i, iface_plateau_rmsd_A=plateau_i,
                    iface_nonstationary=nonstat_i,
                    rmsd_final100ns_mean_A=float(tm[t >= TOTAL_NS - 100].mean()),
                    n_tm_core_ca=int(len(tm_idx)))
    except Exception as exc:
        return dict(system_id=sid, replica=rep, status="error",
                    reason=f"{type(exc).__name__}: {exc}")


def main(workers=32):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "timeseries").mkdir(exist_ok=True)
    tasks = build_tasks()
    assert len(tasks) == 621
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(process_replica_series, t): t for t in tasks}
        done = 0
        for fut in as_completed(futs):
            rows.append(fut.result())
            done += 1
            if done % 25 == 0:
                print(f"{done}/621", flush=True)

    per_rep = pd.DataFrame(rows).sort_values(["system_id", "replica"])
    per_rep.to_csv(OUTDIR / "equilibration_per_replica.csv", index=False)

    ok = per_rep[per_rep.status.eq("ok")].copy()
    ok["t_eq_ns"] = pd.to_numeric(ok["t_eq_ns"])
    per_sys = (ok.groupby(["system_id", "g_protein_family"])
               .agg(t_eq_rep_max_ns=("t_eq_ns", "max"),
                    t_eq_rep_median_ns=("t_eq_ns", "median"),
                    n_nonstationary=("nonstationary", "sum"),
                    plateau_rmsd_mean_A=("plateau_rmsd_A", "mean"))
               .reset_index()
               .rename(columns={"t_eq_rep_max_ns": "recommended_discard_ns"}))
    per_sys.to_csv(OUTDIR / "equilibration_per_system.csv", index=False)

    teq = ok["t_eq_ns"].dropna()
    summary = dict(
        n_replica_tasks=len(per_rep),
        n_evaluated=int(ok.shape[0]),
        n_no_selection=int(per_rep.status.eq("no_selection").sum()),
        n_errors=int(per_rep.status.eq("error").sum()),
        n_nonstationary=int(ok["nonstationary"].sum()),
        n_iface_nonstationary=int(ok["iface_nonstationary"].sum()),
        t_eq_median_ns=float(teq.median()),
        t_eq_mean_ns=float(teq.mean()),
        t_eq_p75_ns=float(teq.quantile(0.75)),
        t_eq_p90_ns=float(teq.quantile(0.9)),
        t_eq_max_ns=float(teq.max()),
        frac_t_eq_le_50ns=float((teq <= 50).mean()),
        frac_t_eq_le_100ns=float((teq <= 100).mean()),
        frac_t_eq_le_250ns=float((teq <= 250).mean()),
        plateau_rmsd_median_A=float(pd.to_numeric(ok["plateau_rmsd_A"]).median()),
    )
    (OUTDIR / "equilibration_cohort_summary.json").write_text(
        json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    errs = per_rep[per_rep.status.eq("error")]
    if len(errs):
        print(errs[["system_id", "replica", "reason"]].to_string(index=False))


if __name__ == "__main__":
    main(workers=int(sys.argv[1]) if len(sys.argv) > 1 else 32)
