#!/usr/bin/env python3
"""Calculate harmonized quantitative validation metrics for 621 replicas.

Each final selected trajectory contributes exactly 1,001 nearest observations at
evenly spaced relative times across its measured production span.  Per-replica
summaries are therefore comparable between the 1,001-point baseline audits and
the 10,000-frame corrected/reviewed audits.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from build_final624_replica_qc_v9 import corrected_sources


ROOT = Path(__file__).resolve().parent
TR = ROOT / "trajectory_readiness"
QC = ROOT / "v9_figure_inputs/scidata_T10_final624_replica_qc_v9.csv"
OUT = ROOT / "v9_figure_inputs/scidata_T11_final624_harmonized_validation_v9.csv"
N_SAMPLES = 1001


def nearest_uniform_indices(time_ps):
    time = np.asarray(time_ps, dtype=float)
    assert len(time) >= N_SAMPLES and np.isfinite(time).all()
    order = np.argsort(time, kind="stable")
    sorted_time = time[order]
    targets = np.linspace(sorted_time[0], sorted_time[-1], N_SAMPLES)
    right = np.searchsorted(sorted_time, targets, side="left")
    right = np.clip(right, 0, len(sorted_time) - 1)
    left = np.clip(right - 1, 0, len(sorted_time) - 1)
    choose_left = np.abs(sorted_time[left] - targets) <= np.abs(sorted_time[right] - targets)
    nearest = np.where(choose_left, left, right)
    assert len(np.unique(nearest)) == N_SAMPLES
    return order[nearest]


def path_maps():
    component = pd.read_csv(TR / "component_audit_summary_nc.csv")
    baseline = {
        (r.system_id, int(r.replica)): Path(r.timeseries_path)
        for r in component.itertuples(index=False) if isinstance(r.timeseries_path, str)
    }
    corrected = {}
    for path, allowed in corrected_sources():
        audit = pd.read_csv(path)
        audit = audit[audit["verdict"].eq("GOOD")]
        if allowed is not None:
            audit = audit[audit["system_id"].isin(allowed)]
        for row in audit.itertuples(index=False):
            key = (row.system_id, int(row.replica))
            assert key not in corrected
            corrected[key] = Path(row.timeseries_path)
    assert len(corrected) == 59
    return baseline, corrected


def summarize(values):
    values = np.asarray(values, dtype=float)
    assert len(values) == N_SAMPLES and np.isfinite(values).all()
    q25, median, q75, p95 = np.quantile(values, [.25, .5, .75, .95])
    return median, q25, q75, p95


def build():
    qc = pd.read_csv(QC)
    assert len(qc) == 621 and qc.groupby("system_id").size().eq(3).all()
    baseline, corrected = path_maps()
    rows = []
    for n, row in enumerate(qc.itertuples(index=False), 1):
        key = (row.system_id, int(row.replica))
        if row.audit_route == "corrected / reviewed evidence":
            path = corrected[key]
        else:
            path = baseline.get(key, TR / "component_timeseries" / f"{row.system_id}_rep{int(row.replica)}_components.csv.gz")
        assert path.exists(), f"missing component time series: {key}: {path}"
        ts = pd.read_csv(path)
        galpha_column = "g_alpha_ca_rmsd_A" if "g_alpha_ca_rmsd_A" in ts else "galpha_ca_rmsd_A"
        required = ["time_ps", "coordinates_finite", "receptor_ca_rmsd_A", galpha_column, "receptor_galpha_min_ca_A"]
        assert all(c in ts for c in required), f"missing metric column: {key}"
        assert ts["coordinates_finite"].astype(bool).all(), f"nonfinite coordinate flag: {key}"
        idx = nearest_uniform_indices(pd.to_numeric(ts["time_ps"], errors="coerce"))
        sample = ts.iloc[idx]
        receptor = summarize(sample["receptor_ca_rmsd_A"])
        galpha = summarize(sample[galpha_column])
        interface = summarize(sample["receptor_galpha_min_ca_A"])
        rows.append({
            "system_id": row.system_id,
            "replica": int(row.replica),
            "g_protein_family": row.g_protein_family,
            "audit_route": row.audit_route,
            "timeseries_path": str(path),
            "source_observations": len(ts),
            "harmonized_observations": N_SAMPLES,
            "relative_start_ns": 0.0,
            "relative_end_ns": 500.0,
            "receptor_rmsd_median_A": receptor[0],
            "receptor_rmsd_q25_A": receptor[1],
            "receptor_rmsd_q75_A": receptor[2],
            "receptor_rmsd_p95_A": receptor[3],
            "galpha_rmsd_median_A": galpha[0],
            "galpha_rmsd_q25_A": galpha[1],
            "galpha_rmsd_q75_A": galpha[2],
            "galpha_rmsd_p95_A": galpha[3],
            "interface_min_ca_median_A": interface[0],
            "interface_min_ca_q25_A": interface[1],
            "interface_min_ca_q75_A": interface[2],
            "interface_min_ca_p95_A": interface[3],
            "interface_contact_fraction_le8A": (sample["receptor_galpha_min_ca_A"].to_numpy(float) <= 8.0).mean(),
        })
        if n % 100 == 0:
            print(f"  summarized {n}/621 replicas")
    result = pd.DataFrame(rows)
    assert len(result) == 621 and not result.duplicated(["system_id", "replica"]).any()
    assert result["harmonized_observations"].eq(N_SAMPLES).all()
    assert result.select_dtypes(include="number").notna().all().all()
    assert result["interface_contact_fraction_le8A"].between(0, 1).all()
    OUT.parent.mkdir(exist_ok=True)
    result.to_csv(OUT, index=False)
    print("Harmonized-validation assertions passed:")
    print("  621 replicas; exactly 1,001 evenly spaced observations per replica")
    print("  finite receptor RMSD, G-alpha RMSD, and receptor--G-alpha distance summaries")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    build()
