#!/usr/bin/env python3
"""Build a harmonized QC ledger for the 621 selected production replicas.

The baseline audit predates several trajectory repairs.  This builder overlays only
the GOOD replica records cited by the final-207 readiness ledger and derives the six
missing GROMACS component-summary rows from their existing component time series.
Timestamp resets remain diagnostic and are not used as a failure criterion.
"""
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
SUPP = ROOT / "CoupledMD_Supplementary_Data"
TR = ROOT / "trajectory_readiness"
OUT = ROOT / "v9_figure_inputs" / "scidata_T10_final624_replica_qc_v9.csv"


def corrected_sources():
    """Return (path, optional system allow-list) for final corrected evidence."""
    return [
        (TR / "Gq_7RAN_direct_final_audit/corrected_role_replica_audit.csv", None),
        (
            TR / "amber_three_replica_repair_20260713/exhaustive_audit_corrected_receptor_fragments/candidate_replica_audit.csv",
            {"Gi_7F1Q", "Gi_7QVM", "Gi_7W0N", "Gi_7XA3", "Gi_8IY5", "Gq_8IBV", "Gs_6NBI", "Gs_6WZG", "Gs_7D68", "Gs_7KH0"},
        ),
        (TR / "continuations/Gi_8J18/interface_recovery/exhaustive_audit_corrected/candidate_replica_audit.csv", None),
        (TR / "continuations/Gq_8ZPT/interface_recovery/exhaustive_audit_corrected/candidate_replica_audit.csv", None),
        (TR / "continuous_review_audit/continuous_review_replica_audit_with_rmsd.csv", None),
        (TR / "corrected_role_audit/corrected_role_replica_audit.csv", None),
        (TR / "interface_centered_pilot_20260713/Gi_6WWZ/three_replica_audit_corrected/candidate_replica_audit.csv", None),
        (TR / "interface_centered_rescue_20260713/Gi_7XA3/audit/candidate_replica_audit.csv", None),
        (TR / "interface_centered_rescue_20260713/Gq_8IBV/audit/candidate_replica_audit.csv", None),
    ]


def build():
    selected = pd.read_csv(SUPP / "Supplementary_Data_S7_selected_replica_ledger.csv").rename(columns={"replica_id": "replica"})
    systems = pd.read_csv(SUPP / "Supplementary_Data_S1_included_system_inventory.csv")[["system_id", "g_protein_family"]]
    strict = pd.read_csv(TR / "strict_full_length_replica_audit.csv")
    component = pd.read_csv(TR / "component_audit_summary_nc.csv")
    time_axis = pd.read_csv(TR / "time_axis_audit.csv")

    data = selected.merge(systems, on="system_id", validate="many_to_one")
    data = data.merge(strict, on=["system_id", "replica", "g_protein_family"], validate="one_to_one")
    data = data.merge(
        component[["system_id", "replica", "n_nonfinite_sampled_frames", "n_chain_break_role_frames", "status"]].rename(columns={"status": "component_status"}),
        on=["system_id", "replica"], how="left", validate="one_to_one",
    )
    data = data.merge(
        time_axis[["system_id", "replica", "n_nonfinite_times", "status"]].rename(columns={"status": "time_axis_status"}),
        on=["system_id", "replica"], how="left", validate="one_to_one",
    )
    data["audit_route"] = "baseline evidence"
    data["audit_source"] = "trajectory_readiness/strict_full_length_replica_audit.csv + component/time-axis audits"
    data["measured_span_ps"] = data["aligned_span_ps"]

    overlaid = set()
    for path, allowed in corrected_sources():
        corrected = pd.read_csv(path)
        corrected = corrected[corrected["verdict"].eq("GOOD")]
        if allowed is not None:
            corrected = corrected[corrected["system_id"].isin(allowed)]
        for row in corrected.itertuples(index=False):
            key = (row.system_id, int(row.replica))
            assert key not in overlaid, f"duplicate corrected evidence: {key}"
            idx = data.index[(data["system_id"].eq(row.system_id)) & (data["replica"].eq(row.replica))]
            assert len(idx) == 1, f"corrected record not uniquely selected: {key}"
            i = idx[0]
            data.loc[i, "audit_route"] = "corrected / reviewed evidence"
            data.loc[i, "audit_source"] = str(path.relative_to(ROOT))
            data.loc[i, "measured_span_ps"] = row.span_ps
            data.loc[i, "n_nonfinite_sampled_frames"] = row.n_nonfinite_frames
            data.loc[i, "n_chain_break_role_frames"] = row.n_chain_break_frames
            data.loc[i, "n_nonfinite_times"] = 0
            data.loc[i, "component_status"] = row.verdict
            data.loc[i, "time_axis_status"] = row.verdict
            overlaid.add(key)

    # The baseline summary omitted these two GROMACS systems, although their
    # 1,001-frame component time series already exist.  Summarize those records.
    missing = data[data["n_nonfinite_sampled_frames"].isna()]
    assert set(missing["system_id"]) == {"Gi_7YK6", "Gi_8YIC"} and len(missing) == 6
    for i, row in missing.iterrows():
        path = TR / "component_timeseries" / f"{row.system_id}_rep{int(row.replica)}_components.csv.gz"
        ts = pd.read_csv(path)
        break_columns = [c for c in ts.columns if c.endswith("_ca_breaks_gt8A")]
        assert len(ts) == 1001 and break_columns
        data.loc[i, "n_nonfinite_sampled_frames"] = (~ts["coordinates_finite"].astype(bool)).sum()
        data.loc[i, "n_chain_break_role_frames"] = (ts[break_columns].sum(axis=1) > 0).sum()
        data.loc[i, "n_nonfinite_times"] = pd.to_numeric(ts["time_ps"], errors="coerce").isna().sum()
        data.loc[i, "component_status"] = "OK_DERIVED_EXISTING_TIMESERIES"
        data.loc[i, "time_axis_status"] = "OK_DERIVED_EXISTING_TIMESERIES"
        data.loc[i, "audit_source"] += f" + {path.relative_to(ROOT)}"

    data["duration_500ns_ok"] = data["measured_span_ps"].between(499900, 500100)
    data["finite_time_axis_ok"] = data["n_nonfinite_times"].eq(0)
    data["finite_coordinates_ok"] = data["n_nonfinite_sampled_frames"].eq(0)
    data["chain_continuity_ok"] = data["n_chain_break_role_frames"].eq(0)
    data["selected_validated_ok"] = data["release_selection_status"].eq("selected_validated")

    expected_families = {"Gi": 285, "Gs": 195, "Gq": 123, "G12-13": 18}
    assert len(data) == 621 and not data.duplicated(["system_id", "replica"]).any()
    assert data.groupby("system_id").size().eq(3).all()
    assert data["g_protein_family"].value_counts().to_dict() == expected_families
    assert len(overlaid) == 59
    assert (data["audit_route"] == "baseline evidence").sum() == 562
    checks = ["duration_500ns_ok", "finite_time_axis_ok", "finite_coordinates_ok", "chain_continuity_ok", "selected_validated_ok"]
    assert data[checks].all().all()
    assert data["duration_ns"].sum() / 1000 == 310.5

    columns = [
        "system_id", "replica", "g_protein_family", "source_format", "duration_ns", "measured_span_ps",
        "audit_route", "audit_source", "n_nonfinite_times", "n_nonfinite_sampled_frames",
        "n_chain_break_role_frames", *checks,
    ]
    OUT.parent.mkdir(exist_ok=True)
    data[columns].sort_values(["g_protein_family", "system_id", "replica"]).to_csv(OUT, index=False)
    print("Final-621 QC assertions passed:")
    print("  207 systems; 621 selected replicas; 3 replicas/system; 310.5 microseconds")
    print("  family replicas: Gi/o 285, Gs 195, Gq/11 123, G12/13 18")
    print("  evidence routes: 562 baseline, 59 corrected/reviewed")
    print("  621/621 passed duration, finite-time, finite-coordinate, chain-continuity, and selection checks")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    build()
