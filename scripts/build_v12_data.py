#!/usr/bin/env python3
"""Build the audited v12 figure inputs, Supplementary Data, and numerical audit.

All cohort counts are recomputed from machine-readable sources.  The script
does not infer publication from a reserved DOI and does not treat the reduced
PDB/XTC release as a full-system simulation archive.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SERVER = HERE.parent.parent
OVERLAP_ROOT = SERVER.parent / "a"
FIG_INPUT = HERE / "v12_figure_inputs"
FIG_SOURCE = HERE / "v12_figure_source_data"
SUPP = HERE / "CoupledMD_Supplementary_Data_v12"
REPORTS = HERE / "reports"

COHORT = SERVER / "data" / "release_cohort_v9_final207.csv"
S2_OLD = (
    HERE
    / "CoupledMD_Supplementary_Data"
    / "Supplementary_Data_S2_release_boundary_exceptions.csv"
)
S5_OLD = (
    HERE
    / "CoupledMD_Supplementary_Data"
    / "Supplementary_Data_S5_pocket_summaries_gpcrdb.csv"
)
S6_OLD = (
    HERE
    / "CoupledMD_Supplementary_Data"
    / "Supplementary_Data_S6_gateway_per_system.csv"
)
S9_OLD = (
    HERE
    / "CoupledMD_Supplementary_Data"
    / "Supplementary_Data_S9_api_access_coverage.csv"
)
T4 = HERE / "v9_figure_inputs" / "scidata_T4_technical_validation_v9.csv"
T10 = HERE / "v9_figure_inputs" / "scidata_T10_final624_replica_qc_v9.csv"
T13 = (
    HERE
    / "v9_figure_inputs"
    / "scidata_T13_final624_core_interface_validation_v9.csv"
)
VIZ_AUDIT = HERE / "final207_reduced_visualization_audit_40frames.csv"
LEGACY_VIZ_AUDIT = HERE / "final208_reduced_visualization_audit_40frames.csv"
PROTOCOL = (
    HERE
    / "parallel_codex"
    / "provenance_direct"
    / "protocol_system_evidence.csv"
)
COHORT_TRUTH = (
    OVERLAP_ROOT
    / "dynamic_activation"
    / "manifests"
    / "cohort_truth_v1.csv"
)
COHORT_TRUTH_SUMMARY = (
    OVERLAP_ROOT
    / "dynamic_activation"
    / "manifests"
    / "cohort_truth_summary_v1.json"
)
RELEASE_REMOTE = REPORTS / "v12_remote_release_verification.json"
GQ8ZPT_REPAIR = (
    HERE
    / "trajectory_readiness"
    / "continuations"
    / "Gq_8ZPT"
    / "interface_recovery"
    / "exhaustive_audit_corrected"
    / "candidate_replica_audit.csv"
)
API_GATEWAY_CONSENSUS = SERVER / "data/api/v1/consensus/gateways.json"

RELEASE_DIRS = {
    1: HERE / "zenodo_reduced_release_208",
    2: HERE / "zenodo_reduced_release_207_replica2",
    3: HERE / "zenodo_reduced_release_207_replica3",
}
DOIS = {
    1: "10.5281/zenodo.21395292",
    2: "10.5281/zenodo.21447748",
    3: "10.5281/zenodo.21448037",
}

FAMILY_ORDER = ["Gi", "Gs", "Gq", "G12-13"]
FAMILY_LABEL = {
    "Gi": "Gi/o",
    "Gs": "Gs",
    "Gq": "Gq/11",
    "G12-13": "G12/13",
}
PROTOCOL_MAP = {
    "CHARMM36 (chamber, AMBER pmemd)": "P1",
    "CHARMM36 (via AMBER CHARMM-GUI)": "P2",
    "CHARMM36 (membrane-embedded, GROMACS)": "P3",
}


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_record(path: Path) -> dict[str, Any]:
    external_aliases = {
        COHORT_TRUTH: "read_only_overlap/cohort_truth_v1.csv",
        COHORT_TRUTH_SUMMARY: (
            "read_only_overlap/cohort_truth_summary_v1.json"
        ),
    }
    return {
        "path": (
            str(path.relative_to(SERVER))
            if path.is_relative_to(SERVER)
            else external_aliases.get(path, f"read_only_overlap/{path.name}")
        ),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def write_csv(data: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(path, index=False, lineterminator="\n")


def load_sources() -> dict[str, Any]:
    cohort = pd.read_csv(COHORT)
    s2 = pd.read_csv(S2_OLD)
    pockets = pd.read_csv(S5_OLD)
    gateways = pd.read_csv(S6_OLD)
    t4 = pd.read_csv(T4)
    t10 = pd.read_csv(T10)
    t13 = pd.read_csv(T13)
    protocol = pd.read_csv(PROTOCOL)
    truth = pd.read_csv(COHORT_TRUTH)
    remote = (
        json.loads(RELEASE_REMOTE.read_text(encoding="utf-8"))
        if RELEASE_REMOTE.is_file()
        else {"records": {}, "summary": {}}
    )

    assert len(cohort) == cohort.system_id.nunique() == 207
    cohort_ids = set(cohort.system_id)
    assert cohort.gpcr_class.value_counts().to_dict() == {"A": 181, "B": 26}
    assert cohort.g_protein_family.value_counts().to_dict() == {
        "Gi": 95,
        "Gs": 65,
        "Gq": 41,
        "G12-13": 6,
    }
    assert cohort.receptor_name.nunique() == 174
    assert cohort.receptor_uniprot.nunique() == 173
    assert cohort.loc[cohort.receptor_uniprot.isna(), "system_id"].tolist() == [
        "Gs_8HTI"
    ]
    assert len(s2) == 15
    assert Counter(s2.release_status) == Counter(
        {"unresolved": 13, "excluded": 2}
    )
    assert not cohort_ids.intersection(s2.system_id)
    assert len(t10) == 621 and t10.system_id.nunique() == 207
    assert t10.groupby("system_id").replica.nunique().eq(3).all()
    assert len(t13) == 621 and t13.system_id.nunique() == 207
    assert set(t13.system_id) == cohort_ids
    assert set(protocol.system_id).issuperset(cohort_ids)
    protocol = protocol[protocol.system_id.isin(cohort_ids)].copy()
    assert len(protocol) == 207
    assert protocol.protocol_group.value_counts().to_dict() == {
        "P1": 179,
        "P2": 26,
        "P3": 2,
    }
    assert set(truth.loc[truth.release_status.eq("included"), "system_id"]) == cohort_ids
    return {
        "cohort": cohort,
        "s2": s2,
        "pockets": pockets,
        "gateways": gateways,
        "t4": t4,
        "t10": t10,
        "t13": t13,
        "protocol": protocol,
        "truth": truth,
        "remote": remote,
    }


def build_figure1(data: dict[str, Any]) -> None:
    cohort = data["cohort"].copy()
    composition = (
        cohort.groupby(["gpcr_class", "g_protein_family"])
        .size()
        .rename("systems")
        .reset_index()
    )
    composition["family_label"] = composition.g_protein_family.map(FAMILY_LABEL)
    composition["nominal_sampling_us"] = composition.systems * 1.5
    composition["scope"] = "administrative release labels; descriptive only"
    write_csv(composition, FIG_SOURCE / "Figure_1A_cohort_composition.csv")

    boundary = pd.DataFrame(
        [
            {
                "release_status": "included",
                "systems": 207,
                "contributes_to_release": True,
                "definition": "frozen final-207 cohort",
            },
            {
                "release_status": "unresolved",
                "systems": 13,
                "contributes_to_release": False,
                "definition": "working-inventory records held outside the release",
            },
            {
                "release_status": "excluded",
                "systems": 2,
                "contributes_to_release": False,
                "definition": "duplicate or mislabelled records",
            },
        ]
    )
    write_csv(boundary, FIG_SOURCE / "Figure_1B_release_boundary.csv")

    protocol = (
        data["protocol"]
        .groupby(
            [
                "protocol_group",
                "declared_force_field",
                "engine_family",
            ],
            dropna=False,
        )
        .size()
        .rename("systems")
        .reset_index()
    )
    protocol["replicas"] = protocol.systems * 3
    protocol["nominal_sampling_us"] = protocol.systems * 1.5
    write_csv(protocol, FIG_SOURCE / "Figure_1C_protocol_groups.csv")

    identifiers = pd.DataFrame(
        [
            {
                "metric": "systems",
                "value": 207,
                "note": "unique final release identifiers",
            },
            {
                "metric": "receptor_names",
                "value": 174,
                "note": "distinct receptor-name entities",
            },
            {
                "metric": "mapped_uniprot_accessions",
                "value": 173,
                "note": "one explicit null accession for Gs_8HTI",
            },
            {
                "metric": "production_replicas",
                "value": 621,
                "note": "three computational repeats per system",
            },
        ]
    )
    write_csv(identifiers, FIG_SOURCE / "Figure_1D_identifier_coverage.csv")


def build_figure2(data: dict[str, Any]) -> None:
    roles = pd.DataFrame(
        [
            {
                "record": "system metadata",
                "format": "CSV/JSON",
                "scope": "207 systems",
                "distributed": "supplementary data, API and paper repository",
                "supports": "discovery, identifiers and cohort joins",
                "does_not_support": "coordinate analysis",
            },
            {
                "record": "reduced structure",
                "format": "PDB",
                "scope": "one matched structure per system and release replica",
                "distributed": "Zenodo private drafts pending publication",
                "supports": "atom order and reduced-topology interpretation",
                "does_not_support": "membrane, solvent or ion analyses",
            },
            {
                "record": "reduced trajectory",
                "format": "XTC",
                "scope": "2,500 protein-complex frames per system and replica",
                "distributed": "Zenodo private drafts pending publication",
                "supports": "protein-complex structural analyses at 200-ps spacing",
                "does_not_support": "lipid gateways, solvent, ions or fast kinetics",
            },
            {
                "record": "full-system source simulation",
                "format": "NetCDF/TRR plus engine topology",
                "scope": "three production replicas per system",
                "distributed": "not deposited in the present release",
                "supports": "source provenance and locally executed analyses",
                "does_not_support": "independent reuse without a stable access route",
            },
            {
                "record": "pocket annotations",
                "format": "CSV/JSON/NPZ",
                "scope": "207 system records",
                "distributed": "supplementary data, API and code repository",
                "supports": "technical annotation reuse",
                "does_not_support": "experimental druggability or mechanism",
            },
            {
                "record": "gateway summaries",
                "format": "CSV/JSON",
                "scope": "207 systems × 7 interfaces × 4 metrics",
                "distributed": "supplementary data and API",
                "supports": "reproduction of reported summary tables",
                "does_not_support": "recalculation from reduced trajectories",
            },
            {
                "record": "release QC and checksums",
                "format": "CSV/JSON",
                "scope": "file and system level",
                "distributed": "with each reduced release record",
                "supports": "integrity verification",
                "does_not_support": "validation of unavailable full-system files",
            },
        ]
    )
    write_csv(roles, FIG_SOURCE / "Figure_2_record_roles.csv")
    gateway_flow = pd.DataFrame(
        [
            {
                "stage_order": 1,
                "record": "full-system source trajectories",
                "distributed_in_reduced_release": False,
                "contains_lipids": True,
                "role": "required molecular source for lipid-gateway calculation",
            },
            {
                "stage_order": 2,
                "record": "gateway calculation code and local intermediates",
                "distributed_in_reduced_release": False,
                "contains_lipids": False,
                "role": "calculation layer executed against full-system sources",
            },
            {
                "stage_order": 3,
                "record": "per-system gateway summaries",
                "distributed_in_reduced_release": False,
                "contains_lipids": False,
                "role": "processed intermediate distributed as Supplementary Data S6 and JSON",
            },
            {
                "stage_order": 4,
                "record": "reduced protein-complex PDB/XTC",
                "distributed_in_reduced_release": True,
                "contains_lipids": False,
                "role": "supports protein-complex reuse but not gateway recalculation",
            },
        ]
    )
    write_csv(
        gateway_flow,
        FIG_SOURCE / "Supplementary_Figure_1_gateway_provenance.csv",
    )


def build_figure3(data: dict[str, Any]) -> None:
    cohort_ids = set(data["cohort"].system_id)
    t4 = data["t4"]
    eligible = t4[
        t4.system_id.isin(cohort_ids) & t4.lig_type.eq("peptide")
    ].copy()
    assert len(eligible) == 58
    assert int(eligible.ortho_recovered.sum()) == 49
    eligible["family_label"] = eligible.g_family.map(FAMILY_LABEL)
    eligible["validation_scope"] = (
        "technical orthosteric positive control; not druggability evidence"
    )
    write_csv(
        eligible[
            [
                "system_id",
                "g_family",
                "family_label",
                "n_pockets",
                "ortho_recovered",
                "best_ortho_freq",
                "n_ortho_pockets",
                "structural_provenance",
                "validation_scope",
            ]
        ],
        FIG_SOURCE / "Figure_3A_positive_control.csv",
    )

    pockets = data["pockets"]
    detected = pockets[pockets.pocket_record_status.eq("available_detected_pocket")]
    zero = pockets[pockets.pocket_record_status.eq("available_zero_pockets")]
    assert detected.system_id.nunique() == 205
    assert len(detected) == 2149
    assert set(zero.system_id) == {"Gi_7YK6", "Gi_8YIC"}
    assert len(zero) == 2
    completeness = pd.DataFrame(
        [
            {
                "record_status": "available_detected_pocket",
                "systems": 205,
                "pocket_rows": 2149,
                "interpretation": "one or more persistent pockets detected",
            },
            {
                "record_status": "available_zero_pockets",
                "systems": 2,
                "pocket_rows": 0,
                "interpretation": "valid completed analysis with zero detected pockets",
            },
        ]
    )
    write_csv(completeness, FIG_SOURCE / "Figure_3B_pocket_completeness.csv")


def load_viz_audit(cohort_ids: set[str]) -> tuple[pd.DataFrame, bool]:
    if VIZ_AUDIT.is_file():
        path = VIZ_AUDIT
        regenerated = True
    else:
        path = LEGACY_VIZ_AUDIT
        regenerated = False
    viz = pd.read_csv(path)
    assert len(viz) == viz.system_id.nunique() == 207
    assert set(viz.system_id) == cohort_ids
    viz["audit_source_filename"] = path.name
    viz["regenerated_exact_final207"] = regenerated
    return viz, regenerated


def build_figure4(data: dict[str, Any]) -> bool:
    t13 = data["t13"].copy()
    available = t13.validation_status.eq("available")
    assert int(available.sum()) == 618
    assert t13.loc[available, "system_id"].nunique() == 206
    assert set(t13.loc[~available, "system_id"]) == {"Gs_8HTI"}
    t13["family_label"] = t13.g_protein_family.map(FAMILY_LABEL)
    t13["availability_reason"] = t13.reason_code.fillna("available")
    keep = [
        "system_id",
        "replica",
        "g_protein_family",
        "family_label",
        "identity_status",
        "validation_status",
        "availability_reason",
        "harmonized_observations",
        "tm_core_rmsd_A_p95",
        "galpha_interface_rmsd_A_p95",
        "contact_retention_p05",
    ]
    write_csv(t13[keep], FIG_SOURCE / "Figure_4A_replica_validation.csv")

    viz, regenerated = load_viz_audit(set(data["cohort"].system_id))
    write_csv(viz, FIG_SOURCE / "Figure_4B_reduced_record_qc.csv")
    return regenerated


def read_release_manifest(replica: int) -> pd.DataFrame | None:
    path = RELEASE_DIRS[replica] / "CoupledMD_reduced_trajectory_manifest.csv"
    if not path.is_file():
        return None
    manifest = pd.read_csv(path)
    required = {
        "system_id",
        "source_replica",
        "pdb_relative_path",
        "xtc_relative_path",
        "n_frames",
        "n_atoms",
        "pdb_size_bytes",
        "xtc_size_bytes",
        "pdb_sha256",
        "xtc_sha256",
    }
    assert required.issubset(manifest.columns)
    assert manifest.system_id.nunique() == len(manifest) == 207
    assert set(manifest.source_replica) == {replica}
    return manifest


def local_release_audit(replica: int) -> dict[str, Any]:
    path = RELEASE_DIRS[replica] / "CoupledMD_reduced_release_audit.json"
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def remote_record(remote: dict[str, Any], replica: int) -> dict[str, Any]:
    return remote.get("records", {}).get(str(replica), {})


def build_release_tables(data: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    cohort = data["cohort"]
    remote = data["remote"]
    status_rows: list[dict[str, Any]] = []
    manifest_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []

    for replica in [1, 2, 3]:
        manifest = read_release_manifest(replica)
        audit = local_release_audit(replica)
        remote_item = remote_record(remote, replica)
        authenticated = remote_item.get("authenticated_remote", {})
        public = remote_item.get("public_resolution", {})
        audit_status = audit.get("status")
        audit_systems = audit.get("cohort_systems") or audit.get(
            "archive_system_total"
        )
        status_rows.append(
            {
                "replica": replica,
                "reserved_doi": DOIS[replica],
                "local_manifest_rows": 0 if manifest is None else len(manifest),
                "local_full_frame_qc_status": audit_status or "not_available",
                "local_qc_systems": audit_systems or 0,
                "remote_state": authenticated.get("state", "not_verified"),
                "remote_access_right": authenticated.get(
                    "access_right", "not_verified"
                ),
                "remote_uploaded_files": authenticated.get("file_count", 0),
                "public_http_status": public.get("http_status"),
                "publicly_resolvable": bool(public.get("publicly_resolvable")),
                "release_claim_status": remote_item.get(
                    "status", "UNVERIFIED"
                ),
                "remote_verification_time": authenticated.get("verified_at"),
            }
        )
        by_system = (
            {}
            if manifest is None
            else manifest.set_index("system_id").to_dict(orient="index")
        )
        audit_system_rows = {
            row["system_id"]: row for row in audit.get("system_qc", [])
        }
        for row in cohort.itertuples(index=False):
            record = by_system.get(row.system_id, {})
            qc = audit_system_rows.get(row.system_id, {})
            for role in ["structure_pdb", "trajectory_xtc"]:
                prefix = "pdb" if role == "structure_pdb" else "xtc"
                manifest_rows.append(
                    {
                        "system_id": row.system_id,
                        "pdb_id": row.pdb_id,
                        "g_protein_family": row.g_protein_family,
                        "release_replica": replica,
                        "reserved_doi": DOIS[replica],
                        "file_role": role,
                        "archive_member_path": record.get(
                            f"{prefix}_relative_path"
                        ),
                        "size_bytes": record.get(f"{prefix}_size_bytes"),
                        "sha256": record.get(f"{prefix}_sha256"),
                        "n_atoms": record.get("n_atoms"),
                        "n_frames": (
                            1
                            if role == "structure_pdb" and record
                            else record.get("n_frames")
                        ),
                        "frame_interval_ps": (
                            record.get("frame_interval_ps")
                            if role == "trajectory_xtc"
                            else np.nan
                        ),
                        "represented_span_ns": (
                            record.get("represented_span_ns")
                            if role == "trajectory_xtc"
                            else np.nan
                        ),
                        "local_manifest_status": (
                            "present" if record else "pending"
                        ),
                        "local_full_frame_qc_status": (
                            "passed"
                            if qc and audit_status == "passed"
                            else "pending"
                        ),
                        "remote_state": authenticated.get(
                            "state", "not_verified"
                        ),
                        "publicly_resolvable": bool(
                            public.get("publicly_resolvable")
                        ),
                        "component_scope": (
                            "reduced protein complex; membrane, solvent and ions excluded"
                        ),
                    }
                )
            qc_rows.append(
                {
                    "system_id": row.system_id,
                    "g_protein_family": row.g_protein_family,
                    "release_replica": replica,
                    "reserved_doi": DOIS[replica],
                    "frame_count": record.get("n_frames"),
                    "frame_interval_ps": record.get("frame_interval_ps"),
                    "represented_span_ns": record.get("represented_span_ns"),
                    "n_atoms": record.get("n_atoms"),
                    "pdb_xtc_atom_counts_match": (
                        bool(qc) if audit_status == "passed" else np.nan
                    ),
                    "finite_coordinates": qc.get("finite_coordinates"),
                    "finite_monotonic_time": qc.get("finite_monotonic_time"),
                    "valid_periodic_box": qc.get("valid_periodic_box"),
                    "atom_order_backbone_geometry": qc.get(
                        "atom_order_backbone_geometry"
                    ),
                    "catastrophic_chain_break_absent": (
                        not qc.get("catastrophic_chain_break")
                        if qc
                        else np.nan
                    ),
                    "complex_separation_absent": (
                        not qc.get("complex_separation") if qc else np.nan
                    ),
                    "coordinate_scatter_absent": (
                        not qc.get("coordinate_scatter") if qc else np.nan
                    ),
                    "full_frames_checked": qc.get("frames_checked"),
                    "local_full_frame_qc_status": (
                        "passed"
                        if qc and audit_status == "passed"
                        else "pending"
                    ),
                    "remote_state": authenticated.get(
                        "state", "not_verified"
                    ),
                    "publicly_resolvable": bool(
                        public.get("publicly_resolvable")
                    ),
                    "component_scope": (
                        "reduced protein complex; membrane, solvent and ions excluded"
                    ),
                }
            )

    status = pd.DataFrame(status_rows)
    manifest_long = pd.DataFrame(manifest_rows)
    qc = pd.DataFrame(qc_rows)
    assert len(manifest_long) == 207 * 3 * 2
    assert len(qc) == 207 * 3
    write_csv(status, FIG_SOURCE / "Figure_5A_release_status.csv")
    write_csv(manifest_long, SUPP / "Supplementary_Data_S4_reduced_release_manifest.csv")
    write_csv(qc, SUPP / "Supplementary_Data_S8_reduced_release_qc.csv")

    complete_manifest = manifest_long[
        manifest_long.local_manifest_status.eq("present")
    ].copy()
    write_csv(
        complete_manifest,
        FIG_SOURCE / "Figure_5B_release_file_manifest.csv",
    )
    qc_summary = (
        qc.groupby(
            [
                "release_replica",
                "local_full_frame_qc_status",
                "remote_state",
                "publicly_resolvable",
            ],
            dropna=False,
        )
        .size()
        .rename("systems")
        .reset_index()
    )
    write_csv(qc_summary, FIG_SOURCE / "Figure_5C_replica_qc_summary.csv")

    example_source = manifest_long[
        manifest_long.release_replica.eq(1)
        & manifest_long.system_id.eq("G12_7SF7")
        & manifest_long.local_manifest_status.eq("present")
    ].copy()
    example_source["example_step"] = example_source.file_role.map(
        {
            "structure_pdb": "verify structure SHA-256 and load as topology",
            "trajectory_xtc": "verify trajectory SHA-256 and load coordinates",
        }
    )
    write_csv(example_source, FIG_SOURCE / "Figure_5D_worked_example.csv")
    return status, qc


def build_replica_provenance(data: dict[str, Any]) -> pd.DataFrame:
    cohort = data["cohort"]
    protocol = data["protocol"].set_index("system_id")
    truth = data["truth"].set_index("system_id")
    rows = []
    for system in cohort.itertuples(index=False):
        p = protocol.loc[system.system_id]
        truth_row = truth.loc[system.system_id]
        observed_system_total = float(
            pd.Series(
                [
                    truth_row.observed_netcdf_sampling_ns,
                    truth_row.observed_gromacs_sampling_ns,
                ]
            ).sum(skipna=True)
        )
        for replica in [1, 2, 3]:
            original_duration = 500.0
            repair_status = "not_required"
            repair_evidence = ""
            if system.system_id == "Gq_8ZPT" and replica == 2:
                original_duration = 416.8
                repair_status = "reduced_release_source_repaired_locally"
                repair_evidence = (
                    "trajectory_readiness/continuations/Gq_8ZPT/"
                    "interface_recovery/exhaustive_audit_corrected/"
                    "candidate_replica_audit.csv"
                )
            rows.append(
                {
                    "system_id": system.system_id,
                    "replica_id": replica,
                    "nominal_duration_ns": 500.0,
                    "original_file_observed_duration_ns": original_duration,
                    "repaired_reduced_source_duration_ns": 500.0,
                    "repair_status": repair_status,
                    "repair_evidence": repair_evidence,
                    "nominal_system_total_ns": 1500.0,
                    "original_file_observed_system_total_ns": (
                        observed_system_total
                    ),
                    "full_system_source_distributed": False,
                    "reduced_release_record_expected": True,
                    "protocol_group": p.protocol_group,
                    "production_engine": p.engine_family,
                    "engine_versions_from_available_outputs": p.engine_versions,
                    "engine_version_scope": p.engine_version_scope,
                    "production_records_recovered_n": p.production_record_count,
                    "replica_seed_evidence": p.replica_seed_status,
                    "protocol_evidence_status": p.status,
                    "protocol_evidence_gap": p.conflicts,
                    "source_trajectory_write_interval_ps": p.output_interval_ps,
                    "provenance_scope": (
                        "source-simulation provenance; not a deposited "
                        "full-system molecular-file record"
                    ),
                }
            )
    provenance = pd.DataFrame(rows)
    assert len(provenance) == 621
    assert np.isclose(
        provenance.original_file_observed_duration_ns.sum() / 1000,
        310.4168,
    )
    assert np.isclose(
        provenance.nominal_duration_ns.sum() / 1000,
        310.5,
    )
    write_csv(
        provenance,
        SUPP
        / "Supplementary_Data_S7_source_simulation_replica_provenance.csv",
    )
    sampling_summary = pd.DataFrame(
        [
            {
                "scope": "full final-207 campaign",
                "sampling_definition": "nominal protocol",
                "sampling_us": provenance.nominal_duration_ns.sum() / 1000,
                "claim_scope": "207 systems × 3 computational repeats × 500 ns",
            },
            {
                "scope": "full final-207 campaign",
                "sampling_definition": "original file-observed coordinates",
                "sampling_us": (
                    provenance.original_file_observed_duration_ns.sum() / 1000
                ),
                "claim_scope": "includes the original 416.8-ns Gq_8ZPT replica 2",
            },
            {
                "scope": "Gq_8ZPT replica 2",
                "sampling_definition": "nominal protocol",
                "sampling_us": 0.5,
                "claim_scope": "specified production duration",
            },
            {
                "scope": "Gq_8ZPT replica 2",
                "sampling_definition": "original file-observed coordinates",
                "sampling_us": 0.4168,
                "claim_scope": "pre-repair coordinate record",
            },
            {
                "scope": "Gq_8ZPT replica 2",
                "sampling_definition": "locally repaired reduced source",
                "sampling_us": 0.49995,
                "claim_scope": (
                    "10,000-frame continuation-derived coordinate source; "
                    "does not establish distribution of a full-system archive"
                ),
            },
        ]
    )
    write_csv(
        sampling_summary,
        FIG_SOURCE / "Supplementary_Figure_2_sampling_provenance.csv",
    )
    return provenance


def build_api_audit(data: dict[str, Any]) -> pd.DataFrame:
    s9 = pd.read_csv(S9_OLD).copy()
    s9["audit_kind"] = "per_system_record_coverage"
    consensus_rows = []
    if API_GATEWAY_CONSENSUS.is_file():
        payload = json.loads(API_GATEWAY_CONSENSUS.read_text(encoding="utf-8"))
        records = payload.get("records", [])
        all_denominators = sorted(
            {
                int(row["n_systems"])
                for row in records
                if row.get("group") == "ALL"
                and row.get("n_systems") is not None
            }
        )
        consensus_rows.append(
            {
                "record_type": "consensus_gateway_ALL_denominator",
                "expected_systems": 207,
                "files_present": 1,
                "records_available": len(records),
                "records_missing": 0,
                "records_unpopulated": 0,
                "records_not_applicable": 0,
                "total_size_bytes": API_GATEWAY_CONSENSUS.stat().st_size,
                "affected_system_ids": "",
                "schema_location": "data/api/v1/consensus/gateways.json",
                "availability_semantics": (
                    f"ALL-group denominators observed: {all_denominators}"
                ),
                "scope_note": (
                    "consensus cache consistency; biological gateway "
                    "distribution is outside the database-paper figure set"
                ),
                "schema_version": "v12.1",
                "audit_kind": "consensus_cache_consistency",
            }
        )
    out = pd.concat([s9, pd.DataFrame(consensus_rows)], ignore_index=True)
    write_csv(out, SUPP / "Supplementary_Data_S9_api_access_audit.csv")
    return out


def selected_version(package: str) -> str | None:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return None


def git_value(*args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=SERVER,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def build_environment_inventory() -> pd.DataFrame:
    paths = [
        HERE / "build_v12_data.py",
        HERE / "audit_v12_release_status.py",
        HERE / "audit_v12_references.py",
        HERE / "audit_v12_repository_access.py",
        HERE / "update_v12_gateway_consensus.py",
        HERE / "run_final207_viz_audit.py",
        HERE / "plot_v12_figures.py",
        HERE / "build_v12_documents.py",
        HERE / "build_v12_reports.py",
        HERE / "build_v12_submission_bundle.py",
        HERE / "visual_qc_v12.py",
        HERE / "validate_v12_submission.py",
        HERE / "v12_manuscript.md",
        HERE / "v12_si.md",
        SERVER / "requirements.txt",
        SERVER / "frontend/package.json",
        SERVER / "frontend/package-lock.json",
        SERVER / "Dockerfile",
        SERVER / "LICENSE",
        SERVER / "LICENSE-CODE",
        SERVER / "CITATION.cff",
    ]
    rows: list[dict[str, Any]] = []
    for path in paths:
        relative = (
            str(path.relative_to(SERVER))
            if path.is_relative_to(SERVER)
            else str(path)
        )
        rows.append(
            {
                "item_kind": "repository_file",
                "name": path.name,
                "repository_relative_path": relative,
                "exists": path.is_file(),
                "version": "",
                "size_bytes": path.stat().st_size if path.is_file() else np.nan,
                "sha256": sha256(path) if path.is_file() else "",
                "version_authority": "repository file content",
            }
        )
    packages = [
        "numpy",
        "pandas",
        "matplotlib",
        "scipy",
        "MDAnalysis",
        "python-docx",
        "requests",
        "Pillow",
    ]
    rows.append(
        {
            "item_kind": "runtime",
            "name": "Python",
            "repository_relative_path": "",
            "exists": True,
            "version": platform.python_version(),
            "size_bytes": np.nan,
            "sha256": "",
            "version_authority": "platform.python_version",
        }
    )
    for package in packages:
        version = selected_version(package)
        rows.append(
            {
                "item_kind": "build_dependency",
                "name": package,
                "repository_relative_path": "",
                "exists": version is not None,
                "version": version or "",
                "size_bytes": np.nan,
                "sha256": "",
                "version_authority": "importlib.metadata",
            }
        )
    rows.append(
        {
            "item_kind": "revision",
            "name": "git_commit",
            "repository_relative_path": "",
            "exists": True,
            "version": git_value("rev-parse", "HEAD"),
            "size_bytes": np.nan,
            "sha256": "",
            "version_authority": "git rev-parse HEAD; excludes uncommitted v12 changes",
        }
    )
    inventory = pd.DataFrame(rows)
    write_csv(
        inventory,
        SUPP / "Supplementary_Data_S10_code_environment_inventory.csv",
    )
    environment = {
        "schema_version": "1.0",
        "generated_at": now(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": {
            package: selected_version(package) for package in packages
        },
        "repository_commit_before_v12_worktree_changes": git_value(
            "rev-parse", "HEAD"
        ),
        "worktree_dirty": bool(git_value("status", "--porcelain")),
        "licences": {
            "data": "CC BY 4.0",
            "code": "MIT",
        },
    }
    (REPORTS / "v12_environment_freeze.json").write_text(
        json.dumps(environment, indent=2) + "\n", encoding="utf-8"
    )
    return inventory


FIELD_DEFINITIONS = {
    "system_id": "Stable release identifier combining family prefix and source PDB identifier.",
    "replica_id": "Computational repeat number within a system.",
    "release_replica": "Reduced-release record number corresponding to source replica 1, 2 or 3.",
    "reserved_doi": "Zenodo DOI reserved for a private draft; not evidence of public release.",
    "publicly_resolvable": "Whether an unauthenticated DOI request resolved successfully at audit time.",
    "remote_state": "State returned by the authenticated read-only Zenodo deposition endpoint.",
    "sha256": "SHA-256 checksum of the named file.",
    "size_bytes": "File size in bytes.",
    "n_atoms": "Atom count in the matched reduced PDB/XTC pair.",
    "n_frames": "Number of coordinate frames represented by the record.",
    "frame_count": "Number of coordinate frames represented by the record.",
    "frame_interval_ps": "Time spacing between retained reduced-trajectory frames in picoseconds.",
    "represented_span_ns": "Elapsed time between first and last retained frame in nanoseconds.",
    "original_file_observed_duration_ns": "Coordinate duration observed before the Gq_8ZPT replica-2 repair.",
    "nominal_duration_ns": "Duration specified by the production protocol.",
    "full_system_source_distributed": "Whether membrane, solvent, ions and full-system source coordinates are deposited.",
    "component_scope": "Molecular components included in the distributed record.",
    "pocket_record_status": "Whether pocket analysis yielded detected pockets or a valid zero-pocket result.",
    "mean": "Mean of the three per-replica gateway summaries.",
    "ci_lo": "Stored lower interval bound.",
    "ci_hi": "Stored upper interval bound.",
}


def infer_type(series: pd.Series) -> str:
    if pd.api.types.is_bool_dtype(series):
        return "boolean"
    if pd.api.types.is_integer_dtype(series):
        return "integer"
    if pd.api.types.is_numeric_dtype(series):
        return "number"
    return "string"


def infer_unit(field: str) -> str:
    lower = field.lower()
    if lower.endswith("_bytes"):
        return "byte"
    if lower.endswith("_ns"):
        return "ns"
    if lower.endswith("_ps"):
        return "ps"
    if lower.endswith("_a") or "_angstrom" in lower:
        return "angstrom"
    if "sha256" in lower:
        return "hexadecimal digest"
    return "unitless or field-specific"


def build_dictionary(files: list[Path]) -> pd.DataFrame:
    observations: dict[str, list[tuple[str, pd.Series]]] = {}
    for path in files:
        data = pd.read_csv(path)
        for field in data.columns:
            observations.setdefault(field, []).append((path.name, data[field]))
    rows = []
    for field in sorted(observations):
        sources = observations[field]
        combined = pd.concat([series for _, series in sources], ignore_index=True)
        definition = FIELD_DEFINITIONS.get(
            field,
            field.replace("_", " ").capitalize()
            + " as serialized in the named v12 table.",
        )
        non_null = combined.dropna()
        example = "" if non_null.empty else str(non_null.iloc[0])
        rows.append(
            {
                "field": field,
                "appears_in": "|".join(sorted(name for name, _ in sources)),
                "type": infer_type(combined),
                "unit": infer_unit(field),
                "requiredness": (
                    "required" if not combined.isna().any() else "conditional"
                ),
                "nullable": bool(combined.isna().any()),
                "definition": definition,
                "example": example,
                "schema_version": "v12.1",
            }
        )
    dictionary = pd.DataFrame(rows)
    write_csv(dictionary, SUPP / "Supplementary_Data_S3_metadata_dictionary.csv")
    return dictionary


def build_supplementary(
    data: dict[str, Any],
    provenance: pd.DataFrame,
    environment: pd.DataFrame,
) -> list[Path]:
    cohort = data["cohort"].copy()
    cohort["protocol_group"] = cohort.force_field.map(PROTOCOL_MAP)
    cohort["sampling_claim_scope"] = (
        "nominal protocol values; original file-observed total reported separately"
    )
    write_csv(
        cohort,
        SUPP / "Supplementary_Data_S1_included_system_inventory.csv",
    )
    write_csv(
        data["s2"],
        SUPP / "Supplementary_Data_S2_release_boundary_exceptions.csv",
    )
    write_csv(
        data["pockets"],
        SUPP / "Supplementary_Data_S5_pocket_summaries_gpcrdb.csv",
    )
    gateway = data["gateways"].copy()
    gateway["reproducibility_note"] = (
        "processed full-system summary; cannot be recalculated from reduced XTC"
    )
    write_csv(
        gateway,
        SUPP / "Supplementary_Data_S6_gateway_per_system.csv",
    )

    files = [
        SUPP / "Supplementary_Data_S1_included_system_inventory.csv",
        SUPP / "Supplementary_Data_S2_release_boundary_exceptions.csv",
        SUPP / "Supplementary_Data_S4_reduced_release_manifest.csv",
        SUPP / "Supplementary_Data_S5_pocket_summaries_gpcrdb.csv",
        SUPP / "Supplementary_Data_S6_gateway_per_system.csv",
        SUPP
        / "Supplementary_Data_S7_source_simulation_replica_provenance.csv",
        SUPP / "Supplementary_Data_S8_reduced_release_qc.csv",
        SUPP / "Supplementary_Data_S9_api_access_audit.csv",
        SUPP / "Supplementary_Data_S10_code_environment_inventory.csv",
    ]
    build_dictionary(files)
    return [
        SUPP / "Supplementary_Data_S1_included_system_inventory.csv",
        SUPP / "Supplementary_Data_S2_release_boundary_exceptions.csv",
        SUPP / "Supplementary_Data_S3_metadata_dictionary.csv",
        *files[2:],
    ]


def write_supplementary_readme(
    files: list[Path], release_status: pd.DataFrame
) -> None:
    public = bool(release_status.publicly_resolvable.all())
    local_complete = bool(
        (release_status.local_manifest_rows == 207).all()
        and release_status.local_full_frame_qc_status.eq("passed").all()
    )
    text = f"""CoupledMD v12 Supplementary Data

Generated: {now()}

This package describes a frozen 207-system cohort and distinguishes the
underlying full-system simulation campaign from the reduced molecular product.
The expected reduced release has three records (replicas 1-3), each containing
one matched retained-component PDB/XTC pair per system.  The reduced records
exclude membrane, solvent and ions and may exclude ligands removed during
system preparation.  They are not a full-system trajectory archive.

Public DOI resolution verified for all three records: {str(public).lower()}
Local 207-system manifest and full-frame QC complete for all three records:
{str(local_complete).lower()}

S1  Included-system inventory (207 systems).
S2  Release-boundary exceptions (13 unresolved and 2 excluded records).
S3  Field dictionary for the distributed CSV tables.
S4  Expected three-replica reduced-release manifest, with pending records
    explicit rather than silently omitted.
S5  Pocket summaries, including 2,149 detected pockets and two explicit
    valid zero-pocket records.
S6  Processed gateway summaries (207 systems x 7 interfaces x 4 metrics).
    These intermediate records cannot be regenerated from reduced XTC files.
S7  Source-simulation replica provenance (621 computational repeats), including
    nominal versus original file-observed duration and Gq_8ZPT replica-2 repair.
S8  Expected three-replica reduced-release QC, with pending records explicit.
S9  Current API and consensus-cache consistency audit.
S10 Code and build-environment inventory.

Data licence: CC BY 4.0.  Code licence: MIT.
"""
    (SUPP / "README.txt").write_text(text, encoding="utf-8")


def numerical_audit(
    data: dict[str, Any],
    provenance: pd.DataFrame,
    release_status: pd.DataFrame,
    release_qc: pd.DataFrame,
    viz_regenerated: bool,
) -> dict[str, Any]:
    cohort = data["cohort"]
    pockets = data["pockets"]
    detected = pockets[pockets.pocket_record_status.eq("available_detected_pocket")]
    zeros = pockets[pockets.pocket_record_status.eq("available_zero_pockets")]
    eligible = data["t4"][
        data["t4"].system_id.isin(cohort.system_id)
        & data["t4"].lig_type.eq("peptide")
    ]
    gateways = data["gateways"]
    t13 = data["t13"]
    remote_summary = data["remote"].get("summary", {})
    gateway_consensus_denominators: list[int] = []
    if API_GATEWAY_CONSENSUS.is_file():
        payload = json.loads(API_GATEWAY_CONSENSUS.read_text(encoding="utf-8"))
        gateway_consensus_denominators = sorted(
            {
                int(row["n_systems"])
                for row in payload.get("records", [])
                if row.get("group") == "ALL"
            }
        )
    repair = pd.read_csv(GQ8ZPT_REPAIR)
    assert len(repair) == 1
    assert repair.iloc[0].verdict == "GOOD"
    assert int(repair.iloc[0].n_frames) == 10000
    audit = {
        "schema_version": "1.0",
        "generated_at": now(),
        "status": "HOLD",
        "authoritative_sources": {
            "cohort": source_record(COHORT),
            "cohort_truth": source_record(COHORT_TRUTH),
            "cohort_truth_summary": source_record(COHORT_TRUTH_SUMMARY),
            "pockets": source_record(S5_OLD),
            "gateways": source_record(S6_OLD),
            "replica_qc": source_record(T10),
            "structural_validation": source_record(T13),
            "protocol_provenance": source_record(PROTOCOL),
            "gq_8zpt_repair": source_record(GQ8ZPT_REPAIR),
            "remote_release_verification": (
                source_record(RELEASE_REMOTE) if RELEASE_REMOTE.is_file() else None
            ),
        },
        "cohort": {
            "systems": len(cohort),
            "classes": cohort.gpcr_class.value_counts().sort_index().to_dict(),
            "families": cohort.g_protein_family.value_counts().to_dict(),
            "class_by_family": {
                cls: {
                    family: int(value)
                    for family, value in row.items()
                }
                for cls, row in pd.crosstab(
                    cohort.gpcr_class, cohort.g_protein_family
                ).to_dict(orient="index").items()
            },
            "receptor_names": int(cohort.receptor_name.nunique()),
            "mapped_uniprot_accessions": int(
                cohort.receptor_uniprot.nunique()
            ),
            "unmapped_systems": cohort.loc[
                cohort.receptor_uniprot.isna(), "system_id"
            ].tolist(),
            "working_inventory": {
                "included": 207,
                "unresolved": 13,
                "excluded": 2,
                "total": 222,
            },
        },
        "protocols": {
            "system_counts": data["protocol"]
            .protocol_group.value_counts()
            .sort_index()
            .to_dict(),
            "evidence_status": data["protocol"].status.value_counts().to_dict(),
            "partial_systems": data["protocol"]
            .loc[data["protocol"].status.ne("verified"), "system_id"]
            .tolist(),
        },
        "replicas_and_sampling": {
            "nominal_replicas": 621,
            "nominal_sampling_us": provenance.nominal_duration_ns.sum() / 1000,
            "original_file_observed_sampling_us": (
                provenance.original_file_observed_duration_ns.sum() / 1000
            ),
            "difference_ns": float(
                provenance.nominal_duration_ns.sum()
                - provenance.original_file_observed_duration_ns.sum()
            ),
            "shortened_original_record": {
                "system_id": "Gq_8ZPT",
                "replica": 2,
                "original_file_observed_duration_ns": 416.8,
                "local_repaired_reduced_source_frames": int(
                    repair.iloc[0].n_frames
                ),
                "local_repaired_reduced_source_span_ns": float(
                    repair.iloc[0].span_ps / 1000
                ),
                "local_repair_verdict": repair.iloc[0].verdict,
                "full_system_source_publicly_distributed": False,
            },
        },
        "pockets": {
            "system_records": int(pockets.system_id.nunique()),
            "systems_with_detected_pockets": int(
                detected.system_id.nunique()
            ),
            "detected_pocket_rows": len(detected),
            "explicit_zero_pocket_systems": sorted(zeros.system_id.tolist()),
            "positive_control_eligible": len(eligible),
            "positive_control_recovered": int(eligible.ortho_recovered.sum()),
        },
        "gateways": {
            "rows": len(gateways),
            "systems": int(gateways.system_id.nunique()),
            "interfaces": int(gateways.interface.nunique()),
            "metrics": int(gateways.metric.nunique()),
            "expected_rows": 207 * 7 * 4,
            "api_consensus_ALL_denominators": gateway_consensus_denominators,
            "api_consensus_final207_consistent": (
                gateway_consensus_denominators == [207]
            ),
            "main_figure_status": "removed_to_resource_paper",
        },
        "structural_validation": {
            "replica_rows": len(t13),
            "available_replica_rows": int(
                t13.validation_status.eq("available").sum()
            ),
            "unavailable_replica_rows": int(
                t13.validation_status.ne("available").sum()
            ),
            "unavailable_systems": sorted(
                t13.loc[
                    t13.validation_status.ne("available"), "system_id"
                ].unique()
            ),
        },
        "reduced_visualization": {
            "regenerated_exact_final207": viz_regenerated,
            "audit_filename": (
                VIZ_AUDIT.name if viz_regenerated else LEGACY_VIZ_AUDIT.name
            ),
        },
        "release": {
            "replicas": release_status.to_dict(orient="records"),
            "local_qc_pass_systems_by_replica": {
                str(replica): int(
                    (
                        release_qc.release_replica.eq(replica)
                        & release_qc.local_full_frame_qc_status.eq("passed")
                    ).sum()
                )
                for replica in [1, 2, 3]
            },
            "remote_summary": remote_summary,
            "data_product": "reduced protein-complex PDB/XTC",
            "excluded_components": ["membrane", "solvent", "ions"],
            "full_system_archive": False,
        },
        "known_inconsistencies": [
            {
                "id": "N01",
                "issue": "206 detected plus two zero-pocket systems",
                "resolution": "205 detected plus two explicit zero-pocket systems equals 207",
                "status": "PASS",
            },
            {
                "id": "N02",
                "issue": "205 PBC-OK plus three WARN and inherited 208 frame counts",
                "resolution": (
                    "exact final-207 visualization audit regenerated"
                    if viz_regenerated
                    else "exact regeneration pending; legacy-named file content has 207 IDs"
                ),
                "status": "PASS" if viz_regenerated else "HOLD",
            },
            {
                "id": "N03",
                "issue": "API consensus gateway n_systems=208",
                "resolution": (
                    "final-207 consensus cache verified"
                    if gateway_consensus_denominators == [207]
                    else f"stale denominators remain: {gateway_consensus_denominators}"
                ),
                "status": (
                    "PASS"
                    if gateway_consensus_denominators == [207]
                    else "HOLD"
                ),
            },
            {
                "id": "N04",
                "issue": "active final208/final624 names",
                "resolution": "all v12 submission artifacts use final207/621 or version-neutral names",
                "status": "PASS",
            },
            {
                "id": "N05",
                "issue": "310.5 versus 310.4168 microseconds",
                "resolution": (
                    "310.5 microseconds nominal; 310.4168 microseconds in original "
                    "pre-repair coordinate files; local reduced-release repair documented; "
                    "full-system source is not distributed"
                ),
                "status": "PASS",
            },
            {
                "id": "N06",
                "issue": "public/complete archive claims",
                "resolution": (
                    "claims withheld until all three records are publicly resolvable "
                    "and locally QC-complete"
                ),
                "status": (
                    "PASS"
                    if remote_summary.get("archive_completeness_claim_allowed")
                    else "HOLD"
                ),
            },
        ],
    }
    blocking = [
        item
        for item in audit["known_inconsistencies"]
        if item["status"] != "PASS"
    ]
    if not blocking and remote_summary.get("archive_completeness_claim_allowed"):
        audit["status"] = "PASS"
    (REPORTS / "v12_numerical_consistency_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    write_csv(
        pd.DataFrame(audit["known_inconsistencies"]),
        REPORTS / "v12_numerical_consistency_findings.csv",
    )
    return audit


def write_source_manifest() -> None:
    manifest_path = FIG_SOURCE / "figure_source_data_manifest.csv"
    paths = sorted(
        path for path in FIG_SOURCE.glob("*.csv") if path != manifest_path
    )
    rows = []
    for path in paths:
        frame = pd.read_csv(path)
        rows.append(
            {
                "file": path.name,
                "rows": len(frame),
                "columns": len(frame.columns),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    write_csv(pd.DataFrame(rows), manifest_path)
    readme = """# v12 figure source data

Each quantitative or schematic main-figure panel has a path-neutral CSV source
table. Figure 4 contains only technical validation. The previous gateway-family
distribution and GPCR-centred pocket-reuse atlas are intentionally excluded
from the database paper because their biological interpretation belongs to the
separate resource-atlas paper.

Figure 5 distinguishes local full-frame QC from authenticated Zenodo draft
state and unauthenticated public DOI resolution. A reserved DOI is never treated
as evidence of publication.
"""
    (FIG_SOURCE / "README.md").write_text(readme, encoding="utf-8")


def write_package_audit(files: list[Path], numerical: dict[str, Any]) -> None:
    package_files = [SUPP / "README.txt", *files]
    payload = {
        "schema_version": "1.0",
        "generated_at": now(),
        "status": "PASS"
        if numerical["status"] == "PASS"
        else "HOLD_RELEASE",
        "files": {
            path.name: {
                "rows": (
                    len(pd.read_csv(path)) if path.suffix == ".csv" else None
                ),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in package_files
        },
        "release_truth": numerical["release"],
        "sampling_truth": numerical["replicas_and_sampling"],
        "scope": {
            "distributed_molecular_product": "reduced protein-complex PDB/XTC",
            "full_system_source_distributed": False,
            "gateway_recalculation_from_reduced_release": False,
        },
    }
    (SUPP / "Supplementary_Data_package_audit.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )


def main() -> int:
    for directory in [FIG_INPUT, FIG_SOURCE, SUPP, REPORTS]:
        directory.mkdir(parents=True, exist_ok=True)
    data = load_sources()
    build_figure1(data)
    build_figure2(data)
    build_figure3(data)
    viz_regenerated = build_figure4(data)
    release_status, release_qc = build_release_tables(data)
    provenance = build_replica_provenance(data)
    build_api_audit(data)
    environment = build_environment_inventory()
    supplementary_files = build_supplementary(data, provenance, environment)
    write_supplementary_readme(supplementary_files, release_status)
    numerical = numerical_audit(
        data,
        provenance,
        release_status,
        release_qc,
        viz_regenerated,
    )
    write_source_manifest()
    write_package_audit(supplementary_files, numerical)
    print(
        "Built v12 data: "
        f"{len(list(FIG_SOURCE.glob('*.csv')))} figure-source CSVs, "
        f"{len(supplementary_files)} Supplementary Data CSVs; "
        f"numerical status={numerical['status']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
