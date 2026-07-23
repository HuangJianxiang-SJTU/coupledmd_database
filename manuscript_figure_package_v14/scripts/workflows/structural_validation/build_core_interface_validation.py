#!/usr/bin/env python3
"""Build harmonized TM-core/interface validation records for selected replicas.

Every replica is sampled at exactly 1,001 evenly spaced frames.  Frames are
superposed on the sequence-verified receptor TM1--TM7 C-alpha core.  The three
reported quantities are TM-core RMSD, G-alpha interface-region RMSD after that
same fit, and retention of initial receptor--G-alpha C-alpha contacts.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse
import warnings

import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from build_final624_replica_qc_v9 import corrected_sources


ROOT = Path(__file__).resolve().parent
QC = ROOT / "v9_figure_inputs/scidata_T10_final624_replica_qc_v9.csv"
SELECTIONS = ROOT / "v9_figure_inputs/scidata_T12_final208_structural_selection_audit_v9.csv"
MANIFEST = ROOT / "trajectory_readiness/replica_manifest.csv"
OUT = ROOT / "v9_figure_inputs/scidata_T13_final624_core_interface_validation_v9.csv"
PROGRESS = ROOT / "v9_figure_inputs/scidata_T13_final624_core_interface_validation_v9.progress.csv"
N_SAMPLES = 1001


def parse_ordinals(value):
    if not isinstance(value, str) or not value:
        return []
    return [int(x) for x in value.split(";")]


def sample_indices(n_frames):
    if n_frames < N_SAMPLES:
        raise ValueError(f"only {n_frames} frames; need at least {N_SAMPLES}")
    idx = np.rint(np.linspace(0, n_frames - 1, N_SAMPLES)).astype(int)
    if len(np.unique(idx)) != N_SAMPLES:
        raise ValueError("uniform frame selection is not unique")
    return idx


def ca_atom_indices(universe, ordinals):
    result = []
    for ordinal in ordinals:
        atoms = universe.residues[ordinal - 1].atoms.select_atoms("name CA")
        if len(atoms) != 1:
            raise ValueError(f"residue ordinal {ordinal} does not have exactly one CA")
        result.append(int(atoms.indices[0]))
    return np.asarray(result, dtype=int)


def kabsch(mobile, reference):
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    x = mobile - mobile_center
    y = reference - reference_center
    u, _, vt = np.linalg.svd(x.T @ y)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    return rotation, mobile_center, reference_center


def summarize(values, lower=False):
    q = np.quantile(np.asarray(values, dtype=float), [.05, .25, .50, .75, .95])
    return {"p05":q[0], "q25":q[1], "median":q[2], "q75":q[3], "p95":q[4]}


def unavailable(task, reason):
    return {
        "system_id":task["system_id"], "replica":task["replica"],
        "g_protein_family":task["g_protein_family"],
        "trajectory_path":task["trajectory_path"], "topology_path":task["topology_path"],
        "identity_status":task["identity_status"], "validation_status":"unavailable",
        "reason_code":reason, "harmonized_observations":0,
    }


def process_replica(task):
    if task["selection_status"] != "available":
        return unavailable(task, task["selection_reason"])
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            universe = mda.Universe(task["topology_path"], task["trajectory_path"])
        tm_idx = ca_atom_indices(universe, parse_ordinals(task["tm_ordinals"]))
        receptor_idx = ca_atom_indices(universe, parse_ordinals(task["receptor_ordinals"]))
        galpha_idx = ca_atom_indices(universe, parse_ordinals(task["galpha_ordinals"]))
        indices = sample_indices(len(universe.trajectory))

        universe.trajectory[indices[0]]
        reference_tm = universe.atoms[tm_idx].positions.astype(float).copy()
        receptor0 = universe.atoms[receptor_idx].positions.astype(float)
        galpha0 = universe.atoms[galpha_idx].positions.astype(float)
        d0 = cdist(receptor0, galpha0)
        interface_local = np.where((d0 <= 12.0).any(axis=0))[0]
        contact_pairs = np.argwhere(d0 <= 8.0)
        if len(interface_local) < 5 or len(contact_pairs) < 3:
            return unavailable(
                task,
                f"insufficient_replica_initial_interface:n_interface={len(interface_local)};n_contacts={len(contact_pairs)}",
            )
        interface_idx = galpha_idx[interface_local]
        reference_interface = universe.atoms[interface_idx].positions.astype(float).copy()
        contact_receptor_idx = receptor_idx[contact_pairs[:, 0]]
        contact_galpha_idx = galpha_idx[contact_pairs[:, 1]]

        tm_rmsd = np.empty(N_SAMPLES, dtype=float)
        interface_rmsd = np.empty(N_SAMPLES, dtype=float)
        contact_retention = np.empty(N_SAMPLES, dtype=float)
        for j, frame in enumerate(indices):
            universe.trajectory[int(frame)]
            mobile_tm = universe.atoms[tm_idx].positions.astype(float)
            rotation, mobile_center, reference_center = kabsch(mobile_tm, reference_tm)
            fitted_tm = (mobile_tm - mobile_center) @ rotation + reference_center
            mobile_interface = universe.atoms[interface_idx].positions.astype(float)
            fitted_interface = (mobile_interface - mobile_center) @ rotation + reference_center
            tm_rmsd[j] = np.sqrt(np.mean(np.sum((fitted_tm - reference_tm) ** 2, axis=1)))
            interface_rmsd[j] = np.sqrt(np.mean(np.sum((fitted_interface - reference_interface) ** 2, axis=1)))
            delta = universe.atoms[contact_receptor_idx].positions - universe.atoms[contact_galpha_idx].positions
            contact_retention[j] = np.mean(np.linalg.norm(delta, axis=1) <= 12.0)

        if not (np.isfinite(tm_rmsd).all() and np.isfinite(interface_rmsd).all() and np.isfinite(contact_retention).all()):
            return unavailable(task, "nonfinite_derived_metric")
        tm = summarize(tm_rmsd)
        interface = summarize(interface_rmsd)
        retention = summarize(contact_retention)
        row = {
            "system_id":task["system_id"], "replica":task["replica"],
            "g_protein_family":task["g_protein_family"],
            "trajectory_path":task["trajectory_path"], "topology_path":task["topology_path"],
            "identity_status":task["identity_status"], "validation_status":"available", "reason_code":"",
            "source_frames":len(universe.trajectory), "harmonized_observations":N_SAMPLES,
            "n_tm_core_ca":len(tm_idx), "n_galpha_interface_ca":len(interface_idx),
            "n_initial_contacts_8A":len(contact_pairs), "contact_retention_cutoff_A":12.0,
        }
        for prefix, values in [("tm_core_rmsd_A", tm), ("galpha_interface_rmsd_A", interface), ("contact_retention", retention)]:
            row.update({f"{prefix}_{key}":value for key, value in values.items()})
        return row
    except Exception as exc:
        return unavailable(task, f"calculation_error:{type(exc).__name__}:{exc}")


def trajectory_paths(qc):
    manifest = pd.read_csv(MANIFEST)
    baseline = {
        (r.system_id, int(r.replica)): str(Path(r.aligned_path))
        for r in manifest.itertuples(index=False)
    }
    # The two GROMACS systems absent from the NC component summary have audited
    # protein-aligned trajectories in the dedicated TRR audit.  The inherited
    # manifest's ``aligned_path`` values for these records point to raw XTCs.
    trr_audit = ROOT / "trajectory_readiness/component_audit_summary_trr.csv"
    if trr_audit.exists():
        trr = pd.read_csv(trr_audit)
        for row in trr.itertuples(index=False):
            baseline[(row.system_id, int(row.replica))] = str(Path(row.aligned_trajectory_path))
    corrected = {}
    for path, allowed in corrected_sources():
        audit = pd.read_csv(path)
        audit = audit[audit.verdict.eq("GOOD")]
        if allowed is not None:
            audit = audit[audit.system_id.isin(allowed)]
        for row in audit.itertuples(index=False):
            p = Path(row.trajectory_path)
            if not p.is_absolute():
                p = ROOT / p
            corrected[(row.system_id, int(row.replica))] = str(p)
    result = {}
    for row in qc.itertuples(index=False):
        key = (row.system_id, int(row.replica))
        result[key] = corrected[key] if row.audit_route == "corrected / reviewed evidence" else baseline[key]
    return result


def build_tasks():
    qc = pd.read_csv(QC)
    selections = pd.read_csv(SELECTIONS).set_index("system_id")
    paths = trajectory_paths(qc)
    tasks = []
    for row in qc.itertuples(index=False):
        selection = selections.loc[row.system_id]
        tasks.append({
            "system_id":row.system_id, "replica":int(row.replica), "g_protein_family":row.g_protein_family,
            "trajectory_path":paths[(row.system_id, int(row.replica))], "topology_path":str(selection.topology_path),
            "selection_status":selection.selection_status,
            "selection_reason":selection.reason_code if isinstance(selection.reason_code, str) else "no_canonical_receptor_accession",
            "tm_ordinals":selection.tm_core_residue_ordinals,
            "receptor_ordinals":selection.verified_receptor_residue_ordinals,
            "galpha_ordinals":selection.galpha_candidate_residue_ordinals,
            "identity_status":"distinct",
        })
    assert len(tasks) == 621
    return tasks


def main(repair_systems=None):
    tasks = build_tasks()
    if repair_systems:
        repair_systems = set(repair_systems)
        if not OUT.exists():
            raise FileNotFoundError(f"cannot repair before full output exists: {OUT}")
        current = pd.read_csv(OUT)
        selected = [task for task in tasks if task["system_id"] in repair_systems]
        assert len(selected) == 3 * len(repair_systems)
        replacements = pd.DataFrame([process_replica(task) for task in selected])
        result = pd.concat([current[~current.system_id.isin(repair_systems)], replacements], ignore_index=True)
        result = result.sort_values(["g_protein_family", "system_id", "replica"])
        assert len(result) == 621 and not result.duplicated(["system_id", "replica"]).any()
        unexpected = result.validation_status.ne("available") & ~result.system_id.eq("Gs_8HTI")
        if unexpected.any():
            print(result.loc[unexpected, ["system_id", "replica", "reason_code"]].to_string(index=False))
            raise RuntimeError("repair produced unexpected unavailable replicas")
        result.to_csv(OUT, index=False)
        print(f"Recomputed {len(replacements)} replicas for {sorted(repair_systems)}")
        print(f"Updated {OUT}")
        return
    done = {}
    if PROGRESS.exists():
        previous = pd.read_csv(PROGRESS)
        done = {(r.system_id, int(r.replica)):r._asdict() for r in previous.itertuples(index=False)}
        print(f"Resuming from {len(done)} completed replicas")
    pending = [task for task in tasks if (task["system_id"], task["replica"]) not in done]
    rows = list(done.values())
    # One reader per physical trajectory volume avoids the severe random-I/O
    # contention produced by assigning four simultaneous NetCDF reads to the
    # same disk.  Corrected trajectories under ROOT share the data02 volume.
    def volume(task):
        path = task["trajectory_path"]
        for name in ("data01", "data02", "data03", "data04"):
            if f"/MDdata/{name}/" in path:
                return name
        return "data02"
    executors = {name:ProcessPoolExecutor(max_workers=1) for name in sorted({volume(task) for task in pending})}
    try:
        futures = {executors[volume(task)].submit(process_replica, task):task for task in pending}
        for n, future in enumerate(as_completed(futures), 1):
            rows.append(future.result())
            if n % 12 == 0 or n == len(pending):
                pd.DataFrame(rows).sort_values(["system_id", "replica"]).to_csv(PROGRESS, index=False)
                print(f"  completed {len(rows)}/621 replicas", flush=True)
    finally:
        for executor in executors.values():
            executor.shutdown(wait=True, cancel_futures=True)
    result = pd.DataFrame(rows).sort_values(["g_protein_family", "system_id", "replica"])
    assert len(result) == 621 and not result.duplicated(["system_id", "replica"]).any()
    expected_null = result.system_id.eq("Gs_8HTI")
    unexpected = result.validation_status.ne("available") & ~expected_null
    if unexpected.any():
        print("Unexpected unavailable replicas:")
        print(result.loc[unexpected, ["system_id", "replica", "reason_code"]].to_string(index=False))
        raise RuntimeError("harmonized validation has unexpected unavailable replicas")
    assert expected_null.sum() == 3 and result.loc[expected_null, "validation_status"].eq("unavailable").all()
    available = result[~expected_null]
    numeric = [
        "tm_core_rmsd_A_p95", "galpha_interface_rmsd_A_p95",
        "contact_retention_p05", "contact_retention_median",
    ]
    assert available[numeric].notna().all().all()
    assert available["harmonized_observations"].eq(N_SAMPLES).all()
    assert available["contact_retention_median"].between(0, 1).all()
    result.to_csv(OUT, index=False)
    if PROGRESS.exists():
        PROGRESS.unlink()
    print("Core/interface validation assertions passed:")
    print("  621 selected replica records; 618 quantitative records + 3 explicit Gs_8HTI nulls")
    print("  exactly 1,001 evenly spaced frames per quantitative replica")
    print("  sequence-verified TM core; per-replica initial G-alpha interface and contacts")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--repair-system", action="append", default=[])
    args = parser.parse_args()
    main(args.repair_system)
