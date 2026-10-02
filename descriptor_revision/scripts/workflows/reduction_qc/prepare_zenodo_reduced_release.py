#!/usr/bin/env python3
"""Build, fully audit, and package the private CoupledMD reduced release."""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import time
import warnings
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.coordinates.XTC import XTCWriter


HERE = Path(__file__).resolve().parent
SERVER = HERE.parent.parent
COHORT = HERE / "CoupledMD_Supplementary_Data/Supplementary_Data_S1_included_system_inventory.csv"
VIZ = SERVER / "data/viz"
STAGE = HERE / "zenodo_reduced_release_207"
CHECKSUM_REPORT = HERE / "zenodo_reduced_release_207_upload_checksums.json"
GS_SOURCE = Path("/MDdata/data04/jxhuang/gpcr_g/a/a_gs/3SN6/traj1.nc")
GQ_SOURCE = HERE / "trajectory_readiness/Gq_7RAN_updated_repair_direct/Gq_7RAN/rep1_aligned.nc"
GI8J18_SOURCE = HERE / "trajectory_readiness/continuations/Gi_8J18/candidate_stage2/rep1_aligned.nc"
MANIFEST = STAGE / "CoupledMD_reduced_trajectory_manifest.csv"
AUDIT = STAGE / "CoupledMD_reduced_release_audit.json"
README = STAGE / "README.txt"
LICENSE = STAGE / "LICENSE.txt"
TARGET_FRAMES = 2500
FRAME_STRIDE = 4
EXPECTED_FAMILIES = {"Gi": 95, "Gs": 65, "Gq": 41, "G12-13": 6}
UNRESOLVED = {
    "G12_8H8J", "Gi_7JVR", "Gi_7V68", "Gi_7VUG", "Gi_8J22", "Gi_8X16",
    "Gq_7F9Z", "Gq_7RYC", "Gq_7XJL", "Gq_8DPF", "Gs_7VUH", "Gs_8GY7", "Gs_8HNK",
}
EXCLUDED = {"Gq_7E9W", "Gq_7DWC"}
REBUILT = {"Gs_3SN6": GS_SOURCE, "Gq_7RAN": GQ_SOURCE, "Gi_8J18": GI8J18_SOURCE}
TRIMMED = {"Gi_7YK6", "Gi_8YIC"}
ARCHIVES = {
    "Gi": "CoupledMD_reduced_trajectories_Gi-o.zip",
    "Gs": "CoupledMD_reduced_trajectories_Gs.zip",
    "Gq": "CoupledMD_reduced_trajectories_Gq-11.zip",
    "G12-13": "CoupledMD_reduced_trajectories_G12-13.zip",
}
FIELDS = [
    "system_id", "pdb_id", "g_protein_family", "source_replica",
    "pdb_relative_path", "xtc_relative_path", "n_frames", "frame_interval_ps",
    "represented_span_ns", "n_atoms", "pdb_size_bytes", "xtc_size_bytes",
    "pdb_sha256", "xtc_sha256",
]

# Per-replica evidence used by build_replica_qc.py.  For replicas 2 and 3 these
# corrected trajectories override the default aligned_output when available.
REPLICA_READINESS = HERE / "trajectory_readiness" / "replica_readiness_final.csv"
CORRECTED_SOURCES = [
    (HERE / "trajectory_readiness" / "Gq_7RAN_direct_final_audit" / "corrected_role_replica_audit.csv", None),
    (
        HERE / "trajectory_readiness" / "amber_three_replica_repair_20260713" / "exhaustive_audit_corrected_receptor_fragments" / "candidate_replica_audit.csv",
        {"Gi_7F1Q", "Gi_7QVM", "Gi_7W0N", "Gi_7XA3", "Gi_8IY5", "Gq_8IBV", "Gs_6NBI", "Gs_6WZG", "Gs_7D68", "Gs_7KH0"},
    ),
    (HERE / "trajectory_readiness" / "continuations" / "Gi_8J18" / "interface_recovery" / "exhaustive_audit_corrected" / "candidate_replica_audit.csv", None),
    (HERE / "trajectory_readiness" / "continuations" / "Gq_8ZPT" / "interface_recovery" / "exhaustive_audit_corrected" / "candidate_replica_audit.csv", None),
    (HERE / "trajectory_readiness" / "continuous_review_audit" / "continuous_review_replica_audit_with_rmsd.csv", None),
    (HERE / "trajectory_readiness" / "corrected_role_audit" / "corrected_role_replica_audit.csv", None),
    (HERE / "trajectory_readiness" / "interface_centered_pilot_20260713" / "Gi_6WWZ" / "three_replica_audit_corrected" / "candidate_replica_audit.csv", None),
    (HERE / "trajectory_readiness" / "interface_centered_rescue_20260713" / "Gi_7XA3" / "audit" / "candidate_replica_audit.csv", None),
    (HERE / "trajectory_readiness" / "interface_centered_rescue_20260713" / "Gq_8IBV" / "audit" / "candidate_replica_audit.csv", None),
]

# Runtime globals set in main() from command-line arguments.
REPLICA = 1
DRAFT_ID = 21395292
SOURCE_LABELS: dict[str, str] = {}

warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_cohort() -> list[dict[str, str]]:
    with COHORT.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    ids = [row["system_id"] for row in rows]
    assert len(rows) == len(set(ids)) == 207
    assert Counter(row["g_protein_family"] for row in rows) == Counter(EXPECTED_FAMILIES)
    assert not set(ids) & (UNRESOLVED | EXCLUDED)
    return sorted(rows, key=lambda row: row["system_id"])


def hashes(path: Path) -> tuple[str, str]:
    sha = hashlib.sha256()
    md5 = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            sha.update(block); md5.update(block)
    return sha.hexdigest(), md5.hexdigest()


def sha256(path: Path) -> str:
    return hashes(path)[0]


def universe(pdb: Path, trajectory: Path) -> mda.Universe:
    return mda.Universe(str(pdb), str(trajectory))


def remove_offset_caches(directory: Path) -> None:
    for path in directory.rglob(".*_offsets.npz"):
        path.unlink()
    for path in directory.rglob(".*_offsets.lock"):
        path.unlink()


def corrected_source_map(replica: int, rows: list[dict[str, str]]) -> dict[str, Path]:
    """Return system_id -> corrected trajectory path for this replica."""
    allowed_sids = {row["system_id"] for row in rows}
    sources: dict[str, Path] = {}
    for path, allow_list in CORRECTED_SOURCES:
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("verdict") != "GOOD":
                    continue
                sid = row["system_id"]
                if sid not in allowed_sids:
                    continue
                if int(row["replica"]) != replica:
                    continue
                if allow_list is not None and sid not in allow_list:
                    continue
                if sid in sources:
                    raise ValueError(f"duplicate corrected source for {sid} replica {replica}")
                tp = row["trajectory_path"]
                src = Path(tp)
                if not src.is_absolute():
                    src = HERE / src
                sources[sid] = src
    return sources


def aligned_output_map(replica: int, rows: list[dict[str, str]]) -> dict[str, Path]:
    """Return system_id -> aligned_output path from replica_readiness_final.csv."""
    allowed_sids = {row["system_id"] for row in rows}
    by_key: dict[tuple[str, int], dict[str, str]] = {}
    with REPLICA_READINESS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (row["system_id"], int(row["replica"]))
            by_key[key] = row
    sources: dict[str, Path] = {}
    for row in rows:
        sid = row["system_id"]
        rec = by_key.get((sid, replica))
        if rec is None:
            raise ValueError(f"no replica {replica} readiness record for {sid}")
        sources[sid] = Path(rec["aligned_output"])
    return sources


def source_map_for_replica(replica: int, rows: list[dict[str, str]]) -> dict[str, tuple[Path, str]]:
    """Return system_id -> (source trajectory path, source label)."""
    corrected = corrected_source_map(replica, rows)
    aligned = aligned_output_map(replica, rows)
    combined: dict[str, tuple[Path, str]] = {}
    for row in rows:
        sid = row["system_id"]
        if sid in corrected:
            combined[sid] = (corrected[sid], "corrected")
        else:
            combined[sid] = (aligned[sid], "aligned_output")
    return combined


def choose_indices(source_path: Path, pdb_path: Path, sid: str, source_top: Path | None = None) -> tuple[list[int], str]:
    """Choose frame indices that yield 2,500 representative frames."""
    top = source_top if source_top is not None else pdb_path
    u = universe(top, source_path)
    n = len(u.trajectory)
    if n == TARGET_FRAMES:
        return list(range(TARGET_FRAMES)), "preserve"
    # Determine the timestep from the first two frames only, to avoid re-reading
    # the entire source trajectory before the conversion pass.
    if n > 1:
        dt = float(u.trajectory[1].time) - float(u.trajectory[0].time)
    else:
        dt = 0.0
    # Raw 50-ps trajectories (10,000 or 10,001 frames): stride by 4.
    if n >= TARGET_FRAMES * 4 and 45.0 <= dt <= 55.0:
        return list(range(0, TARGET_FRAMES * 4, FRAME_STRIDE)), "stride4"
    # Pre-strided 200-ps trajectories with one extra frame: keep the first 2,500.
    if n == TARGET_FRAMES + 1 and 190.0 <= dt <= 210.0:
        return list(range(TARGET_FRAMES)), "trim"
    # Legacy replica-1 special cases already handled by explicit rebuild maps.
    raise ValueError(
        f"{sid}: cannot derive {TARGET_FRAMES} frames from {n} frames (dt={dt:.2f} ps)"
    )


def find_source_topology(source_path: Path) -> Path | None:
    """Return a topology file in the source directory that matches the trajectory."""
    candidates = [
        "protein.pdb",
        "now.prmtop",
        "step5_input.parm7",
        "step5_input.psf",
        "system.parm7",
        "system.prmtop",
        "topol.tpr",
        "prot.pdb",
    ]
    for name in candidates:
        top = source_path.parent / name
        if top.exists():
            try:
                u = mda.Universe(str(top), str(source_path))
            except Exception:
                continue
            if u.atoms.n_atoms > 0:
                return top
    return None


def subset_atom_indices(source_u: mda.Universe, target_u: mda.Universe) -> list[int]:
    """Map target atoms to source atoms by residue name, residue id, and atom name."""
    by_key: dict[tuple[str, int, str], list[int]] = {}
    for idx, atom in enumerate(source_u.atoms):
        key = (atom.resname, atom.resid, atom.name)
        by_key.setdefault(key, []).append(idx)
    indices = []
    for atom in target_u.atoms:
        key = (atom.resname, atom.resid, atom.name)
        if key not in by_key or not by_key[key]:
            raise ValueError(f"cannot map target atom {key} into source topology")
        indices.append(by_key[key].pop(0))
    if len(indices) != target_u.atoms.n_atoms:
        raise ValueError("subset mapping length mismatch")
    return indices


def write_selected_subset(pdb: Path, source: Path, source_top: Path, target: Path, indices) -> None:
    """Write a reduced trajectory using only the source atoms that match the target PDB."""
    source_u = mda.Universe(str(source_top), str(source))
    target_u = mda.Universe(str(pdb))
    subset_indices = subset_atom_indices(source_u, target_u)
    subset = source_u.atoms[subset_indices]
    assert subset.n_atoms == target_u.atoms.n_atoms
    chosen = list(indices)
    assert len(chosen) == TARGET_FRAMES and chosen[-1] < len(source_u.trajectory)
    partial = target.with_suffix(".partial.xtc")
    if partial.exists():
        partial.unlink()
    first_time = last_time = None
    with XTCWriter(str(partial), n_atoms=subset.n_atoms) as writer:
        for count, frame_index in enumerate(chosen, 1):
            ts = source_u.trajectory[frame_index]
            assert np.isfinite(ts.positions).all()
            assert ts.dimensions is not None and np.isfinite(ts.dimensions).all()
            writer.write(subset)
            first_time = float(ts.time) if first_time is None else first_time
            last_time = float(ts.time)
            if count % 500 == 0:
                print(f"  wrote {count}/{TARGET_FRAMES} frames", flush=True)
    partial.replace(target)
    check = universe(pdb, target)
    assert len(check.trajectory) == TARGET_FRAMES and check.atoms.n_atoms == target_u.atoms.n_atoms
    assert first_time is not None and last_time is not None and 499.0 <= (last_time - first_time) / 1000 <= 500.1


def image_xtc(pdb: Path, xtc: Path) -> None:
    """Center and image a reduced XTC trajectory using cpptraj.

    Replica 2/3 source trajectories are raw production outputs that may cross
    periodic boundaries; imaging keeps the protein complex whole and centered
    so that downstream QC coordinate-scatter checks pass.
    """
    tmp = xtc.with_suffix(".imaged.xtc")
    if tmp.exists():
        tmp.unlink()
    script = f"""parm {pdb}
trajin {xtc}
center origin mass
image origin center
trajout {tmp}
"""
    try:
        result = subprocess.run(
            ["cpptraj"],
            input=script,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"cpptraj imaging failed for {xtc}: {exc.stdout}") from exc
    if not tmp.exists():
        raise RuntimeError(f"cpptraj did not create {tmp} for {xtc}\n{result.stdout}")
    tmp.replace(xtc)


def inventory_canonical(rows: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for index, row in enumerate(rows, 1):
        sid = row["system_id"]
        pdb = VIZ / sid / "structure.pdb"; xtc = VIZ / sid / "traj.xtc"
        assert pdb.is_file() and xtc.is_file()
        u = universe(pdb, xtc)
        counts[sid] = len(u.trajectory)
        if index % 25 == 0 or index == len(rows):
            print(f"canonical inventory [{index:03d}/{len(rows)}]", flush=True)
    distribution = Counter(counts.values())
    expected = Counter({2500: 203, 2501: 2, 200: 1, 681: 1})
    assert distribution == expected, (distribution, expected)
    assert {sid for sid, n in counts.items() if n == 2501} == TRIMMED
    assert counts["Gs_3SN6"] == 200 and counts["Gq_7RAN"] == 681
    return counts


def write_selected(pdb: Path, source: Path, target: Path, indices) -> None:
    u = universe(pdb, source)
    assert u.atoms.n_atoms == universe(pdb, pdb).atoms.n_atoms
    chosen = list(indices)
    assert len(chosen) == TARGET_FRAMES and chosen[-1] < len(u.trajectory)
    partial = target.with_suffix(".partial.xtc")
    if partial.exists():
        partial.unlink()
    first_time = last_time = None
    with XTCWriter(str(partial), n_atoms=u.atoms.n_atoms) as writer:
        for count, frame_index in enumerate(chosen, 1):
            ts = u.trajectory[frame_index]
            assert np.isfinite(ts.positions).all()
            assert ts.dimensions is not None and np.isfinite(ts.dimensions).all()
            writer.write(u.atoms)
            first_time = float(ts.time) if first_time is None else first_time
            last_time = float(ts.time)
            if count % 500 == 0:
                print(f"  wrote {count}/{TARGET_FRAMES} frames", flush=True)
    partial.replace(target)
    check = universe(pdb, target)
    assert len(check.trajectory) == TARGET_FRAMES and check.atoms.n_atoms == u.atoms.n_atoms
    assert first_time is not None and last_time is not None and 499.0 <= (last_time - first_time) / 1000 <= 500.1


def stage_release(rows: list[dict[str, str]], resume: bool) -> None:
    if STAGE.exists() and not resume:
        raise SystemExit(f"staging directory exists; use --resume: {STAGE}")
    STAGE.mkdir(exist_ok=True)
    counts = inventory_canonical(rows)
    for index, row in enumerate(rows, 1):
        sid = row["system_id"]
        source_dir = VIZ / sid; target_dir = STAGE / sid
        target_dir.mkdir(exist_ok=True)
        pdb_target = target_dir / "structure.pdb"; xtc_target = target_dir / "traj.xtc"
        if not pdb_target.exists():
            shutil.copy2(source_dir / "structure.pdb", pdb_target)
        if xtc_target.exists():
            u = universe(pdb_target, xtc_target)
            if len(u.trajectory) == TARGET_FRAMES:
                print(f"stage [{index:03d}/207] {sid} resume", flush=True)
                continue
            xtc_target.unlink()
        if sid in REBUILT:
            print(f"stage [{index:03d}/207] {sid} rebuild replica 1", flush=True)
            write_selected(pdb_target, REBUILT[sid], xtc_target, range(0, 10000, FRAME_STRIDE))
        elif sid in TRIMMED:
            assert counts[sid] == 2501
            print(f"stage [{index:03d}/207] {sid} trim 2501 to 2500", flush=True)
            write_selected(pdb_target, source_dir / "traj.xtc", xtc_target, range(TARGET_FRAMES))
        else:
            assert counts[sid] == TARGET_FRAMES
            print(f"stage [{index:03d}/207] {sid} preserve canonical", flush=True)
            shutil.copy2(source_dir / "traj.xtc", xtc_target)
    remove_offset_caches(STAGE)
    dirs = {p.name for p in STAGE.iterdir() if p.is_dir()}
    assert dirs == {row["system_id"] for row in rows}
    assert sum(1 for p in STAGE.rglob("structure.pdb")) == 207
    assert sum(1 for p in STAGE.rglob("traj.xtc")) == 207


def _stage_one_replica(args: tuple[int, dict[str, str], Path, str, int]) -> tuple[int, str, str]:
    """Worker that stages one system for a replica. Returns (index, sid, status)."""
    index, row, source_path, source_label, replica = args
    sid = row["system_id"]
    target_dir = STAGE / sid
    target_dir.mkdir(exist_ok=True)
    pdb_target = target_dir / "structure.pdb"
    xtc_target = target_dir / "traj.xtc"
    if not pdb_target.exists():
        shutil.copy2(VIZ / sid / "structure.pdb", pdb_target)
    if xtc_target.exists():
        try:
            u = universe(pdb_target, xtc_target)
            if len(u.trajectory) == TARGET_FRAMES:
                return index, sid, "resume"
        except Exception:
            pass
        xtc_target.unlink()
    try:
        indices, _action = choose_indices(source_path, pdb_target, sid)
        write_selected(pdb_target, source_path, xtc_target, indices)
    except ValueError as exc:
        if "n_atoms" not in str(exc).lower() and "natom" not in str(exc).lower():
            raise
        source_top = find_source_topology(source_path)
        if source_top is None:
            raise ValueError(f"{sid}: atom-count mismatch and no source topology found for {source_path}")
        indices, _action = choose_indices(source_path, pdb_target, sid, source_top=source_top)
        write_selected_subset(pdb_target, source_path, source_top, xtc_target, indices)
    if replica != 1:
        image_xtc(pdb_target, xtc_target)
    return index, sid, source_label


def stage_release_replica(rows: list[dict[str, str]], replica: int, resume: bool) -> None:
    if STAGE.exists() and not resume:
        raise SystemExit(f"staging directory exists; use --resume: {STAGE}")
    STAGE.mkdir(exist_ok=True)
    SOURCE_LABELS.clear()
    sources = source_map_for_replica(replica, rows)
    tasks = []
    for index, row in enumerate(rows, 1):
        sid = row["system_id"]
        source_path, source_label = sources[sid]
        SOURCE_LABELS[sid] = source_label
        if not source_path.exists():
            raise FileNotFoundError(f"{sid} replica {replica}: source missing: {source_path}")
        tasks.append((index, row, source_path, source_label, replica))
    workers = max(1, int(os.environ.get("STAGE_WORKERS", "4")))
    print(f"staging replica {replica} with {workers} workers", flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for index, sid, status in executor.map(_stage_one_replica, tasks):
            print(f"stage [{index:03d}/{len(rows)}] {sid} {status}", flush=True)
    remove_offset_caches(STAGE)
    dirs = {p.name for p in STAGE.iterdir() if p.is_dir()}
    assert dirs == {row["system_id"] for row in rows}
    assert sum(1 for p in STAGE.rglob("structure.pdb")) == 207
    assert sum(1 for p in STAGE.rglob("traj.xtc")) == 207


def chain_groups(u: mda.Universe) -> list[np.ndarray]:
    try:
        labels = np.asarray(u.atoms.chainIDs, dtype=str)
    except Exception:
        labels = np.asarray(u.atoms.segids, dtype=str)
    if not any(label.strip() for label in labels):
        labels = np.asarray(u.atoms.segids, dtype=str)
    groups = [np.where(labels == label)[0] for label in dict.fromkeys(labels) if label.strip()]
    if not groups:
        groups = [np.arange(u.atoms.n_atoms)]
    return groups


def structural_pairs(u: mda.Universe) -> tuple[np.ndarray, np.ndarray]:
    backbone = []
    ca_pairs = []
    previous_ca: dict[str, int] = {}
    for residue in u.residues:
        by_name = {atom.name: atom.index for atom in residue.atoms}
        if "N" in by_name and "CA" in by_name:
            backbone.append((by_name["N"], by_name["CA"]))
        if "CA" in by_name and "C" in by_name:
            backbone.append((by_name["CA"], by_name["C"]))
        try:
            chain = str(residue.atoms.chainIDs[0]).strip()
        except Exception:
            chain = str(residue.atoms.segids[0]).strip()
        current_ca = by_name.get("CA")
        if current_ca is not None and chain in previous_ca:
            ca_pairs.append((previous_ca[chain], current_ca))
        if current_ca is not None:
            previous_ca[chain] = current_ca
    return np.asarray(backbone, dtype=int), np.asarray(ca_pairs, dtype=int)


def qc_one(row: dict[str, str]) -> tuple[dict[str, str], dict[str, object]]:
    sid = row["system_id"]; directory = STAGE / sid
    pdb = directory / "structure.pdb"; xtc = directory / "traj.xtc"
    # MDAnalysis may create private XTC offset/lock caches while concurrent
    # workers inspect neighbouring systems.  Permit only those transient
    # reader files here; remove them globally after all workers have closed.
    unexpected = {
        p.name for p in directory.iterdir()
        if p.name not in {"structure.pdb", "traj.xtc"}
        and not (p.name.startswith(".traj.xtc_offsets") and p.suffix in {".npz", ".lock"})
    }
    assert not unexpected, (sid, sorted(unexpected))
    u_pdb = mda.Universe(str(pdb)); u = universe(pdb, xtc)
    assert u.atoms.n_atoms == u_pdb.atoms.n_atoms
    assert len(u.trajectory) == TARGET_FRAMES
    groups = chain_groups(u)
    anchor = max(groups, key=len)
    backbone, ca_pairs = structural_pairs(u)
    first_pdb = u_pdb.atoms.positions.copy()
    backbone_total = len(backbone)
    if backbone.size:
        base_backbone = np.linalg.norm(first_pdb[backbone[:, 0]] - first_pdb[backbone[:, 1]], axis=1)
        # Use only covalently plausible reference pairs for the atom-order test.
        # A small number of canonical PDBs contain pre-existing local coordinate
        # gaps; those are not evidence that XTC atom order is wrong.
        backbone = backbone[base_backbone < 2.0]
    if ca_pairs.size:
        base_ca = np.linalg.norm(first_pdb[ca_pairs[:, 0]] - first_pdb[ca_pairs[:, 1]], axis=1)
        ca_pairs = ca_pairs[base_ca < 5.0]
    times = []
    failures = Counter()
    max_span_ratio = 0.0; max_backbone_distance = 0.0
    for frame_number, ts in enumerate(u.trajectory):
        pos = ts.positions
        if not np.isfinite(pos).all():
            failures["nonfinite_coordinates"] += 1; continue
        time_ps = float(ts.time)
        if not math.isfinite(time_ps): failures["nonfinite_time"] += 1
        times.append(time_ps)
        box = ts.dimensions
        if box is None or len(box) < 6 or not np.isfinite(box).all() or np.any(box[:3] <= 0) or np.any(box[3:] <= 0) or np.any(box[3:] > 180):
            failures["invalid_box"] += 1; continue
        box3 = np.asarray(box[:3], dtype=float)
        span_ratio = float(np.max(np.ptp(pos, axis=0) / box3))
        max_span_ratio = max(max_span_ratio, span_ratio)
        if span_ratio > 1.6: failures["coordinate_scatter"] += 1
        if backbone.size:
            delta = pos[backbone[:, 1]] - pos[backbone[:, 0]]
            delta -= box3 * np.round(delta / box3)
            distances = np.linalg.norm(delta, axis=1)
            max_backbone_distance = max(max_backbone_distance, float(distances.max()))
            if np.any(distances > 3.0): failures["atom_order_or_backbone_geometry"] += 1
        if ca_pairs.size:
            delta_raw = pos[ca_pairs[:, 1]] - pos[ca_pairs[:, 0]]
            raw = np.linalg.norm(delta_raw, axis=1)
            delta_mi = delta_raw - box3 * np.round(delta_raw / box3)
            minimum = np.linalg.norm(delta_mi, axis=1)
            if np.any((raw > 8.0) & (minimum <= 8.0)): failures["catastrophic_chain_break"] += 1
        anchor_center = pos[anchor].mean(axis=0)
        for group in groups:
            if group is anchor: continue
            delta = pos[group].mean(axis=0) - anchor_center
            delta -= box3 * np.round(delta / box3)
            if np.linalg.norm(delta) > float(box3.min()):
                failures["complex_separation"] += 1; break
        if (frame_number + 1) % 500 == 0:
            pass
    intervals = np.diff(np.asarray(times, dtype=float))
    if len(times) != TARGET_FRAMES: failures["time_count"] += 1
    if intervals.size != TARGET_FRAMES - 1 or not np.isfinite(intervals).all() or np.any(intervals <= 0):
        failures["nonmonotonic_time"] += 1
    assert not failures, (sid, dict(failures))
    interval = float(statistics.median(intervals))
    span_ns = (times[-1] - times[0]) / 1000.0
    assert 190 <= interval <= 210 and 499.0 <= span_ns <= 500.1
    pdb_sha = sha256(pdb); xtc_sha = sha256(xtc)
    if REPLICA == 1:
        release_action = "rebuilt_from_replica_1" if sid in REBUILT else ("trimmed_to_2500" if sid in TRIMMED else "canonical_preserved")
    else:
        release_action = SOURCE_LABELS.get(sid, f"selected_from_replica_{REPLICA}")
    manifest = {
        "system_id": sid, "pdb_id": row["pdb_id"], "g_protein_family": row["g_protein_family"],
        "source_replica": str(REPLICA), "pdb_relative_path": f"{sid}/structure.pdb",
        "xtc_relative_path": f"{sid}/traj.xtc", "n_frames": str(TARGET_FRAMES),
        "frame_interval_ps": f"{interval:.3f}", "represented_span_ns": f"{span_ns:.3f}",
        "n_atoms": str(u.atoms.n_atoms), "pdb_size_bytes": str(pdb.stat().st_size),
        "xtc_size_bytes": str(xtc.stat().st_size), "pdb_sha256": pdb_sha, "xtc_sha256": xtc_sha,
    }
    detail = {
        "system_id": sid, "frames_checked": TARGET_FRAMES, "n_atoms": u.atoms.n_atoms,
        "frame_interval_ps": interval, "represented_span_ns": span_ns,
        "max_coordinate_span_over_box": max_span_ratio,
        "max_backbone_bond_distance_angstrom": max_backbone_distance,
        "reference_backbone_pairs_total": backbone_total,
        "reference_backbone_pairs_checked": len(backbone),
        "reference_backbone_pairs_excluded_as_baseline_gaps": backbone_total - len(backbone),
        "finite_coordinates": True, "finite_monotonic_time": True, "valid_periodic_box": True,
        "atom_order_backbone_geometry": "passed", "catastrophic_chain_break": False,
        "coordinate_scatter": False, "complex_separation": False,
        "source_replica": REPLICA,
        "release_action": release_action,
    }
    return manifest, detail


def time_check_one(row: dict[str, str]) -> dict[str, object]:
    sid = row["system_id"]; directory = STAGE / sid
    u = universe(directory / "structure.pdb", directory / "traj.xtc")
    times = np.fromiter((float(ts.time) for ts in u.trajectory), dtype=float)
    bad = np.where(np.diff(times) <= 0)[0]
    return {
        "system_id": sid, "n_frames": len(times), "nonfinite": int((~np.isfinite(times)).sum()),
        "nonmonotonic_steps": int(len(bad)),
        "first_bad_step": None if not len(bad) else {"frame": int(bad[0]), "time_before_ps": float(times[bad[0]]), "time_after_ps": float(times[bad[0] + 1])},
    }


def time_preflight(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    remove_offset_caches(STAGE)
    workers = max(1, int(os.environ.get("QC_WORKERS", "4")))
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(time_check_one, rows))
    remove_offset_caches(STAGE)
    failures = [row for row in results if row["nonfinite"] or row["nonmonotonic_steps"] or row["n_frames"] != TARGET_FRAMES]
    print(json.dumps({"systems": len(results), "failures": failures}, indent=2, sort_keys=True))
    return failures


def write_docs() -> None:
    README.write_text(
        f"CoupledMD reduced trajectory release (replica {REPLICA})\n\n"
        f"This release contains one reduced protein-complex trajectory per included system, derived deterministically from replica {REPLICA}. "
        "Each trajectory contains 2,500 frames at approximately 200-ps spacing and represents approximately 500 ns.\n\n"
        "Membrane lipids, solvent and mobile ions are excluded. These reduced records are unsuitable for lipid, solvent, ion or fast-timescale analyses. "
        "The underlying three-replica full-system trajectories are not included.\n\n"
        "Each system directory contains structure.pdb and traj.xtc. CoupledMD_reduced_trajectory_manifest.csv records the system metadata, frame sampling, "
        "atom and byte counts, relative paths and SHA-256 checksums. CoupledMD_reduced_release_audit.json summarizes full-frame release QC.\n",
        encoding="utf-8",
    )
    LICENSE.write_text(
        "CoupledMD reduced trajectory data and accompanying documentation are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0):\n"
        "https://creativecommons.org/licenses/by/4.0/\n\n"
        "You may share and adapt the material for any purpose provided appropriate credit is given, a link to the licence is supplied, and changes are indicated.\n\n"
        "Source PDB coordinate records remain subject to the wwPDB Creative Commons CC0 dedication. Software code is distributed separately under the MIT License.\n",
        encoding="utf-8",
    )


def full_qc(rows: list[dict[str, str]]) -> dict[str, object]:
    remove_offset_caches(STAGE)
    expected_dirs = {row["system_id"] for row in rows}
    actual_dirs = {p.name for p in STAGE.iterdir() if p.is_dir()}
    assert actual_dirs == expected_dirs and not actual_dirs & (UNRESOLVED | EXCLUDED)
    manifests = []; details = []
    started = time.time()
    workers = max(1, int(os.environ.get("QC_WORKERS", "4")))
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for index, (row, result) in enumerate(zip(rows, executor.map(qc_one, rows)), 1):
            manifest, detail = result
            manifests.append(manifest); details.append(detail)
            print(f"full QC [{index:03d}/207] {row['system_id']} PASS ({workers} workers)", flush=True)
    with MANIFEST.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader(); writer.writerows(manifests)
    write_docs()
    remove_offset_caches(STAGE)
    assert len(manifests) == len({r["system_id"] for r in manifests}) == 207
    assert Counter(int(r["n_frames"]) for r in manifests) == Counter({2500: 207})
    assert Counter(r["g_protein_family"] for r in manifests) == Counter(EXPECTED_FAMILIES)
    molecular = [p for p in STAGE.rglob("*") if p.is_file() and p.name in {"structure.pdb", "traj.xtc"}]
    assert len(molecular) == 414
    result = {
        "schema_version": "1.0", "status": "passed", "generated_at": now(),
        "cohort_systems": 207, "system_directories": 207, "molecular_files": 414,
        "family_counts": EXPECTED_FAMILIES, "frame_count_distribution": {"2500": 207},
        "source_replica_distribution": {str(REPLICA): 207},
        "release_actions": dict(Counter(d["release_action"] for d in details)),
        "full_frames_checked": 207 * TARGET_FRAMES,
        "checks": {
            "cohort_exact_no_extras": "passed", "excluded_and_unresolved_absent": "passed",
            "pdb_xtc_atom_counts_match": "passed", "finite_coordinates_all_frames": "passed",
            "finite_monotonic_time": "passed", "valid_periodic_box": "passed",
            "atom_order_backbone_geometry": "passed", "catastrophic_chain_break": "passed",
            "coordinate_scatter": "passed", "complex_separation": "passed",
            f"source_replica_{REPLICA}": "passed", "internal_paths_absent": "passed",
        },
        "sampling": {"frames_per_system": 2500, "nominal_interval_ps": 200, "nominal_span_ns": 500},
        "qc_elapsed_seconds": round(time.time() - started, 1),
        "system_qc": details,
        "zenodo": {"draft_id": DRAFT_ID, "upload_status": "not_started", "published": False},
    }
    AUDIT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def load_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 207
    return rows


def build_archives() -> dict[str, dict[str, object]]:
    rows = load_manifest(); by_family = defaultdict(list)
    for row in rows: by_family[row["g_protein_family"]].append(row)
    metadata = {}
    for family, filename in ARCHIVES.items():
        path = STAGE / filename; partial = path.with_suffix(".partial.zip")
        expected = sorted(r[key] for r in by_family[family] for key in ("pdb_relative_path", "xtc_relative_path"))
        if path.exists():
            with zipfile.ZipFile(path, "r", allowZip64=True) as archive:
                if sorted(archive.namelist()) == expected and archive.testzip() is None and all(info.compress_type == zipfile.ZIP_STORED for info in archive.infolist()):
                    sha, md5 = hashes(path)
                    metadata[family] = {"filename": filename, "systems": len(by_family[family]), "files": len(expected), "size_bytes": path.stat().st_size, "sha256": sha, "md5": md5, "zip_test": "passed", "compression": "stored"}
                    print(f"archive {family}: reused existing {filename} ({len(by_family[family])} systems)", flush=True)
                    continue
        if partial.exists(): partial.unlink()
        print(f"archive {family}: {len(by_family[family])} systems", flush=True)
        with zipfile.ZipFile(partial, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as archive:
            for index, row in enumerate(sorted(by_family[family], key=lambda x: x["system_id"]), 1):
                for key in ("pdb_relative_path", "xtc_relative_path"):
                    relative = row[key]; archive.write(STAGE / relative, arcname=relative)
                if index % 20 == 0: print(f"  archived {index}/{len(by_family[family])}", flush=True)
        partial.replace(path)
        with zipfile.ZipFile(path, "r", allowZip64=True) as archive:
            assert sorted(archive.namelist()) == expected
            assert archive.testzip() is None
            assert all(info.compress_type == zipfile.ZIP_STORED for info in archive.infolist())
        sha, md5 = hashes(path)
        metadata[family] = {"filename": filename, "systems": len(by_family[family]), "files": len(expected), "size_bytes": path.stat().st_size, "sha256": sha, "md5": md5, "zip_test": "passed", "compression": "stored"}
        print(f"  verified {filename} {path.stat().st_size} bytes", flush=True)
    return metadata


def finalize_audit(archive_metadata: dict[str, dict[str, object]]) -> dict[str, object]:
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    audit["archives"] = archive_metadata
    audit["archive_system_total"] = sum(int(v["systems"]) for v in archive_metadata.values())
    audit["release_size_bytes"] = sum((STAGE / row["xtc_relative_path"]).stat().st_size + (STAGE / row["pdb_relative_path"]).stat().st_size for row in load_manifest())
    audit["upload_file_count"] = 8
    audit["upload_files_expected"] = [*ARCHIVES.values(), MANIFEST.name, AUDIT.name, README.name, LICENSE.name]
    audit["upload_files_expected"] = sorted(audit["upload_files_expected"])
    audit["zenodo"] = {"draft_id": DRAFT_ID, "upload_status": "not_started_authentication_unavailable", "published": False, "private_draft_required": True}
    AUDIT.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    upload_files = [STAGE / name for name in audit["upload_files_expected"]]
    assert len(upload_files) == 8 and all(path.is_file() for path in upload_files)
    checksums = {}
    for path in upload_files:
        sha, md5 = hashes(path)
        checksums[path.name] = {"size_bytes": path.stat().st_size, "sha256": sha, "md5": md5}
    CHECKSUM_REPORT.write_text(json.dumps({"generated_at": now(), "upload_file_count": 8, "files": checksums}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    forbidden = ("/MDdata/", "_orig", "_corrected", ".partial", "backup")
    for path in [MANIFEST, AUDIT, README, LICENSE]:
        text = path.read_text(encoding="utf-8")
        assert not any(token in text for token in forbidden), (path, forbidden)
    return checksums


def assert_clean_tree(rows: list[dict[str, str]]) -> None:
    system_ids = {row["system_id"] for row in rows}
    dirs = {p.name for p in STAGE.iterdir() if p.is_dir()}
    assert dirs == system_ids
    for sid in system_ids:
        assert {p.name for p in (STAGE / sid).iterdir()} == {"structure.pdb", "traj.xtc"}
    allowed_root = {MANIFEST.name, AUDIT.name, README.name, LICENSE.name, *ARCHIVES.values()}
    root_files = {p.name for p in STAGE.iterdir() if p.is_file()}
    assert root_files == allowed_root, (root_files - allowed_root, allowed_root - root_files)


def main() -> None:
    global REPLICA, DRAFT_ID, STAGE, CHECKSUM_REPORT, MANIFEST, AUDIT, README, LICENSE
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("stage", "timecheck", "qc", "package", "all"))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--replica", type=int, choices=(1, 2, 3), default=1)
    parser.add_argument("--draft-id", type=int, default=None)
    args = parser.parse_args()
    REPLICA = args.replica
    if args.draft_id is not None:
        DRAFT_ID = args.draft_id
    STAGE = HERE / f"zenodo_reduced_release_207_replica{REPLICA}"
    CHECKSUM_REPORT = HERE / f"zenodo_reduced_release_207_replica{REPLICA}_upload_checksums.json"
    MANIFEST = STAGE / "CoupledMD_reduced_trajectory_manifest.csv"
    AUDIT = STAGE / "CoupledMD_reduced_release_audit.json"
    README = STAGE / "README.txt"
    LICENSE = STAGE / "LICENSE.txt"
    rows = read_cohort()
    if args.action in {"stage", "all"}:
        if REPLICA == 1:
            stage_release(rows, args.resume)
        else:
            stage_release_replica(rows, REPLICA, args.resume)
    if args.action == "timecheck": time_preflight(rows)
    if args.action in {"qc", "all"}: full_qc(rows)
    if args.action in {"package", "all"}:
        archives = build_archives(); checksums = finalize_audit(archives); assert_clean_tree(rows)
        print(json.dumps({"status": "passed", "staging": str(STAGE), "upload_file_count": len(checksums), "checksums": checksums}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
