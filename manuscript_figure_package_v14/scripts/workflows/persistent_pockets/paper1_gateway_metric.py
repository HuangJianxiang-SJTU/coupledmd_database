#!/usr/bin/env python3
"""
paper1_gateway_metric.py  --  Lateral-gateway lipid-penetration metric (Paper 1).

A lateral gateway is a membrane-facing opening between two adjacent TM helices
through which bilayer lipid intrudes toward the receptor core. Because this
dataset retains the POPC bilayer in prod*_now.nc, gateways are measured DIRECTLY
from lipid penetration rather than inferred from protein geometry alone -- the
key advantage over GPCRmd-based work.

For each of the 7 membrane-facing inter-helical portals (consecutive TM pairs
around the bundle), per frame:
  - wedged lipid = lipid heavy atoms within DCUT of BOTH helices, inside the TM
    z-band (so headgroups above/below the bilayer core are excluded).
  - penetration depth = how far the deepest wedged lipid intrudes past the helix
    backbone toward the bundle axis (radial), i.e. r_wall - min(r_lipid).
Aggregated over frames and the 3 replicas (replica-level bootstrap CI per the
project convention; never per-frame bootstrap).

Assumes membrane normal = z (CHARMM-GUI builds the bilayer in xy; center+wrap
reimaging does not rotate, so z is preserved). Receptor TM bundle axis = z
through the xy centre-of-mass of all TM Cα.

NOTE: heavy compute -- run only when CPU is free. Validate single-frame first.
"""
import os, sys, glob, csv, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from MDAnalysis import transformations as trans
from MDAnalysis.lib.distances import distance_array

from difflib import SequenceMatcher
from build_system_config import (load_or_fetch_gpcrdb, gpcrdb_segment_pror_ranges,
                                  detect_chains_by_ca_breaks, assign_chains,
                                  resnames_to_seq)
from paper1_pocket_pilot import build_now_ref, AA

# Corrected chain roles (RCSB-backed, 2026-06-25). Preferred over assign_chains().
CHAIN_ROLES_DIR = "/MDdata/data02/jxhuang/gpcr_g/gpcr_g_server/scripts/audit_output/chain_roles"


def _gw_parse_range(rng):
    if not rng or "-" not in str(rng):
        return None
    parts = str(rng).split("-")
    try:
        return int(parts[0]), int(parts[1])
    except (ValueError, IndexError):
        return None


def load_corrected_roles(sid):
    """Corrected chain_roles.json -> {receptor,galpha,...:(s,e)}. None if absent."""
    path = os.path.join(CHAIN_ROLES_DIR, f"{sid}_chain_roles.json")
    if not os.path.exists(path):
        return None
    data = json.load(open(path))
    roles = {}
    rec_min = rec_max = None
    for role, info in data.get("roles", {}).items():
        rng = _gw_parse_range(info.get("resid_range"))
        if rng is None:
            continue
        s, e = rng
        if role.startswith("Receptor"):
            rec_min = s if rec_min is None else min(rec_min, s)
            rec_max = e if rec_max is None else max(rec_max, e)
        elif role == "G_alpha":
            roles["galpha"] = (s, e)
        elif role == "G_beta":
            roles["gbeta"] = (s, e)
        elif role == "G_gamma":
            roles["g_gamma"] = (s, e)
        elif role == "peptide_ligand":
            roles["peptide"] = (s, e)
    if rec_min is not None:
        roles["receptor"] = (rec_min, rec_max)
    return roles

INVENTORY = "/MDdata/data02/jxhuang/gpcr_g/a/inventory_v3_extended.csv"
OUTDIR = "/MDdata/data02/jxhuang/gpcr_g/a/paper1_gateways"
GCACHE = "/MDdata/data02/jxhuang/gpcr_g/a/cache/gpcrdb"

# consecutive TM pairs around the bundle = the membrane-facing portals
PAIRS = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 1)]
DCUT = 5.0          # A, lipid-to-helix heavy-atom contact for "wedged"
DEPTH_THRESH = 0.5  # A, min radial insertion past the LOCAL wall to count a frame
                    # as a genuine penetration event (de-saturates open_fraction).
STRIDE = 50         # ~200 frames / replica
N_BOOT = 1000
LIPIDS = {"OPC", "POPC", "PC", "POPE", "POPS", "CHL1", "CHOL"}
PENCACHE = "/MDdata/data02/jxhuang/gpcr_g/a/cache/gateway_pen"
PROTO = ["Gs_3SN6"]


def detect_lipid_resname(u):
    rn = set(np.unique(u.atoms.resnames))
    return sorted(rn & LIPIDS)


def robust_pror_map(gd, u_prot, rec_s, rec_e):
    """PROR (GPCRdb sequence_number) -> now_ref/protpdb resid via difflib matching
    blocks. Robust to insertions/deletions between GPCRdb numbering and the modeled
    sequence (the fragile 20-residue-window aligner in build_system_config fails on
    G12 and atypical receptors). Returns (offset, {pror: protpdb})."""
    gd_sorted = sorted(gd, key=lambda x: x["sequence_number"])
    pror_nums = [r["sequence_number"] for r in gd_sorted]
    pror_seq = "".join(str(r.get("amino_acid") or "X") for r in gd_sorted)
    res = u_prot.select_atoms(f"resid {rec_s}:{rec_e} and name CA").residues
    prot_seq = resnames_to_seq(res)
    prot_resids = list(map(int, res.resids))
    sm = SequenceMatcher(None, pror_seq, prot_seq, autojunk=False)
    g2p = {}
    for blk in sm.get_matching_blocks():
        for k in range(blk.size):
            g2p[pror_nums[blk.a + k]] = prot_resids[blk.b + k]
    if not g2p:
        raise RuntimeError("no PROR<->protpdb matching blocks")
    offset = int(np.median([p - r for r, p in g2p.items()]))
    return offset, g2p


def build_gateway_config(sid, row):
    """Return dict: now_ref path, receptor (s,e), tm -> list of now_ref resids.

    Uses the low-level helpers (not build_system_config, which requires a peptide
    ligand and fails on small-molecule systems). Only needs the receptor chain,
    the PROR->protpdb offset, and the TM segment ranges.
    """
    pdb_id = row["pdb_id"]
    sysdir = os.path.dirname(row["trajectory_path"])
    now_ref = build_now_ref(sysdir)
    prot_pdb = os.path.join(os.path.dirname(row["topology_path"]), "protein.pdb")
    roles = load_corrected_roles(sid)
    # Corrected chain_roles ranges are in now_ref/viz (continuous) numbering, so
    # build the PROR map against now_ref — the frame the geometry is measured in.
    # Source protein.pdb often uses per-chain reused resids that don't match.
    u_prot = mda.Universe(now_ref) if os.path.exists(now_ref) else mda.Universe(prot_pdb)
    if roles is None or roles.get("receptor") is None:
        roles = assign_chains(detect_chains_by_ca_breaks(u_prot), row["g_protein_family"])
    if roles["receptor"] is None:
        raise RuntimeError(f"Could not identify receptor chain for {sid}")
    rec_s, rec_e = roles["receptor"]
    # If now_ref lacks the receptor range (numbering mismatch), fall back to source.
    if len(u_prot.select_atoms(f"resid {rec_s}:{rec_e} and name CA")) == 0:
        u_prot = mda.Universe(prot_pdb)
    gd = load_or_fetch_gpcrdb(pdb_id, GCACHE)
    if not gd:
        raise RuntimeError(f"No GPCRdb coverage for {pdb_id}; cannot define TM helices")
    segs = gpcrdb_segment_pror_ranges(gd)
    off, g2p = robust_pror_map(gd, u_prot, rec_s, rec_e)
    tm = {}
    for n in range(1, 8):
        key = f"TM{n}"
        if key not in segs:
            continue
        ps, pe = segs[key]
        tm[n] = [int(g2p.get(p, p + off)) for p in range(ps, pe + 1)]
    return dict(now_ref=now_ref, receptor=(int(rec_s), int(rec_e)),
                tm=tm, sysdir=sysdir)


def _reimaged_universe(now_ref, traj, rec_s, rec_e):
    u = mda.Universe(now_ref, traj)
    rec = u.select_atoms(f"resid {rec_s}:{rec_e}")
    u.trajectory.add_transformations(
        trans.center_in_box(rec, center="mass"),
        trans.wrap(u.atoms, compound="residues"))
    return u


def gateway_frame(tm_heavy_pos, tm_ca_pos, lip_pos, present_tms):
    """One frame -> {pair: (n_wedged, penetration)}.
    tm_heavy_pos[n], tm_ca_pos[n] : Nx3 arrays; lip_pos : Mx3 lipid heavy atoms."""
    all_ca = np.vstack([tm_ca_pos[n] for n in present_tms])
    cx, cy = all_ca[:, 0].mean(), all_ca[:, 1].mean()
    z_lo, z_hi = np.percentile(all_ca[:, 2], [25, 75])

    inband = lip_pos[(lip_pos[:, 2] >= z_lo) & (lip_pos[:, 2] <= z_hi)]
    def radial(p):
        return np.sqrt((p[:, 0] - cx) ** 2 + (p[:, 1] - cy) ** 2)
    lr = radial(inband) if len(inband) else np.zeros(0)

    out = {}
    for a, b in PAIRS:
        if a not in present_tms or b not in present_tms or len(inband) == 0:
            out[(a, b)] = (0, 0.0)
            continue
        ha = tm_heavy_pos[a]; hb = tm_heavy_pos[b]
        ha = ha[(ha[:, 2] >= z_lo) & (ha[:, 2] <= z_hi)]
        hb = hb[(hb[:, 2] >= z_lo) & (hb[:, 2] <= z_hi)]
        if len(ha) == 0 or len(hb) == 0:
            out[(a, b)] = (0, 0.0); continue
        da = distance_array(inband, ha).min(axis=1)
        db = distance_array(inband, hb).min(axis=1)
        wedged = (da < DCUT) & (db < DCUT)
        n = int(wedged.sum())
        if n == 0:
            out[(a, b)] = (0, 0.0); continue
        # penetration vs LOCAL helix wall: for each wedged lipid, compare its
        # radial position to that of its nearest pair-Cα (positive = lipid sits
        # radially inside the helix backbone -> intruded through the gateway).
        pair_ca = np.vstack([tm_ca_pos[a], tm_ca_pos[b]])
        pair_ca = pair_ca[(pair_ca[:, 2] >= z_lo) & (pair_ca[:, 2] <= z_hi)]
        if len(pair_ca) == 0:
            out[(a, b)] = (0, 0.0); continue
        wl = inband[wedged]
        r_wl = lr[wedged]
        nn = distance_array(wl, pair_ca).argmin(axis=1)
        r_ca_local = radial(pair_ca)[nn]
        pen = float(np.maximum(0.0, r_ca_local - r_wl).max())
        out[(a, b)] = (n, pen)
    return out


def compute_replica(cfg, traj):
    rec_s, rec_e = cfg["receptor"]
    tm = cfg["tm"]; present = sorted(tm)
    u = _reimaged_universe(cfg["now_ref"], traj, rec_s, rec_e)
    lipid_rn = detect_lipid_resname(u)
    lip_sel = u.select_atoms("resname " + " ".join(lipid_rn) + " and not name H*")
    tm_heavy = {n: u.select_atoms("resid " + " ".join(map(str, tm[n])) + " and not name H*")
                for n in present}
    tm_ca = {n: u.select_atoms("resid " + " ".join(map(str, tm[n])) + " and name CA")
             for n in present}
    acc = {p: {"n": [], "pen": []} for p in PAIRS}
    nf = 0
    for ts in u.trajectory[::STRIDE]:
        res = gateway_frame({n: tm_heavy[n].positions for n in present},
                            {n: tm_ca[n].positions for n in present},
                            lip_sel.positions, present)
        for p in PAIRS:
            acc[p]["n"].append(res[p][0])
            acc[p]["pen"].append(res[p][1])
        nf += 1
    # per-replica summary per pair. open_fraction is now DEPTH-GATED (fraction of
    # frames with a genuine insertion >= DEPTH_THRESH past the local wall), not the
    # saturated "any wedged lipid" test. penetration_p90 is the robust deep-event
    # magnitude (mean is dragged to ~0 by the many empty frames). Per-frame depth
    # arrays are returned for free re-thresholding (cached by compute_system).
    summ = {}
    arrays = {}
    for p in PAIRS:
        n = np.array(acc[p]["n"]); pen = np.array(acc[p]["pen"], float)
        summ[p] = dict(occupancy=float(n.mean()),
                       penetration=float(pen.mean()),
                       penetration_p90=float(np.percentile(pen, 90)) if len(pen) else 0.0,
                       open_fraction=float((pen >= DEPTH_THRESH).mean()))
        arrays[p] = pen.astype(np.float32)
    return summ, nf, arrays


def boot_ci(vals, n_boot=N_BOOT):
    vals = np.array(vals, float)
    if len(vals) < 2:
        return float(vals.mean()), None, None
    rng = np.random.default_rng(0)
    means = [rng.choice(vals, len(vals), replace=True).mean() for _ in range(n_boot)]
    return float(vals.mean()), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def compute_system(sid, row):
    cfg = build_gateway_config(sid, row)
    replicas = sorted(glob.glob(os.path.join(cfg["sysdir"], "prod*_now.nc")))
    per_rep = []
    cache = {}
    for ri, r in enumerate(replicas):
        summ, nf, arrays = compute_replica(cfg, r)
        per_rep.append(summ)
        for p in PAIRS:
            cache[f"rep{ri}_TM{p[0]}-TM{p[1]}"] = arrays[p]
    os.makedirs(PENCACHE, exist_ok=True)
    np.savez(os.path.join(PENCACHE, f"{sid}_pen.npz"), **cache)
    pairs_out = []
    for p in PAIRS:
        for metric in ("occupancy", "penetration", "penetration_p90", "open_fraction"):
            rep_vals = [pr[p][metric] for pr in per_rep]
            m, lo, hi = boot_ci(rep_vals)
            pairs_out.append(dict(system_id=sid, family=row["g_protein_family"],
                                  pair=f"TM{p[0]}-TM{p[1]}", metric=metric,
                                  mean=round(m, 3),
                                  ci_lo=(round(lo, 3) if lo is not None else None),
                                  ci_hi=(round(hi, 3) if hi is not None else None),
                                  n_replicas=len(per_rep),
                                  replica_values=[round(v, 3) for v in rep_vals]))
    return pairs_out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(GCACHE, exist_ok=True)
    rows = {r["system_id"]: r for r in csv.DictReader(open(INVENTORY))}
    todo = sys.argv[1:] if len(sys.argv) > 1 else PROTO
    allout = []
    for sid in todo:
        out = compute_system(sid, rows[sid])
        allout += out
        with open(os.path.join(OUTDIR, f"{sid}_gateways.json"), "w") as f:
            json.dump(out, f, indent=2)
        occ = [o for o in out if o["metric"] == "open_fraction"]
        print(f"\n[{sid}] gateway open_fraction by portal:")
        for o in sorted(occ, key=lambda x: -x["mean"]):
            print(f"  {o['pair']:8s} open={o['mean']:.2f} reps={o['replica_values']}")
    if allout:
        keys = list(allout[0].keys())
        with open(os.path.join(OUTDIR, "gateways_summary.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
            for o in allout:
                o = dict(o); o["replica_values"] = ";".join(map(str, o["replica_values"]))
                w.writerow(o)


if __name__ == "__main__":
    main()
