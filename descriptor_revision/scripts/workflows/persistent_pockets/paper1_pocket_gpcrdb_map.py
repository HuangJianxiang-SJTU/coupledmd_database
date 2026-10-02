#!/usr/bin/env python3
"""
paper1_pocket_gpcrdb_map.py  --  WP1 prep: map atlas pockets onto GPCRdb generic
numbering and label each whole-complex pocket by location.

Each {sid}_pockets.json pocket carries `lining_resids` in now_ref/protein.pdb
numbering, spanning the WHOLE complex (fpocket was run on the full ternary system).
This script, per pocket:
  - splits lining residues into receptor / Galpha / Gbeta / peptide / other by the
    chain ranges from assign_chains,
  - attaches GPCRdb generic numbers + protein_segment to the RECEPTOR-side residues
    (via robust_pror_map -> sequence_number -> display_generic_number),
  - classifies pocket location: receptor / interface / g_protein / other.

Output: {sid}_pockets_gpcrdb.json (enriched copy) and, when run over many systems,
an aggregate generic-number occurrence table for cross-receptor consensus (WP1 next).

Reuses the same low-level helpers as the gateway metric. Does NOT need the full
sweep to finish -- developed/validated on whatever {sid}_pockets.json already exist.
"""
import os, sys, csv, json, glob, warnings
from collections import Counter, defaultdict
warnings.filterwarnings("ignore")
import numpy as np
import MDAnalysis as mda

from build_system_config import (load_or_fetch_gpcrdb, detect_chains_by_ca_breaks,
                                  assign_chains)
from paper1_gateway_metric import robust_pror_map

BASE = "/MDdata/data02/jxhuang/gpcr_g/a"
INVENTORY = os.path.join(BASE, "inventory_paper1_207.csv")
GCACHE = os.path.join(BASE, "cache/gpcrdb")
ROWS = {r["system_id"]: r for r in csv.DictReader(open(INVENTORY))}

# Corrected chain roles (RCSB-backed, 2026-06-25). Preferred over the
# assign_chains() heuristic, which mislabelled Receptor/Gβ for ~219 systems.
CHAIN_ROLES_DIR = "/MDdata/data02/jxhuang/gpcr_g/gpcr_g_server/scripts/audit_output/chain_roles"


def _parse_range(rng):
    """'1-354' -> (1, 354); None if unparseable."""
    if not rng or "-" not in str(rng):
        return None
    parts = str(rng).split("-")
    try:
        return int(parts[0]), int(parts[1])
    except (ValueError, IndexError):
        return None


def load_corrected_roles(sid):
    """Load corrected chain_roles.json -> {receptor,galpha,gbeta,g_gamma,peptide:(s,e)}.

    Receptor_frag* segments are merged into a single receptor (min,max) range so
    the whole receptor span is covered for role classification. Returns None if no
    corrected file exists (caller falls back to assign_chains).
    """
    path = os.path.join(CHAIN_ROLES_DIR, f"{sid}_chain_roles.json")
    if not os.path.exists(path):
        return None
    data = json.load(open(path))
    roles = {}
    rec_min = rec_max = None
    for role, info in data.get("roles", {}).items():
        rng = _parse_range(info.get("resid_range"))
        if rng is None:
            continue
        s, e = rng
        if role.startswith("Receptor"):
            # merge all receptor fragments into one span
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


def receptor_generic_maps(sid, row):
    """Return (chain_roles, protpdb2gen, protpdb2seg, protpdb2pror) for receptor."""
    roles = load_corrected_roles(sid)
    # The corrected chain_roles.json + pocket lining_resids are in now_ref numbering
    # (continuous 1..N), NOT the source protein.pdb (per-chain, reused resids).
    # Build the pror map against now_ref so resids line up with lining_resids.
    sysdir = os.path.dirname(row["trajectory_path"])
    now_ref = os.path.join(BASE, "cache/now_ref",
                           os.path.basename(sysdir.rstrip("/")) + "_now_ref.pdb")
    u = mda.Universe(now_ref) if os.path.exists(now_ref) else \
        mda.Universe(os.path.join(os.path.dirname(row["topology_path"]), "protein.pdb"))
    if roles is None or roles.get("receptor") is None:
        roles = assign_chains(detect_chains_by_ca_breaks(u), row["g_protein_family"])
    if roles["receptor"] is None:
        raise RuntimeError(f"no receptor chain for {sid}")
    rec_s, rec_e = roles["receptor"]
    gd = load_or_fetch_gpcrdb(row["pdb_id"], GCACHE)
    if not gd:
        raise RuntimeError(f"no GPCRdb coverage for {row['pdb_id']}")
    # pror -> generic / segment
    gen = {r["sequence_number"]: r.get("display_generic_number") for r in gd}
    seg = {r["sequence_number"]: r.get("protein_segment") for r in gd}
    _, g2p = robust_pror_map(gd, u, rec_s, rec_e)      # {pror: protpdb}
    p2pror = {v: k for k, v in g2p.items()}            # {protpdb: pror}
    protpdb2gen = {p: gen.get(pr) for p, pr in p2pror.items()}
    protpdb2seg = {p: seg.get(pr) for p, pr in p2pror.items()}
    return roles, protpdb2gen, protpdb2seg, p2pror


def _role_of(resid, roles):
    for name in ("receptor", "galpha", "gbeta", "g_gamma", "peptide"):
        rng = roles.get(name)
        if rng and rng[0] <= resid <= rng[1]:
            return name
    return "other"


def classify_location(role_counts):
    n = sum(role_counts.values())
    if n == 0:
        return "unknown"
    frec = role_counts.get("receptor", 0) / n
    fg = (role_counts.get("galpha", 0) + role_counts.get("gbeta", 0)) / n
    if frec >= 0.8:
        return "receptor"
    if fg >= 0.8:
        return "g_protein"
    if frec >= 0.2 and fg >= 0.2:
        return "interface"
    return "other"


def receptor_axis(npath, roles):
    """Extracellular->intracellular axis in the now_ref REFERENCE frame (the frame the
    pocket centroids live in: every protein frame was CA-superposed onto now_ref).
    Axis points from receptor Cα COM toward Gα Cα COM (= intracellular direction).
    Returns (R, axis_hat, proj_min, proj_max) or None if it cannot be built."""
    if not os.path.exists(npath) or roles.get("galpha") is None:
        return None
    u = mda.Universe(npath)
    rs, re = roles["receptor"]; gs, ge = roles["galpha"]
    rec_ca = u.select_atoms(f"resid {rs}:{re} and name CA").positions
    ga_ca = u.select_atoms(f"resid {gs}:{ge} and name CA").positions
    if len(rec_ca) < 3 or len(ga_ca) < 3:
        return None
    R = rec_ca.mean(0); A = ga_ca.mean(0)
    ax = A - R; n = float(np.linalg.norm(ax))
    if n == 0:
        return None
    ax = ax / n
    proj = (rec_ca - R) @ ax
    return R, ax, float(proj.min()), float(proj.max())


def classify_zone(centroid, location, is_ortho, axis):
    """Membrane-topology zone for the druggable/allosteric post-filter.
    Coupling cavity (intracellular G-protein interface) is separated from
    receptor allosteric/orthosteric pockets so the druggable map is not swamped
    by the G-protein-binding site."""
    if is_ortho:
        return "orthosteric", None
    if location in ("g_protein", "interface"):
        return "coupling_interface", None
    if location != "receptor":
        return "other", None
    if axis is None:
        return "receptor_unresolved", None
    R, ax, pmin, pmax = axis
    proj = float((np.asarray(centroid) - R) @ ax)
    frac = (proj - pmin) / (pmax - pmin) if pmax > pmin else 0.5   # 0 EC .. 1 IC
    if frac < 0.33:
        zone = "extracellular_vestibule"
    elif frac < 0.66:
        zone = "tm_core_allosteric"
    else:
        zone = "intracellular_allosteric"
    return zone, round(frac, 3)


def enrich_system(sid, agg_generic=None):
    row = ROWS[sid]
    pj = None
    for d in ("atlas", "pilot"):
        cand = os.path.join(BASE, f"paper1_pockets/{d}/{sid}_pockets.json")
        if os.path.exists(cand):
            pj = cand; break
    if pj is None:
        raise FileNotFoundError(f"no pockets.json for {sid}")
    data = json.load(open(pj))
    roles, p2gen, p2seg, p2pror = receptor_generic_maps(sid, row)
    sysdir = os.path.dirname(row["trajectory_path"])
    npath = os.path.join(BASE, "cache/now_ref",
                         os.path.basename(sysdir.rstrip("/")) + "_now_ref.pdb")
    axis = receptor_axis(npath, roles)

    n_rec_mapped = 0
    for p in data["pockets"]:
        role_counts = Counter()
        rec_generics, rec_segs = [], []
        for r in p["lining_resids"]:
            role = _role_of(r, roles)
            role_counts[role] += 1
            if role == "receptor":
                g = p2gen.get(r)
                if g:
                    rec_generics.append(g)
                    n_rec_mapped += 1
                s = p2seg.get(r)
                if s:
                    rec_segs.append(s)
        p["role_counts"] = dict(role_counts)
        p["location"] = classify_location(role_counts)
        p["receptor_generic_numbers"] = sorted(set(rec_generics))
        p["receptor_segments"] = sorted(set(rec_segs))
        zone, frac = classify_zone(p["centroid"], p["location"], p["is_orthosteric"], axis)
        p["zone"] = zone
        p["axis_frac"] = frac
        if agg_generic is not None and not p["is_orthosteric"]:
            for g in set(rec_generics):
                agg_generic[g].add(sid)

    out = pj.replace("_pockets.json", "_pockets_gpcrdb.json")
    json.dump(data, open(out, "w"), indent=2)
    locs = Counter(p["location"] for p in data["pockets"])
    zones = Counter(p["zone"] for p in data["pockets"])
    return dict(sid=sid, family=row["g_protein_family"], n_pockets=len(data["pockets"]),
                locations=dict(locs), zones=dict(zones),
                n_receptor_residues_mapped=n_rec_mapped, out=out)


def main():
    targets = sys.argv[1:]
    if not targets:
        # default: every completed system with a pockets.json
        seen = set()
        for d in ("atlas", "pilot"):
            for f in sorted(glob.glob(os.path.join(BASE, f"paper1_pockets/{d}/*_pockets.json"))):
                sid = os.path.basename(f).replace("_pockets.json", "")
                if sid not in seen and sid in ROWS:
                    targets.append(sid); seen.add(sid)
    print(f"[gpcrdb-map] {len(targets)} system(s)")
    agg_generic = defaultdict(set)
    ok = fail = 0
    for sid in targets:
        try:
            r = enrich_system(sid, agg_generic)
            ok += 1
            print(f"  {sid:12s} {r['family']:7s} pockets={r['n_pockets']:3d} "
                  f"loc={r['locations']}\n               zones={r['zones']}")
        except Exception as e:
            fail += 1
            print(f"  [FAIL {sid}] {type(e).__name__}: {e}")
    print(f"\nmapped: ok={ok} fail={fail}")
    if agg_generic:
        rows = sorted(((g, len(s)) for g, s in agg_generic.items()), key=lambda x: -x[1])
        print(f"\nTop allosteric-pocket generic positions (across {ok} systems):")
        for g, n in rows[:25]:
            print(f"  {g:12s} {n} systems")
        with open(os.path.join(BASE, "paper1_pockets/allosteric_generic_occurrence.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["generic_number", "n_systems"])
            w.writerows(rows)


if __name__ == "__main__":
    main()
