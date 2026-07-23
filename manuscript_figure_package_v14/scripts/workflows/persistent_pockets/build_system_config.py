#!/usr/bin/env python3
"""
build_system_config.py

Provides build_system_config(pdb_id, inventory_row, cache_dir) which
automatically generates the per-system config dict used by the
fingerprint pipeline.

Steps performed:
  1. Parse inventory row for paths, G-protein family, ligand metadata.
  2. Load or fetch GPCRdb generic numbering for the receptor.
  3. Run Cα chain-break detection on protein.pdb → Gα / Gβ(+Gγ) /
     receptor / peptide-ligand boundaries.
  4. Sequence-align receptor chain → identify PROR offset (protein.pdb
     sequential resid = PROR resid + offset).
  5. Identify pocket residues within 5 Å of peptide-ligand heavy atoms.
  6. Return populated config dict matching the SYSTEMS[pdb_id] schema.

Raises SystemConfigError (with descriptive message) rather than
silently proceeding when a required input is missing.
"""

import os
import json
import glob
import warnings
import numpy as np
import requests
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

warnings.filterwarnings("ignore")

GPCRDB_URL  = "https://gpcrdb.org/services/residues/extended/{pdb_id}/"
GPCRDB_TIMEOUT = 30    # seconds
CHAIN_BREAK_DIST = 5.0  # Å Cα–Cα gap → new chain
POCKET_CUTOFF    = 5.0  # Å, heavy-atom distance for pocket definition
MIN_PEPTIDE_LEN  = 5    # residues — shorter = not a peptide agonist
MAX_PEPTIDE_LEN  = 50

# Galpha subtypes with known chain lengths (approximate, for chain assignment)
# Gs ~394 aa, Gi ~354 aa, Gq ~359 aa, G12 ~379 aa
GALPHA_MIN_LEN = 300
GALPHA_MAX_LEN = 450
GBETA_MIN_LEN  = 100   # Gβ + Gγ together

# Generic numbers for key positions
GN_350 = "3.50"  # DRY motif
GN_630 = "6.30"  # TM6 cytoplasmic start reference


class SystemConfigError(Exception):
    pass


# ─── GPCRdb ───────────────────────────────────────────────────────────────────
def load_or_fetch_gpcrdb(pdb_id, cache_dir):
    """
    Return list of GPCRdb residue dicts for pdb_id.
    Uses per-pilot cache if present; otherwise queries GPCRdb REST API
    and caches to {cache_dir}/gpcrdb_{pdb_id}.json.
    Returns None if GPCRdb has no coverage for this PDB.
    """
    cache_path = os.path.join(cache_dir, f"gpcrdb_{pdb_id}.json")
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            return json.load(f)

    url = GPCRDB_URL.format(pdb_id=pdb_id.lower())
    try:
        r = requests.get(url, timeout=GPCRDB_TIMEOUT)
        if r.status_code == 404 or not r.text.strip():
            return None
        data = r.json()
        if not data:
            return None
        os.makedirs(cache_dir, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(data, f, indent=2)
        return data
    except Exception as e:
        print(f"  GPCRdb fetch failed for {pdb_id}: {e}")
        return None


def gpcrdb_segment_pror_ranges(gpcrdb_data):
    """
    Returns dict: segment_name → (pror_start, pror_end)
    Only includes TM1-TM7, ICL2, H8.
    """
    segs = {}
    for r in gpcrdb_data:
        seg = r["protein_segment"]
        pror = r["sequence_number"]
        if seg not in segs:
            segs[seg] = [pror, pror]
        else:
            segs[seg][0] = min(segs[seg][0], pror)
            segs[seg][1] = max(segs[seg][1], pror)
    return {s: tuple(v) for s, v in segs.items()}


def gpcrdb_generic_number_to_pror(gpcrdb_data, gn_target):
    """
    Given a generic number string like '3.50', return the PROR resid.
    GPCRdb format: display_generic_number = '3.50x50' — strip 'x...' suffix.
    """
    for r in gpcrdb_data:
        gn = r.get("display_generic_number", "") or ""
        gn_clean = gn.split("x")[0]
        if gn_clean == gn_target:
            return r["sequence_number"]
    return None


# ─── Cα chain-break detection ─────────────────────────────────────────────────
def detect_chains_by_ca_breaks(u, break_dist=CHAIN_BREAK_DIST):
    """
    Returns list of (start_resid, end_resid, n_residues) tuples,
    one per contiguous Cα chain, sorted by start_resid.
    """
    ca = u.select_atoms("name CA")
    if ca.n_atoms == 0:
        raise SystemConfigError("No CA atoms found in protein.pdb")

    chains = []
    chain_start_idx = 0

    for i in range(len(ca) - 1):
        d = np.linalg.norm(ca.positions[i + 1] - ca.positions[i])
        if d > break_dist:
            chains.append((int(ca.resids[chain_start_idx]),
                           int(ca.resids[i]),
                           i - chain_start_idx + 1))
            chain_start_idx = i + 1

    # Last chain
    chains.append((int(ca.resids[chain_start_idx]),
                   int(ca.resids[-1]),
                   len(ca) - chain_start_idx))
    return chains


def assign_chains(chains, g_protein_family):
    """
    Assign roles to detected chains based on length heuristics and order.
    GPCR-G-protein complex layout in protein.pdb (post-cpptraj strip):
      Chain 0: Gα (longest, 300–450 aa)
      Chain 1: Gβ (next, ~340 aa)
      Chain 2: Gγ (short, ~60–80 aa) — may be absent or merged with Gβ
      Chain 3 (or 2): receptor (>200 aa)
      Chain last: peptide ligand (5–50 aa) — may be absent

    Returns dict with keys: galpha, gbeta, receptor, peptide (or None)
    each as (start_resid, end_resid).
    """
    assigned = {"galpha": None, "gbeta": None, "receptor": None, "peptide": None}

    # Sort by start resid (should already be ordered)
    chains_sorted = sorted(chains, key=lambda c: c[0])

    # Gα: first chain >= GALPHA_MIN_LEN residues
    remaining = list(chains_sorted)
    for i, (s, e, n) in enumerate(remaining):
        if GALPHA_MIN_LEN <= n <= GALPHA_MAX_LEN:
            assigned["galpha"] = (s, e)
            remaining.pop(i)
            break

    if assigned["galpha"] is None:
        # Fallback: largest chain
        c = max(chains_sorted, key=lambda x: x[2])
        assigned["galpha"] = (c[0], c[1])
        remaining = [x for x in chains_sorted if x != c]

    # Receptor: LAST chain with >= 200 residues.
    # In protein.pdb (post-strip) the order is always:
    #   Gα → Gβ → Gγ (short) → [peptide (short)] → receptor
    # so the receptor is the last large chain.
    galpha_start = assigned["galpha"][0]
    large_post_galpha = [(s, e, n) for s, e, n in remaining
                         if s > galpha_start and n >= 200]
    if large_post_galpha:
        rec = large_post_galpha[-1]   # last large chain = receptor
        assigned["receptor"] = (rec[0], rec[1])
        remaining = [(s, e, n) for s, e, n in remaining if not (s == rec[0] and e == rec[1])]

    # Gbeta: all chains between Gα end and receptor start with >= GBETA_MIN_LEN residues.
    # Merges Gβ + Gγ into a single range (matches SYSTEMS["gbeta"] convention).
    if assigned["receptor"]:
        rec_start = assigned["receptor"][0]
        mid = [(s, e, n) for s, e, n in remaining
               if galpha_start < s < rec_start and n >= GBETA_MIN_LEN]
        if mid:
            assigned["gbeta"] = (mid[0][0], mid[-1][1])
            remaining = [(s, e, n) for s, e, n in remaining
                         if not any(s == m[0] for m in mid)]

    # Peptide: any remaining chain in the right size range
    pep_candidates = [(s, e, n) for s, e, n in remaining
                      if MIN_PEPTIDE_LEN <= n <= MAX_PEPTIDE_LEN]
    if pep_candidates:
        pep = pep_candidates[0]
        assigned["peptide"] = (pep[0], pep[1])

    return assigned


# ─── Sequence alignment ───────────────────────────────────────────────────────
THREE_TO_ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "HSD":"H","HSE":"H","HSP":"H","HID":"H","HIE":"H","HIP":"H",
    "CYX":"C","CYM":"C","ASH":"D","GLH":"E","LYN":"K",
}

def resnames_to_seq(residues):
    return "".join(THREE_TO_ONE.get(r.resname, "X") for r in residues)


def sliding_window_align(query_seq, target_seq, min_identity=0.85, window=20):
    """
    Find the best-matching position of query_seq in target_seq.
    Returns (best_start_in_target, identity) or (None, 0.0).
    Searches using a window of `window` residues from the query.
    """
    q_win = query_seq[:window]
    best_start, best_id = None, 0.0
    for i in range(len(target_seq) - len(q_win) + 1):
        t_win = target_seq[i:i + len(q_win)]
        matches = sum(a == b for a, b in zip(q_win, t_win))
        ident = matches / len(q_win)
        if ident > best_id:
            best_id = ident
            best_start = i
    return best_start, best_id


def build_gpcrdb_to_protpdb_map(gpcrdb_data, u_prot, rec_start, rec_end, pror_range):
    """
    Align the GPCRdb/PROR sequence onto the protein.pdb receptor chain
    to establish the PROR resid → protpdb resid mapping.

    gpcrdb_data: list of GPCRdb residue dicts (sorted by sequence_number)
    u_prot: MDAnalysis Universe for protein.pdb
    rec_start, rec_end: protpdb resid range of the receptor chain
    pror_range: (pror_min, pror_max) from GPCRdb or inventory

    Returns (pror_offset, gpcrdb_to_protpdb_dict)
    pror_offset: protpdb_resid = pror_resid + pror_offset
    """
    # Get PROR sequence from GPCRdb.
    # amino_acid field is already single-letter code in the REST API.
    pror_seq = "".join(
        str(r.get("amino_acid") or "X")
        for r in sorted(gpcrdb_data, key=lambda x: x["sequence_number"])
    )
    pror_min = min(r["sequence_number"] for r in gpcrdb_data)

    # Get protein.pdb receptor residues sequence
    prot_residues = u_prot.select_atoms(f"resid {rec_start}:{rec_end} and name CA").residues
    prot_seq = resnames_to_seq(prot_residues)
    prot_resids = list(prot_residues.resids)

    # Align first 20 residues of PROR sequence onto protein.pdb
    best_pos, identity = sliding_window_align(pror_seq, prot_seq, window=20)
    if best_pos is None or identity < 0.70:
        raise SystemConfigError(
            f"Sequence alignment failed: best identity {identity:.2f} < 0.70. "
            f"PROR seq[:20]={pror_seq[:20]}  PROT seq={prot_seq[:40]}"
        )

    # offset: protpdb resid of first receptor residue matching PROR start
    pror_offset = prot_resids[best_pos] - pror_min

    # Build full mapping
    g2p = {}
    for r in gpcrdb_data:
        pror = r["sequence_number"]
        protpdb = pror + pror_offset
        if rec_start <= protpdb <= rec_end:
            g2p[pror] = protpdb

    return pror_offset, g2p


# ─── Pocket residues ──────────────────────────────────────────────────────────
def find_pocket_residues(u_full, peptide_range, receptor_range, cutoff=POCKET_CUTOFF):
    """
    Given a loaded MDAnalysis Universe for the full system PDB
    (step5_input.pdb or ref.pdb), return list of (protpdb_resid, resname)
    tuples for receptor residues with any heavy atom within cutoff Å
    of any peptide heavy atom.

    peptide_range, receptor_range: (start_resid, end_resid) in that PDB's numbering.
    Selects only protein (backbone + sidechain) atoms to exclude lipid/water
    that may share resid numbers.
    """
    p1, p2 = peptide_range
    r1, r2 = receptor_range

    pep_heavy = u_full.select_atoms(f"resid {p1}:{p2} and protein and not name H*")
    if pep_heavy.n_atoms == 0:
        # Fallback without protein keyword (for non-standard residue names)
        pep_heavy = u_full.select_atoms(f"resid {p1}:{p2} and not name H*")
    rec_heavy  = u_full.select_atoms(f"resid {r1}:{r2} and protein and not name H*")

    if pep_heavy.n_atoms == 0:
        raise SystemConfigError(
            f"No heavy atoms found in peptide range {p1}–{p2} in full PDB."
        )
    if rec_heavy.n_atoms == 0:
        raise SystemConfigError(
            f"No heavy atoms found in receptor range {r1}–{r2} in full PDB."
        )

    pep_pos = pep_heavy.positions
    pocket = []
    for res in rec_heavy.residues:
        res_heavy = res.atoms.select_atoms("not name H*")
        if res_heavy.n_atoms == 0:
            continue
        dists = distance_array(res_heavy.positions, pep_pos)
        if dists.min() <= cutoff:
            pocket.append((int(res.resid), str(res.resname)))

    return pocket


# ─── Main entry point ─────────────────────────────────────────────────────────
def build_system_config(pdb_id, inventory_row, cache_dir,
                        pocket_cutoff=POCKET_CUTOFF):
    """
    Build and return a config dict for pdb_id.

    Parameters
    ----------
    pdb_id : str
    inventory_row : dict  (one row from inventory.csv as DictReader)
    cache_dir : str       (directory for GPCRdb JSON cache)
    pocket_cutoff : float (Å, default 5.0)

    Returns
    -------
    cfg : dict with keys matching SYSTEMS[pdb_id] schema:
        label, color, top, trajs,
        galpha, gbeta, receptor, alpha5, hook,
        icl2, res350, res630,
        tm_pror, pror_offset,
        pocket_resids,          # list of protpdb resids
        gpcrdb_to_protpdb,      # {pror_resid: protpdb_resid}

    Raises SystemConfigError on any unrecoverable condition.
    """
    cfg = {}
    pdb_id = pdb_id.upper()

    # ── Step 1: paths from inventory ─────────────────────────────────────────
    top_path = inventory_row.get("topology_path", "").strip()
    traj_path_sample = inventory_row.get("trajectory_path", "").strip()
    g_family = inventory_row.get("g_protein_family", "").strip()
    n_replicas = int(inventory_row.get("n_replicas", 1))

    # Derive protein.pdb from topology directory
    sys_dir = os.path.dirname(top_path)
    prot_pdb = os.path.join(sys_dir, "protein.pdb")
    if not os.path.exists(prot_pdb):
        raise SystemConfigError(f"protein.pdb not found in {sys_dir}")

    cfg["top"] = prot_pdb

    # Trajectories: traj1.nc, traj2.nc, ... (protein-only stripped)
    trajs = sorted(glob.glob(os.path.join(sys_dir, "traj*.nc")))
    if not trajs:
        raise SystemConfigError(f"No traj*.nc files in {sys_dir}")
    cfg["trajs"] = trajs[:n_replicas]  # cap to inventory n_replicas

    # G-protein label and color
    GFAM_LABEL = {"Gs": "Gs", "Gi": "Gi", "Gq": "Gq", "G12-13": "G12"}
    GFAM_COLOR = {"Gs": "#d62728", "Gi": "#1f77b4",
                  "Gq": "#2ca02c", "G12-13": "#9467bd"}
    cfg["label"] = GFAM_LABEL.get(g_family, g_family)
    cfg["color"] = GFAM_COLOR.get(g_family, "#7f7f7f")

    # ── Step 2: GPCRdb generic numbering ─────────────────────────────────────
    gpcrdb_data = load_or_fetch_gpcrdb(pdb_id, cache_dir)

    if gpcrdb_data:
        seg_ranges = gpcrdb_segment_pror_ranges(gpcrdb_data)
        pror_min = min(r["sequence_number"] for r in gpcrdb_data)
        pror_max = max(r["sequence_number"] for r in gpcrdb_data)
    else:
        print(f"  [{pdb_id}] No GPCRdb coverage — TM ranges will be estimated.")
        seg_ranges = {}
        pror_min, pror_max = None, None

    # ── Step 3: Cα chain-break detection on protein.pdb ──────────────────────
    u_prot = mda.Universe(prot_pdb)
    chains = detect_chains_by_ca_breaks(u_prot)
    roles  = assign_chains(chains, g_family)

    if roles["galpha"] is None:
        raise SystemConfigError(f"Could not identify Gα chain in {prot_pdb}")
    if roles["receptor"] is None:
        raise SystemConfigError(f"Could not identify receptor chain in {prot_pdb}")

    galpha_s, galpha_e = roles["galpha"]
    gbeta_s,  gbeta_e  = roles["gbeta"] if roles["gbeta"] else (None, None)
    rec_s,    rec_e    = roles["receptor"]

    # ── Step 4: Sequence alignment → PROR offset ─────────────────────────────
    if gpcrdb_data:
        pror_offset, gpcrdb_to_protpdb = build_gpcrdb_to_protpdb_map(
            gpcrdb_data, u_prot, rec_s, rec_e, (pror_min, pror_max)
        )
    else:
        # Fallback: no GPCRdb — offset estimated as rec_s - 1 (PROR starts at 1)
        # TM ranges will be left empty; downstream will warn
        pror_offset = rec_s - 1
        gpcrdb_to_protpdb = {}

    cfg["pror_offset"] = pror_offset
    cfg["gpcrdb_to_protpdb"] = gpcrdb_to_protpdb

    # ── Step 5: Populate structural ranges in protein.pdb numbering ──────────
    cfg["galpha"]   = (galpha_s, galpha_e)
    cfg["gbeta"]    = (gbeta_s, gbeta_e) if gbeta_s else (None, None)
    cfg["receptor"] = (rec_s, rec_e)

    # alpha5 = last 26 residues of Gα; hook = last 4
    cfg["alpha5"] = (galpha_e - 25, galpha_e)
    cfg["hook"]   = (galpha_e - 3,  galpha_e)

    # ICL2 from GPCRdb
    if "ICL2" in seg_ranges:
        icl2_pror_s, icl2_pror_e = seg_ranges["ICL2"]
        icl2_s = gpcrdb_to_protpdb.get(icl2_pror_s) or (icl2_pror_s + pror_offset)
        icl2_e = gpcrdb_to_protpdb.get(icl2_pror_e) or (icl2_pror_e + pror_offset)
        cfg["icl2"] = (icl2_s, icl2_e)
    else:
        cfg["icl2"] = (None, None)

    # 3.50 and 6.30 key residues
    if gpcrdb_data:
        pror_350 = gpcrdb_generic_number_to_pror(gpcrdb_data, GN_350)
        pror_630 = gpcrdb_generic_number_to_pror(gpcrdb_data, GN_630)
        cfg["res350"] = gpcrdb_to_protpdb.get(pror_350) if pror_350 else None
        cfg["res630"] = gpcrdb_to_protpdb.get(pror_630) if pror_630 else None
    else:
        cfg["res350"] = None
        cfg["res630"] = None

    # TM helix PROR ranges (for RMSD alignment)
    tm_keys = ["TM1","TM2","TM3","TM4","TM5","TM6","TM7"]
    cfg["tm_pror"] = []
    for tm in tm_keys:
        if tm in seg_ranges:
            cfg["tm_pror"].append(seg_ranges[tm])

    # ── Step 6: Pocket residues from full-system PDB ──────────────────────────
    # Locate the full system PDB (step5_input.pdb preferred, fall back to ref.pdb)
    full_pdb = os.path.join(sys_dir, "step5_input.pdb")
    if not os.path.exists(full_pdb):
        full_pdb = os.path.join(sys_dir, "ref.pdb")
    if not os.path.exists(full_pdb):
        raise SystemConfigError(
            f"No step5_input.pdb or ref.pdb found in {sys_dir}. "
            f"Cannot identify pocket residues."
        )

    u_full = mda.Universe(full_pdb)

    # Identify peptide in full PDB.
    # Try PROL segment first (CHARMM-GUI convention), then Cα break detection.
    pep_range_full = None
    pep_chain_info = {}

    prol = u_full.select_atoms("segid PROL")
    if prol.n_residues >= MIN_PEPTIDE_LEN:
        pep_resids = sorted(set(prol.resids))
        pep_range_full = (pep_resids[0], pep_resids[-1])
        pep_chain_info = {"method": "PROL_segment",
                          "n_residues": prol.n_residues}
    else:
        # Cα break detection on full PDB; find chains in peptide-length range
        # that are positioned extracellularly (high |Z| opposite to Gα)
        full_chains = detect_chains_by_ca_breaks(u_full)
        pep_candidates = [(s, e, n) for s, e, n in full_chains
                          if MIN_PEPTIDE_LEN <= n <= MAX_PEPTIDE_LEN]
        if not pep_candidates:
            raise SystemConfigError(
                f"No peptide ligand chain found in {full_pdb}. "
                f"PROL segment absent and no chain with {MIN_PEPTIDE_LEN}–"
                f"{MAX_PEPTIDE_LEN} residues detected. "
                f"If this receptor has a small-molecule ligand, use "
                f"small_molecule_pocket() instead."
            )

        # Pick the candidate most distal from Gα (extracellular face)
        # Gα COM Z vs each peptide candidate Z
        ga_com_z = u_full.select_atoms(
            f"resid {galpha_s}:{galpha_e} and name CA"
        ).center_of_mass()[2] if galpha_s else 0.0

        best_pep, best_dist = None, 0.0
        for s, e, n in pep_candidates:
            pep_com_z = u_full.select_atoms(
                f"resid {s}:{e} and name CA"
            ).center_of_mass()[2]
            dist = abs(pep_com_z - ga_com_z)
            if dist > best_dist:
                best_dist = dist
                best_pep = (s, e, n)

        if best_pep is None or best_dist < 20.0:
            raise SystemConfigError(
                f"Peptide candidate found but not extracellular (Gα–peptide "
                f"Z-distance = {best_dist:.1f} Å < 20 Å). "
                f"Check system layout."
            )

        pep_range_full = (best_pep[0], best_pep[1])
        pep_chain_info = {"method": "ca_break_detection",
                          "n_residues": best_pep[2],
                          "galpha_pep_z_dist": round(best_dist, 1)}

    # Receptor range in full PDB: need to identify by sequence match
    # against protein.pdb receptor residues
    rec_seq_prot = resnames_to_seq(
        u_prot.select_atoms(f"resid {rec_s}:{rec_e} and name CA").residues
    )
    full_ca = u_full.select_atoms("name CA")
    full_seq = resnames_to_seq(full_ca.residues)
    full_resids = list(full_ca.residues.resids)

    best_pos, identity = sliding_window_align(rec_seq_prot, full_seq, window=25)
    if best_pos is None or identity < 0.70:
        raise SystemConfigError(
            f"Cannot align receptor sequence from protein.pdb onto full PDB "
            f"{full_pdb} (identity {identity:.2f}). Cannot compute pocket."
        )
    full_rec_s = full_resids[best_pos]
    full_rec_e = full_resids[min(best_pos + (rec_e - rec_s), len(full_resids) - 1)]

    pocket_raw = find_pocket_residues(u_full, pep_range_full,
                                      (full_rec_s, full_rec_e),
                                      cutoff=pocket_cutoff)

    # Map pocket residues from full-PDB numbering to protein.pdb numbering
    # full_pdb_resid → protpdb_resid via shared offset
    full_rec_offset = rec_s - full_rec_s  # protpdb = full_resid + offset
    pocket_protpdb = []
    for full_resid, resname in pocket_raw:
        protpdb_resid = full_resid + full_rec_offset
        # Also compute PROR resid via inverse of gpcrdb_to_protpdb
        pror_resid = protpdb_resid - pror_offset  # approximate (exact if linear)
        pocket_protpdb.append({
            "pror_resid": pror_resid,
            "protpdb_resid": protpdb_resid,
            "resname": resname,
        })

    cfg["pocket_resids"] = sorted(set(p["protpdb_resid"] for p in pocket_protpdb))
    cfg["pocket_records"] = pocket_protpdb
    cfg["peptide_chain_info"] = pep_chain_info

    return cfg


# ─── Convenience: load inventory as dict keyed by pdb_id ─────────────────────
def load_inventory(inventory_path):
    import csv
    rows = {}
    with open(inventory_path) as f:
        for row in csv.DictReader(f):
            rows[row["pdb_id"].strip().upper()] = row
    return rows


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python build_system_config.py <PDB_ID>")
        sys.exit(1)
    pdb_id = sys.argv[1].upper()
    inv = load_inventory("/MDdata/data02/jxhuang/gpcr_g/a/inventory_v3_extended.csv")
    if pdb_id not in inv:
        print(f"PDB ID {pdb_id} not in inventory")
        sys.exit(1)
    cfg = build_system_config(
        pdb_id, inv[pdb_id],
        cache_dir="/MDdata/data02/jxhuang/gpcr_g/a/cache/cckar_pilot"
    )
    for k, v in cfg.items():
        if k not in ("gpcrdb_to_protpdb", "pocket_records", "trajs"):
            print(f"  {k}: {v}")
    print(f"  trajs: {cfg['trajs']}")
    print(f"  pocket_resids ({len(cfg['pocket_resids'])}): {cfg['pocket_resids'][:5]}...")
    print(f"  gpcrdb_to_protpdb: {len(cfg['gpcrdb_to_protpdb'])} entries")
