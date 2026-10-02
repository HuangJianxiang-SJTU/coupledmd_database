#!/usr/bin/env python3
"""
paper1_pocket_pilot.py  --  Validation pilot for the allosteric-pocket atlas.

Per-frame fpocket + grid aggregation (mdpocket binary is broken in this env).
For each pilot system:
  1. Build water-stripped reference topology (ref.pdb minus water) -> matches prod*_now.nc.
  2. Load each replica, align every protein frame onto the reference (CA superpose).
  3. Run fpocket per subsampled frame; collect alpha-sphere centers (in reference frame).
  4. Aggregate -> per-voxel pocket frequency (fraction of frames a 1A voxel is occupied).
  5. Cluster persistent voxels (freq>=THRESH) -> consensus pockets -> lining residues.
  6. Classify orthosteric (overlaps bound ligand) and report validation.

Run:  mamba run -n base python paper1_pocket_pilot.py   (fpocket called by abs path)
"""
import os, sys, glob, csv, json, shutil, subprocess, tempfile, warnings, atexit
import multiprocessing as mp
import numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from MDAnalysis.analysis import align
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN
from build_system_config import detect_chains_by_ca_breaks, assign_chains

FPOCKET = "/MDdata/data01/jxhuang/miniconda3/envs/fpocket/bin/fpocket"
INVENTORY = "/MDdata/data02/jxhuang/gpcr_g/a/inventory_v3_extended.csv"
OUTDIR = "/MDdata/data02/jxhuang/gpcr_g/a/paper1_pockets/pilot"
CACHE = "/MDdata/data02/jxhuang/gpcr_g/a/cache"

# Corrected chain roles (RCSB-backed, 2026-06-25). Preferred over assign_chains().
CHAIN_ROLES_DIR = "/MDdata/data02/jxhuang/gpcr_g/gpcr_g_server/scripts/audit_output/chain_roles"


def _pp_parse_range(rng):
    if not rng or "-" not in str(rng):
        return None
    parts = str(rng).split("-")
    try:
        return int(parts[0]), int(parts[1])
    except (ValueError, IndexError):
        return None


def load_corrected_roles(sid):
    """Corrected chain_roles.json -> {receptor,galpha,gbeta,g_gamma,peptide:(s,e)}.
    Receptor_frag* merged into one span. None if no corrected file."""
    path = os.path.join(CHAIN_ROLES_DIR, f"{sid}_chain_roles.json")
    if not os.path.exists(path):
        return None
    data = json.load(open(path))
    roles = {}
    rec_min = rec_max = None
    for role, info in data.get("roles", {}).items():
        rng = _pp_parse_range(info.get("resid_range"))
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

PILOT = ["Gs_3SN6", "Gi_7SQO", "Gq_6OIJ", "Gs_7D68"]
NPROC = int(os.environ.get("FP_NPROC", "8"))   # polite default on a shared 48-core box
STRIDE = 50            # ~200 frames per 10000-frame replica
FREQ_THRESH = 0.85     # grid point inside a pocket in >=85% of frames -> persistent
                       # (0.5 merges shallow surface pockets into one blob; 0.8 still
                       #  over-merges the orthosteric pocket with its vestibule on
                       #  β2AR. 0.85 separates cleanly. Full grid is saved so
                       #  re-clustering at any isovalue is free -- see reprocess_saved.)
GRID = 1.5             # A, analysis-grid spacing
LIG_CUT = 5.0          # A, pocket-to-ligand for orthosteric

AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR "
         "TRP TYR VAL HSD HSE HSP HID HIE HIP HSE".split())
WATER = {"WAT", "HOH", "TIP3", "SOL", "T3P"}
IONS = {"SOD", "CLA", "NA", "CL", "POT", "K", "MG", "CAL", "ZN", "ZN2", "NA+", "CL-"}


def build_now_ref(sysdir):
    """Write ref.pdb minus water (matches prod*_now.nc) into cache; return path."""
    ref = os.path.join(sysdir, "ref.pdb")
    sid = os.path.basename(sysdir.rstrip("/"))
    out = os.path.join(CACHE, "now_ref", f"{sid}_now_ref.pdb")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not os.path.exists(out):
        with open(ref) as f, open(out, "w") as o:
            for l in f:
                if l.startswith(("ATOM", "HETATM")) and l[17:21].strip() in WATER:
                    continue
                o.write(l)
    return out


def classify_resnames(u):
    rn = set(np.unique(u.atoms.resnames))
    lipid = {r for r in rn if r in ("OPC", "POPC", "PC", "POPE", "POPS", "CHL1", "CHOL")}
    # leftover non-AA, non-lipid, non-ion, non-water -> ligand candidates.
    # Restrict to HETATM records so modified protein residues (CYM, HSD, MSE, ...)
    # that are part of the protein backbone are NOT mistaken for ligands.
    hetatm_rn = set(np.unique(u.select_atoms("record_type HETATM").resnames))
    lig = (rn - AA - lipid - IONS - WATER) & hetatm_rn
    return lipid, lig


# ─── parallel per-frame worker ─────────────────────────────────────────────────
_W = {}   # per-process cache


def _init_worker(now_ref, prot_sel, grid_xyz):
    _W["now_ref"] = now_ref
    _W["prot_sel"] = prot_sel
    _W["grid"] = grid_xyz
    _W["ref_prot"] = mda.Universe(now_ref).select_atoms(prot_sel)
    _W["univ"] = {}
    _W["tmp"] = tempfile.mkdtemp(prefix="fpw_", dir="/tmp")
    atexit.register(lambda: shutil.rmtree(_W["tmp"], ignore_errors=True))


def _frame_worker(task):
    """Return uint32 indices of grid points inside a pocket this frame."""
    rep, fidx = task
    if rep not in _W["univ"]:
        _W["univ"][rep] = mda.Universe(_W["now_ref"], rep)
    u = _W["univ"][rep]
    mob = u.select_atoms(_W["prot_sel"])
    u.trajectory[fidx]
    align.alignto(mob, _W["ref_prot"], select="name CA")
    fr = os.path.join(_W["tmp"], f"fr_{os.getpid()}.pdb")
    mob.write(fr)
    centers, radii = run_fpocket(fr)
    shutil.rmtree(fr[:-4] + "_out", ignore_errors=True)
    if len(centers) == 0:
        return np.zeros(0, dtype=np.uint32)
    d, idx = cKDTree(centers).query(_W["grid"], k=1)
    return np.where(d < radii[idx])[0].astype(np.uint32)


def run_fpocket(pdb_path):
    """Run fpocket; return (centers Nx3, radii N) of alpha spheres."""
    subprocess.run([FPOCKET, "-f", pdb_path], stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, check=False)
    base = pdb_path[:-4] + "_out"
    pqr = os.path.join(base, os.path.basename(pdb_path)[:-4] + "_pockets.pqr")
    pts, rad = [], []
    if os.path.exists(pqr):
        for l in open(pqr):
            if l.startswith("ATOM"):
                pts.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
                rad.append(float(l.split()[-1]))
    if not pts:
        return np.zeros((0, 3)), np.zeros(0)
    return np.array(pts), np.array(rad)


def get_ligand_coords(ref_u, row):
    """Return (lig_xyz, ligand_type, resnames) for orthosteric classification.

    Small-molecule resname first; if none (peptide-agonist systems, e.g. class B),
    fall back to the bound peptide chain detected by Cα-break chain assignment.
    """
    _, lig_rn = classify_resnames(ref_u)
    if lig_rn:
        xyz = ref_u.select_atoms("resname " + " ".join(sorted(lig_rn)) +
                                 " and not name H*").positions
        return xyz, "small_molecule", sorted(lig_rn)
    try:
        # Prefer corrected chain_roles (RCSB-backed). If a corrected file exists,
        # trust it: a peptide_ligand entry is the peptide, otherwise there is no
        # peptide ligand (do NOT fall back to assign_chains, which mislabels Gγ
        # as a peptide). Only use the assign_chains fallback when no corrected
        # file is available at all.
        roles = load_corrected_roles(row.get("system_id", ""))
        if roles is None:
            prot_pdb = os.path.join(os.path.dirname(row["topology_path"]), "protein.pdb")
            roles = assign_chains(detect_chains_by_ca_breaks(mda.Universe(prot_pdb)),
                                  row["g_protein_family"])
        if roles.get("peptide"):
            s, e = roles["peptide"]
            xyz = ref_u.select_atoms(f"resid {s}:{e} and not name H*").positions
            return xyz, "peptide", [f"peptide:{s}-{e}"]
    except Exception:
        pass
    return np.zeros((0, 3)), "none", []


def call_pockets(grid_xyz, freq, ref_heavy, lig_xyz, thresh=FREQ_THRESH, eps=GRID * 1.9):
    """Cluster persistent grid points into pockets; classify orthosteric.
    Returns (pockets list sorted by size, n_persistent)."""
    keep = freq >= thresh
    persistent = grid_xyz[keep]
    pfreq = freq[keep]
    pockets = []
    if len(persistent) >= 3:
        labels = DBSCAN(eps=eps, min_samples=3).fit_predict(persistent)
        rh = ref_heavy.positions
        for cl in sorted(set(labels) - {-1}):
            pts = persistent[labels == cl]
            d2lig = (cdist(pts, lig_xyz).min() if len(lig_xyz) else np.nan)
            ortho = bool(len(lig_xyz) and d2lig <= LIG_CUT)
            near = ref_heavy[cdist(rh, pts).min(axis=1) < 4.0]
            resids = sorted(set(int(r) for r in near.resids))
            pockets.append(dict(pocket_id=int(cl), n_voxels=int(len(pts)),
                                mean_freq=round(float(pfreq[labels == cl].mean()), 3),
                                max_freq=round(float(pfreq[labels == cl].max()), 3),
                                centroid=[round(float(x), 1) for x in pts.mean(0)],
                                dist_to_ligand=(round(float(d2lig), 2)
                                                if not np.isnan(d2lig) else None),
                                is_orthosteric=ortho, n_lining=len(resids),
                                lining_resids=resids))
    pockets.sort(key=lambda p: -p["n_voxels"])
    return pockets, int(len(persistent))


def process_system(sid, row):
    sysdir = os.path.dirname(row["trajectory_path"])
    now_ref = build_now_ref(sysdir)
    ref_u = mda.Universe(now_ref)
    lipid_rn, _ = classify_resnames(ref_u)
    prot_sel = "resname " + " ".join(sorted(AA))
    ref_prot = ref_u.select_atoms(prot_sel)
    lig_xyz, lig_type, lig_names = get_ligand_coords(ref_u, row)

    replicas = sorted(glob.glob(os.path.join(sysdir, "prod*_now.nc")))
    # static analysis grid over reference protein bounding box
    pmin = ref_prot.positions.min(0) - 6.0
    pmax = ref_prot.positions.max(0) + 6.0
    axes = [np.arange(pmin[i], pmax[i], GRID) for i in range(3)]
    gx, gy, gz = np.meshgrid(*axes, indexing="ij")
    grid_xyz = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
    count = np.zeros(len(grid_xyz), dtype=np.int32)

    # frame tasks across all replicas
    tasks = []
    for rep in replicas:
        nf = mda.Universe(now_ref, rep).trajectory.n_frames
        tasks += [(rep, fi) for fi in range(0, nf, STRIDE)]
    n_frames = len(tasks)

    with mp.Pool(NPROC, initializer=_init_worker,
                 initargs=(now_ref, prot_sel, grid_xyz)) as pool:
        for cov_idx in pool.imap_unordered(_frame_worker, tasks, chunksize=2):
            if len(cov_idx):
                count[cov_idx] += 1

    freq = count / max(n_frames, 1)
    np.savez(os.path.join(OUTDIR, f"{sid}_pocketgrid.npz"),
             grid_xyz=grid_xyz.astype(np.float32), freq=freq.astype(np.float32))

    ref_heavy = ref_prot.select_atoms("not name H*")
    pockets, npers = call_pockets(grid_xyz, freq, ref_heavy, lig_xyz)
    return _write_summary(sid, row, pockets, npers, n_frames, len(replicas),
                          sorted(lipid_rn), lig_type, lig_names, len(lig_xyz))


def _write_summary(sid, row, pockets, npers, n_frames, n_rep,
                   lipid_rn, lig_type, lig_names, n_lig_atoms):
    ortho_recovered = any(p["is_orthosteric"] for p in pockets)
    summary = dict(system_id=sid, uniprot=row["receptor_uniprot"],
                   family=row["g_protein_family"], n_frames=n_frames,
                   n_replicas=n_rep, lipid_resnames=lipid_rn,
                   ligand_type=lig_type, ligand_resnames=lig_names,
                   n_ligand_atoms=int(n_lig_atoms),
                   n_persistent_voxels=npers, n_pockets=len(pockets),
                   orthosteric_recovered=ortho_recovered, pockets=pockets)
    with open(os.path.join(OUTDIR, f"{sid}_pockets.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[{sid}] frames={n_frames} persistent_voxels={npers} "
          f"pockets={len(pockets)} ortho_recovered={ortho_recovered} "
          f"ligand={lig_type}{lig_names}({n_lig_atoms} atoms)", flush=True)
    return summary


def reprocess_saved(sid, row, thresh=FREQ_THRESH):
    """Re-derive pockets from the saved frequency grid (no fpocket re-run)."""
    now_ref = build_now_ref(os.path.dirname(row["trajectory_path"]))
    ref_u = mda.Universe(now_ref)
    lipid_rn, _ = classify_resnames(ref_u)
    ref_heavy = ref_u.select_atoms("resname " + " ".join(sorted(AA))).select_atoms("not name H*")
    lig_xyz, lig_type, lig_names = get_ligand_coords(ref_u, row)
    d = np.load(os.path.join(OUTDIR, f"{sid}_pocketgrid.npz"))
    pockets, npers = call_pockets(d["grid_xyz"], d["freq"], ref_heavy, lig_xyz, thresh)
    return _write_summary(sid, row, pockets, npers, "from_grid", 3,
                          sorted(lipid_rn), lig_type, lig_names, len(lig_xyz))


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    rows = {r["system_id"]: r for r in csv.DictReader(open(INVENTORY))}
    args = sys.argv[1:]
    post = "--post" in args                      # re-derive from saved grids, no recompute
    args = [a for a in args if a != "--post"]
    todo = args if args else PILOT
    fn = reprocess_saved if post else process_system
    results = [fn(sid, rows[sid]) for sid in todo]
    with open(os.path.join(OUTDIR, "pilot_summary.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("\n=== PILOT SUMMARY ===")
    for r in results:
        print(f"  {r['system_id']:10s} {r['family']:6s} pockets={r['n_pockets']:3d} "
              f"ortho_recovered={r['orthosteric_recovered']}")


if __name__ == "__main__":
    main()
