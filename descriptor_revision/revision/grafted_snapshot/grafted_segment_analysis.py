#!/usr/bin/env python3
"""Recover AlphaFold-grafted residue ranges per CoupledMD system and compute
per-grafted-segment Calpha RMSD over the released replica-1 reduced
trajectories (structure.pdb + traj.xtc, 2500 frames, 200 ps/frame).

Rationale: no record of which residue ranges were grafted exists anywhere in
the pipeline. Grafted/modelled residues are recovered by aligning the
sequence of each simulated protein chain (from the released structure.pdb)
against the experimentally resolved residues of the source RCSB PDB entry:
simulated residues with no experimental counterpart (sequence gaps, construct
deletions, or chains with no match) are the AlphaFold-grafted/modelled
residues. Chains are split from the released PDB by consecutive-Calpha
distance (> 4.2 A = peptide-bond break).

RMSD protocol: for frames after the first 50 ns (frame >= 250, stride 10 ->
2 ns sampling), each frame is least-squares fitted on the Calpha of the
experimentally resolved (non-grafted) core of the SAME chain (global core as
fallback), so inter-domain/inter-chain motion does not inflate segment RMSD.
Per-segment RMSD is then computed against the segment's frame-0
(post-equilibration, production-start) conformation. A self-fitted RMSD
(segment internal rigidity) is also reported.
"""
import csv, os, sys, time, json, warnings, difflib, urllib.request
import numpy as np
warnings.filterwarnings("ignore")

BASE = "/MDdata/data02/jxhuang/gpcr_g"
REL = f"{BASE}/gpcr_g_server/manuscript_figures/scidata_figures/zenodo_reduced_release_208"
MASTER = f"{BASE}/gpcr_g_server/data/systems_master.csv"
OUTDIR = f"{BASE}/coupledmd_database/26-10-revision"
CACHE = os.path.join(OUTDIR, "rcsb_pdb_cache")
os.makedirs(CACHE, exist_ok=True)

CN_BOND_A = 2.0       # C(i)-N(i+1) distance above this => no peptide bond
CA_BREAK_A = 4.2      # fallback CA-CA break threshold when C/N missing
STRIDE = 10           # analyse every 10th frame (2 ns)
START_FRAME = 250       # first 50 ns excluded (250 frames x 200 ps)
LAST_NS = 100           # window for "last 100 ns" mean

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
       'GLY':'G','HIS':'H','HSD':'H','HSE':'H','HSP':'H','ILE':'I','LEU':'L',
       'LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}

def parse_exp_pdb(path):
    chains = {}
    for line in open(path, errors="ignore"):
        if line[:6].strip() == "ATOM" and line[12:16].strip() == "CA":
            ch = line[21].strip() or "_"
            key = line[22:27].strip()
            rn = line[17:20].strip()
            lst = chains.setdefault(ch, [])
            if not lst or lst[-1][0] != key:
                lst.append((key, rn))
    return chains

def parse_exp_cif(path):
    import gemmi
    st = gemmi.read_structure(path)
    chains = {}
    for model in st:
        for ch in model:
            lst = []
            for res in ch:
                if res.find_atom("CA", "\0") is not None:
                    lst.append((str(res.seqid.num) + res.seqid.icode.strip(), res.name))
            if lst:
                chains[ch.name] = lst
    return chains

def fetch_pdb(pdb_id):
    p = os.path.join(CACHE, f"{pdb_id}.pdb")
    if os.path.exists(p) and os.path.getsize(p) > 0:
        return p, "pdb"
    try:
        with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60) as r:
            data = r.read()
        if not data.startswith(b"HEADER"):
            raise ValueError("not pdb format")
        open(p, "wb").write(data)
        return p, "pdb"
    except Exception:
        c = os.path.join(CACHE, f"{pdb_id}.cif")
        if not os.path.exists(c):
            with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=60) as r:
                open(c, "wb").write(r.read())
        return c, "cif"

def split_chains(prot):
    """Split a protein AtomGroup into chains by peptide-bond breaks:
    C(i)-N(i+1) distance > 2.0 A (fallback: CA-CA > 4.2 A).
    Returns list of [(resid, resname), ...] per chain plus break indices."""
    residues = list(prot.residues)
    pos = {name: np.array([r.atoms.select_atoms(f"name {name}").positions[0]
                           if len(r.atoms.select_atoms(f"name {name}")) else None
                           for r in residues], dtype=object)
           for name in ("C", "N", "CA")}
    breaks = []
    for i in range(len(residues) - 1):
        c, n = pos["C"][i], pos["N"][i + 1]
        if c is not None and n is not None:
            if np.linalg.norm(c - n) > CN_BOND_A:
                breaks.append(i + 1)
        else:
            ca1, ca2 = pos["CA"][i], pos["CA"][i + 1]
            if ca1 is None or ca2 is None or np.linalg.norm(ca1 - ca2) > CA_BREAK_A:
                breaks.append(i + 1)
    chains, prev = [], 0
    for b in breaks + [len(residues)]:
        chains.append([(str(r.resid), r.resname) for r in residues[prev:b]])
        prev = b
    return chains

def map_grafted(sim_res, exp_chains):
    """sim_res: [(resid, resname)]. Returns (best_chain, grafted_local_idx, n_mut, ratio)."""
    sim_seq = "".join(AA3.get(rn, "X") for _, rn in sim_res)
    best, best_ratio = None, 0.0
    for ch, lst in exp_chains.items():
        exp_seq = "".join(AA3.get(rn, "X") for _, rn in lst)
        r = difflib.SequenceMatcher(None, sim_seq, exp_seq, autojunk=False).ratio()
        if r > best_ratio:
            best_ratio, best = r, ch
    if best is None or best_ratio < 0.25:
        return None, list(range(len(sim_res))), 0, best_ratio
    exp_seq = "".join(AA3.get(rn, "X") for _, rn in exp_chains[best])
    sm = difflib.SequenceMatcher(None, sim_seq, exp_seq, autojunk=False)
    mapped = np.zeros(len(sim_res), bool)
    nmut = 0
    for m in sm.get_matching_blocks():
        for k in range(m.size):
            mapped[m.a + k] = True
            if sim_res[m.a + k][1] != exp_chains[best][m.b + k][1]:
                nmut += 1
    grafted = [i for i in range(len(sim_res)) if not mapped[i]]
    if grafted:  # merge grafted runs separated by <=3 mapped residues
        merged, run = [], [grafted[0]]
        for i in grafted[1:]:
            if i - run[-1] <= 4:
                run.extend(range(run[-1] + 1, i + 1))
            else:
                merged.append(run); run = [i]
        merged.append(run)
        grafted = sorted(set(sum(merged, [])))
    return best, grafted, nmut, best_ratio

def kabsch_rmsd(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return np.sqrt(((Pc @ R.T - Qc) ** 2).sum(-1).mean())

def process_system(row):
    import MDAnalysis as mda
    sysid, pdb_id, fam = row["system_id"], row["pdb_id"], row["g_protein_family"]
    out = {"system_id": sysid, "pdb_id": pdb_id, "family": fam,
           "gpcr_class": master.get(sysid, {}).get("gpcr_class", ""),
           "chains": [], "segments": [], "error": None}
    try:
        u = mda.Universe(os.path.join(REL, sysid, "structure.pdb"),
                         os.path.join(REL, sysid, "traj.xtc"))
        prot_all = u.select_atoms("protein")
        prot = prot_all.select_atoms("name CA")
        if prot.n_atoms < 10:
            out["error"] = f"too few CA atoms ({prot.n_atoms})"; return out
        chains = split_chains(prot_all)
        n_chain_res = sum(len(c) for c in chains)
        if n_chain_res != prot.n_atoms:
            out["error"] = f"residue/CA mismatch {n_chain_res}/{prot.n_atoms}"; return out
        ca0 = prot.positions.copy()
        ep, kind = fetch_pdb(pdb_id)
        exp = parse_exp_pdb(ep) if kind == "pdb" else parse_exp_cif(ep)
        if not exp:
            out["error"] = "no experimental chains parsed"; return out
        graft_mask = np.zeros(prot.n_atoms, bool)
        seg_records = []
        out["n_singleton_graft_residues"] = 0
        idx0 = 0
        for ci, reslist in enumerate(chains):
            chain, gidx, nmut, ratio = map_grafted(reslist, exp)
            out["chains"].append({"chain_index": ci, "length": len(reslist),
                                  "matched_exp_chain": chain or "NONE",
                                  "match_ratio": round(ratio, 3),
                                  "n_grafted_residues": len(gidx)})
            for i in gidx:
                graft_mask[idx0 + i] = True
            runs = []
            for i in gidx:
                if runs and i == runs[-1][-1] + 1:
                    runs[-1].append(i)
                else:
                    runs.append([i])
            for run in runs:
                if len(run) < 2:
                    out["n_singleton_graft_residues"] += len(run)
                    continue  # isolated unmatched residues: alignment noise, not segments
                loc = "loop"
                if run[0] == 0: loc = "N-term"
                if run[-1] == len(reslist) - 1:
                    loc = "C-term" if run[0] != 0 else "whole-chain"
                elif len(run) >= 40:
                    loc = "large-insert"
                seg_records.append(dict(
                    chain_index=ci, matched_chain=chain or "NONE",
                    match_ratio=round(ratio, 3), n_mutations=nmut, location=loc,
                    sim_resid_start=reslist[run[0]][0], sim_resid_end=reslist[run[-1]][0],
                    length=len(run), ca_start=idx0 + run[0], ca_end=idx0 + run[-1],
                ))
            idx0 += len(reslist)
        out["n_residues"] = prot.n_atoms
        out["n_grafted_residues"] = int(graft_mask.sum())
        out["n_chains"] = len(chains)
        if not seg_records:
            return out
        ref = ca0
        # per-chain resolved core for fitting (fall back to global core)
        global_core = np.where(~graft_mask)[0]
        chain_bounds = []
        i0 = 0
        for ci, reslist in enumerate(chains):
            chain_bounds.append((i0, i0 + len(reslist))); i0 += len(reslist)
        def core_for(ci):
            a, b = chain_bounds[ci]
            local = np.array([i for i in range(a, b) if not graft_mask[i]])
            return local if len(local) >= 10 else global_core
        seg_fit_core = [core_for(s["chain_index"]) for s in seg_records]
        seg_idx = [np.arange(s["ca_start"], s["ca_end"] + 1) for s in seg_records]
        seg_ref = [ref[ix] for ix in seg_idx]
        seg_ref_core = [ref[c] for c in seg_fit_core]
        nfr = len(u.trajectory)
        frames = list(range(START_FRAME, nfr, STRIDE))
        last100_start = max(0, nfr - int(LAST_NS * 1000 / u.trajectory.dt))
        nseg = len(seg_records)
        rmsd_all = np.zeros((nseg, len(frames)))
        selffit_all = np.zeros_like(rmsd_all)
        frame_list = []
        for fi, ts in enumerate(u.trajectory[START_FRAME:nfr:STRIDE]):
            frame_list.append(ts.frame)
            mob = prot.positions
            for si, ix in enumerate(seg_idx):
                ci_ = seg_fit_core[si]
                Pc = mob[ci_] - mob[ci_].mean(0)
                Qc = seg_ref_core[si] - seg_ref_core[si].mean(0)
                U_, S_, Vt_ = np.linalg.svd(Pc.T @ Qc)
                d = np.sign(np.linalg.det(Vt_.T @ U_.T))
                R = Vt_.T @ np.diag([1, 1, d]) @ U_.T
                seg_mob_fit = (mob[ix] - mob[ci_].mean(0)) @ R.T + seg_ref_core[si].mean(0)
                rmsd_all[si, fi] = np.sqrt(((seg_mob_fit - seg_ref[si]) ** 2).sum(-1).mean())
                selffit_all[si, fi] = kabsch_rmsd(mob[ix], seg_ref[si])
        for si, s in enumerate(seg_records):
            r, sf = rmsd_all[si], selffit_all[si]
            tail_mask = np.array(frame_list) >= last100_start
            s.update(dict(
                rmsd_mean=round(float(r.mean()), 3), rmsd_max=round(float(r.max()), 3),
                rmsd_last100ns_mean=round(float(r[tail_mask].mean()), 3),
                selffit_rmsd_mean=round(float(sf.mean()), 3),
                selffit_rmsd_max=round(float(sf.max()), 3),
                n_frames=r.size,
            ))
            del s["ca_start"]; del s["ca_end"]
        out["segments"] = seg_records
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
    return out

if __name__ == "__main__":
    global master
    master = {r["system_id"]: r for r in csv.DictReader(open(MASTER))}
    manifest = list(csv.DictReader(open(os.path.join(REL, "CoupledMD_reduced_trajectory_manifest.csv"))))
    only = sys.argv[1:]
    if only:
        manifest = [r for r in manifest if r["system_id"] in only]
    from multiprocessing import Pool
    t0 = time.time()
    with Pool(8) as pool:
        results = pool.map(process_system, manifest)
    print(f"done {len(results)} systems in {time.time()-t0:.0f}s", file=sys.stderr)
    fields = ["system_id", "pdb_id", "family", "gpcr_class", "chain_index",
              "matched_chain", "match_ratio", "n_mutations", "location",
              "sim_resid_start", "sim_resid_end", "length",
              "rmsd_mean", "rmsd_max", "rmsd_last100ns_mean",
              "selffit_rmsd_mean", "selffit_rmsd_max", "n_frames"]
    with open(os.path.join(OUTDIR, "grafted_segment_rmsd_summary_v2.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for r in results:
            for s in r["segments"]:
                w.writerow({**{k: r[k] for k in ("system_id", "pdb_id", "family", "gpcr_class")}, **s})
    with open(os.path.join(OUTDIR, "grafted_segment_run_log_v2.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["system_id", "error", "n_residues", "n_chains",
                    "n_grafted_residues", "n_segments"])
        for r in results:
            w.writerow([r["system_id"], r["error"], r.get("n_residues"),
                        r.get("n_chains"), r.get("n_grafted_residues"),
                        len(r["segments"])])
    with open(os.path.join(OUTDIR, "grafted_segments_recovered_v2.json"), "w") as f:
        json.dump([{k: r[k] for k in ("system_id", "pdb_id", "family", "chains",
                                      "segments", "error")} for r in results], f, indent=1)
    errs = [(r["system_id"], r["error"]) for r in results if r["error"]]
    print(f"errors: {len(errs)}", file=sys.stderr)
    for s, e in errs:
        print(" ", s, e, file=sys.stderr)
