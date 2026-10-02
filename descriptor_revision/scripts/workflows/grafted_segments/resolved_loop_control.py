#!/usr/bin/env python3
"""Control: RMSD of RESOLVED loop segments, measured exactly like grafted
segments (segment excluded from fit; fit on remaining same-chain resolved
core; RMSD vs frame 0 after 50 ns, stride 10)."""
import csv, os, sys, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, "/MDdata/data02/jxhuang/gpcr_g/coupledmd_database/26-10-revision")
from grafted_segment_analysis import (REL, OUTDIR, split_chains, fetch_pdb,
                                      parse_exp_pdb, parse_exp_cif, map_grafted,
                                      kabsch_rmsd, START_FRAME, STRIDE, LAST_NS)

rng = np.random.default_rng(7)

def process(row):
    import MDAnalysis as mda
    sysid, pdb_id = row["system_id"], row["pdb_id"]
    out = []
    try:
        u = mda.Universe(os.path.join(REL, sysid, "structure.pdb"),
                         os.path.join(REL, sysid, "traj.xtc"))
        prot_all = u.select_atoms("protein")
        prot = prot_all.select_atoms("name CA")
        chains = split_chains(prot_all)
        ep, kind = fetch_pdb(pdb_id)
        exp = parse_exp_pdb(ep) if kind == "pdb" else parse_exp_cif(ep)
        ref = prot.positions.copy()
        graft_mask = np.zeros(prot.n_atoms, bool)
        idx0 = 0
        for reslist in chains:
            _, gidx, _, ratio = map_grafted(reslist, exp)
            for i in gidx:
                graft_mask[idx0 + i] = True
            idx0 += len(reslist)
        # candidate resolved internal windows len 6-15
        idx0 = 0
        cands = []
        for ci, reslist in enumerate(chains):
            n = len(reslist)
            resolved = [i for i in range(n) if not graft_mask[idx0 + i]]
            if len(resolved) < 40:
                idx0 += n; continue
            for _ in range(60):
                L = int(rng.integers(6, 16))
                s = int(rng.integers(5, n - 5 - L))
                win = list(range(s, s + L))
                if all(not graft_mask[idx0 + i] for i in win):
                    cands.append((ci, idx0 + s, idx0 + s + L - 1))
            idx0 += n
        if not cands:
            return out
        picks = [cands[i] for i in rng.choice(len(cands), size=min(3, len(cands)), replace=False)]
        nfr = len(u.trajectory)
        frames = list(range(START_FRAME, nfr, STRIDE))
        last100_start = max(0, nfr - int(LAST_NS * 1000 / u.trajectory.dt))
        res = []
        for ci, a, b in picks:
            res.append(dict(ci=ci, ix=np.arange(a, b + 1), ref=ref[a:b + 1],
                            core=np.array([i for i in range(sum(len(c) for c in chains[:ci]),
                                                           sum(len(c) for c in chains[:ci]) + len(chains[ci]))
                                           if not graft_mask[i] and not (a <= i <= b)]),
                            acc=[]))
        for r in res:
            r["ref_core"] = ref[r["core"]]
        frame_list = []
        for fi, ts in enumerate(u.trajectory[START_FRAME:nfr:STRIDE]):
            frame_list.append(ts.frame)
            mob = prot.positions
            for r in res:
                c, rc = r["core"], r["ref_core"]
                Pc = mob[c] - mob[c].mean(0); Qc = rc - rc.mean(0)
                U_, S_, Vt_ = np.linalg.svd(Pc.T @ Qc)
                d = np.sign(np.linalg.det(Vt_.T @ U_.T))
                R = Vt_.T @ np.diag([1, 1, d]) @ U_.T
                fit = (mob[r["ix"]] - mob[c].mean(0)) @ R.T + rc.mean(0)
                r["acc"].append(np.sqrt(((fit - r["ref"]) ** 2).sum(-1).mean()))
        frame_list = np.array(frame_list)
        for r in res:
            acc = np.array(r["acc"])
            out.append(dict(system_id=sysid, pdb_id=pdb_id,
                            family=row["g_protein_family"], length=len(r["ix"]),
                            rmsd_mean=round(float(acc.mean()), 3),
                            rmsd_max=round(float(acc.max()), 3),
                            rmsd_last100ns_mean=round(float(acc[frame_list >= last100_start].mean()), 3)))
    except Exception as e:
        print(sysid, type(e).__name__, e, file=sys.stderr)
    return out

if __name__ == "__main__":
    manifest = list(csv.DictReader(open(os.path.join(REL, "CoupledMD_reduced_trajectory_manifest.csv"))))
    from multiprocessing import Pool
    with Pool(8) as pool:
        results = pool.map(process, manifest)
    flat = [r for sub in results for r in sub]
    with open(os.path.join(OUTDIR, "resolved_loop_control_rmsd.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(flat[0].keys())); w.writeheader()
        w.writerows(flat)
    m = np.array([r["rmsd_mean"] for r in flat])
    print(f"n={len(flat)} median={np.median(m):.2f} <=2A:{100*(m<=2).mean():.0f}% <=3A:{100*(m<=3).mean():.0f}%", file=sys.stderr)
