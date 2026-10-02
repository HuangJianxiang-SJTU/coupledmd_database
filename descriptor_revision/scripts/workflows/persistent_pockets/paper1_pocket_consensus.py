#!/usr/bin/env python3
"""
paper1_pocket_consensus.py  --  WP1 scaffold: location/segment post-filter +
cross-receptor consensus pocket clustering.

Reads every {sid}_pockets_gpcrdb.json produced by paper1_pocket_gpcrdb_map.py and:

  1. POST-FILTER by membrane-topology zone. The whole-complex fpocket inventory is
     dominated by the intracellular G-protein coupling cavity (zone=coupling_interface);
     those are not druggable allosteric sites, so by default they are excluded. The
     druggable-allosteric set = receptor pockets in the extracellular vestibule, the
     TM core, or the intracellular (non-coupling) face. Orthosteric pockets are kept
     in a separate track for the known-site recovery audit.

  2. CONSENSUS CLUSTERING across receptors. Each surviving pocket is represented by its
     set of GPCRdb generic numbers (which ARE comparable across receptors: '3.50x50' is
     the same structural position in any receptor). Pockets are clustered by Jaccard
     distance of their generic-number sets (average-linkage agglomerative, distance
     threshold JTHRESH). Each cluster = a consensus pocket family, summarized by its
     core generic positions, the receptors/families it spans, and its zone.

This is a SCAFFOLD: it runs on whatever systems are mapped so far (41 at last run) and
re-runs unchanged atlas-wide when the sweep finishes. Cluster identities firm up as N
grows; treat current clusters as provisional.

Usage:
    python paper1_pocket_consensus.py              # druggable-allosteric track
    python paper1_pocket_consensus.py orthosteric  # orthosteric track (recovery audit)
"""
import os, sys, csv, json, glob
from collections import Counter, defaultdict
import numpy as np

try:
    from sklearn.cluster import AgglomerativeClustering
    HAVE_SK = True
except Exception:
    HAVE_SK = False

BASE = "/MDdata/data02/jxhuang/gpcr_g/a"
INVENTORY = os.path.join(BASE, "inventory_paper1_207.csv")
ROWS = {r["system_id"]: r for r in csv.DictReader(open(INVENTORY))}
FINAL_COHORT = "/MDdata/data02/jxhuang/gpcr_g/gpcr_g_server/data/release_cohort_v9_final207.csv"
FINAL_IDS = {r["system_id"] for r in csv.DictReader(open(FINAL_COHORT))}
assert len(FINAL_IDS) == 207

DRUGGABLE_ZONES = {"extracellular_vestibule", "tm_core_allosteric", "intracellular_allosteric"}
MIN_GENERICS = 4        # need >=4 mapped positions for a stable cross-receptor fingerprint
JTHRESH = 0.72          # Jaccard distance cut for "same" consensus pocket
MIN_MEMBERS = 2         # a consensus family needs >=2 pockets


def collect_pockets(track):
    """Return list of pocket dicts surviving the post-filter for the given track."""
    files = sorted(glob.glob(os.path.join(BASE, "paper1_pockets/atlas/*_pockets_gpcrdb.json")) +
                   glob.glob(os.path.join(BASE, "paper1_pockets/pilot/*_pockets_gpcrdb.json")))
    pockets, seen = [], set()
    for f in files:
        d = json.load(open(f))
        sid = d["system_id"]
        if sid not in FINAL_IDS:
            continue
        if sid in seen:
            continue
        seen.add(sid)
        row = ROWS.get(sid, {})
        for p in d["pockets"]:
            keep = (p["is_orthosteric"] if track == "orthosteric"
                    else (not p["is_orthosteric"] and p.get("zone") in DRUGGABLE_ZONES))
            if not keep:
                continue
            gens = frozenset(p.get("receptor_generic_numbers", []))
            if len(gens) < MIN_GENERICS:
                continue
            pockets.append(dict(
                sid=sid, uniprot=d.get("uniprot", ""), family=d.get("family", ""),
                receptor=row.get("receptor_name", ""),
                pocket_id=p["pocket_id"], zone=p.get("zone"),
                mean_freq=p.get("mean_freq"), n_voxels=p.get("n_voxels"),
                segments=p.get("receptor_segments", []), generics=gens))
    return pockets


def jaccard_matrix(sets):
    n = len(sets)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = sets[i], sets[j]
            u = len(a | b)
            d = 1.0 - (len(a & b) / u if u else 0.0)
            D[i, j] = D[j, i] = d
    return D


def cluster(pockets):
    if len(pockets) < 2:
        return [0] * len(pockets)
    D = jaccard_matrix([p["generics"] for p in pockets])
    if HAVE_SK:
        cl = AgglomerativeClustering(n_clusters=None, metric="precomputed",
                                     linkage="average", distance_threshold=JTHRESH)
        return cl.fit_predict(D).tolist()
    # greedy fallback (no sklearn): single-link merge under JTHRESH
    labels = list(range(len(pockets)))
    for i in range(len(pockets)):
        for j in range(i + 1, len(pockets)):
            if D[i, j] < JTHRESH:
                old, new = labels[j], labels[i]
                labels = [new if x == old else x for x in labels]
    return labels


def summarize(pockets, labels):
    clusters = defaultdict(list)
    for p, lab in zip(pockets, labels):
        clusters[lab].append(p)
    out = []
    for lab, members in clusters.items():
        if len(members) < MIN_MEMBERS:
            continue
        gen_count = Counter()
        for m in members:
            gen_count.update(m["generics"])
        nm = len(members)
        core = sorted([g for g, c in gen_count.items() if c >= 0.5 * nm],
                      key=lambda g: -gen_count[g])
        out.append(dict(
            consensus_id=int(lab),
            n_pockets=nm,
            n_systems=len(set(m["sid"] for m in members)),
            n_receptors=len(set(m["uniprot"] for m in members)),
            families=dict(Counter(m["family"] for m in members)),
            zones=dict(Counter(m["zone"] for m in members)),
            core_generic_numbers=core,
            all_generic_numbers=sorted(gen_count, key=lambda g: -gen_count[g]),
            mean_freq=round(float(np.mean([m["mean_freq"] for m in members])), 3),
            member_systems=sorted(set(m["sid"] for m in members)),
        ))
    out.sort(key=lambda c: (-c["n_receptors"], -c["n_pockets"]))
    return out


def main():
    track = sys.argv[1] if len(sys.argv) > 1 else "druggable"
    pockets = collect_pockets(track)
    print(f"[{track}] {len(pockets)} pockets pass the post-filter "
          f"(zones={DRUGGABLE_ZONES if track!='orthosteric' else 'orthosteric'}, "
          f">= {MIN_GENERICS} generics)")
    if not pockets:
        print("  nothing to cluster yet"); return
    labels = cluster(pockets)
    cons = summarize(pockets, labels)
    print(f"  {len(set(labels))} raw clusters -> {len(cons)} consensus families "
          f"(>= {MIN_MEMBERS} members)\n")
    print(f"  Top consensus pocket families (by # distinct receptors):")
    for c in cons[:15]:
        seg_hint = "/".join(sorted({g.split('.')[0] for g in c["core_generic_numbers"]}))
        print(f"   C{c['consensus_id']:<3d} recs={c['n_receptors']:2d} pockets={c['n_pockets']:2d} "
              f"fam={c['families']} TM~{seg_hint}")
        print(f"        core: {', '.join(c['core_generic_numbers'][:10])}")

    outj = os.path.join(BASE, f"paper1_pockets/consensus_{track}.json")
    json.dump(cons, open(outj, "w"), indent=2)
    outc = os.path.join(BASE, f"paper1_pockets/consensus_{track}.csv")
    with open(outc, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["consensus_id", "n_pockets", "n_systems", "n_receptors",
                    "families", "zones", "mean_freq", "core_generic_numbers", "member_systems"])
        for c in cons:
            w.writerow([c["consensus_id"], c["n_pockets"], c["n_systems"], c["n_receptors"],
                        json.dumps(c["families"]), json.dumps(c["zones"]), c["mean_freq"],
                        ";".join(c["core_generic_numbers"]), ";".join(c["member_systems"])])
    print(f"\n  wrote {outj}\n        {outc}")


if __name__ == "__main__":
    main()
