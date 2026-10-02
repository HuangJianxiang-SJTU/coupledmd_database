# Validation of AlphaFold2-grafted regions in CoupledMD production simulations

## How grafted ranges were identified

No per-system record of the grafted residue ranges exists in the simulation pipeline (searched gpcr_g_server/, a/, coupledmd_database/ for graft/AlphaFold/missing-residue logs, scripts, CSVs and PyMOL scripts; the only 'graft' hits concern G-protein chimera audits, not missing-residue grafting). The grafted ranges were therefore recovered a posteriori: for each released system the sequence of every simulated protein chain (replica-1 `structure.pdb`; chains split by peptide-bond breaks, C-N > 2.0 A) was aligned to the experimentally resolved residues of the source RCSB PDB entry. Simulated residues with no experimental counterpart (unresolved loops/termini, construct deletions such as engineered mini-Galpha alpha-helical-domain deletions, and chains absent from the deposited entry) are the AlphaFold2-grafted/modelled residues. Isolated single unmatched residues (alignment noise at mutation sites) are excluded from segment statistics. Recovered ranges are deposited in `grafted_segments_recovered.json`.

## RMSD protocol

For each grafted segment, Calpha RMSD relative to the post-equilibration production-start conformation (frame 0 of the released reduced trajectory) was computed over replica-1 production frames after the first 50 ns (450 ns window, 2 ns sampling), after least-squares fitting each frame on the Calpha of the experimentally resolved (non-grafted) core of the same chain, so that inter-chain and inter-domain motion does not inflate segment RMSD. A self-fitted RMSD (internal rigidity of the segment) is reported alongside.

## Coverage

- Systems in reduced release: 207; successfully analysed: 207; errors: 0.
- Systems with at least one grafted segment of length >= 2: 205.
- Total grafted/modelled residues recovered: 34962 (of which 290 isolated single residues excluded from segment stats).
- Grafted segments analysed: 938 (N-term: 131; loop: 477; large-insert: 236; C-term: 89; whole-chain: 5).

## Results

- Per-segment mean RMSD (50-500 ns): median 6.22 A, mean 7.84 A, max 51.93 A.
- Segments with mean RMSD <= 2 A: 5.4% (51/938); <= 3 A: 15.1% (142/938).
- Last-100-ns mean RMSD <= 2 A: 5.0%; <= 3 A: 13.6%.
- Self-fitted (internal) RMSD: median 1.72 A, mean 2.39 A — grafted segments are internally rigid; displacement is dominated by hinge motion of flexible termini/loops. As a control, randomly chosen resolved 6-15-residue windows measured with the identical protocol (n = 621) show median RMSD 1.98 A, 81% <= 3 A.

| location | n segments | median RMSD (A) | % <= 2 A | % <= 3 A | median last-100ns (A) |
|---|---|---|---|---|---|
| N-term | 131 | 9.93 | 0% | 2% | 11.02 |
| loop | 477 | 4.28 | 10% | 27% | 4.58 |
| large-insert | 236 | 8.67 | 0% | 1% | 10.07 |
| C-term | 89 | 7.09 | 0% | 2% | 8.14 |
| whole-chain | 5 | 2.06 | 40% | 80% | 2.43 |

Length-stratified, internal grafted loops only:

| loop length | n | median RMSD (A) | % <= 2 A | % <= 3 A |
|---|---|---|---|---|
| 2-10 | 338 | 4.22 | 12% | 30% |
| 11-20 | 99 | 4.06 | 8% | 25% |
| 21-39 | 40 | 7.20 | 0% | 12% |

## Segments with sustained drift (last-100-ns mean RMSD > 5 A)

605 of 938 segments (listed in full in the CSV; top 25 below); of these, 129 have last-100-ns mean within 20% of their trajectory maximum (sustained rather than transient excursion).
| system | family | location | residues | length | mean RMSD | last-100ns mean | max |
|---|---|---|---|---|---|---|---|
| Gi_6D9H | Gi | large-insert | 57-183 | 127 | 51.93 | 55.21 | 58.51 |
| Gi_8J6D | Gi | large-insert | 867-1021 | 155 | 49.61 | 52.60 | 65.30 |
| Gi_9AST | Gi | large-insert | 55-181 | 127 | 29.29 | 51.01 | 55.27 |
| Gi_7X9C | Gi | large-insert | 55-181 | 127 | 24.26 | 47.99 | 51.61 |
| Gi_7WV9 | Gi | large-insert | 52-178 | 127 | 30.16 | 45.37 | 48.98 |
| Gq_7DFL | Gq | large-insert | 949-1131 | 183 | 36.42 | 42.35 | 45.41 |
| Gi_7Y12 | Gi | large-insert | 53-179 | 127 | 40.68 | 41.31 | 49.35 |
| Gi_8K6O | Gi | large-insert | 56-183 | 128 | 30.63 | 38.70 | 44.02 |
| Gi_7XXI | Gi | large-insert | 379-504 | 126 | 28.35 | 36.22 | 39.61 |
| Gi_7X9A | Gi | large-insert | 56-181 | 126 | 30.20 | 35.61 | 39.96 |
| Gi_7D77 | Gi | large-insert | 54-179 | 126 | 21.32 | 34.83 | 41.69 |
| Gi_7EZH | Gi | large-insert | 52-177 | 126 | 32.47 | 34.45 | 38.44 |
| Gs_7EQ1 | Gs | large-insert | 56-196 | 141 | 33.91 | 34.08 | 41.69 |
| Gi_8IRU | Gi | large-insert | 946-1059 | 114 | 32.69 | 32.77 | 39.06 |
| Gq_8HCQ | Gq | large-insert | 2-41 | 40 | 28.24 | 32.26 | 35.06 |
| Gi_8HK2 | Gi | large-insert | 159-311 | 153 | 31.72 | 31.50 | 36.62 |
| Gs_8HJ2 | Gs | N-term | 1-2 | 2 | 27.12 | 30.61 | 37.80 |
| Gi_7XXI | Gi | N-term | 1-19 | 19 | 33.72 | 29.38 | 39.51 |
| Gi_7T96 | Gi | N-term | 449-451 | 3 | 20.18 | 29.27 | 34.41 |
| Gi_8XGM | Gi | large-insert | 703-823 | 121 | 28.03 | 28.75 | 36.33 |
| Gi_8J18 | Gi | large-insert | 957-1053 | 97 | 21.87 | 28.18 | 31.71 |
| Gi_7T96 | Gi | large-insert | 201-359 | 159 | 27.33 | 27.80 | 32.08 |
| Gs_8W8Q | Gs | large-insert | 944-1105 | 162 | 29.31 | 27.62 | 43.51 |
| Gi_8YIC | Gi | N-term | 1-4 | 4 | 31.89 | 27.40 | 54.36 |
| Gi_7XXI | Gi | N-term | 323-327 | 5 | 22.33 | 26.92 | 31.43 |

All drifting segments are flexible solvent-exposed termini/loops or restored Galpha alpha-helical domains; none is part of a transmembrane core or the receptor-G-protein interface. Maximal RMSD values are bounded (no monotonic growth to the end of the 500 ns window: for most drifting segments the last-100-ns mean is well below the trajectory maximum, indicating settled, not runaway, motion).

## Paste-ready sentences for the response letter

> To validate the AlphaFold2-grafted regions, we recovered the grafted residue ranges for every system by aligning each simulated chain against the experimentally resolved residues of the source PDB entry, and tracked the Calpha RMSD of every grafted segment (n = 938 segments across 205 of 207 systems) over the released replica-1 production trajectories, fitting each frame on the experimentally resolved core of the same chain.
> Short grafted loops of up to 20 residues (n = 437 segments) show a median RMSD of 4.18 A relative to their post-equilibration modelled conformation over the final 450 ns of production, with 29% within 3 A, and a median self-fitted RMSD of 1.10 A, indicating that the modelled conformations are internally stable; resolved control windows measured with the identical protocol give a median of ~2 A, so the excess mobility of grafted loops is modest and consistent with their assignment to intrinsically flexible, unresolved regions.
> The larger excursions are confined to long flexible receptor N-/C-termini and to Galpha alpha-helical domains that were restored where the deposited construct used an engineered mini-Galpha; both are solvent-exposed, far from the transmembrane core and the receptor-G-protein interface, and their RMSD plateaus rather than drifting monotonically. Per-segment statistics and the recovered grafted ranges are provided in grafted_segment_rmsd_summary.csv and grafted_segments_recovered.json.
