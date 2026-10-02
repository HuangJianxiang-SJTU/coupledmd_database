# Supplementary Data S7

**system_id**: Frozen system identifier.

**replica**: Computational repeat, 1–3.

**status**: ok: assessed; no_selection: Gs_8HTI has no mapped TM core; not assessed, not nonstationary.

**reason**: Reason for unavailable assessment, where applicable.

**t_eq_ns**: Operational stationarity onset: 10-ns TM-core Cα RMSD windows within ±0.5 Å of final-100-ns reference and absolute fitted slope ≤0.1 Å/100 ns. Blank if nonstationary or not assessed. Not proof of ensemble convergence.

**plateau_rmsd_A**: TM-core reference RMSD (Å) from the final-100-ns window means; this name does not establish a plateau for flagged records.

**nonstationary_flag**: True: fails the operational test; False: passes; blank: not assessed.

**recommended_discard_ns**: Maximum t_eq among stationary replicas only. Does not make flagged replicas suitable for analysis. Exclude or assess flagged replicas separately. Blank for Gs_8HTI.

**n_nonstationary**: Number of flagged replicas in this system; blank for unassessed Gs_8HTI. Three systems have two flagged replicas, so their discard is based on one passing replica.

**limitations**: RMSD is one observable. Passing late in a 500-ns run may leave very little usable sampling. Supplied S2/S5/S6 describe their original full-window analyses and have not been recomputed after discard.