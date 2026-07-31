# Peptide-stripped sensitivity comparison (v2)

Systems: 66 (22/22/22 Gi/Gq/Gs; 49 Class A / 17 Class B).
Baseline: paired local unstripped control (identical frames and fpocket build); primary final207 values reported as context.

Scope: direct peptide-atom occlusion of fpocket only; does not remove ligand-conditioned receptor conformations.

## Headline metrics

- Median fingerprint Jaccard distance, stripped vs control: **0.212** (replica noise 0.562; control-vs-primary recompute floor 0.000)
- Orthosteric recovery flips (vs control): **0/66**
- Reproducible druggable pockets (66 systems): control **129** vs stripped **130** (primary 125)
- Median persistent-voxel relative change: **0.017**
- R6 Gq support: 5/5 retained; stripped Gq risk ratio nan (control nan)

## Material-change triggers

| rule | triggered |
|---|---|
| M1_fingerprint_shift_exceeds_noise | False |
| M2_orthosteric_flips_gt_10pct | False |
| M3_core_support_loss | False |
| M4_r6_lost_or_direction_reversed | False |

**Material change: False**

## Stable-core support (66 affected systems; control baseline)

The material-change gate uses the frozen preregistered cluster set. The current combined-repair primary stable core is projected below as a non-gating claim-impact check.

| cluster | control support | retained | lost | material |
|---|---|---|---|---|
| R0 | 6 | 6 | 0 | False |
| R1 | 14 | 14 | 0 | False |
| R2 | 7 | 7 | 0 | False |
| R3 | 1 | 1 | 0 | False |
| R4 | 0 | 0 | 0 | False |
| R5 | 5 | 5 | 0 | False |
| R6 | 5 | 5 | 0 | False |
| R7 | 2 | 2 | 0 | False |
| R8 | 3 | 3 | 0 | False |
| R9 | 0 | 0 | 0 | False |
| R11 | 0 | 0 | 0 | False |

## Current-primary stable-core projection (non-gating)

| cluster | control support | retained | lost |
|---|---|---|---|
| R0 | 6 | 6 | 0 |
| R1 | 14 | 14 | 0 |
| R2 | 7 | 7 | 0 |
| R3 | 1 | 1 | 0 |
| R5 | 5 | 5 | 0 |
| R6 | 5 | 5 | 0 |
| R7 | 2 | 2 | 0 |
| R8 | 3 | 3 | 0 |
| R9 | 0 | 0 | 0 |
| R11 | 0 | 0 | 0 |
| R13 | 0 | 0 | 0 |
