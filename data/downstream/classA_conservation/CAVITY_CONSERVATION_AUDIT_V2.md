# Class-A cavity conservation overlay (v2)

Status: **final207 replica-reproducible analysis**.

Eligible for manuscript interpretation:
**yes**.

## Scope

Canonical UniProt amino acids were aligned through the best matching
structure-specific GPCRdb record for each of
151 deduplicated Class-A receptor
accessions. Conservation is calculated after giving each of
62 40%-TM-sequence homology groups equal total
weight. Positions require coverage in at least 15 H40 groups.
Each stable-core position is compared with 20
positions matched by GPCR topology category (individual TM, H8, or loops) and
closest H40 coverage; 10,000 matched null draws are used.

Cluster selection: n_receptors >= 5 and receptor-split recovery probability >= 0.80 at best-core Jaccard >= 0.50.

## Result

- Selected Class-A stable clusters: **11**
  (A0, A1, A2, A3, A5, A6, A7, A8, A9, A10, A11).
- Unique stable-core positions: **127**.
- Stable-core union mean H40 conservation:
  **0.365**.
- Matched-null mean and 95% interval:
  **0.360**
  (0.330–
  0.392).
- Union empirical two-sided p:
  **0.7649**.
- Conservation-enriched cluster-level BH survivors:
  **A6, A11**.
- Conservation-depleted cluster-level BH survivors:
  **A0**.
- Cluster recurrence versus mean core conservation:
  **Spearman ρ=-0.17, p=0.626**.

## Interpretation rule

An enrichment would mean recurrent cavities preferentially line positions more
conserved than topology- and coverage-matched Class-A background. It would not
establish functional conservation or ligandability. A null would mean that
recurrence is not simply a consequence of unusually conserved lining
positions. These final-replica values are eligible for manuscript interpretation after figure and table inspection.

## Outputs

- `source_data/classA_generic_position_conservation_v2.csv`
- `source_data/cavity_classA_conservation_matched_null_v2.csv`
- `source_data/cavity_classA_core_position_conservation_v2.csv`
- `diagnostic_cavity_conservation_v2.pdf/png`
- `cavity_conservation_summary_v2.json`
- generator: `/MDdata/data02/jxhuang/gpcr_g/a/paper1_cavity_conservation_v2.py`
