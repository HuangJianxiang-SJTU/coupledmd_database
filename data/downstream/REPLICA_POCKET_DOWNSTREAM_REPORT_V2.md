# Replica-filtered pocket downstream analysis

Status: **final207**

Systems accepted by the reducer: **207 / 207**.
These values are eligible for updating the manuscript.

## Completed analysis path

- Replica-reproducible receptor-intrinsic system pockets: 426.
- Cross-receptor replica consensus clusters: 54.
- Class-A-only replica consensus clusters: 48.
- Replica-reproducible orthosteric system pockets: 11.
- Orthosteric replica consensus clusters: 3.
- Stable recurrent clusters under the prespecified final rule: ['R0', 'R1', 'R2', 'R3', 'R5', 'R6', 'R7', 'R8', 'R9', 'R11', 'R13'].
- Stable Class-A-only recurrent clusters: ['A0', 'A1', 'A2', 'A3', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11'].
- Prespecified non-redundant known-site families with matched-null q < 0.05: ['IC345', 'IC56', 'NaP', 'OBP', 'TM345/TM35 mid-membrane family', 'UP124'].
- Legacy C7 outcome-blind anatomical candidate: {'replica_cluster': 'R6', 'core_jaccard_to_legacy_C7': 0.75, 'anatomical_replication_eligible': True}.
- Complete same-receptor context contrasts: 17.
- Primary pocket-fingerprint variance fractions: {'GPCR class': 0.02494632752721161, 'receptor identity within class': 0.4879817660542078, 'administrative G-protein family after receptor': 0.0066785446929724535, 'system/context after receptor and family': 0.054741927012904396, 'residual replica variation within system': 0.4256514347127037}.
- Within-receptor family-label permutation p: 0.34865134865134867.

## Interpretation gate

The final script repeats receptor/H40 recurrence, receptor-disjoint split-half
stability, clustering-parameter sensitivity, matched-null site concordance,
the full cluster-by-family max-statistic screen, the prespecified legacy-C7
sensitivity analyses, the three-tier pocket-distance comparison, and a
sequential multivariate pocket-fingerprint variance partition. Final
figure and manuscript replacement remains gated on a 207/207 reducer result.
Curated Class-A known-site definitions are tested primarily against a
Class-A-only reconstructed atlas; the all-cohort result is a sensitivity
analysis, not evidence that Class A and B generic positions are structurally
interchangeable.
