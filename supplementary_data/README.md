# Replica-aware GPCR cavity atlas — Supplementary Data

This directory contains the machine-readable data files supporting the
manuscript's supplementary information (SI v11). All files have been
renamed to remove internal `_v2` and `_v9` version suffixes; the content
matches the SI v11 artifact-identifier and machine-readable-source
references.

---

## Structure

### [S01_system_inventory/](./S01_system_inventory/) — SI Table S1

| File | Description |
|------|-------------|
| `Supplementary_Data_S1_system_inventory.csv` | Released cohort inventory: 207 systems, receptor accessions, administrative G-protein families, ligand presence, production protocol, and receptor-identity notes |
| `Supplementary_Data_S1_system_inventory_audit.json` | Provenance and correction audit for the inventory |

### [S02_fingerprint_variance_partition/](./S02_fingerprint_variance_partition/) — SI Table S15

| File | Description |
|------|-------------|
| `Supplementary_Data_S2_fingerprint_variance_partition.csv` | Complete five-term sequential Hellinger-fingerprint variance partition: GPCR class, receptor identity within class, administrative G-protein family after receptor, system/context, and residual replica variation. 500 receptor-block bootstraps with 95% CI |

### [S03_resampling_multiplicity_design/](./S03_resampling_multiplicity_design/) — SI Table S15

| File | Description |
|------|-------------|
| `Supplementary_Data_S3_resampling_and_multiplicity_design.csv` | Every resampling and permutation design used in the manuscript: iterations, inferential unit, p-value convention, and multiplicity correction for each analysis family |

### [S04_known_site_concordance/](./S04_known_site_concordance/) — SI Table S3

| File | Description |
|------|-------------|
| `Supplementary_Data_S4_known_site_concordance.csv` | Known-site concordance for both all-cohort and Class-A-only atlas scopes: Jaccard overlap, recall, precision, matched-null 95th percentile, empirical p, and Benjamini-Hochberg q across prespecified site families |

### [S05_interface_coherence/](./S05_interface_coherence/) — SI Table S13

| File | Description |
|------|-------------|
| `Supplementary_Data_S5_interface_receptor_level_coherence.csv` | Prespecified interface-descriptor coherence checks: Spearman rho, naive p, and BH q for five descriptor pairs (e.g. α5 tilt vs hook RMSD, contacts vs insertion depth) |

### [S06_alpha5_hook_audit/](./S06_alpha5_hook_audit/) — SI Tables S10–S11

| File | Description |
|------|-------------|
| `Supplementary_Data_S6_alpha5_hook_audit.csv` | Per-system Gα α5 hook sequence audit: deposited, simulated, and canonical family hooks compared for all 207 systems. Records match status, construct evidence, and conflicts |

### [S07_gateway/](./S07_gateway/) — SI Tables S4–S8

| File | Description |
|------|-------------|
| `Supplementary_Data_S7a_gateway_replica_long.csv` | Per-replica, per-system, per-portal open-fraction values (long format, 621 replicas × 7 portals) |
| `Supplementary_Data_S7b_gateway_replica_icc.csv` | Portal-specific ICC(1,1) with system-bootstrap 95% CI |
| `Supplementary_Data_S7c_gateway_threshold_sensitivity.csv` | Depth-threshold sensitivity (0.3, 0.5, 0.7, 1.0 Å) for portal means, medians, and ranks |
| `Supplementary_Data_S7d_gateway_stationarity.csv` | First-vs-last 100-ns stationarity tests: mean difference, receptor-block bootstrap CI, Wilcoxon p, BH q |
| `Supplementary_Data_S7e_gateway_static_associations.csv` | Prepared-geometry (symmetric mean nearest Cα) versus dynamic open fraction: Spearman ρ, block-permutation p, BH q |

### [S08_R6_robustness/](./S08_R6_robustness/) — SI Table S10

| File | Description |
|------|-------------|
| `Supplementary_Data_S8_R6_robustness.csv` | Legacy C7 and final R6 association estimates across all prespecified variants: risk ratio, bootstrap 95% CI, Fisher exact p, one-system-per-receptor resampling, leave-one-out |
| `Supplementary_Data_S8b_pooled_C7_comparator.csv` | Pooled-diagnostic C7 comparison table for the superseded pooled atlas |

### [S09_stable_cluster/](./S09_stable_cluster/) — SI Table S2

| File | Description |
|------|-------------|
| `Supplementary_Data_S9a_stable_cluster_recurrence.csv` | Stable-core recurrence for all-cohort and Class-A-only atlases: systems, receptors, H40 groups, core positions, split-half recovery, characteristic region |
| `Supplementary_Data_S9b_split_half_stability.csv` | Split-half stability metrics: mean/median best-core Jaccard, 90% CI, recovery probability per cluster |

### [S10_conservation_matched_null/](./S10_conservation_matched_null/) — SI Table S3

| File | Description |
|------|-------------|
| `Supplementary_Data_S10_conservation_matched_null.csv` | Class-A conservation overlay: cluster-level mean H40 conservation, matched-null mean and CI, empirical two-sided p, and BH q for 11 stable Class-A clusters |

### [S11_peptide_sensitivity/](./S11_peptide_sensitivity/) — SI Table S14 (M1–M4)

| File | Description |
|------|-------------|
| `comparison_summary.md` | Peptide-stripped sensitivity narrative: four preregistered material-change rules (M1–M4) and their outcomes |
| `comparison_summary.json` | Machine-readable summary of the peptide-sensitivity comparison |

### [S14_statistical_inventory/](./S14_statistical_inventory/) — SI Table S14

Artifact-identifier files for each analysis family in the statistical
test inventory:

| File | Analysis |
|------|----------|
| `gateway_replica_icc.csv` | Gateway ICC(1,1) |
| `gateway_stationarity_tests.csv` | First/last 100-ns stationarity |
| `gateway_threshold_sensitivity.csv` | Depth-threshold sensitivity |
| `gateway_variance_components.csv` | Mixed-model variance partition |
| `gateway_context_noise_controlled_contrasts.csv` | Context contrasts vs sampling null |
| `gprotein_alpha5_hook_sequence_audit.csv` | Gα α5 hook audit |
| `replica_cluster_parameter_sensitivity.csv` | Clustering-parameter sensitivity |
| `contrast_continuous_confound_associations.csv` | Continuous confounds |
| `contrast_binary_confound_associations.csv` | Binary confounds |
| `contrast_gateway_pocket_cross_metric.csv` | Gateway–pocket cross-metric |
| `interface_descriptor_outcome_correlations.csv` | Descriptor–outcome associations |
| `interface_descriptor_coherence_checks.csv` | Descriptor coherence checks |
| `replica_cluster_split_half_stability.csv` | Split-half stability |
| `replica_known_site_null_calibrated_concordance.csv` | Known-site concordance |
| `replica_cluster_family_maxstat_permutation.csv` | Cluster × family screen |
| `legacy_C7_replica_sensitivity.csv` | C7/R6 sensitivity |
| `replica_pocket_context_receptor_inference.csv` | Pocket context vs replica noise |
| `replica_pocket_fingerprint_variance_partition.csv` | Fingerprint variance partition |

---

## Metadata

- **SHA256SUMS.json** — Byte-level checksums for every file in this bundle.

These files are the corrected machine-readable inputs or direct producer
outputs used for the manuscript, SI, and figures. The processed
per-replica gateway tables reproduce the reported inferential summaries
and figures. Reduced receptor/G-protein coordinate deposits omit lipids,
solvent, and ions and therefore do not by themselves permit
coordinate-level recalculation of lipid-facing gateway measures;
full-system trajectory reviewer access is described in the manuscript
data-availability statement.
