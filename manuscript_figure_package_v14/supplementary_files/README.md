CoupledMD Supplementary Data

Generated: 2026-07-23

This package accompanies the CoupledMD Data Descriptor describing a frozen
207-system cohort of active GPCR-G-protein complex molecular dynamics simulations.

S1  System inventory (207 rows).
    Per-system master record: receptor name, UniProt accession, GPCR class,
    G-protein family, explicit retained-small-molecule flag, legacy
    ligand-or-peptide flag, computational-repeat count and production protocol.
    The retained_small_molecule field reproduces the manuscript classification
    (True = 85; False = 122).

S2  Pocket annotations (2,151 rows).
    Detected-pocket records for 205 systems plus 2 explicit validated
    zero-pocket records, providing complete pocket annotation for all 207
    systems. Schema version 1.1 explicitly identifies the 58 peptide controls
    (peptide_control), peptide-proximal persistent clusters
    (is_peptide_proximal), per-control recovery (peptide_site_recovered), and
    the best per-control peptide-site frequency (best_peptide_site_freq).
    The legacy is_orthosteric field is retained for compatibility.

S3  File manifest (1,242 rows).
    SHA-256 checksum, atom count and frame count for every deposited PDB and
    XTC file across the three-replica release (207 systems x 2 files x 3
    replicas = 1,242 file records).

S4  Quality-control results (621 rows).
    Outcome of nine full-frame QC checks for each of the 621 released
    trajectory records: atom-count match, finite coordinates, monotonic
    time stamps, valid box dimensions, backbone geometry atom order,
    no chain breaks, no complex separation, no coordinate scatter and
    full frames checked — plus an overall pass/fail verdict.

S5  Structural validation metrics (621 rows; 618 validated records plus three
    null-metric Gs_8HTI rows).
    Receptor TM-core C-alpha RMSD P95, G-alpha interface-region C-alpha
    RMSD P95 and initial interface-contact retention P05 are populated for 618
    validated replica records. The three Gs_8HTI rows are retained for complete
    release-record coverage but contain null metrics because the absent UniProt
    mapping precludes definition of the TM-core residue set.

S6  Gateway intermediates (5,796 rows).
    Per-system per-interface summary values (occupancy, penetration,
    penetration_p90, open_fraction) for 7 TM helix pairs across all 207
    systems. Each row reports all three replica values, their mean and their
    observed minimum and maximum. With only three replicas, no confidence
    interval is assigned. These processed outputs were computed from the
    full-system source trajectories and are provided so that reported
    manuscript summaries can be reproduced without access to those trajectories.

Data licence: CC BY 4.0.  Code licence: MIT.
