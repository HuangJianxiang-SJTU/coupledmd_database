CoupledMD Supplementary Data

Release scope
-------------
The included cohort comprises 207 GPCR-G-protein systems. The source simulation campaign contains three 500-ns full-system production replicas per system (621 replicas; 310.5 microseconds total). Those full-system source trajectories are not distributed as molecular files.

The reduced molecular release contains one matched PDB/XTC pair per system derived from source replica 1. The reduced files contain proteins plus selected ligands retained during system preparation where applicable. Membrane, solvent, ions, and ligands intentionally removed during system preparation are excluded. Each released XTC contains 2,500 frames at approximately 200-ps intervals and represents 499.8 ns from the first to last retained frame. These files are intended for visualization and conformational analysis of retained components, not membrane, solvent, ion, removed-ligand, or fast-timescale analyses.

Tables
------
S1 lists the 207 included systems and source-campaign metadata.
S2 lists 13 unresolved systems and two excluded records (Gq_7E9W, a non-GPCR duplicate/mislabel, and Gq_7DWC, a duplicate source identity of Gq_8DWC) outside the release cohort.
S3 defines every field serialized in S1-S10.
S4 is the 414-file reduced molecular-release manifest, including archive membership, byte sizes, and SHA-256 checksums.
S5 contains all detected mapped-pocket records plus two explicit valid zero-pocket system records.
S6 contains system/interface/metric gateway summaries with all three source-replica values.
S7 is the 621-row source-simulation replica provenance ledger. It distinguishes replica 1, which generated the reduced molecular files, from replicas 2 and 3, which support analyses but are not distributed, and records path-free protocol-group, engine-version, output-log coverage and seed-evidence scope.
S8 contains one full-frame reduced molecular-file QC record per included system.
S9 reports release-scoped portal/API record coverage; it is not an archive-coverage table.
S10 inventories the current code and build environment. Immutable repository release identifiers remain a final release gate until the public code revision is frozen.

Supplementary_Data_package_audit.json records table counts, package-file hashes, and validation results. DOI publication and immutable code identifiers must be synchronized at final release; no unpublished DOI is claimed here.
