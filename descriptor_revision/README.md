# CoupledMD Scientific Data descriptor workflows

This directory is the code package for “CoupledMD: molecular dynamics trajectories of active GPCR–G-protein complexes”, distinct from the GPCR-atlas analyses elsewhere in the repository. The dataset has 207 systems, 621 repeats and Supplementary Data S1–S8. Dataset DOI: https://doi.org/10.5281/zenodo.21395291; version used: https://doi.org/10.5281/zenodo.23106855.

`python verify_tables.py` checks the deposited-table boundary without opening trajectories. The original analysis code is in `scripts/workflows/`: reduction_qc, structural_validation, persistent_pockets, example_analyses, stationarity and grafted_segments. Environment files and the claim-to-implementation map are in that directory. `source_data/` contains figure inputs; `revision/` contains the S7/S8 methods and frozen outputs; `figures/` contains the current revised figures.

These are the authored scientific workflows, not a zero-configuration reanalysis service. Some preserve institutional absolute paths and repaired-trajectory routing; configure paths before running them. Structural selection/routing CSVs and their correction-audit dependencies are packaged with structural_validation. Its ordinal selections refer to the original source topologies: do not substitute a reduced topology without remapping selections. The historical “624/208” filenames predate the frozen 207-system release and are provenance names, not current cohort counts. Use supplementary_data/Supplementary_Data_S1_system_inventory.csv as the release allow-list.

Persistent-pocket and lipid-gateway recalculation require source trajectories requested through the portal. The released reduced records alone do not reproduce the original annotation frame grid. Grafted-region code is preserved from the validated analysis snapshot; comments in its historical header overstate removal of interdomain motion. The authoritative limitations are in supplementary_data/Supplementary_Data_S8_README.md. Its sequence comparison infers candidate segments and does not prove template provenance, experimental accuracy or absence of drift. The disabled alternative graft script is excluded.

No trajectories were recomputed while assembling this package. Code is MIT licensed; data tables use CC BY 4.0. See LICENSE-CODE and the dataset license.


## Final revision figures and data (2026-10-04)

[Figure inventory](figures/README.md) lists all six main figures and four supplementary figures. The PNG files are the exact embedded image bytes of the finalized manuscript/SI; PDFs are provided for convenient viewing. Their provenance and hashes are in figures/manifest.json. Captions are in figures/captions.md.

The S1–S8 CSVs are unchanged from the validated deposited release (207 / 2,151 / 1,242 / 621 / 621 / 5,796 / 621 / 938 rows). S7 contains TM-core recommendations only. The packaged `revision/equilibration_out/equilibration_per_replica.csv` also exposes the interface assessment: 289 of 618 assessed repeats (46.8%) fail the interface criterion (1.0 Å plateau tolerance and 0.1 Å per 100 ns slope limit). TM-core tolerance is 0.5 Å. There is no cohort-wide interface discard recommendation; this is an operational stationarity test, not proof of convergence. Three Gs_8HTI repeats lack a TM-core selection but remain among the 621 released records.

No trajectory calculation was rerun for this update. The full grafting population remains primary (938 segments); sensitivity subsets remain separate. All archived trajectory layers are water-free. The two GROMACS systems have identical encoded starting states across their three computational repeats; independent sampling is not assumed.

## Replot from packaged tables

From `descriptor_revision/`, run `python scripts/figure2.py`, `python scripts/figure6.py`, or `python scripts/figure_s4.py`. These entry points use local packaged tables and write to `figures/regenerated/`; they do not open trajectories. Existing scripts for Figures 1, 3–5 and S1/S3 also use packaged source tables. Figure S2 is a portal screenshot, not a computational plot. Historical source-data filenames use an earlier panel numbering scheme; the final figure inventory is authoritative.

Replotted images reproduce the packaged statistics; their fonts, cropping and sizing may differ from the author-formatted embedded figures. Regeneration does not overwrite the authoritative final figure files.
