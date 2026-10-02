# Analysis and release workflows

These files support the code-availability statement in the CoupledMD manuscript.
The server-derived scripts are preserved with their original filenames and
scientific parameters. Some retain institutional source-data paths; configure
those paths for another installation. Workflows that recalculate persistent
pockets require the full-system source trajectories described under Data
Availability.

## Manuscript claim → implementation

- Trajectory reduction and nine-point released-record QC:
  `reduction_qc/prepare_zenodo_reduced_release.py`. Its `stage_release*`
  functions generate the reduced records, and `qc_one`/`full_qc` implement the
  full-frame release checks. `build_replica_qc.py` assembles the replica ledger.
- Structural-validation metrics:
  `structural_validation/build_core_interface_validation.py` calculates
  receptor TM-core RMSD, Gα interface-region RMSD, and initial-interface-contact
  retention. The two `build_final624_*` files preserve its selected-replica and
  harmonization dependencies.
- Persistent-pocket annotation:
  `persistent_pockets/paper1_pocket_pilot.py` performs per-frame fpocket
  detection, grid aggregation, persistence thresholding, clustering, and lining
  residue assignment. `paper1_pocket_gpcrdb_map.py` and
  `paper1_pocket_consensus.py` perform receptor mapping and persistent-cluster
  aggregation. `build_system_config.py` and `paper1_gateway_metric.py` are direct
  dependencies.
- Worked ligand-contact and binding-pocket analyses:
  `example_analyses/example_ligand_contact_binding_pocket.py` runs directly on
  a released PDB/XTC pair and writes the contact matrix, per-residue contact
  persistence, and fitted pocket Cα RMSD time series.
- Publication figures: the independently runnable `figure1.py`–`figure5.py`,
  `figure_s1.py`, and `figure_s3.py` scripts in the parent directory.

## Environment

Create the Conda environment with:

```bash
conda env create -f environment.yml
conda activate coupledmd-reproducibility
```

`requirements.txt` provides a pip-compatible Python dependency list. fpocket
4.2.3 is listed in `environment.yml` because it is an external executable rather
than a Python package.

## Worked-example usage

```bash
python example_analyses/example_ligand_contact_binding_pocket.py \
  path/to/system.pdb path/to/system.xtc \
  --ligand-selection "resname S1P and not name H*" \
  --pocket-selection "protein and name CA and resid 62 88 91 684 865" \
  --output-dir example_output
```

Use the complete pocket-residue selection reported for the worked system when
reproducing the manuscript analysis.
