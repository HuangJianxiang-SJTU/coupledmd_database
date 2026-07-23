# Reduced molecular dynamics trajectories of 207 active GPCR–G-protein complexes

Jianxiang Huang^1,2,3^, Xin Qiao^3^, Shaoyong Lu^1,2,3,\*^

^1^ Artificial Intelligence Clinical Research Center for Drug Discovery, Shanghai Key Laboratory of Flexible Medical Robotics, Institute of Medical Robotics, Tongren Hospital, Shanghai Jiao Tong University School of Medicine, Shanghai 200336, China

^2^ Key Laboratory of Protection, Development and Utilization of Medicinal Resources in Liupanshan Area, Ministry of Education; Peptide & Protein Drug Research Center, School of Pharmacy, Ningxia Medical University, Yinchuan 750004, China

^3^ Department of Pharmacology, School of Medicine, Shanghai Jiao Tong University, Shanghai 200025, China

\*Correspondence: Shaoyong Lu (lushaoyong@sjtu.edu.cn)

---

## Abstract

Reusable molecular-dynamics data for G protein-coupled receptor (GPCR)–G-protein complexes require consistent identifiers, molecular representations, provenance and replica-level quality control. We describe a frozen cohort of 207 active or active-like complexes spanning 174 receptor names, two GPCR classes and four G-protein families. The source-simulation design comprises three 500-ns computational repeats per system (621 replicas; 310.5 µs nominal sampling). The molecular release is organized as three replica-specific collections of matched reduced PDB/XTC records with 2,500 frames at 200-ps spacing. These files contain the selected protein-complex representation and exclude membrane lipids, solvent and mobile ions; they are not a full-system primary-trajectory archive. System, replica, file, checksum, protocol and quality-control records accompany the coordinates. Processed pocket and lipid-gateway summaries are supplied for reuse, with source data for every figure. Orthosteric recovery is used only as a technical positive control. The descriptor separates nominal from original file-observed sampling, documents a repaired reduced source for one shortened replica and defines analyses that require unavailable full-system coordinates.

## Background & Summary

G protein-coupled receptors (GPCRs) translate extracellular signals into intracellular responses through conformational coupling to heterotrimeric G proteins. Structures of active GPCR–G-protein complexes provide starting coordinates for molecular simulation, but static structures do not describe the time-dependent variation of receptor and transducer interfaces.^1,2^ Molecular dynamics (MD) can provide these coordinate ensembles, yet reuse across receptors remains difficult when identifiers, constructs, preparation protocols, atom selections, trajectory formats and validation records differ.

Community resources address complementary parts of this problem. GPCRdb and GproteinDb provide curated classifications, generic residue identifiers and structural annotations that permit comparisons across receptors and G-protein subtypes.^3,4^ GPCRmd provides community access to GPCR simulations and established a practical model for trajectory sharing.^5^ General-purpose MD collections such as mdCATH demonstrate the value of a frozen cohort, machine-readable provenance and executable reuse examples.^6^ A reusable collection of receptor–G-protein assemblies additionally requires explicit chain roles, construct provenance, G-protein-family harmonization and a clear distinction between full membrane simulations and molecularly reduced coordinate products.

The present Data Descriptor documents the database-owned parts of the CoupledMD resource: cohort construction, simulation protocols and provenance, identifier harmonization, three-repeat records, technical quality control, checksums, programmatic access and bounded validation of processed annotations. Biological interpretation of recurrent pockets, pocket rankings, partner-associated pocket changes, lipid-gateway patterns and G-protein selectivity is outside this paper and is reserved for a related resource-atlas manuscript. Learned models, predicted selectivity and experimental validation are outside both this descriptor and the present data release.

The frozen cohort contains 207 systems selected from a 222-record working inventory (Fig. 1). It includes 181 Class A and 26 Class B systems; the G-protein-family counts are 95 Gi/o, 65 Gs, 41 Gq/11 and 6 G12/13. The cohort contains 174 receptor names and 173 mapped UniProt accessions. `Gs_8HTI` is retained with an explicit null receptor accession. Thirteen working-inventory records remain unresolved and two are excluded: `Gq_7E9W` is a non-GPCR duplicate or mislabel of `Gq_8E9W`, and `Gq_7DWC` duplicates the source identity assigned to `Gq_8DWC` (MRGPRX1, Q96LB2). Neither unresolved nor excluded identifiers contribute to release denominators.

Each included system has a design of three 500-ns production repeats. These are computational repeats generated from independently assigned production seeds or velocities where evidence is available; they are not biological replicates. The nominal campaign size is therefore 621 trajectories and 310.5 µs. The original coordinate files contain 310.4168 µs because the original `Gq_8ZPT` replica 2 record spans 416.8 ns. A continuation-derived, locally audited reduced source for that replica contains 10,000 frames over 499.95 ns. Nominal, original file-observed and repaired reduced-source quantities are reported separately rather than being combined into one sampling claim.

The molecular product described here is a reduced release, not the underlying full-system simulation archive. Each release replica is organized as 207 matched PDB/XTC pairs, with 2,500 frames at approximately 200-ps spacing and a represented span of 499.8 ns. Membrane lipids, water and mobile ions are excluded. The retained molecular components are defined by each released PDB; ligands removed during preparation or reduction cannot be assumed to be present. Consequently, the reduced trajectories support protein-complex structural reuse but cannot reproduce lipid-gateway, solvent, ion or full-system energetic analyses. Processed intermediate tables are supplied for reported gateway summaries. A stable route to the full-system coordinates remains a submission-readiness requirement for any claim that depends on recalculation from those coordinates.

![Figure 1](publication_figures_v12/figure1_cohort_scope.png)

**Figure 1 | Frozen cohort, release boundary and source-protocol scope.** (**A**) Exact system counts in the GPCR-class × G-protein-family matrix. (**B**) The 222-record working inventory contains 207 included, 13 unresolved and two excluded records. (**C**) Systems and computational repeats assigned to the three source-simulation protocol groups. (**D**) Counts of release systems, computational repeats, receptor-name entities and mapped UniProt accessions. Family composition is descriptive and is not interpreted as physiological coupling prevalence.

## Methods

### Cohort construction and identifiers

Candidate active or active-like GPCR–G-protein assemblies were collected from the Protein Data Bank (PDB).^7^ Receptor chains were cross-referenced against UniProt, GPCRdb and GproteinDb records.^3,4,8^ Gα annotations were harmonized into Gi/o, Gs, Gq/11 or G12/13. Each system has a stable identifier of the form `<family>_<PDBID>`. Structural provenance is tracked separately from activation-state annotation: 197 included systems use experimental source coordinates and 10 have an `engineered_uncertain` provenance flag.

The release boundary was frozen programmatically against the 222-record working inventory. Inclusion required an assignable receptor–G-protein assembly, matched source trajectory and topology evidence, a three-repeat release design and membership in the final identifier table. Supplementary Data S1 lists all 207 included systems. Supplementary Data S2 lists the 13 unresolved and two excluded records with machine-readable reasons. The final list is used as an allow-list by table, figure and API validation.

### System preparation

Preparation records describe membrane-embedded complexes in a homogeneous POPC bilayer constructed through CHARMM-GUI workflows where recorded.^9,10^ Missing atoms or residues were modelled during the original preparation where required. Systems were solvated with TIP3P water, neutralized with Na^+^/Cl^-^ and prepared at a reported ionic strength of 0.15 M. Retained non-protein molecules were parameterized with CGenFF where applicable.^11^ Construct modifications, unresolved residues, retained ligands and engineered G-protein components vary among source structures and are retained as limitations and metadata rather than normalized away.

### Production protocols and provenance

Three protocol groups occur in the final cohort (Table 1). Group P1 contains 179 AMBER PMEMD systems generated from chamber-compatible CHARMM36-family inputs; P2 contains 26 AMBER systems prepared from CHARMM-GUI input sets; P3 contains two GROMACS systems. The CHARMM36-family,^17^ AMBER^18^ and GROMACS^19^ labels are retained at the resolution supported by source records. The source evidence supports AMBER for 205 systems and GROMACS for two systems. Engine-version evidence is heterogeneous because production records were generated over time. Version strings recovered from available records are serialized in Supplementary Data S7 rather than replaced by a single global version.

**Table 1. Source-simulation protocol groups.**

| Protocol | Systems | Computational repeats | Engine family | Declared force-field label |
|---|---:|---:|---|---|
| P1 | 179 | 537 | AMBER PMEMD | CHARMM36 via chamber |
| P2 | 26 | 78 | AMBER PMEMD | CHARMM36 via CHARMM-GUI |
| P3 | 2 | 6 | GROMACS | membrane-embedded CHARMM36 |
| **Total** | **207** | **621** | — | — |

Available input and output records specify a 2-fs integration step, particle-mesh Ewald electrostatics,^12^ constraints on bonds involving hydrogen and NPT production at 310 K and 1 bar. AMBER records use Langevin temperature control and engine-specific pressure coupling; GROMACS records use Nosé–Hoover and Parrinello–Rahman coupling.^13,14^ Engine-specific cutoffs and switching functions are retained in the protocol records. The production design is 500 ns per repeat with coordinates normally written every 50 ps. Replica seed evidence and recovered input/output paths are recorded in Supplementary Data S7. Direct protocol evidence is classified as verified for 202 final-cohort systems and partial for five (`Gi_8EFO`, `Gi_8HK2`, `Gq_7RAN`, `Gq_7XWO` and `Gq_7XXH`); partial records are not presented as fully recovered provenance.

### Nominal and file-observed sampling

Nominal sampling is calculated from the protocol design: 207 systems × 3 repeats × 500 ns = 310.5 µs. File-observed sampling was independently summed from the cohort-truth manifest and equals 310.4168 µs in the original coordinates. The 83.2-ns difference is confined to the original `Gq_8ZPT` replica 2 coordinate record. For preparation of the reduced release, a continuation-derived coordinate source for this replica was audited over all 10,000 frames: timestamps are monotonic, coordinates and boxes are finite, the represented span is 499.95 ns and the audit verdict is `GOOD`. This repair supports the reduced `Gq_8ZPT` replica 2 record. It does not establish that a repaired full-system source trajectory has been deposited.

### Reduced PDB/XTC generation

The three release collections are generated independently from source replicas 1, 2 and 3. For each system, the canonical reduced `structure.pdb` defines the retained atom set and atom order. Source trajectories are matched to that PDB directly or by residue-name, residue-identifier and atom-name mapping when the source contains additional atoms. Ten-thousand-frame, 50-ps sources are strided by four; a 2,501-frame source is trimmed to 2,500 frames. Every release XTC contains 2,500 frames separated by 200 ps and spans 499.8 ns from first to last retained frame.

Each replica-specific release is partitioned by G-protein family. A manifest records `system_id`, PDB identifier, family, source replica, relative PDB and XTC paths, frame interval, represented span, atom count, byte size and SHA-256 checksum. A JSON audit records full-frame coordinate, time, unit-cell, atom-order and component-continuity checks. No absolute local path is permitted in a release manifest.

### Metadata and access records

System, replica and file entities are kept separate and joined by `system_id`, `replica_id` or `release_replica` (Fig. 2). System records contain PDB, receptor, UniProt, class, G-protein family, Gα subtype, construct and protocol fields. Replica records contain nominal and observed duration, seed evidence, source engine and repair status. File records contain role, relative path, size, checksum, atom count, frame count, frame interval and molecular-component scope. Empty strings or JSON `null` represent unavailable values; zero is not used as a missing-value surrogate. The combination of stable identifiers, open formats, checksums and versioned metadata implements the applicable FAIR data principles.^21^

![Figure 2](publication_figures_v12/figure2_records_reuse_boundary.png)

**Figure 2 | Molecular-record boundary and reuse model.** (**A**) Full-system source simulations contain the protein complex, membrane, solvent and ions; reduction produces matched PDB/XTC records that omit the environmental components. (**B**) Supported and unsupported uses of reduced trajectories, processed gateway summaries and full-system source coordinates. `S` denotes a supplied processed summary. (**C**) System-to-replica-to-file joins and checksum-bearing file records. (**D**) Supplementary Data, API, code and repository records form complementary access layers. Public release remains conditional on unauthenticated resolution of the three deposited records.

### Pocket annotation and positive-control definition

Persistent pocket annotations were generated with fpocket.^15^ Approximately 200 frames per replica were aligned to a protein Cα reference, and fpocket was run independently on each selected frame. Alpha-sphere occupancy was accumulated on a 1.5-Å grid. Points occupied in at least 85% of frames were clustered with DBSCAN (ε = 2.85 Å; minimum samples = 3). Protein heavy atoms within 4.0 Å were assigned as lining residues. Pocket residues were mapped to GPCRdb generic positions where a receptor mapping was available.^3,4^

For technical positive-control validation, the analysis used 58 systems with a peptide-defined orthosteric reference. A detected cluster was classed as orthosteric if a retained grid point lay within 5.0 Å of a reference ligand heavy atom. `best_ortho_freq` is the maximum mean grid occupancy among orthosteric clusters. Systems without a ligand-defined reference are not failures and are excluded from this denominator. The test evaluates recovery of a known structural site; it does not validate druggability, biological recurrence or mechanism.

### Gateway intermediate records

Gateway processing used full-system trajectories because lipid coordinates are required. Seven adjacent transmembrane helix pairs were defined from GPCRdb helix assignments. Lipid heavy atoms within 5.0 Å of both helices and within the central transmembrane z band were counted as wedged. Per-system records contain four stored metrics for each helix pair and three replica summaries. Supplementary Data S6 contains 207 × 7 × 4 = 5,796 rows.

The reduced PDB/XTC release contains no lipids and cannot regenerate these gateway records. Accordingly, the database paper provides the processed per-system intermediate table, schema and code path but does not present the family-organized biological gateway distribution. Biological comparison and interpretation of these records are reserved for the related resource-atlas paper.

### Quality-control procedures

Reduced release files were checked over all frames for:

1. exact final-cohort membership and absence of unresolved or excluded identifiers;
2. PDB/XTC atom-count agreement;
3. 2,500 finite coordinate frames with monotonic time;
4. finite, valid periodic boxes;
5. backbone atom-order geometry against the reference;
6. absence of catastrophic chain breaks, coordinate scatter or complex separation; and
7. manifest-relative paths, sizes and checksums.

An independent audit of the portal-oriented reduced records sampled 40 frames per system and checked topology–trajectory agreement, blank element fields, unit-cell consistency, intra-chain breaks, inter-chain splitting and coordinate scatter. This audit was regenerated against the exact final-207 allow-list.

For source-trajectory technical validation, receptor TM-core Cα displacement, Gα interface-region displacement and initial contact retention were evaluated on harmonized per-replica observations. These metrics are available for 618 of 621 replica records from 206 systems. The three unavailable records belong to `Gs_8HTI`, for which the absence of a canonical receptor accession prevents the same mapping procedure.

### Reproducible builds and licensing

The v12 build records the Python version, selected dependency versions, repository base commit, hashes of build scripts and licence files in Supplementary Data S10 and `reports/v12_environment_freeze.json`. Figure scripts consume only the distributed, path-neutral source-data tables. Data tables and reduced molecular records are licensed CC BY 4.0; source code is MIT-licensed. Source PDB coordinates remain subject to the wwPDB terms associated with their records.

## Data Records

### Dataset scope

The final cohort is summarized in Table 2. Counts are regenerated from `data/release_cohort_v9_final207.csv`; they are not corrected by manual subtraction.

**Table 2. Frozen cohort by GPCR class and G-protein family.**

| GPCR class | G-protein family | Systems | Nominal sampling (µs) |
|---|---|---:|---:|
| Class A | Gi/o | 92 | 138.0 |
| Class A | Gs | 46 | 69.0 |
| Class A | Gq/11 | 39 | 58.5 |
| Class A | G12/13 | 4 | 6.0 |
| Class B | Gi/o | 3 | 4.5 |
| Class B | Gs | 19 | 28.5 |
| Class B | Gq/11 | 2 | 3.0 |
| Class B | G12/13 | 2 | 3.0 |
| **Total** | **All** | **207** | **310.5 nominal** |

The 310.5-µs total is the protocol design. The original file-observed total is 310.4168 µs and is reported in Supplementary Data S7. Neither number is described as the size of a deposited full-system archive.

### Reduced molecular records

The intended three-record release contains one record per source replica. Within each record, four family archives contain one directory per included system:

```text
<family archive>/
  <system_id>/
    structure.pdb
    traj.xtc
```

Each record also contains:

```text
CoupledMD_reduced_trajectory_manifest.csv
CoupledMD_reduced_release_audit.json
README.txt
LICENSE.txt
```

`structure.pdb` is the topology and atom-order reference for `traj.xtc`. The XTC has 2,500 frames at approximately 200-ps spacing. The files exclude membrane, solvent and mobile ions. They should not be called full-system or primary trajectory records.

### Supplementary Data

The v12 Supplementary Data package contains:

| File | Granularity | Content |
|---|---:|---|
| S1 | 207 systems | Included-system inventory |
| S2 | 15 records | Thirteen unresolved and two excluded records |
| S3 | one row per field | Metadata dictionary |
| S4 | 1,242 expected files | Three-replica reduced PDB/XTC manifest, with pending entries explicit |
| S5 | 2,151 rows | 2,149 detected pockets plus two explicit zero-pocket records |
| S6 | 5,796 rows | Processed gateway intermediates |
| S7 | 621 repeats | Source-simulation provenance and sampling definitions |
| S8 | 621 expected records | Three-replica reduced-release QC, with pending entries explicit |
| S9 | record types | API and consensus-cache consistency audit |
| S10 | files and packages | Code and environment inventory |

Supplementary Data S4 and S8 represent the expected 207 × 3 release grid and keep unavailable or pending release records explicit. They must be regenerated after the replica-release worker completes. The package audit records file hashes and current release truth.

### Figure source data

Every quantitative or schematic panel is backed by a CSV in `v12_figure_source_data/`, accompanied by a SHA-256 manifest. The old GPCR-centred pocket-reuse atlas is absent. Figure 4 source data contain only technical validation. Gateway intermediates needed to reproduce manuscript summaries are supplied in Supplementary Data S6, while the unavailable full-system molecular source is explicitly identified.

### Portal, API and code records

The public code and web repository is `https://github.com/HuangJianxiang-SJTU/coupledmd`. The manuscript, figure and source-data repository is `https://github.com/HuangJianxiang-SJTU/coupledmd_database`. The portal uses a browser-based molecular viewer for reduced coordinate inspection.^20^ The REST API schema and static response files are versioned with the code. The current gateway consensus cache has been regenerated from the authoritative final-207 summary: the `ALL` denominator is 207, with Class A 181, Class B 26, Gi/o 95, Gs 65, Gq/11 41 and G12/13 6.

## Technical Validation

### Numerical and identifier consistency

A machine-readable numerical audit checks all occurrences and source denominators involving 207/208 systems, 621/624 replicas, 310.5/310.4168/312 µs, 181/182 Class A systems, 41/42 Gq systems, protocol counts, pocket totals, zero-pocket systems and validation denominators. The active v12 artifacts contain 207 final identifiers and no active `final208` or `final624` label. Historical filenames are treated only as provenance inputs and are not copied into the submission package.

The final boundary is 207 included + 13 unresolved + 2 excluded = 222 records. The 208→207 change is traced to exclusion of `Gq_7DWC`; it is not a hidden deletion or arithmetic adjustment. The API consensus cache formerly retained stale 208-system denominators and has been replaced with an asserted final-207 cache.

### Pocket-record completeness and positive control

Persistent pockets were detected in 205 systems, giving 2,149 detected-pocket rows. `Gi_7YK6` and `Gi_8YIC` are completed analyses with zero detected persistent pockets. Thus the pocket layer contains 205 detected + 2 explicit zero-pocket system records = 207, not 206 + 2.

Forty-nine of 58 eligible peptide-reference systems recovered an orthosteric cluster at the stated frequency threshold (Fig. 3). The family-specific technical counts are 16/18 Gi/o, 15/18 Gs and 18/22 Gq/11; no G12/13 system is eligible. These denominators are reported to make the positive control auditable, not to compare G-protein families biologically.

![Figure 3](publication_figures_v12/figure3_annotation_validation.png)

**Figure 3 | Technical validation and completeness of the pocket records.** (**A**) Best orthosteric-cluster frequency for 58 peptide-reference systems. Crosses denote nine systems without a recovered orthosteric cluster at the 0.85 threshold. (**B**) Two hundred and five systems have one or more detected persistent pockets and two have valid zero-pocket results, yielding 207 completed system records. The 2,149 count refers to detected-pocket rows, not systems. This figure is a technical positive control and does not support a druggability or mechanistic claim.

### Replica-level structural metrics

The structural metric table has 621 expected replica rows. Metrics are available for 618 rows across 206 systems and unavailable for the three `Gs_8HTI` repeats. Each available row contains 1,001 harmonized observations and finite values for receptor TM-core RMSD P95, Gα interface-region RMSD P95 and initial-contact retention P05. Figure 4 displays distributions by computational repeat rather than by G-protein family, avoiding biological comparison in this Data Descriptor.

### Reduced-record validation

The exact final-207 independent audit classifies all 207 portal-oriented reduced records as `OK`. Each sampled record has matching PDB/XTC atom counts, populated element fields, a unit cell, no box mismatch, no sampled intra-chain break, no coordinate-scatter frame and no inter-chain split. Frame counts are 2,500 for 203 systems, 2,501 for two systems, 681 for `Gq_7RAN` and 200 for `Gs_3SN6`. These irregular legacy portal records are not confused with the replica-specific Zenodo release, whose build target is exactly 2,500 frames per system.

Full-frame release QC has been completed locally for replicas 1 and 3 at the time of the v12 build. Replica 2 remains under parallel processing and is represented as pending in the current v12 tables. Release completion and public resolution are therefore not claimed.

![Figure 4](publication_figures_v12/figure4_replica_technical_validation.png)

**Figure 4 | Replica-level structural validation and reduced-record quality control.** (**A–C**) Distributions across computational repeats for receptor TM-core Cα RMSD P95, Gα interface-region Cα RMSD P95 and initial-interface-contact retention P05. The available denominator is 618/621 records. (**D**) Completeness of structural metrics and exact final-207 reduced-record checks. The frame-count distribution refers to the independently audited portal-oriented records. Panels are technical and are not organized to support family-level biological interpretation.

### Release integrity and reproducible reuse

The release gate distinguishes four separate facts: a complete local manifest, passed local full-frame QC, the authenticated state of the private repository draft and unauthenticated public DOI resolution. A reserved DOI alone satisfies none of the last three. Figure 5 is regenerated from the machine-readable remote-status audit and therefore acts as a release gate rather than a publication claim.

Available replica manifests show closely matching technical distributions of XTC size and retained atom count because each replica uses the same per-system PDB atom definition. The worked example uses `G12_7SF7` replica 1. It verifies both file digests, loads the PDB as the XTC topology with MDAnalysis^16^ and checks the expected atom and frame counts.

![Figure 5](publication_figures_v12/figure5_release_integrity_reuse.png)

**Figure 5 | Three-replica release integrity and checksum-first reuse.** (**A**) Current machine-verified status for local manifests, full-frame QC, authenticated private-draft state and public DOI resolution. A reserved identifier is not treated as publication. (**B,C**) Per-system reduced-XTC sizes and retained atom counts for locally complete replica manifests; a blank replica position denotes a pending manifest, not zero-sized data. (**D**) Worked checksum and MDAnalysis loading sequence for one matched PDB/XTC pair. Panel A must be regenerated immediately before submission.

## Usage Notes

### Appropriate uses of the reduced records

The PDB/XTC pairs are suitable for protein-complex visualization, topology-aware structural measurements on retained atoms, method prototyping and cross-receptor joins that do not require removed molecular components. Users should load `structure.pdb` as the topology for `traj.xtc`, verify both SHA-256 values before analysis and retain `system_id` and `release_replica` in all derived tables.

The three production repeats are correlated computational repeats of one prepared system. They should not be treated as biological replicates or as three independent receptor structures. When constructing statistical or machine-learning splits, users should group by receptor identity—or a stricter biological grouping—rather than splitting replicas across training and test sets.

### Checksum-first worked example

After extracting one family archive, a user can select the `G12_7SF7` pair from replica 1. The complete expected digests and relative paths are in `Figure_5D_worked_example.csv` and Supplementary Data S4.

```bash
sha256sum G12_7SF7/structure.pdb
sha256sum G12_7SF7/traj.xtc
```

The returned digests must match the manifest exactly. The pair can then be loaded:

```python
from pathlib import Path
import hashlib
import pandas as pd
import MDAnalysis as mda

root = Path("CoupledMD_reduced_trajectories_G12-13/G12_7SF7")
manifest = pd.read_csv("CoupledMD_reduced_trajectory_manifest.csv")
row = manifest.loc[manifest["system_id"].eq("G12_7SF7")].iloc[0]

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

pdb_path = root / "structure.pdb"
xtc_path = root / "traj.xtc"
assert digest(pdb_path) == row["pdb_sha256"]
assert digest(xtc_path) == row["xtc_sha256"]

u = mda.Universe(pdb_path, xtc_path)
assert u.atoms.n_atoms == int(row["n_atoms"])
assert len(u.trajectory) == int(row["n_frames"]) == 2500
```

### Unsupported uses and limitations

- The reduced XTC files do not contain membrane, solvent or mobile ions. They cannot reproduce lipid-gateway, hydration, ion-binding or membrane-property analyses.
- Ligand retention varies with the canonical reduced PDB. Users must inspect the released atom inventory and must not assume that the source ligand is present.
- A 200-ps frame interval is inappropriate for fast kinetic or water-mediated measurements.
- Full-system source trajectories, engine topologies and preparation inputs are not part of the present reduced release. Until a stable access route exists, independent recalculation of full-system-only analyses remains unavailable.
- The cohort reflects available active or active-like structures, not physiological coupling prevalence. Related structures, constructs and receptor identities are not statistically independent.
- Ten systems carry engineered or uncertain structural provenance. `Gs_8HTI` lacks a canonical receptor accession.
- The simulations use a homogeneous POPC environment and three protocol groups. Protocol, construct, source PDB, ligand context and structural provenance should be treated as possible batch variables.
- Nominal sampling and file-observed coordinate sampling are different quantities. Cite the definition used.
- Pocket and gateway records are computational annotations. Their presence is not experimental evidence of druggability, selectivity or mechanism.

## Data Availability

The reduced release is assigned three reserved Zenodo identifiers, one per source replica: `10.5281/zenodo.21395292` (replica 1), `10.5281/zenodo.21447748` (replica 2) and `10.5281/zenodo.21448037` (replica 3). On the v12 remote audit date, authenticated read-only checks classified all three as restricted, unsubmitted private drafts, and unauthenticated DOI requests did not resolve publicly. Replica 2 was still being generated. These identifiers are therefore reported as reserved identifiers only; they are not asserted to be published data records, and archive completeness is not claimed.

An anonymous reviewer-access route has not been machine-verified in this build. Before submission, the three records must either be made public or supplied through tested anonymous reviewer links, their files must match the local manifests and audits, and Figure 5, Supplementary Data S4/S8, the checksum reports and this section must be regenerated.

The present release contains reduced protein-complex PDB/XTC files and excludes the underlying full-system trajectories. Processed source data for all manuscript figures are included in the paper-content repository. Supplementary Data S6 and the accompanying code provide the intermediate gateway records needed to reproduce summaries, but recalculation from molecular coordinates requires full-system trajectories that are not distributed here. A stable full-system access route is therefore a blocker for any manuscript claim that requires independent molecular recalculation of lipid-dependent results.

## Code Availability

The portal, API, processing and validation code is available at `https://github.com/HuangJianxiang-SJTU/coupledmd` under the MIT licence. Manuscript source, figure scripts, figure source data, Supplementary Data and validation reports are maintained at `https://github.com/HuangJianxiang-SJTU/coupledmd_database`. The v12 environment table records the base repository commit and file hashes used for this build. A v12 commit and immutable release tag must be pushed to the paper-content repository before submission; the current text does not invent a tag or commit for unpublished changes.

## Related manuscripts

A related resource-atlas manuscript analyses recurrent pockets, pocket clusters, lipid-gateway patterns, partner-associated changes and G-protein-family selectivity. The present Data Descriptor owns cohort construction, simulation and identifier provenance, reduced records, technical quality control, access and reuse instructions, and limited positive-control validation. Source-data overlap is disclosed in Supplementary Table S2, and biological conclusions from the resource-atlas study are not repeated here. A planned modelling study involving learned architectures, predictions or experimental validation does not contribute claims or results to this descriptor.

## Author Contributions

J.H. and S.L. conceived the resource. J.H. curated the inventory, performed the simulations, developed the processing and annotation workflows, implemented the portal and API, performed technical validation, and generated the figures and tables. X.Q. contributed to data curation and validation. J.H. and S.L. wrote the manuscript. All authors reviewed and approved the manuscript.

## Competing Interests

The authors declare no competing interests.

## Acknowledgements

The authors acknowledge computational resources provided by Shanghai Jiao Tong University.

## Funding

This work was supported by the Noncommunicable Chronic Diseases–National Science and Technology Major Project (2024ZD0531200) and the Innovative Research Team of High-Level Local Universities in Shanghai.

## References

1. Weis, W. I. & Kobilka, B. K. The molecular basis of G protein-coupled receptor activation. *Annu. Rev. Biochem.* **87**, 897–919 (2018). https://doi.org/10.1146/annurev-biochem-060614-033910
2. Hauser, A. S. *et al.* Pharmacogenomics of GPCR drug targets. *Cell* **172**, 41–54.e19 (2018). https://doi.org/10.1016/j.cell.2017.11.033
3. Taracena Herrera, L. P. *et al.* GPCRdb in 2025: adding odorant receptors, data mapper, structure similarity search and models of physiological ligand complexes. *Nucleic Acids Res.* **53**, D425–D435 (2025). https://doi.org/10.1093/nar/gkae1065
4. Pándy-Szekeres, G. *et al.* GproteinDb in 2024: new G protein–GPCR couplings, AlphaFold2-multimer models and interface interactions. *Nucleic Acids Res.* **52**, D466–D475 (2024). https://doi.org/10.1093/nar/gkad1089
5. Rodríguez-Espigares, I. *et al.* GPCRmd uncovers the dynamics of the 3D-GPCRome. *Nat. Methods* **17**, 777–787 (2020). https://doi.org/10.1038/s41592-020-0884-y
6. Mirarchi, A., Giorgino, T. & De Fabritiis, G. mdCATH: a large-scale MD dataset for data-driven computational biophysics. *Sci. Data* **11**, 1299 (2024). https://doi.org/10.1038/s41597-024-04140-z
7. Berman, H. M. *et al.* The Protein Data Bank. *Nucleic Acids Res.* **28**, 235–242 (2000). https://doi.org/10.1093/nar/28.1.235
8. UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res.* **53**, D609–D617 (2025). https://doi.org/10.1093/nar/gkae1010
9. Jo, S., Kim, T., Iyer, V. G. & Im, W. CHARMM-GUI: a web-based graphical user interface for CHARMM. *J. Comput. Chem.* **29**, 1859–1865 (2008). https://doi.org/10.1002/jcc.20945
10. Lee, J. *et al.* CHARMM-GUI Membrane Builder for complex biological membrane simulations with glycolipids and lipoglycans. *J. Chem. Theory Comput.* **15**, 775–786 (2019). https://doi.org/10.1021/acs.jctc.8b01066
11. Vanommeslaeghe, K. *et al.* CHARMM General Force Field: a force field for drug-like molecules compatible with the CHARMM all-atom additive biological force fields. *J. Comput. Chem.* **31**, 671–690 (2010). https://doi.org/10.1002/jcc.21367
12. Darden, T., York, D. & Pedersen, L. Particle mesh Ewald: an N·log(N) method for Ewald sums in large systems. *J. Chem. Phys.* **98**, 10089–10092 (1993). https://doi.org/10.1063/1.464397
13. Nosé, S. A unified formulation of the constant temperature molecular dynamics methods. *J. Chem. Phys.* **81**, 511–519 (1984). https://doi.org/10.1063/1.447334
14. Parrinello, M. & Rahman, A. Polymorphic transitions in single crystals: a new molecular dynamics method. *J. Appl. Phys.* **52**, 7182–7190 (1981). https://doi.org/10.1063/1.328693
15. Le Guilloux, V., Schmidtke, P. & Tuffery, P. Fpocket: an open source platform for ligand pocket detection. *BMC Bioinformatics* **10**, 168 (2009). https://doi.org/10.1186/1471-2105-10-168
16. Michaud-Agrawal, N. *et al.* MDAnalysis: a toolkit for the analysis of molecular dynamics simulations. *J. Comput. Chem.* **32**, 2319–2327 (2011). https://doi.org/10.1002/jcc.21787
17. Huang, J. *et al.* CHARMM36m: an improved force field for folded and intrinsically disordered proteins. *Nat. Methods* **14**, 71–73 (2017). https://doi.org/10.1038/nmeth.4067
18. Salomon-Ferrer, R. *et al.* Routine microsecond molecular dynamics simulations with AMBER on GPUs. 2. Explicit solvent particle mesh Ewald. *J. Chem. Theory Comput.* **9**, 3878–3888 (2013). https://doi.org/10.1021/ct400314y
19. Abraham, M. J. *et al.* GROMACS: high performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX* **1–2**, 19–25 (2015). https://doi.org/10.1016/j.softx.2015.06.001
20. Rose, A. S. *et al.* NGL viewer: web-based molecular graphics for large complexes. *Bioinformatics* **34**, 3755–3758 (2018). https://doi.org/10.1093/bioinformatics/bty419
21. Wilkinson, M. D. *et al.* The FAIR Guiding Principles for scientific data management and stewardship. *Sci. Data* **3**, 160018 (2016). https://doi.org/10.1038/sdata.2016.18
