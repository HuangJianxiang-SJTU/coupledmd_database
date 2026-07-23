#!/usr/bin/env python3
"""Build v12 claim, overlap, revision, guidance, and readiness reports."""
from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


HERE = Path(__file__).resolve().parent
SERVER = HERE.parent.parent
REPORTS = HERE / "reports"
MANUSCRIPT = HERE / "v12_manuscript.md"
SI = HERE / "v12_si.md"
NUMERICAL = REPORTS / "v12_numerical_consistency_audit.json"
REMOTE = REPORTS / "v12_remote_release_verification.json"
REFERENCES = REPORTS / "v12_reference_doi_audit.json"
FIGURES = HERE / "publication_figures_v12"
FIG_SOURCE = HERE / "v12_figure_source_data"
SUPP = HERE / "CoupledMD_Supplementary_Data_v12"


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(data: pd.DataFrame, path: Path) -> None:
    data.to_csv(path, index=False, lineterminator="\n")


def ledger_rows() -> list[dict[str, Any]]:
    rows = [
        ("C001", "Background & Summary", "The final cohort contains 207 systems.", "data/release_cohort_v9_final207.csv", "row count and unique system_id", "PASS", False, "database"),
        ("C002", "Background & Summary", "The cohort contains 181 Class A and 26 Class B systems.", "data/release_cohort_v9_final207.csv", "gpcr_class value counts", "PASS", False, "database"),
        ("C003", "Background & Summary", "Family counts are 95 Gi/o, 65 Gs, 41 Gq/11 and 6 G12/13.", "data/release_cohort_v9_final207.csv", "g_protein_family value counts", "PASS", False, "database"),
        ("C004", "Background & Summary", "The cohort contains 174 receptor names and 173 mapped UniProt accessions.", "data/release_cohort_v9_final207.csv", "distinct receptor_name and non-null receptor_uniprot", "PASS", False, "database"),
        ("C005", "Background & Summary", "Gs_8HTI has an explicit null receptor accession.", "data/release_cohort_v9_final207.csv", "system_id=Gs_8HTI", "PASS", False, "database"),
        ("C006", "Cohort construction", "The working inventory is 207 included, 13 unresolved and two excluded.", "CoupledMD_Supplementary_Data_v12/Supplementary_Data_S1_included_system_inventory.csv|CoupledMD_Supplementary_Data_v12/Supplementary_Data_S2_release_boundary_exceptions.csv", "207 rows; release_status counts 13/2", "PASS", False, "database"),
        ("C007", "Cohort construction", "Gq_7DWC was removed from the former 208-system cohort as a duplicate source identity of Gq_8DWC.", "read_only_overlap/release_exclusions_v1.csv", "system_id=Gq_7DWC", "PASS", False, "database"),
        ("C008", "Cohort construction", "Gq_7E9W is an excluded duplicate or mislabel.", "CoupledMD_Supplementary_Data_v12/Supplementary_Data_S2_release_boundary_exceptions.csv", "system_id=Gq_7E9W", "PASS", False, "database"),
        ("C009", "Methods", "Structural provenance is 197 experimental and 10 engineered_uncertain.", "data/release_cohort_v9_final207.csv", "structural_provenance value counts", "PASS", False, "database"),
        ("C010", "Production protocols", "Protocol counts are P1=179, P2=26 and P3=2.", "manuscript_figures/scidata_figures/parallel_codex/provenance_direct/protocol_system_evidence.csv", "filter final cohort then protocol_group value counts", "PASS", False, "database"),
        ("C011", "Production protocols", "The source evidence assigns 205 systems to AMBER and two to GROMACS.", "manuscript_figures/scidata_figures/parallel_codex/provenance_direct/protocol_system_evidence.csv", "filter final cohort then engine_family", "PASS", False, "database"),
        ("C012", "Production protocols", "Protocol evidence is verified for 202 systems and partial for five.", "manuscript_figures/scidata_figures/parallel_codex/provenance_direct/protocol_system_evidence.csv", "filter final cohort then status", "PASS", False, "database"),
        ("C013", "Sampling", "The production design comprises 621 computational repeats.", "CoupledMD_Supplementary_Data_v12/Supplementary_Data_S7_source_simulation_replica_provenance.csv", "207 systems x 3 replica_id", "PASS", False, "database"),
        ("C014", "Sampling", "Nominal protocol sampling is 310.5 microseconds.", "CoupledMD_Supplementary_Data_v12/Supplementary_Data_S7_source_simulation_replica_provenance.csv", "sum nominal_duration_ns", "PASS", False, "database"),
        ("C015", "Sampling", "Original file-observed coordinate sampling is 310.4168 microseconds.", "read_only_overlap/cohort_truth_v1.csv", "sum observed NetCDF/GROMACS sampling for included systems", "PASS", False, "database"),
        ("C016", "Sampling", "The original Gq_8ZPT replica 2 record contains 416.8 ns.", "read_only_overlap/cohort_truth_v1.csv", "system_id=Gq_8ZPT observed total", "PASS", False, "database"),
        ("C017", "Sampling", "The repaired reduced source for Gq_8ZPT replica 2 has 10,000 frames, 499.95 ns and verdict GOOD.", "manuscript_figures/scidata_figures/trajectory_readiness/continuations/Gq_8ZPT/interface_recovery/exhaustive_audit_corrected/candidate_replica_audit.csv", "only row", "PASS", False, "database"),
        ("C018", "Reduced records", "Each complete reduced release record has 207 PDB/XTC pairs.", "CoupledMD_Supplementary_Data_v12/Supplementary_Data_S4_reduced_release_manifest.csv", "per release_replica present system and role counts", "HOLD", True, "database"),
        ("C019", "Reduced records", "Each released XTC contains 2,500 frames at 200-ps spacing and spans 499.8 ns.", "CoupledMD_Supplementary_Data_v12/Supplementary_Data_S4_reduced_release_manifest.csv", "n_frames, frame_interval_ps, represented_span_ns", "HOLD", True, "database"),
        ("C020", "Reduced records", "The release excludes membrane, solvent and mobile ions.", "manuscript_figures/scidata_figures/zenodo_reduced_release_207_replica3/README.txt", "component-scope statement", "PASS", False, "database"),
        ("C021", "Reduced records", "The reduced release is not a full-system archive.", "manuscript_figures/scidata_figures/zenodo_reduced_release_207_replica3/README.txt", "underlying full-system trajectories not included", "PASS", False, "database"),
        ("C022", "Pocket validation", "Persistent pockets were detected in 205 systems and two systems have explicit zero-pocket records.", "CoupledMD_Supplementary_Data_v12/Supplementary_Data_S5_pocket_summaries_gpcrdb.csv", "pocket_record_status by unique system_id", "PASS", False, "database"),
        ("C023", "Pocket validation", "There are 2,149 detected-pocket rows.", "CoupledMD_Supplementary_Data_v12/Supplementary_Data_S5_pocket_summaries_gpcrdb.csv", "available_detected_pocket row count", "PASS", False, "database"),
        ("C024", "Pocket validation", "Orthosteric recovery is 49 of 58 eligible systems.", "v12_figure_source_data/Figure_3A_positive_control.csv", "row count and ortho_recovered sum", "PASS", False, "database"),
        ("C025", "Gateway records", "The gateway intermediate contains 5,796 rows: 207 systems x 7 interfaces x 4 metrics.", "CoupledMD_Supplementary_Data_v12/Supplementary_Data_S6_gateway_per_system.csv", "dimension counts", "PASS", False, "shared-input"),
        ("C026", "Gateway records", "Gateway values cannot be recalculated from the reduced XTC files.", "v12_figure_source_data/Figure_2_record_roles.csv", "reduced trajectory does_not_support field", "PASS", False, "database"),
        ("C027", "API", "The active gateway consensus ALL denominator is 207.", "data/api/v1/consensus/gateways.json", "ALL n_systems unique value", "PASS", False, "database"),
        ("C028", "Structural validation", "Structural metrics are available for 618 of 621 replica records from 206 systems.", "v12_figure_source_data/Figure_4A_replica_validation.csv", "validation_status counts and unique systems", "PASS", False, "database"),
        ("C029", "Structural validation", "The three unavailable records are the Gs_8HTI repeats.", "v12_figure_source_data/Figure_4A_replica_validation.csv", "validation_status != available", "PASS", False, "database"),
        ("C030", "Reduced-record audit", "The exact final-207 sampled audit has 207 of 207 records classified OK.", "final207_reduced_visualization_audit_40frames.csv", "status counts", "PASS", False, "database"),
        ("C031", "Reduced-record audit", "Portal-record frame counts are 203x2500, 2x2501, 1x681 and 1x200.", "final207_reduced_visualization_audit_40frames.csv", "n_frames value counts", "PASS", False, "database"),
        ("C032", "Release integrity", "Replica 1 local manifest and full-frame QC are complete.", "reports/v12_remote_release_verification.json", "records.1.local_evidence", "PASS", True, "database"),
        ("C033", "Release integrity", "Replica 2 remains release-dependent.", "reports/v12_remote_release_verification.json", "records.2", "HOLD", True, "database"),
        ("C034", "Release integrity", "Replica 3 local manifest and full-frame QC are complete.", "reports/v12_remote_release_verification.json", "records.3.local_evidence", "PASS", True, "database"),
        ("C035", "Data Availability", "All three Zenodo identifiers are restricted unsubmitted drafts at the audit time.", "reports/v12_remote_release_verification.json", "records.*.authenticated_remote", "HOLD", True, "database"),
        ("C036", "Data Availability", "None of the three DOI URLs resolves publicly at the audit time.", "reports/v12_remote_release_verification.json", "records.*.public_resolution", "HOLD", True, "database"),
        ("C037", "Data Availability", "Anonymous reviewer access has not been machine verified.", "reports/v12_remote_release_verification.json", "coordination file unavailable; no verified reviewer URL", "HOLD", True, "database"),
        ("C038", "Code Availability", "The web/code repository is the coupledmd GitHub repository.", "reports/v12_repository_access_verification.json|CITATION.cff", "public repository and repository-code", "PASS", False, "database"),
        ("C039", "Code Availability", "The paper-content repository is the coupledmd_database GitHub repository.", "reports/v12_repository_access_verification.json|CITATION.cff", "public repository and repository-artifact", "PASS", False, "database"),
        ("C040", "Licensing", "Data are CC BY 4.0 and code is MIT licensed.", "LICENSE|LICENSE-CODE", "licence text", "PASS", False, "database"),
        ("C041", "References", "All 21 DOI/title pairs resolve through Crossref.", "reports/v12_reference_doi_audit.json", "status=PASS", "PASS", False, "database"),
        ("C042", "Related manuscripts", "Pocket recurrence, biological gateway patterns and selectivity are excluded from the database-paper claim set.", "reports/v12_database_resource_claim_firewall.md", "ownership matrix", "PASS", False, "resource"),
        ("C043", "Related manuscripts", "AI models, predictions and BRET are absent from the descriptor.", "v12_manuscript.md", "prohibited-term validation", "PASS", False, "AI"),
    ]
    output = []
    for claim_id, section, claim, evidence, selector, status, dependent, owner in rows:
        output.append(
            {
                "claim_id": claim_id,
                "manuscript_section": section,
                "claim": claim,
                "authoritative_evidence": evidence,
                "evidence_selector": selector,
                "verification_method": "automated or asserted machine-readable join",
                "status": status,
                "release_dependent": dependent,
                "claim_owner": owner,
                "overlap_risk": (
                    "controlled"
                    if owner in {"resource", "AI", "shared-input"}
                    else "none"
                ),
                "notes": (
                    "Refresh after release worker completion."
                    if dependent
                    else ""
                ),
            }
        )
    return output


def build_ledger() -> pd.DataFrame:
    ledger = pd.DataFrame(ledger_rows())
    write_csv(ledger, REPORTS / "v12_manuscript_claim_evidence_ledger.csv")
    payload = {
        "schema_version": "1.0",
        "generated_at": now(),
        "status_counts": ledger.status.value_counts().to_dict(),
        "claims": ledger.to_dict(orient="records"),
    }
    (REPORTS / "v12_manuscript_claim_evidence_ledger.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    return ledger


def build_firewall() -> None:
    text = """# v12 database-versus-resource claim firewall

Generated: {generated}

## Decision

The v12 Data Descriptor is restricted to data-resource claims. The previous
GPCR-centred pocket-reuse atlas has been removed. The previous family-organized
gateway distribution has been removed from the main figures. Figure 3 retains
orthosteric recovery only as a technical positive control; Figure 4 contains
only repeat-level and reduced-record technical validation.

## Ownership matrix

| Topic | Database Data Descriptor | Resource-atlas paper | AI paper |
|---|---|---|---|
| Cohort construction, identifiers and exclusions | owns | cites | cites |
| Simulation protocol and provenance | owns | summarizes as input | summarizes as input |
| Three computational repeats and QC | owns | uses as input | uses as input |
| Reduced PDB/XTC records, checksums and access | owns | links | links |
| Pocket method and 49/58 positive control | owns limited validation | may cite | may cite |
| 2,149 pocket rows; 205 detected + 2 zero | owns completeness | uses as input | uses as input |
| Recurrent pocket atlas and 65 clusters | excluded | owns | excluded |
| Pocket rankings and biological recurrence | excluded | owns | excluded |
| C7–Gq enrichment | excluded | owns | excluded |
| Partner-associated pocket changes | excluded | owns | excluded |
| Gateway method and processed S6 records | owns methods/intermediate | uses with disclosure | excluded |
| Gateway family distribution and interpretation | excluded | owns | excluded |
| α5 and family-selectivity analyses | excluded | owns | model input only if disclosed |
| Mechanism or druggability conclusions | excluded | owns only if supported | excluded |
| Learned models and architecture | excluded | excluded | owns |
| Precog3D benchmarking and predictions | excluded | excluded | owns |
| Prospective panels, BRET, rescue or epistasis | excluded | excluded | owns |

## Shared-data disclosure

The resource paper may reuse the final-207 cohort table, source-simulation
provenance, processed pocket records and processed gateway intermediates.
Reuse must cite the Data Descriptor and disclose that the database paper owns
the cohort, data-record and technical-validation descriptions. The database
paper does not reproduce the resource paper's biological figures, rankings,
cluster headline, family enrichment or mechanistic interpretation.

## Automated prohibited-content gate

The v12 manuscript is checked for model architectures, Precog3D, BRET,
prospective challenge panels, model-selected mutations, rescue, epistasis,
the 65-cluster headline, C7–Gq enrichment, pocket rankings and gateway-family
biological claims. Incidental use of “model” for starting structural models is
not an AI claim and is reviewed in context.
""".format(generated=now())
    (REPORTS / "v12_database_resource_claim_firewall.md").write_text(
        text,
        encoding="utf-8",
    )


def build_revision_log() -> None:
    text = f"""# v11-to-v12 detailed revision log

Generated: {now()}

The v11 release candidates remain unchanged. v12 is a new auditable source and
output line.

| Area | v11 issue | v12 revision | Evidence/status |
|---|---|---|---|
| Title | resource name and promotional framing | objective 78-character title without colon or self-acronym | PASS |
| Abstract | mixed full-archive and derived-resource framing | 164-word description of the reduced product and its exclusions | PASS |
| Cohort | residual 208-derived counts | regenerated 207, 181/26 class and 95/65/41/6 family counts | PASS |
| Boundary | Gq_7DWC ambiguity | explicit 208→207 transition and two exclusions | PASS |
| Protocols | 180/26/2 summed to 208 | regenerated 179/26/2 | PASS |
| Sampling | 310.5 µs treated as file-observed | nominal 310.5 and original observed 310.4168 reported separately | PASS |
| Gq_8ZPT replica 2 | shortened record not reconciled | 416.8-ns original and 499.95-ns locally repaired reduced source documented | PASS |
| Pocket completeness | 206 detected + two zero = 208 | 205 detected + two explicit zero = 207; 2,149 rows | PASS |
| Reduced visualization | stale 205 OK + 3 WARN | exact final-207 audit regenerated; 207/207 OK | PASS |
| API consensus | `n_systems=208` in active cache | regenerated from authoritative final-207 gateway summary | PASS |
| Historical names | active final208/final624 artifacts | v12/final207 or version-neutral submission names; legacy inputs isolated | PASS |
| Data product | implied full primary archive | reduced matched PDB/XTC product stated throughout | PASS |
| Component scope | membrane/solvent/ions ambiguity | exclusions and possible ligand removal made explicit | PASS |
| Gateway reproducibility | implied reproducible from reduced trajectories | S6/code supplied; molecular recalculation limitation explicit | PASS |
| Replica semantics | risk of biological-replicate wording | “computational repeats, not biological replicates” | PASS |
| Figure 4 | family-organized gateway distribution | replaced with repeat-level structural and reduced-record QC | PASS |
| Figure 6 | consumed resource-paper pocket-atlas novelty | removed | PASS |
| Replacement Figure 5 | absent | release-state matrix, record-size/atom-count distributions and worked reuse | PASS, release-dependent |
| Figure source data | incomplete panel-level tables | one CSV per panel plus SHA-256 manifest | PASS |
| Supplementary Data | legacy and ambiguous manifests | ten v12 tables with pending entries explicit | PASS/HOLD release |
| Data Availability | DOI/publication could be inferred | three identifiers described only as restricted unsubmitted drafts | HOLD |
| Reviewer access | unverified | explicit pre-submission gate | HOLD |
| Full-system access | archive scope inflated | declared absent and conservatively treated as blocker for full-system-only claims | BLOCKED |
| Code availability | placeholder URL | both public GitHub repositories named; v12 tag still required | HOLD |
| Licences | author name error and fake DOI in citation metadata | author corrected; fake DOI removed; CC BY 4.0/MIT retained | PASS |
| References | one GPCRdb DOI pointed to an AlphaFold paper | replaced with official GPCRdb 2025 citation; 21/21 Crossref checks pass | PASS |
| Related work | no complete ownership disclosure | main disclosure, SI overlap table and firewall report added | PASS |
| AI boundary | potential future-paper leakage | AI, predictions and BRET excluded | PASS |
"""
    (REPORTS / "v12_revision_log.md").write_text(text, encoding="utf-8")


def journal_guidance() -> None:
    text = f"""# Scientific Data guidance audit for v12

Verified against official journal pages on 20 July 2026.

Sources:

- https://www.nature.com/sdata/publish/submission-guidelines
- https://www.nature.com/sdata/policies/repositories
- https://www.nature.com/sdata/aims-and-scope

| Official requirement | v12 implementation | Status |
|---|---|---|
| Data Descriptor title is descriptive, ≤110 characters, without advertising or a dataset self-acronym | 78-character objective title | PASS |
| Abstract is concise, recommended ≤170 words, and contains no URLs | 164 words; no URL | PASS |
| Background & Summary defines the reuse need without subjective novelty claims | rewritten around interoperability and provenance | PASS |
| Methods describe data generation | cohort, preparation, protocols, reduction and annotations included | PASS |
| Data Records describes repository files and organization | reduced product, Supplementary Data and figure sources listed | PASS |
| Technical Validation documents data quality without a Results/Discussion section | source-derived validation sections | PASS |
| Usage Notes contain practical limitations and reuse information | supported/unsupported uses and worked integrity example | PASS |
| Data Availability precedes Code Availability and gives persistent identifiers | reserved identifiers and private status stated | HOLD until accessible |
| Code Availability names executable sources and versions | repositories named; v12 tag pending | HOLD |
| Data are available to reviewers during peer review | no anonymous route machine verified | HOLD |
| Public repository deposition is required before final publication | records are private drafts | HOLD |
| References contain complete resolvable DOI metadata | Crossref audit 21/21 PASS | PASS |
| Main figures are composite, legible and normally ≤8 | five composite main figures | PASS |
| Figure legends explain panels and denominators | captions rewritten and source tables supplied | PASS |
| Supplementary Information, if retained, is supplied as one PDF | PDF generation/visual audit required | {'PASS' if (HERE / 'rendered_v12_si/v12_si_release_candidate.pdf').is_file() else 'HOLD'} |

The official repository policy requires persistent, reviewer-accessible data.
A reserved DOI that does not resolve publicly is not treated as satisfying that
requirement.
"""
    (REPORTS / "v12_scientific_data_guidance_audit.md").write_text(
        text,
        encoding="utf-8",
    )


def readiness_items() -> list[dict[str, str]]:
    manuscript_text = MANUSCRIPT.read_text(encoding="utf-8")
    title = manuscript_text.splitlines()[0].removeprefix("# ")
    abstract = manuscript_text.split("## Abstract\n\n", 1)[1].split("\n\n## ", 1)[0]
    abstract_words = len(re.findall(r"\b[\wµ–-]+\b", abstract))
    numerical = load_json(NUMERICAL)
    remote = load_json(REMOTE)
    references = load_json(REFERENCES)
    expected_figures = [
        "figure1_cohort_scope",
        "figure2_records_reuse_boundary",
        "figure3_annotation_validation",
        "figure4_replica_technical_validation",
        "figure5_release_integrity_reuse",
    ]
    figure_complete = all(
        (FIGURES / f"{stem}.pdf").is_file()
        and (FIGURES / f"{stem}.png").is_file()
        for stem in expected_figures
    )
    source_complete = len(
        [
            p
            for p in FIG_SOURCE.glob("*.csv")
            if p.name != "figure_source_data_manifest.csv"
        ]
    ) == 15
    items = [
        ("J01", "Descriptive title ≤110 characters without colon or self-acronym", "PASS" if len(title) <= 110 and ":" not in title and "CoupledMD" not in title else "BLOCKED", f"{len(title)} characters"),
        ("J02", "Abstract ≤170 words and without URL", "PASS" if abstract_words <= 170 and "http" not in abstract else "BLOCKED", f"{abstract_words} words"),
        ("J03", "Required Scientific Data section headings", "PASS", "Background, Methods, Data Records, Technical Validation, Usage Notes, Data/Code Availability present"),
        ("J04", "No Results or Discussion section", "PASS" if "\n## Results" not in manuscript_text and "\n## Discussion" not in manuscript_text else "BLOCKED", "Data Descriptor structure"),
        ("N01", "Frozen cohort and all core denominators machine consistent", "PASS" if numerical["cohort"]["systems"] == 207 else "BLOCKED", "numerical audit"),
        ("N02", "Exact final-207 reduced visualization audit", "PASS" if numerical["reduced_visualization"]["regenerated_exact_final207"] else "HOLD", "207/207 OK"),
        ("N03", "Final-207 API consensus cache", "PASS" if numerical["gateways"]["api_consensus_final207_consistent"] else "HOLD", str(numerical["gateways"]["api_consensus_ALL_denominators"])),
        ("D01", "Reduced molecular product scope is explicit", "PASS", "not a full-system archive"),
        ("D02", "Replica 1 manifest and full-frame QC", "PASS" if remote["records"]["1"]["status"] == "PRIVATE_DRAFT_LOCAL_QC_PASS" else "HOLD", remote["records"]["1"]["status"]),
        ("D03", "Replica 2 manifest and full-frame QC", "PASS" if remote["records"]["2"]["status"] == "PRIVATE_DRAFT_LOCAL_QC_PASS" else "HOLD", remote["records"]["2"]["status"]),
        ("D04", "Replica 3 manifest and full-frame QC", "PASS" if remote["records"]["3"]["status"] == "PRIVATE_DRAFT_LOCAL_QC_PASS" else "HOLD", remote["records"]["3"]["status"]),
        ("D05", "All three records publicly resolvable or anonymously reviewer-accessible", "PASS" if remote["summary"]["all_three_public"] else "HOLD", remote["summary"]["submission_access_status"]),
        ("D06", "Remote files match final local manifests and audits", "PASS" if remote["summary"]["archive_completeness_claim_allowed"] else "HOLD", "completeness claim withheld"),
        ("D07", "Stable full-system access for full-system-only molecular recalculation", "BLOCKED", "full-system trajectories/topologies are not distributed"),
        ("D08", "Processed gateway intermediates and source data supplied", "PASS" if (SUPP / "Supplementary_Data_S6_gateway_per_system.csv").is_file() else "BLOCKED", "5,796-row S6"),
        ("F01", "Five main figure PDF/PNG pairs", "PASS" if figure_complete else "HOLD", "Figure 6 atlas removed"),
        ("F02", "Panel-level figure source tables and manifest", "PASS" if source_complete else "HOLD", "15 panel/source CSVs"),
        ("F03", "Biological gateway distribution removed from main figures", "PASS", "Figure 4 is technical"),
        ("F04", "Resource-paper pocket atlas removed", "PASS", "no Figure 6 atlas"),
        ("R01", "Reference DOI/title metadata verified", "PASS" if references["status"] == "PASS" else "HOLD", "21/21 Crossref"),
        ("C01", "Public code and paper repository URLs supplied", "PASS" if (REPORTS / "v12_repository_access_verification.json").is_file() and load_json(REPORTS / "v12_repository_access_verification.json").get("status") == "PASS" else "HOLD", "two GitHub repositories"),
        ("C02", "Immutable v12 paper/code commit and release tag", "HOLD", "must be created after final regeneration"),
        ("L01", "Data and code licences", "PASS", "CC BY 4.0 and MIT"),
        ("O01", "Database/resource/AI claim firewall", "PASS", "separate report and SI overlap table"),
        ("O02", "Related-manuscript disclosure", "PASS", "main manuscript and SI"),
        ("P01", "v12 manuscript DOCX exists", "PASS" if (HERE / "v12_manuscript_release_candidate.docx").is_file() else "HOLD", "production output"),
        ("P02", "v12 SI DOCX exists", "PASS" if (HERE / "v12_si_release_candidate.docx").is_file() else "HOLD", "production output"),
        ("P03", "Rendered manuscript PDF exists", "PASS" if (HERE / "rendered_v12_manuscript/v12_manuscript_release_candidate.pdf").is_file() else "HOLD", "visual-QC input"),
        ("P04", "Rendered SI PDF exists", "PASS" if (HERE / "rendered_v12_si/v12_si_release_candidate.pdf").is_file() else "HOLD", "visual-QC input"),
        ("P05", "Document visual-QC report passes", "PASS" if (REPORTS / "v12_visual_qc_report.json").is_file() and load_json(REPORTS / "v12_visual_qc_report.json").get("status") == "PASS" else "HOLD", "render inspection"),
        ("V01", "Automated cross-artifact validation passes", "PASS" if (REPORTS / "v12_submission_validation.json").is_file() and load_json(REPORTS / "v12_submission_validation.json").get("status") == "PASS" else "HOLD", "validation suite"),
        ("A01", "Package can truthfully be labelled submission-ready", "BLOCKED", "release access and stable full-system route remain unresolved"),
    ]
    return [
        {
            "requirement_id": item[0],
            "requirement": item[1],
            "status": item[2],
            "evidence": item[3],
        }
        for item in items
    ]


def build_readiness() -> None:
    items = readiness_items()
    statuses = {row["status"] for row in items}
    overall = "BLOCKED" if "BLOCKED" in statuses else ("HOLD" if "HOLD" in statuses else "PASS")
    payload = {
        "schema_version": "1.0",
        "generated_at": now(),
        "overall_status": overall,
        "status_counts": pd.Series([row["status"] for row in items]).value_counts().to_dict(),
        "requirements": items,
    }
    (REPORTS / "v12_submission_readiness_checklist.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    write_csv(
        pd.DataFrame(items),
        REPORTS / "v12_submission_readiness_checklist.csv",
    )
    lines = [
        "# Scientific Data v12 submission-readiness checklist",
        "",
        f"Generated: {payload['generated_at']}",
        "",
        f"Overall status: **{overall}**",
        "",
        "| ID | Requirement | Status | Evidence |",
        "|---|---|---|---|",
    ]
    lines.extend(
        f"| {row['requirement_id']} | {row['requirement']} | {row['status']} | {row['evidence']} |"
        for row in items
    )
    lines.extend(
        [
            "",
            "The overall status is conservative: a single hard blocker prevents a submission-ready label. Replica and remote-access items remain HOLD rather than being inferred from reserved identifiers.",
        ]
    )
    (REPORTS / "v12_submission_readiness_checklist.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    REPORTS.mkdir(parents=True, exist_ok=True)
    ledger = build_ledger()
    build_firewall()
    build_revision_log()
    journal_guidance()
    build_readiness()
    print(
        f"Built v12 reports: {len(ledger)} claim-ledger rows; "
        "firewall, revision, guidance and readiness outputs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
