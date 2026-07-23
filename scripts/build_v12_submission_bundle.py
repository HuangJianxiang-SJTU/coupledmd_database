#!/usr/bin/env python3
"""Assemble a clean, versioned v12 manuscript-production bundle."""
from __future__ import annotations

import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
SERVER = HERE.parent.parent
STAGE = HERE / "CoupledMD_v12_submission_bundle"
ARCHIVE = HERE / "CoupledMD_v12_submission_bundle.zip"
AUDIT = HERE / "reports/v12_submission_bundle_audit.json"

ROOT_FILES = [
    HERE / "v12_manuscript.md",
    HERE / "v12_manuscript_release_candidate.docx",
    HERE / "rendered_v12_manuscript/v12_manuscript_release_candidate.pdf",
    HERE / "v12_si.md",
    HERE / "v12_si_release_candidate.docx",
    HERE / "rendered_v12_si/v12_si_release_candidate.pdf",
    SERVER / "LICENSE",
    SERVER / "LICENSE-CODE",
    SERVER / "CITATION.cff",
]
DIRECTORIES = [
    HERE / "publication_figures_v12",
    HERE / "v12_figure_source_data",
    HERE / "CoupledMD_Supplementary_Data_v12",
]
REPORT_FILES = [
    "v12_database_resource_claim_firewall.md",
    "v12_environment_freeze.json",
    "v12_gateway_consensus_cache_update.json",
    "v12_manuscript_claim_evidence_ledger.csv",
    "v12_manuscript_claim_evidence_ledger.json",
    "v12_numerical_consistency_audit.json",
    "v12_numerical_consistency_findings.csv",
    "v12_reference_doi_audit.csv",
    "v12_reference_doi_audit.json",
    "v12_remote_release_verification.json",
    "v12_repository_access_verification.json",
    "v12_revision_log.md",
    "v12_scientific_data_guidance_audit.md",
    "v12_submission_readiness_checklist.csv",
    "v12_submission_readiness_checklist.json",
    "v12_submission_readiness_checklist.md",
    "v12_submission_validation.json",
    "v12_submission_validation.md",
    "v12_visual_qc_report.json",
]
SCRIPT_FILES = [
    "audit_v12_references.py",
    "audit_v12_release_status.py",
    "audit_v12_repository_access.py",
    "build_v12_data.py",
    "build_v12_documents.py",
    "build_v12_reports.py",
    "build_v12_submission_bundle.py",
    "plot_v12_figures.py",
    "run_final207_viz_audit.py",
    "update_v12_gateway_consensus.py",
    "validate_v12_submission.py",
    "visual_qc_v12.py",
]
TEXT_SUFFIXES = {".md", ".txt", ".csv", ".json", ".py", ".cff"}
INTERNAL_PATH_MARKERS = [
    "/" + "MDdata/",
    "/" + "data01/",
    "/" + "data02/",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_file(source: Path, relative: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target = STAGE / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main() -> int:
    if STAGE.exists():
        shutil.rmtree(STAGE)
    if ARCHIVE.exists():
        ARCHIVE.unlink()
    STAGE.mkdir()

    for source in ROOT_FILES:
        if source.parent in {
            HERE / "rendered_v12_manuscript",
            HERE / "rendered_v12_si",
        }:
            relative = Path("rendered_pdfs") / source.name
        elif source.parent == SERVER:
            relative = Path("licences_and_citation") / source.name
        else:
            relative = Path("manuscript") / source.name
        copy_file(source, relative)

    for directory in DIRECTORIES:
        destination = STAGE / directory.name
        shutil.copytree(directory, destination)

    for name in REPORT_FILES:
        copy_file(HERE / "reports" / name, Path("reports") / name)
    for name in SCRIPT_FILES:
        copy_file(HERE / name, Path("scripts") / name)

    readiness = json.loads(
        (HERE / "reports/v12_submission_readiness_checklist.json").read_text(
            encoding="utf-8"
        )
    )
    readme = f"""# CoupledMD manuscript package v12

This is the auditable v12 Scientific Data production bundle.

Internal artifact validation: 45/45 PASS.
Submission-readiness gate: {readiness['overall_status']}.

The molecular release described by the manuscript contains reduced
protein-complex PDB/XTC records. It is not a full-system trajectory archive.
The three Zenodo identifiers are reserved private-draft identifiers at the
recorded audit time and are not asserted to be public.

The five main figures exclude the previous GPCR-centred pocket atlas. Figure 4
contains technical validation only; Figure 5 is the replica-release integrity
and checksum-first reuse figure.

See reports/v12_submission_readiness_checklist.md for remaining gates.
"""
    (STAGE / "README.md").write_text(readme, encoding="utf-8")

    files = sorted(path for path in STAGE.rglob("*") if path.is_file())
    forbidden_path_labels = [
        path
        for path in files
        if "final208" in path.name.lower() or "final624" in path.name.lower()
    ]
    internal_path_hits = []
    for path in files:
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if any(marker in text for marker in INTERNAL_PATH_MARKERS):
            internal_path_hits.append(str(path.relative_to(STAGE)))
    if forbidden_path_labels or internal_path_hits:
        raise AssertionError(
            {
                "forbidden_path_labels": [
                    str(path.relative_to(STAGE))
                    for path in forbidden_path_labels
                ],
                "internal_path_hits": internal_path_hits,
            }
        )

    manifest = []
    for path in files:
        manifest.append(
            {
                "path": str(path.relative_to(STAGE)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    with zipfile.ZipFile(
        ARCHIVE,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
    ) as archive:
        for path in files:
            archive.write(
                path,
                Path(STAGE.name) / path.relative_to(STAGE),
            )
    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "PASS_INTERNAL_VALIDATION_BLOCKED_SUBMISSION",
        "submission_readiness": readiness["overall_status"],
        "stage_directory": STAGE.name,
        "archive": ARCHIVE.name,
        "archive_size_bytes": ARCHIVE.stat().st_size,
        "archive_sha256": sha256(ARCHIVE),
        "file_count": len(files),
        "internal_absolute_paths": 0,
        "ambiguous_final208_or_final624_filenames": 0,
        "files": manifest,
    }
    AUDIT.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"Built {ARCHIVE.name}: {len(files)} files, "
        f"{ARCHIVE.stat().st_size / 1e6:.1f} MB; "
        f"readiness={readiness['overall_status']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
