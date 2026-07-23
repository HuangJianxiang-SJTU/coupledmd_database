#!/usr/bin/env python3
"""Read-only verification of the three CoupledMD Zenodo drafts and local release audits.

This script never creates, edits, publishes, or deletes a Zenodo deposition.  It
uses the authenticated deposition endpoint only to determine the current draft
state and uploaded-file inventory.  No access token is serialized or printed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests


HERE = Path(__file__).resolve().parent
REPORT = HERE / "reports" / "v12_remote_release_verification.json"
COORDINATION = HERE / "reports" / "zenodo_final207_three_replica_status_v1.json"

RECORDS = {
    1: {
        "draft_id": 21395292,
        "doi": "10.5281/zenodo.21395292",
        "local_dir": HERE / "zenodo_reduced_release_208",
        "local_label": "replica1_local_stage",
    },
    2: {
        "draft_id": 21447748,
        "doi": "10.5281/zenodo.21447748",
        "local_dir": HERE / "zenodo_reduced_release_207_replica2",
        "local_label": "replica2_local_stage",
    },
    3: {
        "draft_id": 21448037,
        "doi": "10.5281/zenodo.21448037",
        "local_dir": HERE / "zenodo_reduced_release_207_replica3",
        "local_label": "replica3_local_stage",
    },
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_token() -> str:
    token = os.environ.get("ZENODO_TOKEN")
    if token:
        return token.strip()
    token_path = Path.home() / ".zenodo_token"
    if not token_path.is_file():
        raise FileNotFoundError("ZENODO_TOKEN is unset and ~/.zenodo_token is absent")
    return token_path.read_text(encoding="utf-8").strip()


def local_evidence(directory: Path, label: str) -> dict[str, Any]:
    audit_path = directory / "CoupledMD_reduced_release_audit.json"
    manifest_path = directory / "CoupledMD_reduced_trajectory_manifest.csv"
    evidence: dict[str, Any] = {
        "local_location_label": label,
        "legacy_directory_name_withheld": True,
        "directory_exists": directory.is_dir(),
        "audit_exists": audit_path.is_file(),
        "manifest_exists": manifest_path.is_file(),
        "remotely_verified": False,
    }
    if audit_path.is_file():
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        evidence["audit"] = {
            "status": audit.get("status"),
            "cohort_systems": audit.get("cohort_systems")
            or audit.get("archive_system_total"),
            "molecular_files": audit.get("molecular_files"),
            "full_frames_checked": audit.get("full_frames_checked"),
            "frame_count_distribution": audit.get("frame_count_distribution"),
            "source_replica_distribution": audit.get("source_replica_distribution"),
            "checks": audit.get("checks"),
            "generated_at": audit.get("generated_at"),
            "zenodo": audit.get("zenodo"),
            "sha256": sha256(audit_path),
        }
    if manifest_path.is_file():
        with manifest_path.open(encoding="utf-8") as handle:
            row_count = max(0, sum(1 for _ in handle) - 1)
        evidence["manifest"] = {
            "rows": row_count,
            "sha256": sha256(manifest_path),
            "size_bytes": manifest_path.stat().st_size,
        }
    return evidence


def normalize_checksum(value: str | None) -> str | None:
    if not value:
        return None
    return value.removeprefix("md5:")


def authenticated_remote_record(
    session: requests.Session, token: str, draft_id: int
) -> dict[str, Any]:
    url = f"https://zenodo.org/api/deposit/depositions/{draft_id}"
    response = session.get(
        url, headers={"Authorization": f"Bearer {token}"}, timeout=60
    )
    response.raise_for_status()
    data = response.json()
    metadata = data.get("metadata", {})
    files = []
    for item in sorted(data.get("files", []), key=lambda value: value["filename"]):
        files.append(
            {
                "name": item["filename"],
                "size_bytes": item.get("filesize") or item.get("size"),
                "md5": normalize_checksum(item.get("checksum")),
            }
        )
    return {
        "remotely_verified": True,
        "verified_at": utc_now(),
        "endpoint": url,
        "state": data.get("state"),
        "submitted": bool(data.get("submitted")),
        "title": metadata.get("title"),
        "doi": metadata.get("doi"),
        "prereserve_doi": metadata.get("prereserve_doi"),
        "access_right": metadata.get("access_right"),
        "license": metadata.get("license"),
        "file_count": len(files),
        "files": files,
    }


def public_resolution(session: requests.Session, doi: str) -> dict[str, Any]:
    url = f"https://doi.org/{doi}"
    response = session.get(url, allow_redirects=True, timeout=60)
    return {
        "remotely_verified": True,
        "verified_at": utc_now(),
        "url": url,
        "http_status": response.status_code,
        "final_url": response.url,
        "publicly_resolvable": response.status_code == 200,
    }


def coordination_evidence() -> dict[str, Any]:
    if not COORDINATION.is_file():
        return {
            "available": False,
            "path": str(COORDINATION.relative_to(HERE)),
            "used_for_public_release_claims": False,
        }
    data = json.loads(COORDINATION.read_text(encoding="utf-8"))
    return {
        "available": True,
        "path": str(COORDINATION.relative_to(HERE)),
        "sha256": sha256(COORDINATION),
        "content": data,
        "used_for_public_release_claims": True,
    }


def derive_status(record: dict[str, Any]) -> str:
    remote = record["authenticated_remote"]
    public = record["public_resolution"]
    local = record["local_evidence"]
    if remote.get("submitted") and public.get("publicly_resolvable"):
        return "PUBLIC"
    if (
        local.get("audit", {}).get("status") == "passed"
        and remote.get("state") == "unsubmitted"
    ):
        return "PRIVATE_DRAFT_LOCAL_QC_PASS"
    if remote.get("state") == "unsubmitted":
        return "PRIVATE_DRAFT_INCOMPLETE"
    return "UNVERIFIED"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPORT)
    args = parser.parse_args()

    token = read_token()
    session = requests.Session()
    records: dict[str, Any] = {}
    for replica, config in RECORDS.items():
        record = {
            "replica": replica,
            "draft_id": config["draft_id"],
            "doi": config["doi"],
            "local_evidence": local_evidence(
                config["local_dir"], config["local_label"]
            ),
            "authenticated_remote": authenticated_remote_record(
                session, token, config["draft_id"]
            ),
            "public_resolution": public_resolution(session, config["doi"]),
        }
        record["status"] = derive_status(record)
        records[str(replica)] = record

    coordination = coordination_evidence()
    all_public = all(
        record["status"] == "PUBLIC" for record in records.values()
    )
    all_local_qc = all(
        record["local_evidence"].get("audit", {}).get("status") == "passed"
        for record in records.values()
    )
    payload = {
        "schema_version": "1.0",
        "generated_at": utc_now(),
        "mode": "read_only_remote_verification",
        "records": records,
        "coordination_file": coordination,
        "summary": {
            "all_three_public": all_public,
            "all_three_local_qc_pass": all_local_qc,
            "public_release_claim_allowed": all_public,
            "archive_completeness_claim_allowed": all_public and all_local_qc,
            "submission_access_status": (
                "PASS" if all_public else "HOLD_PRIVATE_OR_INCOMPLETE"
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2))
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
