#!/usr/bin/env python3
"""Verify public GitHub repository identity and current default-branch commits."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import requests


HERE = Path(__file__).resolve().parent
REPORT = HERE / "reports/v12_repository_access_verification.json"
REPOSITORIES = [
    "HuangJianxiang-SJTU/coupledmd",
    "HuangJianxiang-SJTU/coupledmd_database",
]


def main() -> int:
    session = requests.Session()
    session.headers.update(
        {
            "Accept": "application/vnd.github+json",
            "User-Agent": "CoupledMD-v12-repository-audit/1.0",
        }
    )
    records = []
    for name in REPOSITORIES:
        repository_response = session.get(
            f"https://api.github.com/repos/{name}",
            timeout=30,
        )
        repository_response.raise_for_status()
        repository = repository_response.json()
        branch = repository["default_branch"]
        commit_response = session.get(
            f"https://api.github.com/repos/{name}/commits/{branch}",
            timeout=30,
        )
        commit_response.raise_for_status()
        commit = commit_response.json()
        records.append(
            {
                "repository": name,
                "html_url": repository["html_url"],
                "private": repository["private"],
                "archived": repository["archived"],
                "default_branch": branch,
                "default_branch_commit": commit["sha"],
                "commit_html_url": commit["html_url"],
                "pushed_at": repository["pushed_at"],
                "license_spdx": (
                    repository.get("license", {}) or {}
                ).get("spdx_id"),
                "remotely_verified": True,
                "status": (
                    "PASS_PUBLIC"
                    if not repository["private"]
                    and not repository["archived"]
                    else "REVIEW"
                ),
            }
        )
    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "registry": "GitHub REST API",
        "status": (
            "PASS"
            if all(row["status"] == "PASS_PUBLIC" for row in records)
            else "REVIEW"
        ),
        "records": records,
        "paper_repository_v12_commit_present": False,
        "note": (
            "Public repository identity is verified. The v12 manuscript "
            "changes remain local and require a dedicated commit/tag."
        ),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"GitHub repository audit: {len(records)} public repositories; "
        f"status={payload['status']}"
    )
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
