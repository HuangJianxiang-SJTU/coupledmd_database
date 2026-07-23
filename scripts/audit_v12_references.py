#!/usr/bin/env python3
"""Validate manuscript DOI metadata against the Crossref REST API."""
from __future__ import annotations

import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests


HERE = Path(__file__).resolve().parent
MANUSCRIPT = HERE / "v12_manuscript.md"
REPORTS = HERE / "reports"


def normalize(value: str) -> str:
    value = re.sub(r"[*_`]", "", value)
    value = re.sub(r"<[^>]+>", " ", value)
    value = re.sub(r"[^a-z0-9]+", " ", value.lower())
    return " ".join(value.split())


def manuscript_references() -> list[dict]:
    text = MANUSCRIPT.read_text(encoding="utf-8").split("## References", 1)[1]
    rows = []
    pattern = re.compile(
        r"^(\d+)\.\s+(.+?)\s+\((\d{4})\)\."
        r"\s+https://doi\.org/(\S+)$"
    )
    for line in text.splitlines():
        if not re.match(r"^\d+\.", line):
            continue
        match = pattern.match(line)
        if match is None:
            raise ValueError(f"cannot parse reference: {line}")
        rows.append(
            {
                "reference_number": int(match.group(1)),
                "manuscript_reference": match.group(2),
                "manuscript_year": int(match.group(3)),
                "doi": match.group(4),
            }
        )
    return rows


def main() -> int:
    REPORTS.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = (
        "CoupledMD-v12-reference-audit/1.0 "
        "(mailto:lushaoyong@sjtu.edu.cn)"
    )
    output = []
    for row in manuscript_references():
        endpoint = f"https://api.crossref.org/works/{quote(row['doi'], safe='')}"
        response = session.get(endpoint, timeout=30)
        item = dict(row)
        item["http_status"] = response.status_code
        item["registry"] = "Crossref REST API"
        if response.status_code == 200:
            message = response.json()["message"]
            crossref_title = (message.get("title") or [""])[0]
            years = (
                message.get("published-print", {}).get("date-parts")
                or message.get("published-online", {}).get("date-parts")
                or message.get("issued", {}).get("date-parts")
                or [[]]
            )
            crossref_year = years[0][0] if years and years[0] else None
            title_in_reference = normalize(crossref_title) in normalize(
                row["manuscript_reference"]
            )
            item.update(
                {
                    "crossref_title": crossref_title,
                    "crossref_year": crossref_year,
                    "title_in_manuscript_reference": title_in_reference,
                    "doi_normalized": str(message.get("DOI", "")).lower(),
                    "status": (
                        "PASS"
                        if title_in_reference
                        and str(message.get("DOI", "")).lower()
                        == row["doi"].lower()
                        else "REVIEW"
                    ),
                }
            )
        else:
            item.update(
                {
                    "crossref_title": "",
                    "crossref_year": None,
                    "title_in_manuscript_reference": False,
                    "doi_normalized": "",
                    "status": "FAIL",
                }
            )
        output.append(item)
        time.sleep(0.05)
    status = "PASS" if all(row["status"] == "PASS" for row in output) else "REVIEW"
    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": str(MANUSCRIPT.name),
        "registry": "https://api.crossref.org",
        "status": status,
        "references": output,
    }
    (REPORTS / "v12_reference_doi_audit.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    pd.DataFrame(output).to_csv(
        REPORTS / "v12_reference_doi_audit.csv",
        index=False,
        lineterminator="\n",
    )
    print(
        f"Crossref audit: {len(output)} references; "
        f"{sum(row['status'] == 'PASS' for row in output)} PASS; "
        f"status={status}"
    )
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
