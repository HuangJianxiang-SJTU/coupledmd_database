#!/usr/bin/env python3
"""Replace the stale 208-system gateway consensus cache with final-207 data."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
SERVER = HERE.parent.parent
OVERLAP_ROOT = SERVER.parent / "a"
SOURCE = (
    OVERLAP_ROOT
    / "paper1_gateways"
    / "gateway_atlas_summary_final207.csv"
)
TARGET = SERVER / "data/api/v1/consensus/gateways.json"
REPORT = HERE / "reports/v12_gateway_consensus_cache_update.json"

EXPECTED_DENOMINATORS = {
    "ALL": 207,
    "ClassA": 181,
    "ClassB": 26,
    "Gi": 95,
    "Gs": 65,
    "Gq": 41,
    "G12-13": 6,
}
EXPECTED_COLUMNS = ["pair", "metric", "group", "n_systems", "mean", "ci_lo", "ci_hi"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def denominators(records: list[dict]) -> dict[str, list[int]]:
    output: dict[str, set[int]] = {}
    for record in records:
        output.setdefault(str(record["group"]), set()).add(
            int(record["n_systems"])
        )
    return {
        group: sorted(values)
        for group, values in sorted(output.items())
    }


def main() -> int:
    source = pd.read_csv(SOURCE)
    assert source.columns.tolist() == EXPECTED_COLUMNS
    assert len(source) == 7 * 7 * 4 == 196
    assert source.pair.nunique() == 7
    assert source.metric.nunique() == 4
    assert source.group.nunique() == 7
    assert not source.isna().any().any()
    for group, expected in EXPECTED_DENOMINATORS.items():
        values = source.loc[source.group.eq(group), "n_systems"].unique()
        assert values.tolist() == [expected], (group, values)
        assert len(source[source.group.eq(group)]) == 28

    before = json.loads(TARGET.read_text(encoding="utf-8"))
    before_hash = sha256(TARGET)
    generated_at = datetime.now(timezone.utc).isoformat()
    records = source.to_dict(orient="records")
    payload = {
        "_schema_version": "v12.1",
        "_generated_at": generated_at,
        "_source": (
            "paper1_gateways/gateway_atlas_summary_final207.csv; "
            "frozen final-207 processed full-system summary"
        ),
        "_scope_note": (
            "Processed intermediate records; values cannot be recalculated "
            "from the reduced protein-complex XTC release."
        ),
        "n_records": len(records),
        "columns": EXPECTED_COLUMNS,
        "records": records,
    }
    temporary = TARGET.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(TARGET)

    after = json.loads(TARGET.read_text(encoding="utf-8"))
    assert after["n_records"] == 196
    after_denominators = denominators(after["records"])
    assert after_denominators == {
        group: [value]
        for group, value in sorted(EXPECTED_DENOMINATORS.items())
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "generated_at": generated_at,
                "status": "PASS",
                "source": (
                    "read_only_overlap/"
                    "paper1_gateways/gateway_atlas_summary_final207.csv"
                ),
                "source_sha256": sha256(SOURCE),
                "target": str(TARGET.relative_to(SERVER)),
                "before_sha256": before_hash,
                "before_denominators": denominators(before["records"]),
                "after_sha256": sha256(TARGET),
                "after_denominators": after_denominators,
                "records": len(after["records"]),
                "biological_distribution_in_database_main_figures": False,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        "Updated gateway consensus cache: "
        f"ALL={after_denominators['ALL'][0]}, records={len(records)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
