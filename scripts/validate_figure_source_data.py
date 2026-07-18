#!/usr/bin/env python3
"""Independently audit exported source data for main Figures 5 and 6."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import re


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = (
    SCRIPT_DIR / "source_data"
    if (SCRIPT_DIR / "source_data").is_dir()
    else SCRIPT_DIR.parent / "source_data"
)
OUTPUT = ROOT / "independent_audit.json"
PRIVATE_PATTERN = re.compile(r"/MDdata/|(?:^|[^A-Za-z0-9])(?:v9|final208)(?:[^A-Za-z0-9]|$)")


def read_csv(name: str) -> list[dict[str, str]]:
    with (ROOT / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def is_true(value: str) -> bool:
    return value.lower() == "true"


def main() -> None:
    global ROOT, OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    ROOT = args.source_dir
    OUTPUT = args.output or ROOT / "independent_audit.json"
    findings: list[str] = []

    manifest = read_csv("figure_source_data_manifest.csv")
    for row in manifest:
        path = ROOT / row["file"]
        if not path.is_file():
            findings.append(f"manifest file missing: {row['file']}")
            continue
        rows = sum(1 for _ in path.open(encoding="utf-8")) - 1
        if rows != int(row["rows"]):
            findings.append(f"row-count mismatch: {row['file']}")
        if digest(path) != row["sha256"]:
            findings.append(f"checksum mismatch: {row['file']}")

    for path in ROOT.iterdir():
        if path.suffix.lower() not in {".csv", ".md"}:
            continue
        text = path.read_text(encoding="utf-8")
        if PRIVATE_PATTERN.search(text):
            findings.append(f"private or internal-version token in {path.name}")

    composition = read_csv("Figure_1_system_composition_source.csv")
    ligand_context = read_csv("Figure_1_ligand_context_source.csv")
    if len(composition) != 207 or len({row["system_id"] for row in composition}) != 207:
        findings.append("Figure 1 composition is not a 207-system unique cohort")
    if Counter(row["g_protein_family"] for row in composition) != Counter(
        {"Gi": 95, "Gs": 65, "Gq": 41, "G12-13": 6}
    ):
        findings.append("Figure 1 family counts differ from 95/65/41/6")
    if Counter(row["gpcr_class"] for row in composition) != Counter({"A": 181, "B": 26}):
        findings.append("Figure 1 GPCR-class counts differ from 181/26")
    if len({row["receptor_name"] for row in composition}) != 174:
        findings.append("Figure 1 does not contain 174 distinct receptor names")
    mapped_accessions = {row["receptor_uniprot"] for row in composition if row["receptor_uniprot"]}
    if len(mapped_accessions) != 173:
        findings.append("Figure 1 does not contain 173 mapped UniProt accessions")
    if len(ligand_context) != 207 or Counter(
        row["ligand_context"] for row in ligand_context
    ) != Counter({"none": 149, "peptide": 58}):
        findings.append("Figure 1 ligand context differs from 149 none/58 peptide")

    release_facts = {
        row["quantity"]: (float(row["value"]), row["unit"])
        for row in read_csv("Figure_2_release_architecture_source.csv")
    }
    expected_release_facts = {
        "included systems": (207.0, "systems"),
        "source simulation replicas": (621.0, "replicas"),
        "released molecular files": (414.0, "files"),
        "released family archives": (4.0, "ZIP archives"),
        "source replica for released pair": (1.0, "replica"),
        "frames per released trajectory": (2500.0, "frames"),
    }
    for key, expected in expected_release_facts.items():
        if release_facts.get(key) != expected:
            findings.append(f"Figure 2 release fact mismatch: {key}")

    positive = read_csv("Figure_3_positive_control_source.csv")
    availability = read_csv("Figure_3_pocket_availability_source.csv")
    if len(positive) != 58 or sum(is_true(row["ortho_recovered"]) for row in positive) != 49:
        findings.append("Figure 3 positive control differs from 49/58")
    if Counter(row["g_protein_family"] for row in positive) != Counter(
        {"Gi": 18, "Gs": 18, "Gq": 22}
    ):
        findings.append("Figure 3 positive-control family counts differ from 18/18/22")
    if len(availability) != 207 or sum(int(row["detected_pockets"]) for row in availability) != 2149:
        findings.append("Figure 3 availability does not represent 2,149 pockets across 207 systems")
    if sum(int(row["detected_pockets"]) == 0 for row in availability) != 2:
        findings.append("Figure 3 availability does not contain two zero-pocket systems")

    gateway = read_csv("Figure_4_gateway_open_fraction_source.csv")
    gateway_keys = {
        (row["system_id"], row["interface"])
        for row in gateway
    }
    if len(gateway) != 1449 or len(gateway_keys) != 1449:
        findings.append("Figure 4 does not contain 1,449 unique system-interface rows")
    gateway_by_system: dict[str, int] = Counter(row["system_id"] for row in gateway)
    if len(gateway_by_system) != 207 or set(gateway_by_system.values()) != {7}:
        findings.append("Figure 4 does not contain seven interfaces for every system")
    for row in gateway:
        values = [float(row[key]) for key in ("mean", "ci_lo", "ci_hi")]
        if any(not 0 <= value <= 1 for value in values):
            findings.append("Figure 4 contains an open fraction outside [0, 1]")
            break

    profiles = read_csv("Supplementary_Figure_S1_receptor_profile_source.csv")
    if len(profiles) != 174 or Counter(
        row["uniprot_mapping_status"] for row in profiles
    ) != Counter({"mapped": 173, "unmapped": 1}):
        findings.append("Supplementary Figure S1 receptor profiles are incomplete")
    gateway_method = read_csv("Supplementary_Figure_S2_gateway_method_source.csv")
    method_names = {row["parameter_or_metric"] for row in gateway_method}
    required_method_names = {
        "helix-pair distance cutoff",
        "open-frame penetration threshold",
        "adjacent TM interfaces",
        "source replicas",
        "bootstrap resamples",
        "penetration",
        "penetration_p90",
        "open_fraction",
        "occupancy",
    }
    if len(gateway_method) != 12 or not required_method_names.issubset(method_names):
        findings.append("Supplementary Figure S2 method dictionary is incomplete")

    ledger = read_csv("Figure_5_replica_validation_ledger.csv")
    panel5 = read_csv("Figure_5_panel_source_data.csv")
    excluded5 = read_csv("Figure_5_exclusion_audit.csv")
    ledger_keys = {(row["system_id"], row["replica"]) for row in ledger}
    panel5_keys = {(row["system_id"], row["replica"]) for row in panel5}
    excluded5_keys = {(row["system_id"], row["replica"]) for row in excluded5}
    if len(ledger) != 621 or len(ledger_keys) != 621:
        findings.append("Figure 5 ledger is not 621 unique replica records")
    replicas_by_system: dict[str, set[str]] = defaultdict(set)
    for row in ledger:
        replicas_by_system[row["system_id"]].add(row["replica"])
    if len(replicas_by_system) != 207 or any(
        replicas != {"1", "2", "3"} for replicas in replicas_by_system.values()
    ):
        findings.append("Figure 5 ledger does not contain replicas 1-3 for 207 systems")
    expected_panel5 = {
        (row["system_id"], row["replica"])
        for row in ledger
        if is_true(row["included_in_figure5"])
    }
    if panel5_keys != expected_panel5 or len(panel5) != 618:
        findings.append("Figure 5 plotted subset does not match ledger inclusion flags")
    if excluded5_keys != ledger_keys - expected_panel5 or len(excluded5) != 3:
        findings.append("Figure 5 exclusion audit does not match ledger exclusions")
    if any(not row["figure5_exclusion_reason"] for row in excluded5):
        findings.append("Figure 5 exclusion reason is missing")
    family_counts = Counter(row["g_protein_family"] for row in panel5)
    if family_counts != Counter({"Gi": 285, "Gs": 192, "Gq": 123, "G12-13": 18}):
        findings.append("Figure 5 family counts differ from 285/192/123/18")
    if len({row["system_id"] for row in panel5}) != 206:
        findings.append("Figure 5 plotted subset does not contain 206 systems")
    metric_names = [
        "tm_core_rmsd_A_p95",
        "galpha_interface_rmsd_A_p95",
        "contact_retention_p05",
    ]
    for row in panel5:
        if row["harmonized_observations"] not in {"1001", "1001.0"}:
            findings.append("Figure 5 plotted row has a non-1001 observation count")
            break
        values = [float(row[name]) for name in metric_names]
        if not all(math.isfinite(value) for value in values):
            findings.append("Figure 5 plotted row has a nonfinite metric")
            break
        if not 0 <= values[2] <= 1:
            findings.append("Figure 5 contact retention lies outside [0, 1]")
            break

    mapping = read_csv("Figure_6_pocket_record_mapping.csv")
    panel_b = read_csv("Figure_6_panel_B_class_counts.csv")
    panel_c = read_csv("Figure_6_panel_C_system_summary.csv")
    panel_d = read_csv("Figure_6_panel_D_position_counts.csv")
    if len(mapping) != 2151:
        findings.append("Figure 6 mapping does not contain 2,151 status rows")
    detected = [row for row in mapping if row["pocket_id"]]
    zero = [row for row in mapping if not row["pocket_id"]]
    detected_keys = {(row["system_id"], row["pocket_id"]) for row in detected}
    if len(detected) != 2149 or len(detected_keys) != 2149:
        findings.append("Figure 6 mapping does not contain 2,149 unique detected pockets")
    if len(zero) != 2 or {row["system_id"] for row in zero} != {"Gi_7YK6", "Gi_8YIC"}:
        findings.append("Figure 6 explicit zero-pocket records are incorrect")
    availability_by_system = {
        row["system_id"]: int(row["detected_pockets"]) for row in availability
    }
    mapping_counts = Counter(row["system_id"] for row in detected)
    if any(mapping_counts[system_id] != count for system_id, count in availability_by_system.items()):
        findings.append("Figure 3 pocket counts disagree with the Figure 6 record mapping")
    included6 = [row for row in mapping if is_true(row["included_in_figure6_panels_b_c"])]
    if len(included6) != 669 or len({row["system_id"] for row in included6}) != 201:
        findings.append("Figure 6 panels B-C do not map to 669 pockets from 201 systems")
    if any(
        not row["figure6_panels_b_c_exclusion_reason"]
        for row in mapping
        if not is_true(row["included_in_figure6_panels_b_c"])
    ):
        findings.append("Figure 6 panels B-C have an unexplained excluded record")

    expected_classes = Counter(row["anatomical_class"] for row in included6)
    exported_classes = {
        row["anatomical_class"]: int(row["pocket_records"]) for row in panel_b
    }
    if exported_classes != dict(expected_classes):
        findings.append("Figure 6 panel B class counts do not match the record mapping")
    if exported_classes != {
        "orthosteric": 82,
        "extracellular_vestibule": 135,
        "tm_core_allosteric": 287,
        "intracellular_allosteric": 165,
    }:
        findings.append("Figure 6 panel B class counts differ from 82/135/287/165")

    by_system: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in included6:
        by_system[row["system_id"]].append(row)
    panel_c_by_system = {row["system_id"]: row for row in panel_c}
    if set(panel_c_by_system) != set(by_system):
        findings.append("Figure 6 panel C system set differs from the included mapping")
    else:
        for system_id, rows in by_system.items():
            summary = panel_c_by_system[system_id]
            expected_mean = sum(float(row["mean_freq"]) for row in rows) / len(rows)
            if int(summary["gpcr_centered_pockets"]) != len(rows) or not math.isclose(
                float(summary["mean_gpcr_pocket_occupancy"]),
                expected_mean,
                rel_tol=1e-12,
                abs_tol=1e-12,
            ):
                findings.append(f"Figure 6 panel C mismatch for {system_id}")
                break

    position_counts: Counter[str] = Counter()
    for row in mapping:
        if row["anatomical_class"] != "orthosteric":
            continue
        positions = row["receptor_generic_numbers"]
        if positions:
            position_counts.update(positions.split(";"))
    exported_positions = {
        row["gpcrdb_generic_position"]: int(row["orthosteric_pocket_records"])
        for row in panel_d
    }
    if len(panel_d) != 15:
        findings.append("Figure 6 panel D does not contain 15 positions")
    if any(position_counts[position] != count for position, count in exported_positions.items()):
        findings.append("Figure 6 panel D position count differs from the mapping")
    omitted = [count for position, count in position_counts.items() if position not in exported_positions]
    if omitted and min(exported_positions.values()) < max(omitted):
        findings.append("Figure 6 panel D omits a position with a higher count")

    result = {
        "status": "passed" if not findings else "failed",
        "findings": findings,
        "figures_1_4_and_supplementary": {
            "cohort_rows": len(composition),
            "positive_control_rows": len(positive),
            "gateway_rows": len(gateway),
            "receptor_profiles": len(profiles),
            "gateway_method_rows": len(gateway_method),
        },
        "figure5": {
            "ledger_rows": len(ledger),
            "plotted_rows": len(panel5),
            "excluded_rows": len(excluded5),
        },
        "figure6": {
            "status_rows": len(mapping),
            "detected_pockets": len(detected),
            "panels_b_c_pockets": len(included6),
            "panels_b_c_systems": len(by_system),
            "panel_d_positions": len(panel_d),
        },
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if findings:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
