#!/usr/bin/env python3
"""Reproduce ligand-contact and binding-pocket RMSD analyses from a PDB/XTC pair."""
from __future__ import annotations

import argparse
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
import pandas as pd


def kabsch_rmsd(mobile: np.ndarray, reference: np.ndarray) -> float:
    """Return RMSD after optimal rigid-body superposition."""
    x = mobile - mobile.mean(axis=0)
    y = reference - reference.mean(axis=0)
    u, _, vt = np.linalg.svd(x.T @ y)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    fitted = x @ rotation
    return float(np.sqrt(np.mean(np.sum((fitted - y) ** 2, axis=1))))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("topology", type=Path, help="Released PDB topology")
    parser.add_argument("trajectory", type=Path, help="Released XTC trajectory")
    parser.add_argument("--ligand-selection", required=True,
                        help='MDAnalysis selection, e.g. "resname S1P and not name H*"')
    parser.add_argument("--pocket-selection", required=True,
                        help='Pocket C-alpha selection, e.g. "protein and name CA and resid 62 88 91"')
    parser.add_argument("--protein-selection", default="protein and not name H*",
                        help="Atoms used for ligand-contact analysis")
    parser.add_argument("--contact-cutoff", type=float, default=4.0,
                        help="Heavy-atom contact cutoff in Å (default: 4.0)")
    parser.add_argument("--frame-interval-ps", type=float, default=200.0,
                        help="Released-record frame interval in ps (default: 200)")
    parser.add_argument("--output-dir", type=Path, default=Path("example_analysis_output"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    universe = mda.Universe(args.topology, args.trajectory)
    ligand = universe.select_atoms(args.ligand_selection)
    protein = universe.select_atoms(args.protein_selection)
    pocket = universe.select_atoms(args.pocket_selection)
    if not len(ligand):
        raise ValueError("ligand selection is empty")
    if not len(protein):
        raise ValueError("protein selection is empty")
    if len(pocket) < 3:
        raise ValueError("pocket selection must contain at least three atoms")

    universe.trajectory[0]
    pocket_reference = pocket.positions.copy()
    residues = protein.residues
    residue_atoms = [res.atoms.intersection(protein) for res in residues]
    contacts = np.zeros((len(residues), len(universe.trajectory)), dtype=np.uint8)
    pocket_rmsd = np.zeros(len(universe.trajectory), dtype=float)

    for frame_index, ts in enumerate(universe.trajectory):
        ligand_xyz = ligand.positions
        for residue_index, atoms in enumerate(residue_atoms):
            if len(atoms) and np.min(distance_array(atoms.positions, ligand_xyz,
                                                    box=ts.dimensions)) <= args.contact_cutoff:
                contacts[residue_index, frame_index] = 1
        pocket_rmsd[frame_index] = kabsch_rmsd(pocket.positions, pocket_reference)

    time_ns = np.arange(len(universe.trajectory)) * args.frame_interval_ps / 1000.0
    labels = [f"{res.resname}{res.resid}" for res in residues]
    contact_table = pd.DataFrame(contacts.T, columns=labels)
    contact_table.insert(0, "time_ns", time_ns)
    contact_table.to_csv(args.output_dir / "ligand_contact_matrix.csv", index=False)
    pd.DataFrame({
        "residue": labels,
        "contact_persistence": contacts.mean(axis=1),
    }).sort_values("contact_persistence", ascending=False).to_csv(
        args.output_dir / "ligand_contact_persistence.csv", index=False
    )
    pd.DataFrame({
        "frame": np.arange(len(universe.trajectory)),
        "time_ns": time_ns,
        "binding_pocket_ca_rmsd_A": pocket_rmsd,
    }).to_csv(args.output_dir / "binding_pocket_rmsd.csv", index=False)
    print(f"Wrote analysis tables to {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
