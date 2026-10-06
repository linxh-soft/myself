#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import itertools
import math
import os
import numpy as np


def read_xyz(filename):
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if not lines:
        raise ValueError(f"{filename} is empty.")

    try:
        n_atoms = int(lines[0].strip())
    except Exception:
        raise ValueError(f"{filename}: first line must be the number of atoms.")

    if len(lines) < n_atoms + 2:
        raise ValueError(f"{filename}: incomplete XYZ file.")

    atoms = []
    coords = []

    for i, line in enumerate(lines[2:2 + n_atoms], start=1):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"{filename}: bad XYZ line for atom {i}: {line.rstrip()}")
        atoms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return atoms, np.asarray(coords, dtype=float)


def kabsch_rmsd(P, Q):
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)

    if P.shape != Q.shape:
        raise ValueError("Coordinate arrays have different shapes.")

    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)

    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)

    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt

    P_fit = Pc @ R
    diff = P_fit - Qc
    return math.sqrt(np.mean(np.sum(diff * diff, axis=1)))


def distance(a, b):
    return float(np.linalg.norm(a - b))


def bond_angle(a, center, b):
    v1 = a - center
    v2 = b - center

    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < 1e-12 or n2 < 1e-12:
        return None

    c = float(np.dot(v1, v2) / (n1 * n2))
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


def plane_normal(a, center, b):
    v1 = a - center
    v2 = b - center

    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < 1e-12 or n2 < 1e-12:
        return None

    cross = np.cross(v1, v2)
    cross_norm = np.linalg.norm(cross)
    quality = cross_norm / (n1 * n2)

    if quality < 1e-6:
        return None

    return cross / cross_norm


def plane_plane_angle(coords, metal, pair1, pair2):
    i, j = pair1
    k, l = pair2

    n1 = plane_normal(coords[i], coords[metal], coords[j])
    n2 = plane_normal(coords[k], coords[metal], coords[l])

    if n1 is None or n2 is None:
        return None

    c = abs(float(np.dot(n1, n2)))
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


def fmt(x, ndigits=6):
    if x is None:
        return ""
    return f"{x:.{ndigits}f}"


def short_name(path):
    return os.path.basename(path)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare multiple XYZ structures. The FIRST XYZ is the reference. "
            "Outputs RMSD, metal-ligand distances, all ligand-metal-ligand angles, "
            "all disjoint plane-plane angles, and pairwise RMSD matrices."
        )
    )

    parser.add_argument(
        "xyz",
        nargs="+",
        help="XYZ files. The first file is used as the reference."
    )

    parser.add_argument(
        "--metal",
        type=int,
        default=None,
        help="1-based metal atom number. If omitted, a unique Fe atom is auto-detected."
    )

    parser.add_argument(
        "--coord",
        type=int,
        nargs="+",
        required=True,
        help="1-based coordinating atom numbers, e.g. --coord 2 3 4 5 66 73"
    )

    parser.add_argument(
        "--prefix",
        default="geometry_multi",
        help="Prefix for output CSV files (default: geometry_multi)"
    )

    args = parser.parse_args()

    if len(args.xyz) < 2:
        raise ValueError("Please provide at least two XYZ files.")

    structures = []
    for filename in args.xyz:
        atoms, coords = read_xyz(filename)
        structures.append({
            "file": filename,
            "name": short_name(filename),
            "atoms": atoms,
            "coords": coords,
        })

    ref_atoms = structures[0]["atoms"]
    n_atoms = len(ref_atoms)

    for s in structures[1:]:
        if len(s["atoms"]) != n_atoms:
            raise ValueError(
                f"{s['file']}: number of atoms differs from reference {structures[0]['file']}."
            )
        if s["atoms"] != ref_atoms:
            raise ValueError(
                f"{s['file']}: atom types/order differ from reference. "
                "Atom i must correspond to atom i in every XYZ."
            )

    if args.metal is None:
        fe_indices = [i for i, el in enumerate(ref_atoms) if el.upper() == "FE"]
        if len(fe_indices) == 1:
            metal = fe_indices[0]
        elif len(fe_indices) == 0:
            raise ValueError("No Fe atom found. Please specify --metal.")
        else:
            raise ValueError("More than one Fe atom found. Please specify --metal.")
    else:
        metal = args.metal - 1

    if not (0 <= metal < n_atoms):
        raise ValueError("--metal atom number is out of range.")

    coord = [x - 1 for x in args.coord]

    if len(set(coord)) != len(coord):
        raise ValueError("Duplicate atom numbers found in --coord.")

    for i in coord:
        if not (0 <= i < n_atoms):
            raise ValueError(f"Coordinating atom {i + 1} is out of range.")
        if i == metal:
            raise ValueError("Metal atom cannot also be listed in --coord.")

    metal_label = f"{ref_atoms[metal]}{metal + 1}"

    def atom_label(i):
        return f"{ref_atoms[i]}{i + 1}"

    names = [s["name"] for s in structures]
    ref_name = names[0]
    ref_coords = structures[0]["coords"]

    heavy_mask = np.array([el.upper() != "H" for el in ref_atoms], dtype=bool)
    sphere_indices = [metal] + coord
    ligand_pairs = list(itertools.combinations(coord, 2))

    plane_pair_combinations = []
    for pair1, pair2 in itertools.combinations(ligand_pairs, 2):
        if set(pair1).isdisjoint(pair2):
            plane_pair_combinations.append((pair1, pair2))

    metric_rows = []

    rmsd_all_values = []
    rmsd_heavy_values = []
    rmsd_sphere_values = []

    for s in structures:
        c = s["coords"]
        rmsd_all_values.append(kabsch_rmsd(c, ref_coords))
        rmsd_heavy_values.append(kabsch_rmsd(c[heavy_mask], ref_coords[heavy_mask]))
        rmsd_sphere_values.append(
            kabsch_rmsd(c[sphere_indices], ref_coords[sphere_indices])
        )

    metric_rows.append({
        "Category": "RMSD vs reference",
        "Measurement": "All atoms",
        "Unit": "Å",
        "Values": rmsd_all_values,
        "Note": f"Reference = {ref_name}; best-fit Kabsch RMSD",
    })

    metric_rows.append({
        "Category": "RMSD vs reference",
        "Measurement": "Heavy atoms",
        "Unit": "Å",
        "Values": rmsd_heavy_values,
        "Note": f"Reference = {ref_name}; hydrogens excluded",
    })

    metric_rows.append({
        "Category": "RMSD vs reference",
        "Measurement": f"{metal_label} + {len(coord)} coordinating atoms",
        "Unit": "Å",
        "Values": rmsd_sphere_values,
        "Note": f"Reference = {ref_name}; coordination-sphere RMSD",
    })

    for i in coord:
        values = [
            distance(s["coords"][metal], s["coords"][i])
            for s in structures
        ]
        metric_rows.append({
            "Category": "M-X distance",
            "Measurement": f"{metal_label}-{atom_label(i)}",
            "Unit": "Å",
            "Values": values,
            "Note": "",
        })

    angle_values_by_pair = {}
    for i, j in ligand_pairs:
        values = [
            bond_angle(s["coords"][i], s["coords"][metal], s["coords"][j])
            for s in structures
        ]
        angle_values_by_pair[(i, j)] = values

        metric_rows.append({
            "Category": "L-M-L angle",
            "Measurement": f"{atom_label(i)}-{metal_label}-{atom_label(j)}",
            "Unit": "deg",
            "Values": values,
            "Note": "",
        })

    for pair1, pair2 in plane_pair_combinations:
        i, j = pair1
        k, l = pair2

        values = [
            plane_plane_angle(s["coords"], metal, pair1, pair2)
            for s in structures
        ]

        notes = []
        for idx, s in enumerate(structures):
            for pair in (pair1, pair2):
                key = tuple(sorted(pair))
                a = angle_values_by_pair[key][idx]
                if a is not None and (a < 10.0 or a > 170.0):
                    notes.append(f"{s['name']}: near-linear defining angle")

        if any(v is None for v in values):
            notes.append("undefined plane present because defining atoms are collinear")

        plane1_name = f"{atom_label(i)}-{metal_label}-{atom_label(j)}"
        plane2_name = f"{atom_label(k)}-{metal_label}-{atom_label(l)}"

        metric_rows.append({
            "Category": "Plane-plane angle",
            "Measurement": f"({plane1_name}) vs ({plane2_name})",
            "Unit": "deg",
            "Values": values,
            "Note": "; ".join(sorted(set(notes))),
        })

    values_csv = f"{args.prefix}_values.csv"
    with open(values_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Category", "Measurement", "Unit"] + names + ["Note"])
        for row in metric_rows:
            writer.writerow(
                [row["Category"], row["Measurement"], row["Unit"]]
                + [fmt(v) for v in row["Values"]]
                + [row["Note"]]
            )

    delta_csv = f"{args.prefix}_deltas_vs_{os.path.splitext(ref_name)[0]}.csv"
    with open(delta_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["Category", "Measurement", "Unit"]
            + [f"{name} - {ref_name}" for name in names[1:]]
            + ["Note"]
        )

        for row in metric_rows:
            ref_value = row["Values"][0]
            deltas = []

            for v in row["Values"][1:]:
                if ref_value is None or v is None:
                    deltas.append("")
                else:
                    deltas.append(fmt(v - ref_value))

            writer.writerow(
                [row["Category"], row["Measurement"], row["Unit"]]
                + deltas
                + [row["Note"]]
            )

    pairwise_files = []

    selections = [
        ("all_atoms", lambda c: c),
        ("heavy_atoms", lambda c: c[heavy_mask]),
        ("coordination_sphere", lambda c: c[sphere_indices]),
    ]

    for tag, selector in selections:
        matrix_csv = f"{args.prefix}_pairwise_rmsd_{tag}.csv"
        pairwise_files.append(matrix_csv)

        matrix = []
        for s1 in structures:
            row = []
            P = selector(s1["coords"])
            for s2 in structures:
                Q = selector(s2["coords"])
                row.append(kabsch_rmsd(P, Q))
            matrix.append(row)

        with open(matrix_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow([""] + names)
            for name, row in zip(names, matrix):
                writer.writerow([name] + [fmt(v) for v in row])

    print()
    print("=" * 80)
    print("MULTI-STRUCTURE GEOMETRY COMPARISON")
    print("=" * 80)
    print(f"Reference          : {ref_name}")
    print(f"Number structures  : {len(structures)}")
    print(f"Metal              : {metal_label}")
    print("Coord atoms        : " + " ".join(atom_label(i) for i in coord))
    print(f"M-X distances      : {len(coord)}")
    print(f"L-M-L angles       : {len(ligand_pairs)}")
    print(f"Plane-plane angles : {len(plane_pair_combinations)}")
    print()

    width_name = max(len("Structure"), max(len(n) for n in names))
    print(
        f"{'Structure':<{width_name}}  "
        f"{'All RMSD':>10}  "
        f"{'Heavy RMSD':>11}  "
        f"{'Coord RMSD':>10}"
    )
    print("-" * (width_name + 37))

    for idx, name in enumerate(names):
        print(
            f"{name:<{width_name}}  "
            f"{rmsd_all_values[idx]:>10.4f}  "
            f"{rmsd_heavy_values[idx]:>11.4f}  "
            f"{rmsd_sphere_values[idx]:>10.4f}"
        )

    print()
    print("Files written:")
    print(f"  {values_csv}")
    print(f"  {delta_csv}")
    for filename in pairwise_files:
        print(f"  {filename}")


if __name__ == "__main__":
    main()






