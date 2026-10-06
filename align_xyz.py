#!/usr/bin/env python3

import numpy as np
from pathlib import Path

REF_FILE = "sy_cat.xyz"
FIT_IDX = [0, 1, 2, 3, 4]   # Fe + N4


def read_xyz(filename):
    with open(filename) as f:
        lines = f.readlines()

    nat = int(lines[0])
    atoms = []
    coords = []

    for line in lines[2:2+nat]:
        p = line.split()
        atoms.append(p[0])
        coords.append([
            float(p[1]),
            float(p[2]),
            float(p[3])
        ])

    return atoms, np.array(coords)


def write_xyz(filename, atoms, coords, rmsd):
    with open(filename, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"Aligned to {REF_FILE}; FeN4 RMSD = {rmsd:.5f} A\n")

        for atom, xyz in zip(atoms, coords):
            f.write(
                f"{atom:<3s}"
                f"{xyz[0]:16.8f}"
                f"{xyz[1]:16.8f}"
                f"{xyz[2]:16.8f}\n"
            )


def kabsch(mobile_fit, ref_fit, mobile_all):

    mob_center = mobile_fit.mean(axis=0)
    ref_center = ref_fit.mean(axis=0)

    P = mobile_fit - mob_center
    Q = ref_fit - ref_center

    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)

    R = U @ Vt

    # 禁止镜像
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = U @ Vt

    aligned_all = (mobile_all - mob_center) @ R + ref_center

    aligned_fit = aligned_all[FIT_IDX]

    rmsd = np.sqrt(
        np.mean(
            np.sum((aligned_fit - ref_fit)**2, axis=1)
        )
    )

    return aligned_all, rmsd


ref_atoms, ref_coords = read_xyz(REF_FILE)
ref_fit = ref_coords[FIT_IDX]

for xyzfile in sorted(Path(".").glob("*.xyz")):

    if xyzfile.name == REF_FILE:
        continue

    if xyzfile.name.startswith("aligned_"):
        continue

    atoms, coords = read_xyz(xyzfile)

    mobile_fit = coords[FIT_IDX]

    aligned, rmsd = kabsch(
        mobile_fit,
        ref_fit,
        coords
    )

    outfile = "aligned_" + xyzfile.name
    write_xyz(outfile, atoms, aligned, rmsd)

    print(
        f"{xyzfile.name:30s} -> {outfile:35s}"
        f" FeN4 RMSD = {rmsd:.5f} A"
    )


