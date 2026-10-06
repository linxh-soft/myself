#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from pathlib import Path


# ============================================================
# XYZ -> GJF 时使用的 Gaussian header
# 这里你自己改
#
# {name} 会自动替换成文件名
# ============================================================

GJF_HEADER = """%chk={name}.chk
%nprocshared=16
%mem=32GB
#p YOUR_ROUTE_SECTION_HERE

{name}

"""


# ============================================================
# Periodic table
# 用于 gjf 中坐标第一列如果是原子序数的情况
# ============================================================

ELEMENTS = [
    "",
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
]


# ============================================================
# XYZ reader
# ============================================================

def read_xyz(filename):

    with open(filename, "r") as f:
        lines = f.readlines()

    if len(lines) < 2:
        raise ValueError(f"{filename}: invalid XYZ file")

    try:
        natoms = int(lines[0].strip())
    except ValueError:
        raise ValueError(
            f"{filename}: first line must be atom number"
        )

    coords = []

    for line in lines[2:2 + natoms]:

        p = line.split()

        if len(p) < 4:
            raise ValueError(
                f"{filename}: bad coordinate line:\n{line}"
            )

        element = p[0]
        x = float(p[1])
        y = float(p[2])
        z = float(p[3])

        coords.append(
            (element, x, y, z)
        )

    if len(coords) != natoms:
        raise ValueError(
            f"{filename}: expected {natoms} atoms, "
            f"but found {len(coords)}"
        )

    return coords


# ============================================================
# XYZ -> GJF
# ============================================================

def xyz_to_gjf(filename, charge, mult):

    path = Path(filename)

    coords = read_xyz(path)

    name = path.stem

    outfile = path.with_suffix(".gjf")

    header = GJF_HEADER.format(
        name=name
    )

    with open(outfile, "w") as f:

        f.write(header)

        # charge multiplicity
        f.write(f"{charge} {mult}\n")

        # coordinates
        for element, x, y, z in coords:

            f.write(
                f"{element:<3s}"
                f"{x:16.8f}"
                f"{y:16.8f}"
                f"{z:16.8f}\n"
            )

        # Gaussian 文件末尾留空行
        f.write("\n\n\n")

    print(
        f"{path.name}  ->  {outfile.name}"
    )


# ============================================================
# Gaussian coordinate parser
# ============================================================

def normalize_element(token):

    # 比如 "C" / "Fe"
    if re.fullmatch(
        r"[A-Za-z]{1,3}",
        token
    ):

        return (
            token[0].upper()
            +
            token[1:].lower()
        )

    # 比如原子序数 "6"
    if re.fullmatch(
        r"\d+",
        token
    ):

        z = int(token)

        if 1 <= z < len(ELEMENTS):
            return ELEMENTS[z]

    raise ValueError(
        f"Cannot recognize element: {token}"
    )


def read_gjf_geometry(filename):

    with open(filename, "r") as f:
        lines = f.readlines()

    # Gaussian charge/multiplicity:
    # 0 1
    # -1 2
    charge_mult_pattern = re.compile(
        r"^\s*[+-]?\d+\s+\d+\s*$"
    )

    cm_index = None

    for i, line in enumerate(lines):

        if charge_mult_pattern.match(line):

            cm_index = i
            break

    if cm_index is None:

        raise ValueError(
            f"{filename}: cannot find charge/multiplicity line"
        )

    coords = []

    for line in lines[cm_index + 1:]:

        # 第一空行结束 geometry
        if not line.strip():

            if coords:
                break
            else:
                continue

        p = line.split()

        if len(p) < 4:

            if coords:
                break

            continue

        try:
            element = normalize_element(
                p[0]
            )

            x = float(p[1])
            y = float(p[2])
            z = float(p[3])

        except (ValueError, IndexError):

            if coords:
                break

            continue

        coords.append(
            (element, x, y, z)
        )

    if not coords:

        raise ValueError(
            f"{filename}: no Cartesian coordinates found"
        )

    return coords


# ============================================================
# GJF -> XYZ
# ============================================================

def gjf_to_xyz(filename):

    path = Path(filename)

    coords = read_gjf_geometry(path)

    outfile = path.with_suffix(".xyz")

    with open(outfile, "w") as f:

        f.write(
            f"{len(coords)}\n"
        )

        f.write(
            f"Converted from {path.name}\n"
        )

        for element, x, y, z in coords:

            f.write(
                f"{element:<3s}"
                f"{x:16.8f}"
                f"{y:16.8f}"
                f"{z:16.8f}\n"
            )

    print(
        f"{path.name}  ->  {outfile.name}"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description="Batch XYZ <-> Gaussian GJF converter"
    )

    subparsers = parser.add_subparsers(
        dest="mode",
        required=True
    )

    # xyz2gjf
    p_xyz = subparsers.add_parser(
        "xyz2gjf"
    )

    p_xyz.add_argument(
        "files",
        nargs="+",
        help="XYZ files"
    )

    p_xyz.add_argument(
        "-c",
        "--charge",
        type=int,
        default=0,
        help="Charge, default = 0"
    )

    p_xyz.add_argument(
        "-m",
        "--mult",
        type=int,
        default=1,
        help="Multiplicity, default = 1"
    )

    # gjf2xyz
    p_gjf = subparsers.add_parser(
        "gjf2xyz"
    )

    p_gjf.add_argument(
        "files",
        nargs="+",
        help="GJF files"
    )

    args = parser.parse_args()

    if args.mode == "xyz2gjf":

        for filename in args.files:

            try:
                xyz_to_gjf(
                    filename,
                    args.charge,
                    args.mult
                )

            except Exception as e:

                print(
                    f"ERROR: {filename}: {e}"
                )

    elif args.mode == "gjf2xyz":

        for filename in args.files:

            try:
                gjf_to_xyz(
                    filename
                )

            except Exception as e:

                print(
                    f"ERROR: {filename}: {e}"
                )


if __name__ == "__main__":
    main()




