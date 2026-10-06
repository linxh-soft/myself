#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import shutil


def get_basis(filename):

    name = filename.lower()

    # gen2: SDD / 6-31g*
    if "gen2" in name:
        return """Fe 0
SDD
****
C Cl H N 0
6-31g*
****

Fe 0
SDD
"""

    # gen1: LANL2DZ / 6-31g*
    elif "gen1" in name:
        return """Fe 0
lanl2dz
****
C Cl H N 0
6-31g*
****

Fe 0
lanl2dz
"""

    # gen without number: def2-TZVP / def2-SVP
    elif "gen" in name:
        return """Fe 0
def2tzvp
****
C Cl H N 0
def2svp
****

Fe 0
def2tzvp
"""

    else:
        return ""


def find_charge_mult(lines):

    for i, line in enumerate(lines):
        parts = line.split()

        if len(parts) == 2:
            try:
                int(parts[0])
                int(parts[1])
                return i
            except:
                continue

    return None


for gjf in Path(".").glob("*.gjf"):

    print("Processing:", gjf)

    # backup
    shutil.copy(
        gjf,
        str(gjf) + ".bak"
    )

    with open(gjf, "r", encoding="utf-8") as f:
        lines = f.readlines()


    # remove old chk
    if lines and lines[0].startswith("%chk"):
        lines.pop(0)


    # locate charge/multiplicity
    cm = find_charge_mult(lines)

    if cm is None:
        print("Cannot find charge/multiplicity:", gjf)
        continue


    # keep header + charge/multiplicity
    output = lines[:cm+1]


    coord_finished = False

    for line in lines[cm+1:]:

        output.append(line)

        # first blank line after coordinates
        if line.strip() == "":
            coord_finished = True
            break


    if not coord_finished:
        print("Cannot locate coordinate end:", gjf)
        continue


    # remove extra blank lines
    while output and output[-1].strip() == "":
        output.pop()


    # add chk
    output.insert(
        0,
        f"%chk={gjf.stem}.chk\n"
    )


    # add basis
    basis = get_basis(gjf.name)

    if basis:

        output.append("\n")
        output.append(basis)


    # remove trailing blanks
    while output and output[-1].strip() == "":
        output.pop()


    # three blank lines at end
    output.extend([
        "\n",
        "\n",
        "\n"
    ])


    with open(gjf, "w", encoding="utf-8") as f:
        f.writelines(output)


print("Done.")



