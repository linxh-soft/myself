#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path


def get_basis(filename):

    name = filename.lower()

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

    elif "gen" in name:
        return """Fe 0
def2tzvp
****
C Cl H N 0
def2svp
6-31g*
****

Fe 0
def2tzvp
"""

    else:
        return ""


for gjf in Path(".").glob("*.gjf"):

    with open(gjf, "r", encoding="utf-8") as f:
        lines = f.readlines()


    # 删除198-217行
    lines = [
        line for i, line in enumerate(lines, start=1)
        if not (198 <= i <= 217)
    ]


    # 加入chk
    chk = f"%chk={gjf.stem}.chk\n"

    if lines and lines[0].startswith("%chk"):
        lines[0] = chk
    else:
        lines.insert(0, chk)

    # 找 charge multiplicity
    coord_start = None

    for i,line in enumerate(lines):
        if len(line.split()) == 2:
            try:
                int(line.split()[0])
                int(line.split()[1])
                coord_start = i
                break
            except:
                pass

    # 保留到坐标结束
    new_lines = lines[:coord_start+1]

    for line in lines[coord_start+1:]:
        new_lines.append(line)
        if line.strip() == "":
            break


    # 添加gen/genecp内容
    basis = get_basis(gjf.name)

    if basis:
        lines.append("\n")
        lines.append(basis)


    # 删除尾部空行
    while lines and lines[-1].strip() == "":
        lines.pop()


    # 最后三个空行
    lines.append("\n")
    lines.append("\n")
    lines.append("\n")


    with open(gjf, "w", encoding="utf-8") as f:
        f.writelines(lines)


    print(f"Processed: {gjf}")

print("Done.")



