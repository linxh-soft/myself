#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import shutil


finish_dir = Path("finished")
finish_dir.mkdir(exist_ok=True)


def get_gjf_name(output_file):
    """
    A.123.log -> A.gjf
    A.123.out -> A.gjf
    A.log     -> A.gjf
    """

    stem = output_file.stem   # 去掉 .log/.out

    # 如果有数字后缀，例如 A.123
    parts = stem.split(".")

    if len(parts) > 1 and parts[-1].isdigit():
        base = ".".join(parts[:-1])
    else:
        base = stem

    return Path(base + ".gjf")


def add_scf_xqc(gjf):

    if not gjf.exists():
        print("Missing gjf:", gjf)
        return

    with open(gjf, "r", encoding="utf-8") as f:
        lines = f.readlines()


    for i, line in enumerate(lines):

        if line.lstrip().startswith("#"):

            if "scf=xqc" not in line.lower():

                lines[i] = line.rstrip() + " scf=xqc\n"

                with open(gjf, "w", encoding="utf-8") as f:
                    f.writelines(lines)

                print("Added scf=xqc:", gjf)

            else:
                print("Already has scf=xqc:", gjf)

            return


    print("Cannot find route line:", gjf)



def process_file(output):

    with open(output, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()


    last20 = "".join(lines[-20:])


    gjf = get_gjf_name(output)


    # SCF失败
    if "Convergence failure" in last20:

        print("\nSCF failure:", output)

        add_scf_xqc(gjf)


    # 正常结束
    elif "Normal termination" in last20:

        print("\nNormal termination:", output)

        files = [
            output,
            gjf
        ]

        for file in files:

            if file.exists():

                shutil.move(
                    str(file),
                    str(finish_dir / file.name)
                )

                print("Moved:", file)



# 同时处理 log 和 out
for output in list(Path(".").glob("*.log")) + list(Path(".").glob("*.out")):

    process_file(output)


print("\nDone.")


