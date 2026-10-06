#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re
import csv


terms = [
    "Bond Energy",
    "Orbital Energy",
    "Electrostatic Energy",
    "Pauli Energy",
    "Delta E^0(XC)",
    "Delta Dispersion"
]


# 用于匹配数字
num = r"[-+]?\d*\.\d+"


results = []

for outfile in Path(".").rglob("*.out"):

    text = outfile.read_text(errors="ignore")

    data = {
        "File": str(outfile)
    }

    for term in terms:
        pattern = (
            re.escape(term)
            + r"\s+("
            + num
            + r")\s+("
            + num
            + r")"
        )

        match = re.search(pattern, text)

        if match:
            data[term + "_Eh"] = float(match.group(1))
            data[term + "_kcal"] = float(match.group(2))
        else:
            data[term + "_Eh"] = None
            data[term + "_kcal"] = None

    results.append(data)


# 输出CSV
output = "EDA_energy_table.csv"

headers = ["File"]

for term in terms:
    headers.append(term + "_Eh")
    headers.append(term + "_kcal")


with open(output, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=headers)
    writer.writeheader()
    writer.writerows(results)


print(f"Done! Extracted {len(results)} files.")
print(f"Saved as: {output}")


