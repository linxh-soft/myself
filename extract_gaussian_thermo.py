#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract Gaussian thermochemistry / frequency / spin information from .log/.out files.

Output columns:
File, ZPE, TCH, TCG, TS_term, H-Gas, G-Gas, SP,
First Frequency, Second Frequency, S**2, S, ddG, ddE, Normal termination

Definitions:
  ZPE     = Zero-point correction
  TCH     = Thermal correction to Enthalpy
  TCG     = Thermal correction to Gibbs Free Energy
  TS_term = TCH - TCG = T*S
  H-Gas   = Sum of electronic and thermal Enthalpies
  G-Gas   = Sum of electronic and thermal Free Energies
  SP      = last "SCF Done" electronic energy
  S**2    = last <S^2> before annihilation
  S       = effective spin quantum number from S(S+1)=<S^2>
  ddG     = relative G-Gas within an R/S filename pair, kcal/mol
  ddE     = relative SP within an R/S filename pair, kcal/mol

Usage:
  python extract_gaussian_thermo.py
  python extract_gaussian_thermo.py -r
  python extract_gaussian_thermo.py -o thermo.csv
  python extract_gaussian_thermo.py --no-pair

The script recognizes both .log and .out.
"""

import argparse
import csv
import math
import re
from pathlib import Path

HARTREE_TO_KCAL = 627.509474

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"

def as_float(s):
    if s is None:
        return None
    return float(s.replace("D", "E").replace("d", "e"))

def last_match_float(pattern, text, flags=0, group=1):
    matches = list(re.finditer(pattern, text, flags))
    if not matches:
        return None
    return as_float(matches[-1].group(group))

def get_last_frequency_block(text):
    """
    Return frequencies from the last Gaussian harmonic-frequency section.
    If a thermochemistry/ZPE section exists, isolate the harmonic-frequency
    section immediately preceding the last ZPE block.
    """
    zpe_pos = text.rfind("Zero-point correction=")
    end = zpe_pos if zpe_pos >= 0 else len(text)
    prefix = text[:end]

    # Gaussian usually prints this header before the frequency list
    harm_pos = prefix.rfind("Harmonic frequencies")
    if harm_pos >= 0:
        region = text[harm_pos:end]
    else:
        # Fallback: use all frequency lines before the last ZPE
        region = prefix

    freqs = []
    for m in re.finditer(r"Frequencies\s+--\s+([^\r\n]+)", region):
        for tok in re.findall(NUM, m.group(1)):
            try:
                freqs.append(as_float(tok))
            except Exception:
                pass
    return freqs

def strip_jobid(stem):
    # e.g. con_R.7876459.log -> stem may be con_R.7876459
    return re.sub(r"\.\d+$", "", stem)

def rs_group_key(filename):
    """
    Convert a filename into an R/S pairing key.
    Examples:
      con_R -> con_X
      unsy-Concerted_4TS4-R_Sty-FCOnly -> unsy-Concerted_4TS4-X_Sty-FCOnly
    """
    stem = strip_jobid(Path(filename).stem)
    # Replace a standalone R or S token bounded by _ or - or string boundaries
    key = re.sub(
        r"(^|[_-])([RS])(?=([_-]|$))",
        lambda m: m.group(1) + "X",
        stem,
        count=1,
    )
    return key if key != stem else None

def parse_gaussian(path):
    text = path.read_text(errors="ignore")

    zpe = last_match_float(
        rf"Zero-point correction=\s*({NUM})",
        text
    )
    tch = last_match_float(
        rf"Thermal correction to Enthalpy=\s*({NUM})",
        text
    )
    tcg = last_match_float(
        rf"Thermal correction to Gibbs Free Energy=\s*({NUM})",
        text
    )
    hgas = last_match_float(
        rf"Sum of electronic and thermal Enthalpies=\s*({NUM})",
        text
    )
    ggas = last_match_float(
        rf"Sum of electronic and thermal Free Energies=\s*({NUM})",
        text
    )

    # Last SCF electronic energy
    sp = last_match_float(
        rf"SCF Done:\s+E\([^)]+\)\s*=\s*({NUM})",
        text
    )

    # Last <S^2> before/after annihilation
    s2_before = None
    s2_after = None
    spin_matches = list(re.finditer(
        rf"S\*\*2 before annihilation\s+({NUM})\s*,\s*after\s+({NUM})",
        text
    ))
    if spin_matches:
        s2_before = as_float(spin_matches[-1].group(1))
        s2_after = as_float(spin_matches[-1].group(2))
    else:
        # Fallback for occasional alternative formatting
        s2_before = last_match_float(
            rf"<S\*\*2>\s*=\s*({NUM})",
            text
        )

    s_eff = None
    if s2_before is not None and s2_before >= -0.25:
        s_eff = (-1.0 + math.sqrt(1.0 + 4.0 * s2_before)) / 2.0

    freqs = get_last_frequency_block(text)
    first_freq = freqs[0] if len(freqs) >= 1 else None
    second_freq = freqs[1] if len(freqs) >= 2 else None

    ts_term = None
    if tch is not None and tcg is not None:
        ts_term = tch - tcg

    return {
        "File": path.stem,
        "Path": str(path),
        "ZPE": zpe,
        "TCH": tch,
        "TCG": tcg,
        "TS_term": ts_term,
        "H-Gas": hgas,
        "G-Gas": ggas,
        "SP": sp,
        "First Frequency": first_freq,
        "Second Frequency": second_freq,
        "S**2": s2_before,
        "S": s_eff,
        "S**2_after": s2_after,
        "ddG": None,
        "ddE": None,
        "Normal termination": "Yes" if "Normal termination of Gaussian" in text else "No",
    }

def add_rs_relative_energies(rows):
    """
    For files that differ only by an R/S token, set the lower-energy member to 0.
    ddG uses G-Gas; ddE uses SP. Units: kcal/mol.
    """
    groups = {}
    for i, row in enumerate(rows):
        key = rs_group_key(row["File"])
        if key is not None:
            groups.setdefault(key, []).append(i)

    for key, idxs in groups.items():
        if len(idxs) < 2:
            continue

        gvals = [rows[i]["G-Gas"] for i in idxs if rows[i]["G-Gas"] is not None]
        evals = [rows[i]["SP"] for i in idxs if rows[i]["SP"] is not None]

        gmin = min(gvals) if len(gvals) >= 2 else None
        emin = min(evals) if len(evals) >= 2 else None

        for i in idxs:
            if gmin is not None and rows[i]["G-Gas"] is not None:
                rows[i]["ddG"] = (rows[i]["G-Gas"] - gmin) * HARTREE_TO_KCAL
            if emin is not None and rows[i]["SP"] is not None:
                rows[i]["ddE"] = (rows[i]["SP"] - emin) * HARTREE_TO_KCAL

def fmt(v, ndp=6):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.{ndp}f}"
    return str(v)

def main():
    ap = argparse.ArgumentParser(
        description="Extract Gaussian thermochemistry from .log/.out files."
    )
    ap.add_argument(
        "-r", "--recursive", action="store_true",
        help="Search subdirectories recursively"
    )
    ap.add_argument(
        "-o", "--output", default="gaussian_thermo_summary.csv",
        help="Output CSV filename (default: gaussian_thermo_summary.csv)"
    )
    ap.add_argument(
        "--no-pair", action="store_true",
        help="Do not calculate R/S-pair ddG and ddE"
    )
    args = ap.parse_args()

    root = Path(".")
    if args.recursive:
        files = sorted(
            [p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in {".log", ".out"}]
        )
    else:
        files = sorted(
            [p for p in root.iterdir() if p.is_file() and p.suffix.lower() in {".log", ".out"}]
        )

    if not files:
        print("No .log or .out files found.")
        return

    rows = []
    for f in files:
        try:
            row = parse_gaussian(f)
            rows.append(row)
            print(f"Parsed: {f}")
        except Exception as e:
            print(f"FAILED: {f} -> {e}")

    if not args.no_pair:
        add_rs_relative_energies(rows)

    fields = [
        "File",
        "ZPE",
        "TCH",
        "TCG",
        "TS_term",
        "H-Gas",
        "G-Gas",
        "SP",
        "First Frequency",
        "Second Frequency",
        "S**2",
        "S",
        "ddG",
        "ddE",
        "Normal termination",
        "Path",
    ]

    with open(args.output, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            # Hartree-like quantities
            for k in ["ZPE", "TCH", "TCG", "TS_term", "H-Gas", "G-Gas", "SP"]:
                out[k] = fmt(out.get(k), 9)
            # Frequency / spin
            for k in ["First Frequency", "Second Frequency", "S**2", "S"]:
                out[k] = fmt(out.get(k), 4)
            # Relative energies in kcal/mol
            for k in ["ddG", "ddE"]:
                out[k] = fmt(out.get(k), 2)
            writer.writerow({k: out.get(k, "") for k in fields})

    print()
    print(f"Done. Wrote {len(rows)} rows to: {args.output}")
    print("ddG/ddE are in kcal/mol and are paired automatically by R/S in filenames.")
    print("TS_term = TCH - TCG = T*S (Hartree).")

if __name__ == "__main__":
    main()
