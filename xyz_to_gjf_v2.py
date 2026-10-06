#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
from pathlib import Path

ATOMIC_SYMBOLS = {
    1:"H", 2:"He", 3:"Li", 4:"Be", 5:"B", 6:"C", 7:"N", 8:"O", 9:"F", 10:"Ne",
    11:"Na", 12:"Mg", 13:"Al", 14:"Si", 15:"P", 16:"S", 17:"Cl", 18:"Ar",
    19:"K", 20:"Ca", 21:"Sc", 22:"Ti", 23:"V", 24:"Cr", 25:"Mn", 26:"Fe",
    27:"Co", 28:"Ni", 29:"Cu", 30:"Zn", 31:"Ga", 32:"Ge", 33:"As", 34:"Se",
    35:"Br", 36:"Kr", 37:"Rb", 38:"Sr", 39:"Y", 40:"Zr", 41:"Nb", 42:"Mo",
    43:"Tc", 44:"Ru", 45:"Rh", 46:"Pd", 47:"Ag", 48:"Cd", 49:"In", 50:"Sn",
    51:"Sb", 52:"Te", 53:"I", 54:"Xe", 55:"Cs", 56:"Ba", 57:"La", 58:"Ce",
    59:"Pr", 60:"Nd", 61:"Pm", 62:"Sm", 63:"Eu", 64:"Gd", 65:"Tb", 66:"Dy",
    67:"Ho", 68:"Er", 69:"Tm", 70:"Yb", 71:"Lu", 72:"Hf", 73:"Ta", 74:"W",
    75:"Re", 76:"Os", 77:"Ir", 78:"Pt", 79:"Au", 80:"Hg", 81:"Tl", 82:"Pb",
    83:"Bi", 84:"Po", 85:"At", 86:"Rn"
}

DEFAULT_ROUTE = "opt b3lyp/6-31g(d)"
DEFAULT_MEM = "32GB"
DEFAULT_NPROC = 32
DEFAULT_CHARGE = 0
DEFAULT_MULTIPLICITY = 1


def sanitize_name(text):
    text = text.strip()
    if not text:
        return ""
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_.-]", "_", text)
    text = re.sub(r"_+", "_", text).strip("_.")
    return text[:100]


def normalize_element(token):
    token = token.strip()
    if re.fullmatch(r"\d+", token):
        z = int(token)
        if z not in ATOMIC_SYMBOLS:
            raise ValueError(f"Unsupported atomic number: {z}")
        return ATOMIC_SYMBOLS[z]

    m = re.match(r"^([A-Za-z]{1,2})", token)
    if not m:
        raise ValueError(f"Cannot interpret element token: {token}")

    sym = m.group(1)
    return sym[0].upper() + sym[1:].lower()


def read_multi_xyz(path):
    path = Path(path)
    lines = path.read_text(errors="replace").splitlines()
    conformers = []
    i = 0
    block_number = 0

    while i < len(lines):
        while i < len(lines) and not lines[i].strip():
            i += 1

        if i >= len(lines):
            break

        try:
            natoms = int(lines[i].strip())
        except ValueError as exc:
            raise ValueError(
                f"Expected atom count at line {i+1}, but found: {lines[i]!r}"
            ) from exc

        if natoms <= 0:
            raise ValueError(f"Invalid atom count {natoms} at line {i+1}")

        if i + 1 >= len(lines):
            raise ValueError(f"Missing XYZ comment line after line {i+1}")

        comment = lines[i + 1].strip()
        start = i + 2
        end = start + natoms

        if end > len(lines):
            raise ValueError(
                f"Conformer starting at line {i+1} is incomplete: expected {natoms} atoms."
            )

        atoms = []
        for line_no, line in enumerate(lines[start:end], start=start + 1):
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"Invalid XYZ coordinate line {line_no}: {line!r}")

            element = normalize_element(parts[0])
            try:
                x, y, z = map(float, parts[1:4])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid coordinates at line {line_no}: {line!r}"
                ) from exc

            atoms.append((element, x, y, z))

        block_number += 1
        conformers.append({
            "index": block_number,
            "comment": comment,
            "atoms": atoms
        })

        i = end

    if not conformers:
        raise ValueError(f"No XYZ structures found in {path}")

    return conformers


def read_genecp_block(path):
    if path is None:
        return None

    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"GenECP file not found: {path}")

    text = path.read_text(errors="replace").strip("\n")
    if not text.strip():
        raise ValueError(f"GenECP file is empty: {path}")

    return text


def write_gjf(
    path, atoms, title, route, charge, multiplicity, mem, nproc,
    chk=True, genecp_block=None
):
    path = Path(path)
    lines = []

    if nproc:
        lines.append(f"%nprocshared={nproc}")
    if mem:
        lines.append(f"%mem={mem}")
    if chk:
        lines.append(f"%chk={path.stem}.chk")

    route = route.strip()
    if route.startswith("#"):
        lines.append(route)
    else:
        lines.append(f"#p {route}")

    lines.extend(["", title, "", f"{charge} {multiplicity}"])

    for element, x, y, z in atoms:
        lines.append(f"{element:<3s} {x:16.8f} {y:16.8f} {z:16.8f}")

    # Gaussian requires a blank line between Cartesian coordinates
    # and the Gen/GenECP basis/ECP section.
    lines.append("")

    if genecp_block:
        lines.extend(genecp_block.rstrip().splitlines())
        lines.append("")

    lines.append("")
    content = "\n".join(lines).rstrip()
    path.write_text(content + "\n\n\n\n")


def xyz_to_gjf(
    xyz_path,
    output_dir=None,
    prefix=None,
    route=DEFAULT_ROUTE,
    charge=DEFAULT_CHARGE,
    multiplicity=DEFAULT_MULTIPLICITY,
    mem=DEFAULT_MEM,
    nproc=DEFAULT_NPROC,
    use_comment_names=False,
    genecp_file=None
):
    xyz_path = Path(xyz_path).resolve()
    conformers = read_multi_xyz(xyz_path)
    genecp_block = read_genecp_block(genecp_file)

    if genecp_block and "genecp" not in route.lower():
        print(
            "Warning: a GenECP block was supplied, but the route section "
            "does not contain 'GenECP'. Make sure the Gaussian route is correct.",
            file=sys.stderr
        )

    if output_dir is None:
        output_dir = xyz_path.parent / f"{xyz_path.stem}_gjf"
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if prefix is None or not prefix.strip():
        prefix = xyz_path.stem
    prefix = sanitize_name(prefix) or "conf"

    digits = max(3, len(str(len(conformers))))
    written = []
    used_names = set()

    for conf in conformers:
        comment_name = sanitize_name(conf["comment"]) if use_comment_names else ""

        if comment_name:
            stem = f"{prefix}_{comment_name}"
        else:
            stem = f"{prefix}_{conf['index']:0{digits}d}"

        base_stem = stem
        counter = 2
        while stem.lower() in used_names or (output_dir / f"{stem}.gjf").exists():
            stem = f"{base_stem}_{counter}"
            counter += 1

        used_names.add(stem.lower())
        out = output_dir / f"{stem}.gjf"

        title = conf["comment"] if conf["comment"] else stem
        write_gjf(
            out,
            conf["atoms"],
            title=title,
            route=route,
            charge=charge,
            multiplicity=multiplicity,
            mem=mem,
            nproc=nproc,
            chk=True,
            genecp_block=genecp_block
        )
        written.append(out)

    return written


def is_charge_mult_line(line):
    parts = line.split()
    if len(parts) < 2:
        return False
    return bool(re.fullmatch(r"[+-]?\d+", parts[0]) and re.fullmatch(r"[+-]?\d+", parts[1]))


def parse_gjf_coordinates(path):
    path = Path(path)
    lines = path.read_text(errors="replace").splitlines()

    route_start = None
    for i, line in enumerate(lines):
        if line.lstrip().startswith("#"):
            route_start = i
            break

    if route_start is None:
        raise ValueError(f"No Gaussian route section found in {path}")

    i = route_start
    while i < len(lines) and lines[i].strip():
        i += 1

    while i < len(lines) and not lines[i].strip():
        i += 1

    title_lines = []
    while i < len(lines) and lines[i].strip():
        title_lines.append(lines[i].strip())
        i += 1

    title = " ".join(title_lines).strip() or path.stem

    while i < len(lines) and not lines[i].strip():
        i += 1

    charge_mult_index = None
    for j in range(i, min(i + 20, len(lines))):
        if is_charge_mult_line(lines[j]):
            charge_mult_index = j
            break

    if charge_mult_index is None:
        raise ValueError(f"Could not find charge/multiplicity line in {path}")

    parts = lines[charge_mult_index].split()
    charge = int(parts[0])
    multiplicity = int(parts[1])

    atoms = []
    j = charge_mult_index + 1

    while j < len(lines):
        line = lines[j].strip()
        if not line:
            break

        parts = line.split()
        if len(parts) < 4:
            break

        try:
            element = normalize_element(parts[0])
        except ValueError:
            break

        xyz = None
        try:
            xyz = tuple(map(float, parts[1:4]))
        except ValueError:
            pass

        if xyz is None and len(parts) >= 5:
            try:
                xyz = tuple(map(float, parts[2:5]))
            except ValueError:
                pass

        if xyz is None:
            break

        atoms.append((element, *xyz))
        j += 1

    if not atoms:
        raise ValueError(
            f"No Cartesian coordinates found in {path}. "
            "This script expects a standard Cartesian Gaussian input."
        )

    return {
        "title": title,
        "charge": charge,
        "multiplicity": multiplicity,
        "atoms": atoms
    }


def write_xyz(path, atoms, comment):
    path = Path(path)
    lines = [str(len(atoms)), comment]
    for element, x, y, z in atoms:
        lines.append(f"{element:<3s} {x:16.8f} {y:16.8f} {z:16.8f}")
    lines.append("")
    path.write_text("\n".join(lines))


def collect_gjf_files(target, recursive=False):
    target = Path(target).resolve()

    if target.is_file():
        if target.suffix.lower() not in {".gjf", ".com"}:
            raise ValueError(f"Not a Gaussian input file: {target}")
        return [target]

    if not target.is_dir():
        raise ValueError(f"Path does not exist: {target}")

    patterns = ("*.gjf", "*.GJF", "*.com", "*.COM")
    files = []

    if recursive:
        for pattern in patterns:
            files.extend(target.rglob(pattern))
    else:
        for pattern in patterns:
            files.extend(target.glob(pattern))

    return sorted(set(p.resolve() for p in files), key=lambda p: str(p).lower())


def gjf_to_xyz(target, output_dir=None, recursive=False, combine=None):
    files = collect_gjf_files(target, recursive=recursive)

    if not files:
        raise ValueError("No .gjf or .com files found.")

    target = Path(target).resolve()

    if combine:
        combine = Path(combine).resolve()
        combine.parent.mkdir(parents=True, exist_ok=True)

        with combine.open("w") as out:
            for gjf in files:
                data = parse_gjf_coordinates(gjf)
                out.write(f"{len(data['atoms'])}\n")
                out.write(
                    f"{gjf.stem} | charge={data['charge']} "
                    f"multiplicity={data['multiplicity']} | {data['title']}\n"
                )
                for element, x, y, z in data["atoms"]:
                    out.write(f"{element:<3s} {x:16.8f} {y:16.8f} {z:16.8f}\n")
        return files, combine

    if output_dir is None:
        if target.is_file():
            output_dir = target.parent
        else:
            output_dir = target / "xyz"

    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    written = []
    for gjf in files:
        data = parse_gjf_coordinates(gjf)

        if recursive and target.is_dir():
            try:
                rel_parent = gjf.parent.relative_to(target)
            except ValueError:
                rel_parent = Path()

            out_dir = output_dir / rel_parent
            out_dir.mkdir(parents=True, exist_ok=True)
        else:
            out_dir = output_dir

        out = out_dir / f"{gjf.stem}.xyz"
        comment = (
            f"{gjf.stem} | charge={data['charge']} "
            f"multiplicity={data['multiplicity']} | {data['title']}"
        )
        write_xyz(out, data["atoms"], comment)
        written.append(out)

    return written, None


def ask(prompt, default=None):
    if default is None:
        return input(prompt + ": ").strip()
    value = input(f"{prompt} [{default}]: ").strip()
    return value if value else str(default)


def ask_yes_no(prompt, default=False):
    suffix = "[y/N]" if not default else "[Y/n]"
    value = input(f"{prompt} {suffix}: ").strip().lower()
    if not value:
        return default
    return value in {"y", "yes", "1", "true"}


def interactive():
    print("=== Gaussian / XYZ Conformer Converter ===")
    print()
    print("1. Multi-conformer XYZ -> multiple GJF files")
    print("2. GJF files in current directory -> individual XYZ files")
    print("3. GJF files in a directory -> individual XYZ files")
    print("4. GJF files in a directory -> one combined multi-XYZ file")
    print("5. Single GJF -> XYZ")
    print("6. Exit")
    print()

    choice = ask("Select", "1")

    if choice == "1":
        xyz_path = ask("Multi-conformer XYZ file")
        output_dir = ask("Output directory", f"{Path(xyz_path).stem}_gjf")
        prefix = ask("Output filename prefix", Path(xyz_path).stem)
        charge = int(ask("Charge", DEFAULT_CHARGE))
        multiplicity = int(ask("Multiplicity", DEFAULT_MULTIPLICITY))

        use_genecp = ask_yes_no("Append a GenECP basis/ECP block?", False)
        genecp_file = None
        if use_genecp:
            genecp_file = ask("GenECP block file")
            route_default = "opt b3lyp/genecp"
        else:
            route_default = DEFAULT_ROUTE

        route = ask("Gaussian route (without #p)", route_default)
        mem = ask("Memory", DEFAULT_MEM)
        nproc = int(ask("Number of cores", DEFAULT_NPROC))
        use_comment_names = ask_yes_no("Use XYZ comment lines in output filenames?", False)

        written = xyz_to_gjf(
            xyz_path,
            output_dir=output_dir,
            prefix=prefix,
            route=route,
            charge=charge,
            multiplicity=multiplicity,
            mem=mem,
            nproc=nproc,
            use_comment_names=use_comment_names,
            genecp_file=genecp_file
        )

        print(f"\nCreated {len(written)} GJF files in:")
        print(Path(output_dir).resolve())

    elif choice in {"2", "3", "4"}:
        if choice == "2":
            target = "."
        else:
            target = ask("Directory", ".")

        recursive = ask_yes_no("Search subdirectories recursively?", False)

        if choice == "4":
            default_name = str(Path(target).resolve() / "combined.xyz")
            combine = ask("Combined XYZ output file", default_name)
            files, combined = gjf_to_xyz(target, recursive=recursive, combine=combine)
            print(f"\nConverted {len(files)} Gaussian input files.")
            print(f"Combined XYZ: {combined}")
        else:
            default_out = str(Path(target).resolve() / "xyz")
            output_dir = ask("Output directory", default_out)
            written, _ = gjf_to_xyz(target, output_dir=output_dir, recursive=recursive)
            print(f"\nCreated {len(written)} XYZ files in:")
            print(Path(output_dir).resolve())

    elif choice == "5":
        gjf = ask("Gaussian input file")
        written, _ = gjf_to_xyz(gjf)
        print(f"\nCreated: {written[0]}")

    elif choice == "6":
        print("Exit.")
        return

    else:
        raise SystemExit(f"Invalid selection: {choice}")


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Convert concatenated multi-conformer XYZ files to Gaussian GJF files, "
            "and convert Gaussian GJF/COM files to XYZ."
        )
    )
    sub = parser.add_subparsers(dest="command")

    p1 = sub.add_parser("xyz2gjf", help="Split a multi-XYZ file into multiple GJF files")
    p1.add_argument("xyz", help="Input multi-conformer XYZ file")
    p1.add_argument("-o", "--output-dir", help="Output directory")
    p1.add_argument("-p", "--prefix", help="Output filename prefix")
    p1.add_argument("--route", default=DEFAULT_ROUTE, help="Gaussian route section without #p")
    p1.add_argument("--charge", type=int, default=DEFAULT_CHARGE)
    p1.add_argument("--multiplicity", type=int, default=DEFAULT_MULTIPLICITY)
    p1.add_argument("--mem", default=DEFAULT_MEM)
    p1.add_argument("--nproc", type=int, default=DEFAULT_NPROC)
    p1.add_argument(
        "--genecp",
        metavar="FILE",
        help=(
            "Text file containing the complete Gaussian GenECP basis/ECP block. "
            "The file is appended verbatim after the Cartesian coordinates of every GJF."
        )
    )
    p1.add_argument(
        "--comment-names",
        action="store_true",
        help="Use sanitized XYZ comment lines in output filenames"
    )

    p2 = sub.add_parser("gjf2xyz", help="Convert GJF/COM file(s) to XYZ")
    p2.add_argument("target", help="A GJF/COM file or a directory")
    p2.add_argument("-o", "--output-dir", help="Directory for individual XYZ files")
    p2.add_argument("-r", "--recursive", action="store_true")
    p2.add_argument(
        "--combine",
        help="Write all structures into one concatenated multi-XYZ file"
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    try:
        if args.command == "xyz2gjf":
            written = xyz_to_gjf(
                args.xyz,
                output_dir=args.output_dir,
                prefix=args.prefix,
                route=args.route,
                charge=args.charge,
                multiplicity=args.multiplicity,
                mem=args.mem,
                nproc=args.nproc,
                use_comment_names=args.comment_names,
                genecp_file=args.genecp
            )
            print(f"Created {len(written)} GJF files.")
            if written:
                print(f"Output directory: {written[0].parent}")

        elif args.command == "gjf2xyz":
            written, combined = gjf_to_xyz(
                args.target,
                output_dir=args.output_dir,
                recursive=args.recursive,
                combine=args.combine
            )
            if combined:
                print(f"Converted {len(written)} GJF/COM files.")
                print(f"Combined XYZ: {combined}")
            else:
                print(f"Created {len(written)} XYZ files.")

        else:
            interactive()

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()









