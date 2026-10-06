#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
substrate_mapping_align_v2.py

Purpose
-------
For multiple catalyst-aligned XYZ files:

1) Chemically map the styrene substrate by connectivity (not by absolute position).
2) Resolve graph-symmetry mappings by ACTUAL rigid-body superposition RMSD.
3) Map attached substrate hydrogens as well.
4) Write FULL XYZ files in which:
       - Fe / ligand / reactive O coordinates are UNCHANGED
       - substrate coordinates are UNCHANGED
       - BUT substrate atom ORDER is changed to match the reference structure
   => use these files for later C-O distance / XY / contact / regression work.

5) Also write substrate-only superposed XYZ files, in reference atom order:
       A. all-heavy best-fit superposition
       B. reactive-center superposition using the styrene
          terminal alkene C / internal alkene C / ipso C anchors
   => use these only for visual comparison of substrate conformation.

IMPORTANT
---------
The "reordered_aligned_xyz" files preserve the original catalyst-relative geometry.
The "substrate_superposed_*" files DO NOT preserve catalyst-relative geometry.

No numpy required.
"""

import math
import csv
import itertools
from pathlib import Path


# ============================================================
# SETTINGS
# ============================================================

FILES = {
    "con_R":       "aligned_TS-44R.xyz",
    "con_S":       "aligned_TS-44S.xyz",
    "con_R_1":       "aligned_TS-44R-1.xyz",
    "con_S_1":       "aligned_TS-44S-1.xyz",
    "step_R":      "aligned_unsy_TS-42R.xyz",
    "step_R_1":     "aligned_unsy_TS-42R-1.xyz",
    "step_S_1":      "aligned_unsy_TS_42S_Q23.xyz",

    "step_S":      "aligned_unsy_TS-42S.xyz",
}

REFERENCE = "step_R"

# 1 Fe
# 2-65 ligand
# 66-67 reactive O
# 65-end substrate
SUB_START =65

BOND_SCALE = 1.25
MAX_MAPPINGS = 10000

# If auto-detection of styrene anchors is wrong, set manually, e.g.
# REACTIVE_ANCHORS_REF = (65, 67, 70)
#
# Meaning:
#   terminal alkene C, internal alkene C, phenyl ipso C
REACTIVE_ANCHORS_REF = (69,70,66)


OUT_REORDERED = "reordered_aligned_xyz"
OUT_SUPER_ALL = "substrate_superposed_allheavy"
OUT_SUPER_REACTIVE = "substrate_superposed_reactive3"


COV = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
}


# ============================================================
# XYZ / GEOMETRY
# ============================================================

def read_xyz(filename):

    with open(filename) as f:
        lines = f.readlines()

    nat = int(lines[0].strip())
    comment = (
        lines[1].rstrip("\n")
        if len(lines) > 1
        else ""
    )

    atoms = []

    for i, line in enumerate(
        lines[2:2 + nat],
        start=1
    ):

        p = line.split()

        if len(p) < 4:
            raise ValueError(
                f"Invalid XYZ line in {filename}: {line}"
            )

        atoms.append({
            "index": i,
            "element": p[0],
            "x": float(p[1]),
            "y": float(p[2]),
            "z": float(p[3]),
        })

    if len(atoms) != nat:
        raise ValueError(
            f"{filename}: expected {nat} atoms, got {len(atoms)}"
        )

    return atoms, comment


def atom_dict(atoms):

    return {
        a["index"]: a
        for a in atoms
    }


def xyz(a):

    return [
        a["x"],
        a["y"],
        a["z"],
    ]


def distance(a, b):

    return math.sqrt(
        (a["x"] - b["x"]) ** 2
        +
        (a["y"] - b["y"]) ** 2
        +
        (a["z"] - b["z"]) ** 2
    )


def point_distance(a, b):

    return math.sqrt(
        sum(
            (a[k] - b[k]) ** 2
            for k in range(3)
        )
    )


def bond_cutoff(e1, e2):

    if e1 not in COV or e2 not in COV:
        raise KeyError(
            f"Missing covalent radius for {e1} or {e2}"
        )

    return (
        BOND_SCALE
        *
        (COV[e1] + COV[e2])
    )


# ============================================================
# SUBSTRATE HEAVY-ATOM GRAPH
# ============================================================

def heavy_substrate(atoms):

    return [
        a for a in atoms
        if a["index"] >= SUB_START
        and a["element"] != "H"
    ]


def substrate_atoms(atoms):

    return [
        a for a in atoms
        if a["index"] >= SUB_START
    ]


def build_graph(atoms):

    heavy = heavy_substrate(atoms)

    graph = {
        a["index"]: set()
        for a in heavy
    }

    adict = {
        a["index"]: a
        for a in heavy
    }

    for i in range(len(heavy)):

        for j in range(i + 1, len(heavy)):

            a = heavy[i]
            b = heavy[j]

            if distance(a, b) <= bond_cutoff(
                a["element"],
                b["element"]
            ):

                graph[a["index"]].add(
                    b["index"]
                )

                graph[b["index"]].add(
                    a["index"]
                )

    return graph, adict


def signature(
    idx,
    graph,
    atoms
):

    neigh_elements = sorted(
        atoms[n]["element"]
        for n in graph[idx]
    )

    return (
        atoms[idx]["element"],
        len(graph[idx]),
        tuple(neigh_elements),
    )


# ============================================================
# GRAPH MAPPING
# ============================================================

def find_mappings(
    ref_graph,
    ref_atoms,
    target_graph,
    target_atoms
):

    candidates = {}

    for r in ref_graph:

        sr = signature(
            r,
            ref_graph,
            ref_atoms
        )

        candidates[r] = [
            t
            for t in target_graph
            if signature(
                t,
                target_graph,
                target_atoms
            ) == sr
        ]

        if not candidates[r]:
            return []

    order = sorted(
        ref_graph,
        key=lambda r: (
            len(candidates[r]),
            -len(ref_graph[r])
        )
    )

    solutions = []

    def backtrack(
        pos,
        mapping,
        used
    ):

        if len(solutions) >= MAX_MAPPINGS:
            return

        if pos == len(order):

            solutions.append(
                mapping.copy()
            )

            return

        r = order[pos]

        for t in candidates[r]:

            if t in used:
                continue

            ok = True

            for r2, t2 in mapping.items():

                ref_edge = (
                    r2 in ref_graph[r]
                )

                tar_edge = (
                    t2 in target_graph[t]
                )

                if ref_edge != tar_edge:
                    ok = False
                    break

            if not ok:
                continue

            mapping[r] = t
            used.add(t)

            backtrack(
                pos + 1,
                mapping,
                used
            )

            used.remove(t)
            del mapping[r]

    backtrack(
        0,
        {},
        set()
    )

    return solutions


# ============================================================
# RIGID-BODY SUPERPOSITION
# Horn quaternion method, pure Python
# ============================================================

def centroid(points):

    n = float(len(points))

    return [
        sum(p[k] for p in points) / n
        for k in range(3)
    ]


def mat_vec(M, v):

    return [
        sum(
            M[i][j] * v[j]
            for j in range(len(v))
        )
        for i in range(len(M))
    ]


def normalize(v):

    n = math.sqrt(
        sum(x * x for x in v)
    )

    if n == 0:
        return [
            1.0,
            0.0,
            0.0,
            0.0
        ]

    return [
        x / n
        for x in v
    ]


def largest_eigenvector_power(
    K,
    iterations=500
):

    bound = (
        max(
            sum(abs(x) for x in row)
            for row in K
        )
        + 1.0
    )

    Ks = [
        row[:]
        for row in K
    ]

    for i in range(4):
        Ks[i][i] += bound

    q = normalize(
        [1.0, 0.3, 0.2, 0.1]
    )

    for _ in range(iterations):

        q_new = normalize(
            mat_vec(Ks, q)
        )

        d1 = sum(
            (q_new[i] - q[i]) ** 2
            for i in range(4)
        )

        d2 = sum(
            (q_new[i] + q[i]) ** 2
            for i in range(4)
        )

        q = q_new

        if min(d1, d2) < 1e-30:
            break

    return q


def horn_rotation(
    target_points,
    ref_points
):

    """
    Return R,t such that:
        R * target + t
    best matches reference.
    """

    if len(target_points) != len(ref_points):
        raise ValueError(
            "Point-set sizes do not match."
        )

    if len(target_points) < 3:
        raise ValueError(
            "At least 3 points are needed."
        )

    ct = centroid(target_points)
    cr = centroid(ref_points)

    X = [
        [
            p[k] - ct[k]
            for k in range(3)
        ]
        for p in target_points
    ]

    Y = [
        [
            p[k] - cr[k]
            for k in range(3)
        ]
        for p in ref_points
    ]

    Sxx = sum(
        x[0] * y[0]
        for x, y in zip(X, Y)
    )

    Sxy = sum(
        x[0] * y[1]
        for x, y in zip(X, Y)
    )

    Sxz = sum(
        x[0] * y[2]
        for x, y in zip(X, Y)
    )

    Syx = sum(
        x[1] * y[0]
        for x, y in zip(X, Y)
    )

    Syy = sum(
        x[1] * y[1]
        for x, y in zip(X, Y)
    )

    Syz = sum(
        x[1] * y[2]
        for x, y in zip(X, Y)
    )

    Szx = sum(
        x[2] * y[0]
        for x, y in zip(X, Y)
    )

    Szy = sum(
        x[2] * y[1]
        for x, y in zip(X, Y)
    )

    Szz = sum(
        x[2] * y[2]
        for x, y in zip(X, Y)
    )

    K = [
        [
            Sxx + Syy + Szz,
            Syz - Szy,
            Szx - Sxz,
            Sxy - Syx
        ],
        [
            Syz - Szy,
            Sxx - Syy - Szz,
            Sxy + Syx,
            Szx + Sxz
        ],
        [
            Szx - Sxz,
            Sxy + Syx,
            -Sxx + Syy - Szz,
            Syz + Szy
        ],
        [
            Sxy - Syx,
            Szx + Sxz,
            Syz + Szy,
            -Sxx - Syy + Szz
        ],
    ]

    w, x, y, z = (
        largest_eigenvector_power(K)
    )

    R = [
        [
            1 - 2 * (y * y + z * z),
            2 * (x * y - z * w),
            2 * (x * z + y * w)
        ],
        [
            2 * (x * y + z * w),
            1 - 2 * (x * x + z * z),
            2 * (y * z - x * w)
        ],
        [
            2 * (x * z - y * w),
            2 * (y * z + x * w),
            1 - 2 * (x * x + y * y)
        ],
    ]

    Rct = mat_vec(
        R,
        ct
    )

    t = [
        cr[k] - Rct[k]
        for k in range(3)
    ]

    return R, t


def transform_point(
    p,
    R,
    t
):

    v = mat_vec(
        R,
        p
    )

    return [
        v[k] + t[k]
        for k in range(3)
    ]


def superposition_rmsd(
    target_points,
    ref_points,
    R,
    t
):

    ss = 0.0

    for p, q in zip(
        target_points,
        ref_points
    ):

        pp = transform_point(
            p,
            R,
            t
        )

        ss += sum(
            (pp[k] - q[k]) ** 2
            for k in range(3)
        )

    return math.sqrt(
        ss / len(target_points)
    )


# ============================================================
# CHOOSE GRAPH MAPPING BY ACTUAL SUPERPOSITION RMSD
# ============================================================

def evaluate_mapping_rigid_rmsd(
    mapping,
    ref_heavy,
    target_heavy
):

    ref_order = sorted(mapping)

    ref_points = [
        xyz(ref_heavy[r])
        for r in ref_order
    ]

    target_points = [
        xyz(target_heavy[mapping[r]])
        for r in ref_order
    ]

    R, t = horn_rotation(
        target_points,
        ref_points
    )

    fit_rmsd = superposition_rmsd(
        target_points,
        ref_points,
        R,
        t
    )

    return fit_rmsd, R, t


def best_mapping(
    ref_atoms_all,
    target_atoms_all
):

    rg, ra = build_graph(
        ref_atoms_all
    )

    tg, ta = build_graph(
        target_atoms_all
    )

    mappings = find_mappings(
        rg,
        ra,
        tg,
        ta
    )

    if not mappings:

        raise RuntimeError(
            "No substrate graph mapping found."
        )

    scored = []

    for m in mappings:

        fit_rmsd, R, t = (
            evaluate_mapping_rigid_rmsd(
                m,
                ra,
                ta
            )
        )

        scored.append(
            (
                fit_rmsd,
                m,
                R,
                t
            )
        )

    scored.sort(
        key=lambda x: x[0]
    )

    return (
        scored[0][1],
        scored[0][0],
        len(mappings),
        scored[0][2],
        scored[0][3]
    )


# ============================================================
# SUBSTRATE HYDROGEN PARENT / MAPPING
# ============================================================

def assign_substrate_hydrogens(
    atoms
):

    subs = substrate_atoms(atoms)

    heavy = [
        a for a in subs
        if a["element"] != "H"
    ]

    hydrogens = [
        a for a in subs
        if a["element"] == "H"
    ]

    result = {
        a["index"]: []
        for a in heavy
    }

    for h in hydrogens:

        candidates = []

        for hv in heavy:

            d = distance(
                h,
                hv
            )

            if d <= bond_cutoff(
                "H",
                hv["element"]
            ):

                candidates.append(
                    (
                        d,
                        hv["index"]
                    )
                )

        if not candidates:

            raise RuntimeError(
                f"Could not assign substrate H{h['index']} "
                "to a heavy-atom parent."
            )

        candidates.sort()

        parent = (
            candidates[0][1]
        )

        result[parent].append(
            h["index"]
        )

    for parent in result:

        result[parent].sort()

    return result


def best_h_permutation(
    ref_h_indices,
    target_h_indices,
    ref_ad,
    target_ad,
    R,
    t
):

    if len(ref_h_indices) != len(
        target_h_indices
    ):

        raise RuntimeError(
            "Attached-H count mismatch: "
            f"reference {ref_h_indices}, "
            f"target {target_h_indices}"
        )

    if not ref_h_indices:
        return {}

    best = None

    for perm in itertools.permutations(
        target_h_indices
    ):

        ss = 0.0

        for rh, th in zip(
            ref_h_indices,
            perm
        ):

            p = transform_point(
                xyz(target_ad[th]),
                R,
                t
            )

            q = xyz(
                ref_ad[rh]
            )

            ss += sum(
                (p[k] - q[k]) ** 2
                for k in range(3)
            )

        if (
            best is None
            or ss < best[0]
        ):

            best = (
                ss,
                perm
            )

    return {
        rh: th
        for rh, th in zip(
            ref_h_indices,
            best[1]
        )
    }


def build_all_substrate_mapping(
    ref_atoms,
    target_atoms,
    heavy_mapping,
    heavy_fit_R,
    heavy_fit_t
):

    ref_ad = atom_dict(
        ref_atoms
    )

    target_ad = atom_dict(
        target_atoms
    )

    ref_hparents = (
        assign_substrate_hydrogens(
            ref_atoms
        )
    )

    target_hparents = (
        assign_substrate_hydrogens(
            target_atoms
        )
    )

    full = dict(
        heavy_mapping
    )

    for ref_heavy_idx, target_heavy_idx \
            in heavy_mapping.items():

        ref_hs = (
            ref_hparents.get(
                ref_heavy_idx,
                []
            )
        )

        target_hs = (
            target_hparents.get(
                target_heavy_idx,
                []
            )
        )

        hmap = best_h_permutation(
            ref_hs,
            target_hs,
            ref_ad,
            target_ad,
            heavy_fit_R,
            heavy_fit_t
        )

        full.update(
            hmap
        )

    ref_sub_indices = sorted(
        a["index"]
        for a in ref_atoms
        if a["index"] >= SUB_START
    )

    missing = [
        i for i in ref_sub_indices
        if i not in full
    ]

    if missing:

        raise RuntimeError(
            "Full substrate mapping is incomplete. "
            f"Missing reference atoms: {missing}"
        )

    return full


# ============================================================
# AUTO-DETECT STYRENE REACTIVE ANCHORS
# terminal alkene C -> internal alkene C -> phenyl ipso C
# ============================================================

def detect_styrene_anchors(
    ref_atoms
):

    graph, ad = build_graph(
        ref_atoms
    )

    candidates = []

    for terminal in graph:

        if (
            ad[terminal]["element"] != "C"
            or len(graph[terminal]) != 1
        ):
            continue

        vinyl = next(
            iter(graph[terminal])
        )

        if (
            ad[vinyl]["element"] != "C"
            or len(graph[vinyl]) != 2
        ):
            continue

        others = [
            x for x in graph[vinyl]
            if x != terminal
        ]

        if len(others) != 1:
            continue

        ipso = others[0]

        if (
            ad[ipso]["element"] == "C"
            and len(graph[ipso]) == 3
        ):

            candidates.append(
                (
                    terminal,
                    vinyl,
                    ipso
                )
            )

    if len(candidates) == 1:
        return candidates[0]

    if len(candidates) == 0:

        raise RuntimeError(
            "Could not auto-detect styrene reactive anchors. "
            "Set REACTIVE_ANCHORS_REF manually, e.g. "
            "(65, 67, 70)."
        )

    raise RuntimeError(
        "Multiple possible styrene anchor triplets found: "
        f"{candidates}. "
        "Set REACTIVE_ANCHORS_REF manually."
    )


# ============================================================
# OUTPUT
# ============================================================

def write_full_reordered_xyz(
    outfile,
    ref_atoms,
    target_atoms,
    full_mapping,
    source_name
):

    """
    Write complete target structure:
      atoms 1..64 stay target original order
      substrate 65..end is reordered into REFERENCE atom order

    COORDINATES ARE NOT MOVED.
    """

    target_ad = atom_dict(
        target_atoms
    )

    ref_sub_indices = sorted(
        a["index"]
        for a in ref_atoms
        if a["index"] >= SUB_START
    )

    output_atoms = []

    # catalyst / O part unchanged
    for a in target_atoms:

        if a["index"] < SUB_START:

            output_atoms.append(
                a
            )

    # substrate in reference chemical order
    for ref_idx in ref_sub_indices:

        target_idx = (
            full_mapping[ref_idx]
        )

        output_atoms.append(
            target_ad[target_idx]
        )

    if len(output_atoms) != len(
        target_atoms
    ):

        raise RuntimeError(
            "Atom-count mismatch while writing reordered XYZ."
        )

    with open(
        outfile,
        "w"
    ) as f:

        f.write(
            f"{len(output_atoms)}\n"
        )

        f.write(
            "FULL structure; substrate reordered to "
            f"{REFERENCE}; coordinates unchanged; "
            f"source={source_name}\n"
        )

        for a in output_atoms:

            f.write(
                f"{a['element']:<3s}"
                f"{a['x']:15.8f}"
                f"{a['y']:15.8f}"
                f"{a['z']:15.8f}\n"
            )


def write_substrate_superposed_xyz(
    outfile,
    ref_atoms,
    target_atoms,
    full_mapping,
    R,
    t,
    source_name,
    fit_label,
    fit_rmsd
):

    target_ad = atom_dict(
        target_atoms
    )

    ref_sub_indices = sorted(
        a["index"]
        for a in ref_atoms
        if a["index"] >= SUB_START
    )

    with open(
        outfile,
        "w"
    ) as f:

        f.write(
            f"{len(ref_sub_indices)}\n"
        )

        f.write(
            f"substrate reordered + {fit_label} "
            f"superposition to {REFERENCE}; "
            f"source={source_name}; "
            f"fit_RMSD={fit_rmsd:.6f} A\n"
        )

        for ref_idx in ref_sub_indices:

            target_idx = (
                full_mapping[ref_idx]
            )

            a = target_ad[
                target_idx
            ]

            p = transform_point(
                xyz(a),
                R,
                t
            )

            f.write(
                f"{a['element']:<3s}"
                f"{p[0]:15.8f}"
                f"{p[1]:15.8f}"
                f"{p[2]:15.8f}\n"
            )


# ============================================================
# MAIN
# ============================================================

def main():

    structures = {}
    comments = {}

    for name, filename in FILES.items():

        if not Path(filename).exists():

            raise FileNotFoundError(
                filename
            )

        atoms, comment = read_xyz(
            filename
        )

        structures[name] = atoms
        comments[name] = comment

    ref_atoms = structures[
        REFERENCE
    ]

    ref_graph, ref_heavy = build_graph(
        ref_atoms
    )

    # Determine reactive-center anchors.
    if REACTIVE_ANCHORS_REF is None:

        reactive_anchors = (
            detect_styrene_anchors(
                ref_atoms
            )
        )

    else:

        reactive_anchors = tuple(
            REACTIVE_ANCHORS_REF
        )

    print()
    print(
        "Reference:",
        REFERENCE
    )

    print(
        "Reactive-center reference anchors:",
        reactive_anchors
    )

    ref_ad = atom_dict(
        ref_atoms
    )

    print(
        "Anchor elements:",
        " ".join(
            f"{ref_ad[i]['element']}{i}"
            for i in reactive_anchors
        )
    )

    print()

    # output dirs
    for d in [
        OUT_REORDERED,
        OUT_SUPER_ALL,
        OUT_SUPER_REACTIVE
    ]:

        Path(d).mkdir(
            parents=True,
            exist_ok=True
        )

    mapping_rows = []
    report_rows = []

    # --------------------------------------------------------
    # Every structure
    # --------------------------------------------------------

    for name, atoms in \
            structures.items():

        target_ad = atom_dict(
            atoms
        )

        # ---- heavy mapping
        if name == REFERENCE:

            heavy_mapping = {
                i: i
                for i in ref_graph
            }

            n_mappings = 1

            R_heavy = [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]

            t_heavy = [
                0.0,
                0.0,
                0.0
            ]

            heavy_rmsd = 0.0

        else:

            (
                heavy_mapping,
                heavy_rmsd,
                n_mappings,
                R_heavy,
                t_heavy
            ) = best_mapping(
                ref_atoms,
                atoms
            )

        # ---- all substrate atom mapping including H
        full_mapping = (
            build_all_substrate_mapping(
                ref_atoms,
                atoms,
                heavy_mapping,
                R_heavy,
                t_heavy
            )
        )

        # ----------------------------------------------------
        # A) Write full structure with substrate atom ORDER fixed
        #    Geometry remains untouched.
        # ----------------------------------------------------

        write_full_reordered_xyz(
            Path(OUT_REORDERED)
            /
            f"{name}_reordered.xyz",
            ref_atoms,
            atoms,
            full_mapping,
            FILES[name]
        )

        # ----------------------------------------------------
        # B) All-heavy substrate superposition
        # ----------------------------------------------------

        ref_order = sorted(
            ref_graph
        )

        ref_points_all = [
            xyz(
                ref_ad[r]
            )
            for r in ref_order
        ]

        target_points_all = [
            xyz(
                target_ad[
                    heavy_mapping[r]
                ]
            )
            for r in ref_order
        ]

        if name == REFERENCE:

            R_all = [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]

            t_all = [
                0.0,
                0.0,
                0.0
            ]

        else:

            R_all, t_all = (
                horn_rotation(
                    target_points_all,
                    ref_points_all
                )
            )

        rmsd_all = (
            superposition_rmsd(
                target_points_all,
                ref_points_all,
                R_all,
                t_all
            )
        )

        write_substrate_superposed_xyz(
            Path(OUT_SUPER_ALL)
            /
            f"{name}_substrate_allheavy.xyz",
            ref_atoms,
            atoms,
            full_mapping,
            R_all,
            t_all,
            FILES[name],
            "all-heavy",
            rmsd_all
        )

        # ----------------------------------------------------
        # C) Reactive-center 3-atom superposition
        # ----------------------------------------------------

        ref_points_reactive = [
            xyz(
                ref_ad[r]
            )
            for r in reactive_anchors
        ]

        target_points_reactive = [
            xyz(
                target_ad[
                    heavy_mapping[r]
                ]
            )
            for r in reactive_anchors
        ]

        if name == REFERENCE:

            R_reactive = [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]

            t_reactive = [
                0.0,
                0.0,
                0.0
            ]

        else:

            R_reactive, t_reactive = (
                horn_rotation(
                    target_points_reactive,
                    ref_points_reactive
                )
            )

        rmsd_reactive = (
            superposition_rmsd(
                target_points_reactive,
                ref_points_reactive,
                R_reactive,
                t_reactive
            )
        )

        write_substrate_superposed_xyz(
            Path(OUT_SUPER_REACTIVE)
            /
            f"{name}_substrate_reactive3.xyz",
            ref_atoms,
            atoms,
            full_mapping,
            R_reactive,
            t_reactive,
            FILES[name],
            "reactive-3-atom",
            rmsd_reactive
        )

        # ----------------------------------------------------
        # Reports
        # ----------------------------------------------------

        report_rows.append({
            "structure":
                name,

            "source_file":
                FILES[name],

            "number_of_heavy_graph_mappings":
                n_mappings,

            "best_allheavy_superposition_RMSD_A":
                f"{rmsd_all:.6f}",

            "reactive3_superposition_RMSD_A":
                f"{rmsd_reactive:.6f}",

            "reactive_anchor_reference_atoms":
                ";".join(
                    str(i)
                    for i in reactive_anchors
                ),

            "reactive_anchor_mapped_atoms":
                ";".join(
                    str(
                        heavy_mapping[i]
                    )
                    for i in reactive_anchors
                ),
        })

        for ref_idx in sorted(
            full_mapping
        ):

            target_idx = (
                full_mapping[ref_idx]
            )

            mapping_rows.append({
                "structure":
                    name,

                "reference_atom":
                    ref_idx,

                "reference_element":
                    ref_ad[ref_idx]["element"],

                "mapped_original_atom":
                    target_idx,

                "mapped_original_element":
                    target_ad[target_idx]["element"],

                "new_reordered_atom":
                    ref_idx,
            })

        print(
            f"{name:<14s} "
            f"graph_maps={n_mappings:<3d} "
            f"all-heavy RMSD={rmsd_all:8.4f} A   "
            f"reactive3 RMSD={rmsd_reactive:8.4f} A"
        )

    # --------------------------------------------------------
    # CSV outputs
    # --------------------------------------------------------

    with open(
        "substrate_atom_mapping_v2.csv",
        "w",
        newline=""
    ) as f:

        fields = [
            "structure",
            "reference_atom",
            "reference_element",
            "mapped_original_atom",
            "mapped_original_element",
            "new_reordered_atom",
        ]

        w = csv.DictWriter(
            f,
            fieldnames=fields
        )

        w.writeheader()
        w.writerows(
            mapping_rows
        )

    with open(
        "substrate_alignment_report_v2.csv",
        "w",
        newline=""
    ) as f:

        fields = [
            "structure",
            "source_file",
            "number_of_heavy_graph_mappings",
            "best_allheavy_superposition_RMSD_A",
            "reactive3_superposition_RMSD_A",
            "reactive_anchor_reference_atoms",
            "reactive_anchor_mapped_atoms",
        ]

        w = csv.DictWriter(
            f,
            fieldnames=fields
        )

        w.writeheader()
        w.writerows(
            report_rows
        )

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)
    print()
    print("1) USE FOR YOUR NEXT DESCRIPTOR / C-O / XY ANALYSIS:")
    print(
        f"   {OUT_REORDERED}/"
        "<structure>_reordered.xyz"
    )
    print(
        "   These keep the original catalyst-relative coordinates, "
        "but substrate atom numbers now follow the reference."
    )
    print()
    print("2) USE ONLY TO VISUALLY CHECK SUBSTRATE OVERLAP:")
    print(
        f"   {OUT_SUPER_ALL}/"
        "*_substrate_allheavy.xyz"
    )
    print(
        f"   {OUT_SUPER_REACTIVE}/"
        "*_substrate_reactive3.xyz"
    )
    print()
    print("3) MAPPING QC:")
    print("   substrate_atom_mapping_v2.csv")
    print("   substrate_alignment_report_v2.csv")
    print()
    print(
        "If reactive-center anchors are not the two alkene carbons + ipso carbon,"
    )
    print(
        "set REACTIVE_ANCHORS_REF manually at the top of the script."
    )


if __name__ == "__main__":
    main()
