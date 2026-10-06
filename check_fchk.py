import sys
import re
import math


def read_array(lines, key, dtype=float):
    for i, line in enumerate(lines):
        if line.startswith(key):
            m = re.search(r'N=\s*(\d+)', line)
            if not m:
                raise RuntimeError(f"Cannot find N= for {key}")

            n = int(m.group(1))
            vals = []
            j = i + 1

            while len(vals) < n:
                vals.extend(lines[j].split())
                j += 1

            return [dtype(x.replace("D", "E")) for x in vals[:n]]

    raise RuntimeError(f"Cannot find {key}")


def read_fchk(filename):
    with open(filename) as f:
        lines = f.readlines()

    atoms = read_array(lines, "Atomic numbers", int)
    coords = read_array(lines, "Current cartesian coordinates", float)

    xyz = [
        tuple(coords[i:i+3])
        for i in range(0, len(coords), 3)
    ]

    return atoms, xyz


def dist(a, b):
    return math.sqrt(
        (a[0]-b[0])**2 +
        (a[1]-b[1])**2 +
        (a[2]-b[2])**2
    )


def map_fragment(name, zfrag, xfrag, zcomp, xcomp, tol=1e-5):

    print(f"\n{name} mapping:")
    print("frag_atom   -> complex_atom    Z       distance(Bohr)")
    print("------------------------------------------------------")

    used = set()
    mapping = []

    for i, (z, xyz) in enumerate(zip(zfrag, xfrag), 1):

        candidates = []

        for j, (zc, xyzc) in enumerate(zip(zcomp, xcomp), 1):

            if j in used:
                continue

            if zc != z:
                continue

            d = dist(xyz, xyzc)

            if d < tol:
                candidates.append((d, j))

        if not candidates:
            print(f"{i:5d}      -> NOT FOUND       Z={z}")
            mapping.append(None)
            continue

        candidates.sort()
        d, j = candidates[0]

        used.add(j)
        mapping.append(j)

        print(
            f"{i:5d}      -> {j:5d}"
            f"             {z:3d}      {d:.3e}"
        )

    print("\nComplex atom sequence represented by this fragment:")
    print(mapping)

    return mapping


def main():

    if len(sys.argv) != 4:
        print(
            "Usage: python check_fchk.py "
            "complex.fchk cat.fchk sub.fchk"
        )
        sys.exit()

    complex_file = sys.argv[1]
    cat_file = sys.argv[2]
    sub_file = sys.argv[3]

    zc, xc = read_fchk(complex_file)
    zcat, xcat = read_fchk(cat_file)
    zsub, xsub = read_fchk(sub_file)

    print(f"Complex atoms : {len(zc)}")
    print(f"CAT atoms     : {len(zcat)}")
    print(f"SUB atoms     : {len(zsub)}")

    catmap = map_fragment(
        "CAT", zcat, xcat, zc, xc
    )

    submap = map_fragment(
        "SUB", zsub, xsub, zc, xc
    )

    print("\n================ SUMMARY ================")

    print("\nCAT complex indices:")
    print(catmap)

    print("\nSUB complex indices:")
    print(submap)

    cat_clean = [x for x in catmap if x is not None]
    sub_clean = [x for x in submap if x is not None]

    if cat_clean == list(range(1, len(cat_clean)+1)) and \
       sub_clean == list(
           range(len(cat_clean)+1,
                 len(cat_clean)+len(sub_clean)+1)
       ):

        print("\n>>> Complex ordering is CAT + SUB.")

    elif sub_clean == list(range(1, len(sub_clean)+1)) and \
         cat_clean == list(
             range(len(sub_clean)+1,
                   len(sub_clean)+len(cat_clean)+1)
         ):

        print("\n>>> Complex ordering is SUB + CAT.")

    else:
        print(
            "\n>>> CAT and SUB atoms are NOT two contiguous blocks "
            "in the complex."
        )


if __name__ == "__main__":
    main()



