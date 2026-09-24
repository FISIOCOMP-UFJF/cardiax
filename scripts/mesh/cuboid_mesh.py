#!/usr/bin/env python3
"""
cuboid_mesh.py - Generate a structured cuboid mesh and write it as XML.

Defaults reproduce a 4000 x 4000 x 4000 cube with dx = dy = dz = 100
(40 layers per axis).

Per-axis spacing starts at --start-dx/dy/dz and grows geometrically by a
ratio. If no --ratio is given, the ratio is solved automatically so that the
layers exactly fill the side length (ratio = 1.0 -> uniform spacing).

Output formats:
  xml    : <mesh><nodes><elements><element_data>  (default, with fibers)
  dolfin : FEniCS/DOLFIN legacy XML  (<dolfin><mesh>...)
  vtu    : VTK XML UnstructuredGrid   (open in ParaView, fibers as CellData)

Examples:
  python cuboid_mesh.py -o cube.xml
  python cuboid_mesh.py --lx 3000 --ly 7000 --lz 20000 --fiber 0,0,1 -o slab.xml
  python cuboid_mesh.py --lx 8000 --ly 4000 --lz 2000 --nx 80 --ny 40 --nz 20 -o box.xml
  python cuboid_mesh.py --num-layers 30 --start-dz 20 -o graded.xml   # z graded to fit 4000
  python cuboid_mesh.py --cell tet --format vtu --center -o cube.vtu
"""

import argparse
import sys

# Local hex corner c -> offsets (bit0 = x, bit1 = y, bit2 = z)
CORNERS = [(c & 1, (c >> 1) & 1, (c >> 2) & 1) for c in range(8)]
# DOLFIN hexahedron ordering is tensor-product (x fastest) = corners 0..7
DOLFIN_HEX = [0, 1, 2, 3, 4, 5, 6, 7]
# VTK_HEXAHEDRON: bottom face CCW, then top face CCW
VTK_HEX = [0, 1, 3, 2, 4, 5, 7, 6]
# Kuhn split of a hex into 6 tets around the 0-7 diagonal (conforming mesh)
KUHN_TETS = [(0, 1, 3, 7), (0, 1, 5, 7), (0, 2, 3, 7),
             (0, 2, 6, 7), (0, 4, 5, 7), (0, 4, 6, 7)]


def series_length(first, n, r):
    if abs(r - 1.0) < 1e-14:
        return first * n
    return first * (r ** n - 1.0) / (r - 1.0)


def solve_ratio(first, n, length, axis):
    """Find growth ratio r so that sum(first * r**i, i < n) == length."""
    if first <= 0 or length <= 0 or n < 1:
        sys.exit(f"[{axis}] start spacing, length and layers must be positive")
    if abs(first * n - length) <= 1e-9 * length:
        return 1.0
    if n == 1:
        sys.exit(f"[{axis}] 1 layer requires start spacing == side length")
    if first >= length:
        sys.exit(f"[{axis}] start spacing {first} >= side length {length}")
    if first * n > length:          # need shrinking cells
        lo, hi = 1e-12, 1.0
    else:                           # need growing cells
        lo, hi = 1.0, 2.0
        while series_length(first, n, hi) < length:
            hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if series_length(first, n, mid) < length:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def axis_coords(origin, first, n, length, ratio, axis):
    if ratio is None:
        ratio = solve_ratio(first, n, length, axis)
        snap = True
    else:
        snap = False
    coords, d = [origin], first
    for _ in range(n):
        coords.append(coords[-1] + d)
        d *= ratio
    if snap:  # kill round-off so the far face sits exactly at origin + length
        coords[-1] = origin + length
    return coords, ratio


def tet_volume_sign(p):
    ax, ay, az = (p[1][i] - p[0][i] for i in range(3))
    bx, by, bz = (p[2][i] - p[0][i] for i in range(3))
    cx, cy, cz = (p[3][i] - p[0][i] for i in range(3))
    return ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)


def build_cells(xs, ys, zs, cell, hex_order):
    nx, ny, nz = len(xs) - 1, len(ys) - 1, len(zs) - 1
    sx, sxy = nx + 1, (nx + 1) * (ny + 1)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                corner = [(i + a) + sx * (j + b) + sxy * (k + c) for a, b, c in CORNERS]
                if cell == "hex":
                    yield [corner[c] for c in hex_order]
                else:
                    for t in KUHN_TETS:
                        pts = [(xs[i + CORNERS[c][0]], ys[j + CORNERS[c][1]],
                                zs[k + CORNERS[c][2]]) for c in t]
                        v = [corner[c] for c in t]
                        if tet_volume_sign(pts) < 0:  # enforce positive orientation
                            v[2], v[3] = v[3], v[2]
                        yield v


def vertices(xs, ys, zs):
    for z in zs:
        for y in ys:
            for x in xs:
                yield x, y, z


def write_dolfin(out, xs, ys, zs, cell, fmt, fiber=None, fiber_type=None):
    nv = len(xs) * len(ys) * len(zs)
    nh = (len(xs) - 1) * (len(ys) - 1) * (len(zs) - 1)
    nc = nh if cell == "hex" else 6 * nh
    tag = "hexahedron" if cell == "hex" else "tetrahedron"

    out.write('<?xml version="1.0"?>\n')
    out.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
    out.write(f'  <mesh celltype="{tag}" dim="3">\n')
    out.write(f'    <vertices size="{nv}">\n')
    for idx, (x, y, z) in enumerate(vertices(xs, ys, zs)):
        out.write(f'      <vertex index="{idx}" x="{x:{fmt}}" y="{y:{fmt}}" z="{z:{fmt}}" />\n')
    out.write('    </vertices>\n')
    out.write(f'    <cells size="{nc}">\n')
    for idx, v in enumerate(build_cells(xs, ys, zs, cell, DOLFIN_HEX)):
        attrs = " ".join(f'v{n}="{vi}"' for n, vi in enumerate(v))
        out.write(f'      <{tag} index="{idx}" {attrs} />\n')
    out.write('    </cells>\n')
    out.write('  </mesh>\n')
    out.write('</dolfin>\n')


def fiber_str(f):
    return ",".join(f"{c:.6f}" for c in f)


def write_xml(out, xs, ys, zs, cell, fmt, fiber=None, fiber_type=None):
    nv = len(xs) * len(ys) * len(zs)
    nh = (len(xs) - 1) * (len(ys) - 1) * (len(zs) - 1)
    nc = nh if cell == "hex" else 6 * nh
    tag = "hexahedron" if cell == "hex" else "tetrahedron"

    out.write('<?xml version="1.0"?>\n')
    out.write(f'<mesh celltype="{tag}" dim="3">\n')
    out.write(f'  <nodes size="{nv}">\n')
    for idx, (x, y, z) in enumerate(vertices(xs, ys, zs)):
        out.write(f'    <node id="{idx}" x="{x:{fmt}}" y="{y:{fmt}}" z="{z:{fmt}}" />\n')
    out.write('  </nodes>\n')
    out.write(f'  <elements size="{nc}">\n')
    for idx, v in enumerate(build_cells(xs, ys, zs, cell, VTK_HEX)):
        attrs = " ".join(f'v{n}="{vi}"' for n, vi in enumerate(v))
        out.write(f'    <element id="{idx}" {attrs}  />\n')
    out.write('  </elements>\n')
    if fiber is not None:
        fs = fiber_str(fiber)
        out.write(f'  <element_data type="{fiber_type}">\n')
        for idx in range(nc):
            out.write(f'    <element id="{idx}">\n'
                      f'        <fiber>{fs}</fiber>\n'
                      f'    </element>\n')
        out.write('  </element_data>\n')
    out.write('</mesh>\n')


def write_vtu(out, xs, ys, zs, cell, fmt, fiber=None, fiber_type=None):
    nv = len(xs) * len(ys) * len(zs)
    nh = (len(xs) - 1) * (len(ys) - 1) * (len(zs) - 1)
    nc = nh if cell == "hex" else 6 * nh
    npc, vtk_type = (8, 12) if cell == "hex" else (4, 10)

    out.write('<?xml version="1.0"?>\n')
    out.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
    out.write('  <UnstructuredGrid>\n')
    out.write(f'    <Piece NumberOfPoints="{nv}" NumberOfCells="{nc}">\n')
    if fiber is not None:
        out.write('      <CellData Vectors="fiber">\n')
        out.write('        <DataArray type="Float64" Name="fiber" '
                  'NumberOfComponents="3" format="ascii">\n')
        line = '          ' + " ".join(f"{c:{fmt}}" for c in fiber) + '\n'
        for _ in range(nc):
            out.write(line)
        out.write('        </DataArray>\n')
        out.write('      </CellData>\n')
    out.write('      <Points>\n')
    out.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
    for x, y, z in vertices(xs, ys, zs):
        out.write(f'          {x:{fmt}} {y:{fmt}} {z:{fmt}}\n')
    out.write('        </DataArray>\n')
    out.write('      </Points>\n')
    out.write('      <Cells>\n')
    out.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
    for v in build_cells(xs, ys, zs, cell, VTK_HEX):
        out.write('          ' + " ".join(map(str, v)) + '\n')
    out.write('        </DataArray>\n')
    out.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
    for s in range(0, nc, 16):
        out.write('          ' + " ".join(str((c + 1) * npc)
                                          for c in range(s, min(s + 16, nc))) + '\n')
    out.write('        </DataArray>\n')
    out.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
    for s in range(0, nc, 32):
        out.write('          ' + " ".join([str(vtk_type)] * (min(s + 32, nc) - s)) + '\n')
    out.write('        </DataArray>\n')
    out.write('      </Cells>\n')
    out.write('    </Piece>\n')
    out.write('  </UnstructuredGrid>\n')
    out.write('</VTKFile>\n')


def parse_vector(text):
    try:
        v = [float(c) for c in text.replace(" ", "").split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError(f"bad vector '{text}', use e.g. 1,0,0")
    if len(v) != 3:
        raise argparse.ArgumentTypeError(f"vector needs 3 components, got '{text}'")
    if sum(c * c for c in v) == 0.0:
        raise argparse.ArgumentTypeError("fiber vector must be non-zero")
    return v


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate a structured cuboid mesh as XML.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    g = p.add_argument_group("geometry")
    g.add_argument("--side-length", type=float, default=4000.0,
                   help="default side length for all axes")
    g.add_argument("--lx", type=float, help="length in x (overrides --side-length)")
    g.add_argument("--ly", type=float, help="length in y (overrides --side-length)")
    g.add_argument("--lz", type=float, help="length in z (overrides --side-length)")
    g.add_argument("--origin", type=float, nargs=3, default=[0.0, 0.0, 0.0],
                   metavar=("X0", "Y0", "Z0"), help="minimum corner of the cuboid")
    g.add_argument("--center", action="store_true",
                   help="center the cuboid on (0,0,0) (ignores --origin)")

    d = p.add_argument_group("discretization")
    d.add_argument("--num-layers", type=int, default=40,
                   help="default number of layers (cells) per axis")
    d.add_argument("--nx", type=int, help="layers in x (overrides --num-layers)")
    d.add_argument("--ny", type=int, help="layers in y (overrides --num-layers)")
    d.add_argument("--nz", type=int, help="layers in z (overrides --num-layers)")
    d.add_argument("--start-dx", type=float, default=100.0, help="first cell size in x")
    d.add_argument("--start-dy", type=float, default=100.0, help="first cell size in y")
    d.add_argument("--start-dz", type=float, default=100.0, help="first cell size in z")
    d.add_argument("--ratio", type=float,
                   help="fixed growth ratio for all axes; the side lengths then "
                        "follow from start size, layers and ratio "
                        "(default: solve ratio to fit side length)")

    o = p.add_argument_group("output")
    o.add_argument("--cell", choices=["hex", "tet"], default="hex",
                   help="hexahedra, or tetrahedra (6 per hex)")
    o.add_argument("--format", choices=["xml", "dolfin", "vtu"], default="xml",
                   help="XML flavour")
    o.add_argument("--precision", type=int, default=6,
                   help="decimals (xml) / significant digits (dolfin, vtu)")

    f = p.add_argument_group("fibers")
    f.add_argument("--fiber", type=parse_vector, metavar="FX,FY,FZ",
                   help="constant fiber direction written for every element, "
                        "e.g. 1,0,0 (omit for no <element_data>)")
    f.add_argument("--fiber-type", default="fiber_transversely_isotropic",
                   help="type attribute of <element_data>")
    f.add_argument("--no-normalize", action="store_true",
                   help="write the fiber vector as given instead of unit length")
    o.add_argument("-o", "--output", default="-", help="output file ('-' = stdout)")
    return p.parse_args()


def main():
    a = parse_args()
    L = [a.lx or a.side_length, a.ly or a.side_length, a.lz or a.side_length]
    N = [a.nx or a.num_layers, a.ny or a.num_layers, a.nz or a.num_layers]
    D = [a.start_dx, a.start_dy, a.start_dz]
    if min(N) < 1:
        sys.exit("number of layers must be >= 1")

    origin = list(a.origin)
    if a.center:
        if a.ratio is not None:  # lengths not known yet; compute them first
            L = [series_length(D[i], N[i], a.ratio) for i in range(3)]
        origin = [-l / 2.0 for l in L]

    axes = []
    for i, name in enumerate("xyz"):
        coords, r = axis_coords(origin[i], D[i], N[i], L[i], a.ratio, name)
        axes.append(coords)
        sizes = [coords[k + 1] - coords[k] for k in range(N[i])]
        print(f"{name}: [{coords[0]:g}, {coords[-1]:g}]  layers={N[i]}  "
              f"ratio={r:.8g}  d_min={min(sizes):g}  d_max={max(sizes):g}",
              file=sys.stderr)

    nv = (N[0] + 1) * (N[1] + 1) * (N[2] + 1)
    nc = N[0] * N[1] * N[2] * (1 if a.cell == "hex" else 6)
    print(f"vertices={nv}  cells={nc} ({a.cell})  format={a.format}", file=sys.stderr)

    fiber = a.fiber
    if fiber is not None:
        if not a.no_normalize:
            n = sum(c * c for c in fiber) ** 0.5
            fiber = [c / n for c in fiber]
        if a.format == "dolfin":
            print("warning: dolfin format has no fiber section; --fiber ignored",
                  file=sys.stderr)
        else:
            print(f"fiber={fiber_str(fiber)}  type={a.fiber_type}", file=sys.stderr)

    fmt = f".{a.precision}f" if a.format == "xml" else f".{a.precision}g"
    writer = {"xml": write_xml, "dolfin": write_dolfin, "vtu": write_vtu}[a.format]
    args = (*axes, a.cell, fmt, fiber, a.fiber_type)
    if a.output == "-":
        writer(sys.stdout, *args)
    else:
        with open(a.output, "w", buffering=1 << 20) as f:
            writer(f, *args)
        print(f"written: {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
