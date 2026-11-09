#!/usr/bin/env python3
"""
Slab for debugging the coupled electromechanics solver (Electromechanic).

Geometry (micrometres, as in Cardiax):
    x in [0, nx*h]    fibre direction, wave travels left -> right
    y in [0, ny*h]
    z in [0, nz*hz]   thickness (hz defaults to h)

Hexahedral elements, same XML layout and local node order (v0..v7) as
coupled_electromech_Tomek.xml:

    <elasticity>         material, ninc, pressure, prescribed displacements
    <mesh>               nodes, elements, fibres (transversely isotropic, along x)
    <electrophysiology>  one stimulus: a strip at x = 0, full width and thickness

Boundary conditions (--bc):
    free (default) - slab that contracts freely ("3-2-1" support)
        u_z = 0 as set by --zbc
        left  mid-edge node(s) (0, Ly/2)  : u_x = u_y = 0  (no rigid translation)
        right mid-edge node(s) (Lx, Ly/2) : u_y = 0        (no rigid rotation)
        With --zbc all these are the whole vertical lines (plane strain);
        otherwise only the bottom nodes (z = 0), so nothing else is constrained.
    clamp
        x = 0 face : u_x = u_y = u_z = 0, plus u_z as set by --zbc
    symmetry (symmetry planes, like the reference cube)
        x = 0 face : u_x = 0,  y = 0 face : u_y = 0, plus u_z as set by --zbc

--zbc (out-of-plane condition):
    all    : u_z = 0 at every node -> plane strain, a 2D-like problem
    bottom : u_z = 0 on the z = 0 face only (free thickening)
    both   : u_z = 0 on the bottom and top faces
    auto   : (default) "all" for --bc free with nz = 1, otherwise "bottom"

Usage:
    python3 make_em_slab.py                          # 50 x 50 x 1, h = 200 um, plane strain
    python3 make_em_slab.py --nx 200 --ny 10         # long thin strip, 40 x 2 mm
    python3 make_em_slab.py --nz 5                   # 3D slab, 1 mm thick
    python3 make_em_slab.py --nz 4 --hz 100          # 4 elements over 0.4 mm
"""

import argparse


def node_id(i, j, k, nx, ny):
    """Node numbering: x fastest, then y, then z."""
    return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i


def element_nodes(i, j, k, nx, ny):
    """Local node order of the reference file:
       v0..v3 on the face x = x_i, v4..v7 on the face x = x_{i+1};
       on each face: (y0,z1), (y0,z0), (y1,z0), (y1,z1)."""
    n = lambda a, b, c: node_id(a, b, c, nx, ny)
    return [n(i,     j,     k + 1), n(i,     j,     k), n(i,     j + 1, k), n(i,     j + 1, k + 1),
            n(i + 1, j,     k + 1), n(i + 1, j,     k), n(i + 1, j + 1, k), n(i + 1, j + 1, k + 1)]


def boundary_conditions(a, nx, ny, nz):
    """Return a sorted list of (node, direction) with prescribed u = 0."""
    zbc = a.zbc
    if zbc == "auto":
        zbc = "all" if (a.bc == "free" and nz == 1) else "bottom"

    bcs = set()   # a set, so no dof is prescribed twice
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nid = node_id(i, j, k, nx, ny)

                # out-of-plane condition
                if zbc == "all" or k == 0 or (zbc == "both" and k == nz):
                    bcs.add((nid, 2))

                # in-plane condition
                if a.bc == "clamp" and i == 0:
                    bcs.update({(nid, 0), (nid, 1), (nid, 2)})
                elif a.bc == "symmetry":
                    if i == 0:
                        bcs.add((nid, 0))
                    if j == 0:
                        bcs.add((nid, 1))

    if a.bc == "free":
        jm = ny // 2
        layers = range(nz + 1) if zbc == "all" else [0]
        for k in layers:
            left, right = node_id(0, jm, k, nx, ny), node_id(nx, jm, k, nx, ny)
            bcs.update({(left, 0), (left, 1), (right, 1)})

    return sorted(bcs, key=lambda b: (b[1], b[0])), zbc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nx", type=int, default=50, help="elements along x (fibre)")
    ap.add_argument("--ny", type=int, default=50, help="elements along y")
    ap.add_argument("--nz", type=int, default=1, help="elements along z (thickness)")
    ap.add_argument("--h", type=float, default=200.0, help="element size in x and y [um]")
    ap.add_argument("--hz", type=float, default=None, help="element size in z [um] (default: h)")
    ap.add_argument("--stim-width", type=float, default=400.0,
                    help="width of the stimulated strip at x = 0 [um]")
    ap.add_argument("--stim-value", type=float, default=1.0,
                    help="stimulus value (NP: positive depolarises; TT2/ToRORd: negative)")
    ap.add_argument("--stim-duration", type=float, default=2.0)
    ap.add_argument("--bc", choices=["free", "clamp", "symmetry"], default="free")
    ap.add_argument("--zbc", choices=["auto", "all", "bottom", "both"], default="auto")
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()

    nx, ny, nz = a.nx, a.ny, a.nz
    h  = a.h
    hz = a.hz if a.hz is not None else h
    Lx, Ly, Lz = nx * h, ny * h, nz * hz
    out = a.out or f"em_slab_{nx}x{ny}x{nz}_{int(h)}um.xml"

    n_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    n_elems = nx * ny * nz

    bcs, zbc = boundary_conditions(a, nx, ny, nz)

    with open(out, "w") as f:
        w = f.write

        # ---- elasticity --------------------------------------------------
        w('<?xml version="1.0"?>\n')
        w('<elasticity type="THREE_DIM">\n')
        w('  <parameters>\n')
        w('  <material>HolzapfelOgden</material>\n')
        w('    <coefficients>0.15, 0.11685, 0.372, 0.410, 6.0, 11.83425, 5.16, 11.3, 0, 333.0</coefficients>\n')
        w('    <ninc>1</ninc>\n')
        w('  </parameters>\n')
        w('  <pressure>\n')
        w('    <node id="1" marker="30" value="0.0"/>\n')
        w('  </pressure>\n')
        w('  <prescribed_displacement>\n')
        for nid, d in bcs:
            w(f'    <node id="{nid}" direction="{d}" value="0.0" />\n')
        w('  </prescribed_displacement>\n')
        w('</elasticity>\n\n')

        # ---- mesh --------------------------------------------------------
        w('<mesh celltype="hexahedron" dim="3">\n')
        w(f'  <nodes size="{n_nodes}">\n')
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    nid = node_id(i, j, k, nx, ny)
                    w(f'    <node id="{nid}" x="{i*h:.6f}" y="{j*h:.6f}" z="{k*hz:.6f}" />\n')
        w('  </nodes>\n')

        w(f'  <elements size="{n_elems}">\n')
        eid = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    v = element_nodes(i, j, k, nx, ny)
                    attrs = " ".join(f'v{m}="{v[m]}"' for m in range(8))
                    w(f'    <element id="{eid}" {attrs}  />\n')
                    eid += 1
        w('  </elements>\n')

        w('  <element_data type="fiber_transversely_isotropic">\n')
        for e in range(n_elems):
            w(f'    <element id="{e}">\n')
            w('        <fiber>1.000000,0.000000,0.000000</fiber>\n')
            w('    </element>\n')
        w('  </element_data>\n')
        w('</mesh>\n')

        # ---- electrophysiology -------------------------------------------
        # Small margins so nodes exactly on the box faces are inside.
        eps = 0.01 * min(h, hz)
        w('<electrophysiology>\n')
        w('  <stimuli number="1">\n')
        w(f'      <stim start="0.0" duration="{a.stim_duration:g}" value="{a.stim_value:g}"'
          f' x0="{-eps:.2f}" x1="{a.stim_width + eps:.2f}"'
          f' y0="{-eps:.2f}" y1="{Ly + eps:.2f}"'
          f' z0="{-eps:.2f}" z1="{Lz + eps:.2f}" />\n')
        w('  </stimuli>\n')
        w('</electrophysiology>\n')

    print(f"{out}: {nx}x{ny}x{nz} hexahedra, h = {h:g} um, hz = {hz:g} um, "
          f"{Lx/1000:g} x {Ly/1000:g} x {Lz/1000:g} mm, "
          f"{n_nodes} nodes, {n_elems} elements, "
          f"bc = {a.bc}, zbc = {zbc}, {len(bcs)} prescribed dofs")


if __name__ == "__main__":
    main()
