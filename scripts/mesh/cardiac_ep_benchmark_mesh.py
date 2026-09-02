#!/usr/bin/env python3
"""
Cardiac monodomain benchmark mesh generator (Niederer et al., 2011,
Phil. Trans. R. Soc. A, doi:10.1098/rsta.2011.0139).

Geometry (Figure 1): a cuboid slab
    3 mm (thickness) x 7 mm (width) x 20 mm (length).
The stimulus is applied in a 1.5 mm cube located at one corner of the slab.

This script uses the gmsh Python API to build a *structured* hexahedral mesh
at a chosen edge length h (the benchmark uses h = 0.5, 0.2, 0.1 mm) and exports
it to the required XML format:

    <mesh celltype="hexahedron" dim="3">
        <nodes size="N"> <node id=.. x=.. y=.. z=../> ... </nodes>
        <elements size="M"> <element id=.. v0..v7=../> ... </elements>
        <element_data type="fiber_isotropic">          <!-- ISOTROPIC, empty fibers -->
            <element id=..> <fiber></fiber> </element> ...
        </element_data>
    </mesh>
    <electrophysiology>
        <stimuli number="1"> <stim .../> </stimuli>
    </electrophysiology>

Coordinates are written in micrometres (µm), matching the reference file.

Usage:
    python3 benchmark_mesh.py                 # generates h=0.5,0.2,0.1 mm
    python3 benchmark_mesh.py --h 0.5         # one resolution
    python3 benchmark_mesh.py --h 0.5 --units mm
"""

import argparse
import os
import gmsh

# ---------------------------------------------------------------------------
# Benchmark geometry (in mm), following the paper's convention:
#   x : thickness (3 mm), y : width (7 mm), z : length (20 mm)
# ---------------------------------------------------------------------------
LX, LY, LZ = 3.0, 7.0, 20.0          # slab dimensions [mm]
STIM_SIZE = 1.5                       # stimulus cube edge [mm]

# Fibre direction for the benchmark is along the long (z) axis, but because the
# request is ISOTROPIC we leave the fibre tag empty and tag element_data as
# fiber_isotropic.


def build_slab(h_mm):
    """Create a structured hex mesh of the slab with edge length h_mm (mm).

    Returns (node_ids, coords, hexes) where:
        node_ids : list[int]              gmsh node tags (renumbered 0..N-1)
        coords   : list[(x,y,z)]          coordinates in mm
        hexes    : list[tuple(8 ints)]    0-based node indices per hexahedron
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("benchmark_slab")

    # Box occupying [0,LX] x [0,LY] x [0,LZ]
    box = gmsh.model.occ.addBox(0, 0, 0, LX, LY, LZ)
    gmsh.model.occ.synchronize()

    # Number of divisions per edge (at least 1). With L in {3,7,20} mm and
    # h in {0.5,0.2,0.1} mm this divides evenly, so nodes land exactly on a
    # dx = h grid (e.g. 0, 500, 1000, ... um for h = 0.5 mm), matching the
    # benchmark and the Cardiax um coordinate convention.
    nx = max(1, round(LX / h_mm))
    ny = max(1, round(LY / h_mm))
    nz = max(1, round(LZ / h_mm))

    # Transfinite + recombine -> structured hexahedra.
    # Assign transfinite node counts (divisions + 1) to every curve based on
    # its dominant direction.
    for dim, tag in gmsh.model.getEntities(1):
        # bounding box of the curve tells us its orientation
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        dxl, dyl, dzl = xmax - xmin, ymax - ymin, zmax - zmin
        if dxl >= dyl and dxl >= dzl:
            n = nx
        elif dyl >= dxl and dyl >= dzl:
            n = ny
        else:
            n = nz
        gmsh.model.mesh.setTransfiniteCurve(tag, n + 1)

    for dim, tag in gmsh.model.getEntities(2):
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(dim, tag)

    gmsh.model.mesh.setTransfiniteVolume(box)
    gmsh.model.mesh.setRecombine(3, box)

    gmsh.model.mesh.generate(3)

    # ---- extract nodes ----
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_coords = node_coords.reshape(-1, 3)
    # map gmsh tag -> contiguous 0-based id
    tag2id = {int(t): i for i, t in enumerate(node_tags)}
    coords = [tuple(node_coords[i]) for i in range(len(node_tags))]

    # ---- extract hexahedra (element type 5 = 8-node hexahedron) ----
    etypes, etags, enodes = gmsh.model.mesh.getElements(3)
    hexes = []
    for et, ns in zip(etypes, enodes):
        if et == 5:  # hexahedron
            ns = ns.reshape(-1, 8)
            for row in ns:
                hexes.append(tuple(tag2id[int(x)] for x in row))

    gmsh.finalize()
    return coords, hexes


def write_xml(path, coords, hexes, scale, precision=6, fiber_dir=(1.0, 0.0, 0.0)):
    """Write the mesh + electrophysiology XML.

    scale     : multiply mm coordinates by this (1000 -> micrometres).
    fiber_dir : constant unit fibre direction written for every element.
                (1,0,0) -> fibres along x (the request).
                Note: the slab's long axis is z (20 mm); the benchmark's
                longitudinal fibre is along that long axis, i.e. (0,0,1).
                If None, fibres are left empty and the mesh is tagged
                isotropic instead.
    """
    N = len(coords)
    M = len(hexes)

    # Stimulus box (a 1.5 mm cube at the origin corner), expressed in the
    # written coordinate system (after scaling).
    s = STIM_SIZE * scale
    x0, y0, z0 = 0.0, 0.0, 0.0
    x1, y1, z1 = s, s, s

    fx = f"%.{precision}f"

    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<mesh celltype="hexahedron" dim="3">\n')

        # nodes
        f.write(f'  <nodes size="{N}">\n')
        for i, (x, y, z) in enumerate(coords):
            f.write(f'    <node id="{i}" x="{fx % (x*scale)}" '
                    f'y="{fx % (y*scale)}" z="{fx % (z*scale)}" />\n')
        f.write('  </nodes>\n')

        # elements
        f.write(f'  <elements size="{M}">\n')
        for i, hx in enumerate(hexes):
            f.write(f'    <element id="{i}" '
                    + " ".join(f'v{j}="{hx[j]}"' for j in range(8))
                    + '  />\n')
        f.write('  </elements>\n')

        # element_data : fibre direction per element.
        #   fiber_dir given   -> fibres along that direction, tagged
        #                        fiber_transversely_isotropic so Cardiax's
        #                        calc_cond_tensor actually uses get_fiber().
        #   fiber_dir is None -> empty fibres, tagged fiber_isotropic.
        if fiber_dir is None:
            f.write('  <element_data type="fiber_isotropic">\n')
            for i in range(M):
                f.write(f'    <element id="{i}">\n')
                f.write('        <fiber></fiber>\n')
                f.write('    </element>\n')
        else:
            fdx, fdy, fdz = fiber_dir
            fiber_str = f"{fx % fdx},{fx % fdy},{fx % fdz}"
            f.write('  <element_data type="fiber_transversely_isotropic">\n')
            for i in range(M):
                f.write(f'    <element id="{i}">\n')
                f.write(f'        <fiber>{fiber_str}</fiber>\n')
                f.write('    </element>\n')
        f.write('  </element_data>\n')
        f.write('</mesh>\n')

        # electrophysiology : single stimulus in the 1.5 mm corner cube.
        #
        # In Cardiax the <stim> "value" is passed straight to the ionic cell
        # model (Monodomain::solve_odes -> cells->advance(...stim_val...)), so
        # it is a per-cell stimulus current in the model's units (uA/uF for
        # human models such as ToR-ORd), NOT the paper's volumetric current
        # density 50000 uA/cm^3. That volumetric figure only becomes a per-cell
        # amplitude after dividing by the surface-to-volume ratio and membrane
        # capacitance, which the cell model handles internally.
        #
        # A typical suprathreshold amplitude for ToR-ORd-type models is around
        # -53 uA/uF for ~1-2 ms. Adjust value/duration to your cell model.
        f.write('<electrophysiology>\n')
        f.write('  <stimuli number="1">\n')
        f.write(f'      <stim start="0.0" duration="2.0" value="-50.0" '
                f'x0="{fx % x0}" x1="{fx % x1}" '
                f'y0="{fx % y0}" y1="{fx % y1}" '
                f'z0="{fx % z0}" z1="{fx % z1}" />\n')
        f.write('  </stimuli>\n')
        f.write('</electrophysiology>\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, nargs="*", default=[0.5, 0.2, 0.1],
                    help="edge length(s) in mm (benchmark uses 0.5, 0.2, 0.1)")
    ap.add_argument("--units", choices=["um", "mm"], default="um",
                    help="output coordinate units (default: um, matches reference)")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--fiber", choices=["x", "y", "z", "iso"], default="x",
                    help="constant fibre direction (default: x, as requested; "
                         "'z' is the benchmark's longitudinal long-axis; "
                         "'iso' writes empty isotropic fibres)")
    args = ap.parse_args()

    fiber_map = {"x": (1.0, 0.0, 0.0),
                 "y": (0.0, 1.0, 0.0),
                 "z": (0.0, 0.0, 1.0),
                 "iso": None}
    fiber_dir = fiber_map[args.fiber]

    scale = 1000.0 if args.units == "um" else 1.0
    os.makedirs(args.outdir, exist_ok=True)

    for h in args.h:
        coords, hexes = build_slab(h)
        tag = f"{h:g}mm".replace(".", "p")
        out = os.path.join(args.outdir, f"benchmark_slab_{tag}.xml")
        write_xml(out, coords, hexes, scale, fiber_dir=fiber_dir)
        print(f"h={h} mm  ->  {out}   nodes={len(coords)}  hexes={len(hexes)}  fiber={args.fiber}")


if __name__ == "__main__":
    main()
