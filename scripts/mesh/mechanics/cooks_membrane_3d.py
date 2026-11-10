"""
Cook's Membrane — 3D Mesh Generator
=====================================
Generates .geo, .msh, and .xml files matching the exact format
of cook_hex0_k1.xml .. cook_hex5_k1.xml.

XML structure produced:
  <elasticity type="THREE_DIM">
    <parameters>
      <material>NeoHookean</material>
      <coefficients>mu, kappa, kappa</coefficients>
      <ninc>N</ninc>
    </parameters>
    <neumann>
      <node id="1" marker="10" t0="0.0" t1="6.25" t2="0.0" />
    </neumann>
    <prescribed_displacement>
      <!-- left face nodes: directions 0 and 1 fixed (clamped) -->
      <!-- ALL nodes: direction 2 fixed (plane strain in Z)     -->
    </prescribed_displacement>
  </elasticity>
  <mesh celltype="hexahedron" dim="3">
    <nodes size="N"> ... </nodes>
    <elements size="M"> ... </elements>
    <boundary celltype="quadrilateral" dim="2">
      <!-- right face quads only, marker=10 -->
    </boundary>
  </mesh>

Hex node ordering convention (matching cook_hex*_k1.xml):
  Nodes are stored in alternating z pairs within each face:
    v0 = (x0,y0, z=h)   v1 = (x1,y1, z=h)
    v2 = (x1,y1, z=0)   v3 = (x0,y0, z=0)   <- bottom face (z=h top, z=0 bottom interleaved)
    v4 = (x2,y2, z=h)   v5 = (x3,y3, z=h)
    v6 = (x3,y3, z=0)   v7 = (x2,y2, z=0)   <- top face

  In practice: GMSH extrudes from z=0 upward, so the z=0 face is the
  "front" surface and z=h is the "back". The reference files store nodes
  column-by-column with z=h BEFORE z=0 for each (x,y) pair.

  The reorder_hex() function remaps GMSH's default hex connectivity
  [z=0 bottom: 0,1,2,3 | z=h top: 4,5,6,7]
  to the reference convention
  [v0,v1,v2,v3,v4,v5,v6,v7] = [top1, top2, bot2, bot1, top3, top4, bot4, bot3]

Geometry:
  Trapezoid extruded by h in Z:
    A = (0,  0)   left-bottom  -> clamped (directions 0 and 1)
    B = (48, 44)  right-bottom
    C = (48, 60)  right-top    -> shear load t1=6.25
    D = (0,  44)  left-top     -> clamped (directions 0 and 1)
  All nodes: direction 2 (Z) fixed -> plane strain

Usage:
  pip install gmsh
  python cooks_membrane_3d.py [--n 4] [--h 1.0] [--ninc 10] [--out_dir .]

  --n     : in-plane divisions (default 4)
  --h     : thickness in Z (default 1.0)
  --ninc  : load increments for solver (default 10)
  --out_dir: output directory (default .)

Refinement guide (matching cook_hex0..hex5):
  --n 3   -> cook_hex0  (3x3,   27 elements)
  --n 6   -> cook_hex1  (6x6,  ~112 elements)
  --n 9   -> cook_hex2
  --n 12  -> cook_hex3
  --n 18  -> cook_hex4
  --n 24  -> cook_hex5

Examples:
    python cooks_membrane_3d.py --n 3   # matches cook_hex0 (27 elements)
    python cooks_membrane_3d.py --n 6   # matches cook_hex1
    python cooks_membrane_3d.py --n 9 --ninc 20
    python cooks_membrane_3d.py --n 4 --out_dir ./results

"""

import argparse
import os
import sys

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Cook's Membrane 3D mesh generator")
parser.add_argument("--n",       type=int,   default=4,   help="In-plane divisions per edge (default 4)")
parser.add_argument("--h",       type=float, default=1.0, help="Thickness in Z (default 1.0)")
parser.add_argument("--ninc",    type=int,   default=10,  help="Load increments for solver (default 10)")
parser.add_argument("--out_dir", type=str,   default=".", help="Output directory (default: .)")
args = parser.parse_args()

N      = args.n
H      = args.h
NINC   = args.ninc
OUTDIR = args.out_dir
os.makedirs(OUTDIR, exist_ok=True)

# Material
MU    = 80.1938
KAPPA = 400942.0

# Shear load (per unit length on right face)
T1 = 6.25

# Physical group markers
MARKER_SHEAR   = 10   # right face -> Neumann + boundary output
MARKER_CLAMPED = 20   # left face  -> prescribed_displacement dir 0,1
MARKER_DOMAIN  = 100  # volume

GEO_FILE = os.path.join(OUTDIR, f"cooks_membrane_3d_n{N}.geo")
MSH_FILE = os.path.join(OUTDIR, f"cooks_membrane_3d_n{N}.msh")
XML_FILE = os.path.join(OUTDIR, f"cooks_membrane_3d_n{N}.xml")


# ---------------------------------------------------------------------------
# 1. Write .geo
# ---------------------------------------------------------------------------
def write_geo(path, n, h):
    lc = 48.0 / n
    geo = f"""// ============================================================
//  Cook's Membrane — 3D GMSH geometry
//  n={n}, h={h}
//
//  Trapezoid corners at Z=0:
//    A=(0,0)   B=(48,44)   C=(48,60)   D=(0,44)
//  Extruded by h={h} in Z (1 layer -> plane strain).
//
//  Physical markers:
//   {MARKER_SHEAR}   shear_load  (right face B-C, Neumann t1={T1})
//   {MARKER_CLAMPED}   clamped     (left face  D-A, u0=u1=0)
//   100  domain      (volume)
// ============================================================

lc = {lc:.6f};
n  = {n};
h  = {h};

Point(1) = {{  0,  0, 0, lc }};   // A
Point(2) = {{ 48, 44, 0, lc }};   // B
Point(3) = {{ 48, 60, 0, lc }};   // C  (tip, reference point)
Point(4) = {{  0, 44, 0, lc }};   // D

Line(1) = {{1, 2}};   // bottom A->B
Line(2) = {{2, 3}};   // right  B->C  (shear)
Line(3) = {{3, 4}};   // top    C->D
Line(4) = {{4, 1}};   // left   D->A  (clamped)

Curve Loop(1)    = {{1, 2, 3, 4}};
Plane Surface(1) = {{1}};

Transfinite Curve{{1}} = n+1;
Transfinite Curve{{2}} = n+1;
Transfinite Curve{{3}} = n+1;
Transfinite Curve{{4}} = n+1;
Transfinite Surface{{1}};
Recombine Surface{{1}};

// Extrude 1 layer in Z (plane strain / thin slab)
out[] = Extrude {{0, 0, h}} {{
    Surface{{1}};
    Layers{{1}};
    Recombine;
}};

// out[0]=back(Z=h), out[1]=vol, out[2]=lat(l1/bottom),
// out[3]=lat(l2/right/shear), out[4]=lat(l3/top), out[5]=lat(l4/left/clamped)

Physical Volume("domain",      {MARKER_DOMAIN}) = {{out[1]}};
Physical Surface("front",       50) = {{1}};
Physical Surface("back",        60) = {{out[0]}};
Physical Surface("bottom_edge", 30) = {{out[2]}};
Physical Surface("shear_load",  {MARKER_SHEAR}) = {{out[3]}};
Physical Surface("top_edge",    40) = {{out[4]}};
Physical Surface("clamped",     {MARKER_CLAMPED}) = {{out[5]}};
"""
    with open(path, "w") as f:
        f.write(geo)
    print(f"[1/3] Wrote {path}")


# ---------------------------------------------------------------------------
# 2. Generate mesh via GMSH Python API
# ---------------------------------------------------------------------------
def generate_mesh(geo_path, msh_path, n, h):
    try:
        import gmsh
    except ImportError:
        print("ERROR: gmsh not found. Install with:  pip install gmsh")
        sys.exit(1)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm",    8)
    gmsh.option.setNumber("Mesh.Algorithm3D",  4)

    gmsh.model.add("cooks_membrane_3d")
    fac = gmsh.model.geo
    lc  = 48.0 / n

    p1 = fac.addPoint(0.0,  0.0,  0.0, lc)
    p2 = fac.addPoint(48.0, 44.0, 0.0, lc)
    p3 = fac.addPoint(48.0, 60.0, 0.0, lc)
    p4 = fac.addPoint(0.0,  44.0, 0.0, lc)

    l1 = fac.addLine(p1, p2)   # bottom A->B
    l2 = fac.addLine(p2, p3)   # right  B->C
    l3 = fac.addLine(p3, p4)   # top    C->D
    l4 = fac.addLine(p4, p1)   # left   D->A

    cl  = fac.addCurveLoop([l1, l2, l3, l4])
    sur = fac.addPlaneSurface([cl])

    for line in [l1, l2, l3, l4]:
        fac.mesh.setTransfiniteCurve(line, n + 1)
    fac.mesh.setTransfiniteSurface(sur)
    fac.mesh.setRecombine(2, sur)
    fac.synchronize()

    out = fac.extrude([(2, sur)], 0, 0, h, numElements=[1], recombine=True)
    fac.synchronize()

    back_tag = out[0][1]
    vol_tag  = out[1][1]
    lat_bot  = out[2][1]   # bottom edge face
    lat_rgt  = out[3][1]   # right face (shear)
    lat_top  = out[4][1]   # top edge face
    lat_lft  = out[5][1]   # left face (clamped)

    gmsh.model.addPhysicalGroup(3, [vol_tag],  MARKER_DOMAIN); gmsh.model.setPhysicalName(3, MARKER_DOMAIN, "domain")
    gmsh.model.addPhysicalGroup(2, [sur],       50);            gmsh.model.setPhysicalName(2, 50,            "front")
    gmsh.model.addPhysicalGroup(2, [back_tag],  60);            gmsh.model.setPhysicalName(2, 60,            "back")
    gmsh.model.addPhysicalGroup(2, [lat_bot],   30);            gmsh.model.setPhysicalName(2, 30,            "bottom_edge")
    gmsh.model.addPhysicalGroup(2, [lat_rgt],  MARKER_SHEAR);  gmsh.model.setPhysicalName(2, MARKER_SHEAR,  "shear_load")
    gmsh.model.addPhysicalGroup(2, [lat_top],   40);            gmsh.model.setPhysicalName(2, 40,            "top_edge")
    gmsh.model.addPhysicalGroup(2, [lat_lft],  MARKER_CLAMPED);gmsh.model.setPhysicalName(2, MARKER_CLAMPED,"clamped")

    gmsh.model.mesh.generate(3)
    gmsh.write(msh_path)
    gmsh.finalize()
    print(f"[2/3] Wrote {msh_path}")


# ---------------------------------------------------------------------------
# 3. Parse .msh (format 2 or 4)  ->  nodes, hexes, quads
# ---------------------------------------------------------------------------
def _block(text, tag):
    s = text.find(f"${tag}\n")
    e = text.find(f"$End{tag}")
    if s == -1 or e == -1: return ""
    return text[s + len(tag) + 2: e].strip()


def parse_msh(msh_path):
    with open(msh_path) as f:
        content = f.read()
    fmt = _block(content, "MeshFormat")
    ver = float(fmt.split()[0]) if fmt else 2.0
    return _parse_msh4(content) if ver >= 4.0 else _parse_msh2(content)


def _parse_msh2(content):
    nodes = {}
    hexes, quads = [], []

    for line in _block(content, "Nodes").splitlines()[1:]:
        p = line.split()
        if len(p) < 4: continue
        nodes[int(p[0])-1] = (float(p[1]), float(p[2]), float(p[3]))

    eid = 0
    for line in _block(content, "Elements").splitlines()[1:]:
        p = line.split()
        if len(p) < 5: continue
        etype  = int(p[1])
        ntags  = int(p[2])
        marker = int(p[3]) if ntags > 0 else 0
        conn   = [int(x)-1 for x in p[3+ntags:]]
        if   etype == 5: hexes.append((eid, conn, marker)); eid += 1
        elif etype == 3: quads.append((eid, conn, marker)); eid += 1
        else: eid += 1

    return nodes, hexes, quads


def _parse_msh4(content):
    nodes = {}
    hexes, quads = [], []

    # ------------------------------------------------------------------
    # Build entity_tag -> physical_marker map from $Entities block.
    # In MSH4, surfaces (dim=2) and volumes (dim=3) list their physical
    # group tags in the header line. We need this to map the entity tag
    # stored per element block back to the user-defined marker (10, 20...).
    # ------------------------------------------------------------------
    entity_to_marker = {}   # (dim, entity_tag) -> physical_marker

    ent_block = _block(content, "Entities")
    if ent_block:
        ent_lines = ent_block.splitlines()
        idx = 0
        counts = ent_lines[idx].split(); idx += 1
        n_pts, n_crv, n_sur, n_vol = int(counts[0]), int(counts[1]), int(counts[2]), int(counts[3])

        # Points (dim=0): tag, x,y,z, nPhysicals, [physicals]
        for _ in range(n_pts):
            p = ent_lines[idx].split(); idx += 1
            tag = int(p[0])
            n_phys = int(p[4])
            if n_phys > 0:
                entity_to_marker[(0, tag)] = int(p[5])

        # Curves (dim=1): tag, xmin,ymin,zmin,xmax,ymax,zmax, nPhysicals, [physicals], nBounding, [bounding]
        for _ in range(n_crv):
            p = ent_lines[idx].split(); idx += 1
            tag = int(p[0])
            n_phys = int(p[7])
            if n_phys > 0:
                entity_to_marker[(1, tag)] = int(p[8])

        # Surfaces (dim=2): tag, xmin..zmax, nPhysicals, [physicals], nBounding, [bounding]
        for _ in range(n_sur):
            p = ent_lines[idx].split(); idx += 1
            tag = int(p[0])
            n_phys = int(p[7])
            if n_phys > 0:
                entity_to_marker[(2, tag)] = int(p[8])

        # Volumes (dim=3): tag, xmin..zmax, nPhysicals, [physicals], nBounding, [bounding]
        for _ in range(n_vol):
            p = ent_lines[idx].split(); idx += 1
            tag = int(p[0])
            n_phys = int(p[7])
            if n_phys > 0:
                entity_to_marker[(3, tag)] = int(p[8])

    # ------------------------------------------------------------------
    # Nodes
    # ------------------------------------------------------------------
    node_lines = _block(content, "Nodes").splitlines()
    idx = 0
    nb  = int(node_lines[idx].split()[0]); idx += 1
    for _ in range(nb):
        n_in = int(node_lines[idx].split()[3]); idx += 1
        ids  = [int(node_lines[idx+k])-1 for k in range(n_in)]; idx += n_in
        for nid in ids:
            p = node_lines[idx].split(); idx += 1
            nodes[nid] = (float(p[0]), float(p[1]), float(p[2]))

    # ------------------------------------------------------------------
    # Elements — resolve entity tag to physical marker
    # ------------------------------------------------------------------
    elem_lines = _block(content, "Elements").splitlines()
    idx = 0
    nb  = int(elem_lines[idx].split()[0]); idx += 1
    eid = 0
    for _ in range(nb):
        hdr    = elem_lines[idx].split(); idx += 1
        dim    = int(hdr[0])
        entity = int(hdr[1])
        etype  = int(hdr[2])
        n_el   = int(hdr[3])
        marker = entity_to_marker.get((dim, entity), 0)
        for _ in range(n_el):
            p    = elem_lines[idx].split(); idx += 1
            conn = [int(x)-1 for x in p[1:]]
            if   etype == 5: hexes.append((eid, conn, marker))
            elif etype == 3: quads.append((eid, conn, marker))
            eid += 1

    return nodes, hexes, quads


# ---------------------------------------------------------------------------
# 4. Reorder hex nodes to match cook_hex*_k1.xml convention
# ---------------------------------------------------------------------------
# GMSH hex connectivity (after extrude from z=0 upward):
#   bottom face (z=0): [0, 1, 2, 3]
#   top    face (z=h): [4, 5, 6, 7]
#
# Reference cook_hex*_k1.xml convention (column-major, alternating z):
#   Each (x,y) column stores z=h node BEFORE z=0 node.
#   For a hex with bottom face nodes b0,b1,b2,b3 (z=0)
#                  and top  face nodes t0,t1,t2,t3 (z=h):
#
#   v0=t0  v1=t1  v2=b1  v3=b0   <- first  quad column (z=h, z=h, z=0, z=0)
#   v4=t3  v5=t2  v6=b2  v7=b3   <- second quad column
#
# This permutation maps GMSH [b0,b1,b2,b3,t0,t1,t2,t3]
#                      to ref [t0,t1,b1,b0,t3,t2,b2,b3]
#
# Permutation index: ref[i] = gmsh[perm[i]]
GMSH_TO_REF = [4, 5, 1, 0, 7, 6, 2, 3]


def reorder_hex(conn):
    """Reorder GMSH hex node list to cook_hex*_k1.xml convention."""
    return [conn[i] for i in GMSH_TO_REF]


def reorder_quad(conn, nodes):
    """
    Reorder boundary quad nodes to match cook_hex*_k1.xml convention.
    The reference boundary quads on the shear face are ordered so that
    the z=0 nodes come first (v0,v1) then z=h nodes (v2,v3), matching
    the pattern seen in cook_hex1_k1.xml boundary elements.
    """
    # Split into z=0 and z=h nodes
    z0 = [v for v in conn if abs(nodes[v][2]) < 1e-10]
    zh = [v for v in conn if abs(nodes[v][2]) >= 1e-10]
    # Reference order: v0,v1 are z=0; v2,v3 are z=h (bottom to top in Y)
    z0_sorted = sorted(z0, key=lambda v: nodes[v][1])
    zh_sorted = sorted(zh, key=lambda v: nodes[v][1], reverse=True)
    return z0_sorted + zh_sorted


# ---------------------------------------------------------------------------
# 5. Write XML — exactly matching cook_hex*_k1.xml format
# ---------------------------------------------------------------------------
def write_xml(xml_path, nodes, hexes, quads,
              mu, kappa, t1, ninc,
              clamped_marker=MARKER_CLAMPED,
              shear_marker=MARKER_SHEAR):

    nodes_sorted = sorted(nodes.items())
    n_nodes = len(nodes_sorted)
    n_hexes = len(hexes)
    all_ids = [nid for nid, _ in nodes_sorted]

    # Node sets by face marker
    clamped_nodes = set()
    shear_quads   = []
    for _, conn, marker in quads:
        if marker == clamped_marker:
            clamped_nodes.update(conn)
        if marker == shear_marker:
            shear_quads.append(conn)

    lines = []

    # Header comment
    lines.append("<!--")
    lines.append("      Nonlinear elasticity test problem")
    lines.append("      Cook's problem")
    lines.append("-->")
    lines.append('<?xml version="1.0"?>')

    # --- elasticity block ---
    lines.append('<elasticity type="THREE_DIM">')
    lines.append('  <parameters> ')
    lines.append('    <material>NeoHookean</material>')
    lines.append(f'    <coefficients>{mu}, {kappa}, {kappa}</coefficients>')
    lines.append(f'    <ninc>{ninc}</ninc>')
    lines.append('  </parameters>')
    lines.append('')

    # Neumann: shear on right face
    lines.append('  <neumann>')
    lines.append(f'    <node id="1" marker="{shear_marker}" t0="0.0" t1="{t1}" t2="0.0" />')
    lines.append('  </neumann>')
    lines.append('')

    # prescribed_displacement — match cook_hex*_k1.xml ordering:
    # first all clamped nodes dir 0 and dir 1, then all nodes dir 2
    lines.append('  <prescribed_displacement>')

    for nid in sorted(clamped_nodes):
        lines.append(f'    <node id="{nid}" direction="0" value="0.000000" />')
        lines.append(f'    <node id="{nid}" direction="1" value="0.000000" />')

    for nid in sorted(all_ids):
        lines.append(f'    <node id="{nid}" direction="2" value="0.000000" />')

    lines.append('  </prescribed_displacement>')
    lines.append('</elasticity>')
    lines.append('')

    # --- mesh block ---
    lines.append('<mesh celltype="hexahedron" dim="3">')

    lines.append(f'  <nodes size="{n_nodes}">')
    for nid, (x, y, z) in nodes_sorted:
        lines.append(f'    <node id="{nid}" x="{x:.6f}" y="{y:.6f}" z="{z:.6f}" />')
    lines.append('  </nodes>')
    lines.append('')

    # Elements — apply hex reordering to match cook_hex*_k1.xml convention
    lines.append(f'  <elements size="{n_hexes}">')
    for i, (_, conn, _marker) in enumerate(hexes):
        reordered = reorder_hex(conn)
        verts = " ".join(f'v{j}="{v}"' for j, v in enumerate(reordered))
        lines.append(f'    <element id="{i}" {verts} />')
    lines.append('  </elements>')
    lines.append('')

    # Boundary — shear face quads with marker=10, reordered to match reference
    lines.append('  <boundary celltype="quadrilateral" dim="2">')
    for bnd_id, conn in enumerate(shear_quads):
        reordered = reorder_quad(conn, nodes)
        verts = " ".join(f'v{j}="{v}"' for j, v in enumerate(reordered))
        lines.append(f'    <element id="{bnd_id}" marker="{shear_marker}" {verts} />')
    lines.append('  </boundary>')
    lines.append('</mesh>')

    with open(xml_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"[3/3] Wrote {xml_path}")
    print(f"      Nodes          : {n_nodes}")
    print(f"      Hex elements   : {n_hexes}")
    print(f"      Clamped nodes  : {len(clamped_nodes)}  (dir 0,1 fixed, marker {clamped_marker})")
    print(f"      Plane strain   : all {n_nodes} nodes fixed in dir 2")
    print(f"      Shear faces    : {len(shear_quads)}  (marker {shear_marker}, t1={t1})")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print(f"\nCook's Membrane 3D mesh generator")
    print(f"  In-plane divisions : {N}x{N}  ({N*N} quads per layer)")
    print(f"  Thickness h        : {H}")
    print(f"  Load increments    : {NINC}")
    print(f"  Material           : NeoHookean  mu={MU}  kappa={KAPPA}")
    print(f"  Shear load t1      : {T1}")
    print(f"  Output dir         : {os.path.abspath(OUTDIR)}\n")

    write_geo(GEO_FILE, N, H)
    generate_mesh(GEO_FILE, MSH_FILE, N, H)

    nodes, hexes, quads = parse_msh(MSH_FILE)
    write_xml(XML_FILE, nodes, hexes, quads,
              mu=MU, kappa=KAPPA, t1=T1, ninc=NINC)

    print("\nDone!")
    for fp in [GEO_FILE, MSH_FILE, XML_FILE]:
        size = os.path.getsize(fp) if os.path.exists(fp) else 0
        print(f"  {fp}  ({size} bytes)")
