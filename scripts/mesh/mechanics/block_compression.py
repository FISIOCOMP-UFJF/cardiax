"""
Nearly Incompressible Block Under Compression
===============================================
Reference: Fig. 4 — Block under compression.

Geometry (quarter model):
  Full block L x L x H = 1 mm x 1 mm x 1 mm (default L=2 so quarter = 1x1x1 cube)
  Quarter: X1 in [0, L/2],  X2 in [0, L/2],  X3 in [0, H]

Boundary Conditions:
  X1=0 plane   (marker 10): symmetry  -> dir 0 = 0
  X2=0 plane   (marker 20): symmetry  -> dir 1 = 0
  X3=0 plane   (marker 30): bottom    -> dir 2 = 0
  X3=H plane   (marker 40): Neumann   -> pressure p0 in -X3 (t2 = -p0)

XML format matches cook_hex*_k1.xml exactly:
  <problem>
    <elasticity type="THREE_DIM">
      <parameters>
        <material>NeoHookean</material>
        <coefficients>mu, kappa, kappa</coefficients>
        <ninc>N</ninc>
      </parameters>
      <neumann>
        <node id="1" marker="40" t0="0.0" t1="0.0" t2="-p0" />
      </neumann>
      <prescribed_displacement>
        <!-- sym X1=0: direction 0 -->
        <!-- sym X2=0: direction 1 -->
        <!-- bottom X3=0: direction 2 -->
      </prescribed_displacement>
    </elasticity>
    <mesh celltype="hexahedron" dim="3">
      <nodes> ... </nodes>
      <elements> ... </elements>
      <element_data type="fiber_isotropic"> </element_data>
      <boundary celltype="quadrilateral" dim="2">
        <!-- top face quads only (marker 40) -->
      </boundary>
    </mesh>
  </problem>

Material (Neo-Hookean):
  mu    = 80.194   N/mm^2
  kappa = 400889.806  N/mm^2

Load:
  p0 = 4 N/mm^2 on top face (-X3 direction)

Usage:
  pip install gmsh
  python block_compression.py [--n 2] [--L 2.0] [--H 1.0] [--p0 4.0] [--ninc 10] [--full] [--out_dir .]

  --n      : divisions per edge (default 2 -> 2x2x2 elements in quarter)
  --L      : full block side length (default 2.0 -> quarter is 1x1x1 cube)
  --H      : block height (default 1.0)
  --p0     : applied pressure N/mm^2 (default 4.0)
  --ninc   : load increments for solver (default 10)
  --full   : mesh full block instead of quarter
  --out_dir: output directory (default .)

Refinement guide:
  --n 2  ->   8 elements (quarter)
  --n 4  ->  64 elements
  --n 8  -> 512 elements
  --n 16 -> 4096 elements
"""

import argparse
import os
import sys

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Nearly incompressible block under compression")
parser.add_argument("--n",       type=int,   default=2,         help="Divisions per edge (default 2)")
parser.add_argument("--L",       type=float, default=2.0,       help="Full block side length mm (default 2.0 -> quarter=1x1x1)")
parser.add_argument("--H",       type=float, default=1.0,       help="Block height mm (default 1.0)")
parser.add_argument("--p0",      type=float, default=4.0,       help="Pressure on top N/mm^2 (default 4.0)")
parser.add_argument("--ninc",    type=int,   default=10,        help="Load increments (default 10)")
parser.add_argument("--full",    action="store_true",           help="Mesh full block (default: quarter symmetry)")
parser.add_argument("--out_dir", type=str,   default=".",       help="Output directory (default: .)")
args = parser.parse_args()

N      = args.n
L      = args.L
H      = args.H
P0     = args.p0
NINC   = args.ninc
FULL   = args.full
OUTDIR = args.out_dir
os.makedirs(OUTDIR, exist_ok=True)

# Quarter dimensions (L/2 x L/2 x H)
# With default L=2, H=1: Lq=1=H -> perfect cube
Lq = L if FULL else L / 2.0
Hq = H

# Material
MU    = 80.194
KAPPA = 400889.806

# Physical markers
MARKER_SYM_X1 = 10   # X1=0 symmetry -> dir 0 = 0
MARKER_SYM_X2 = 20   # X2=0 symmetry -> dir 1 = 0
MARKER_BOTTOM  = 30   # X3=0 bottom   -> dir 2 = 0
MARKER_TOP     = 40   # X3=H top      -> Neumann pressure
MARKER_FREE_X1 = 50   # X1=Lq free
MARKER_FREE_X2 = 60   # X2=Lq free
MARKER_DOMAIN  = 100  # volume

suffix = "full" if FULL else "quarter"
GEO_FILE = os.path.join(OUTDIR, f"block_compression_{suffix}_n{N}.geo")
MSH_FILE = os.path.join(OUTDIR, f"block_compression_{suffix}_n{N}.msh")
XML_FILE = os.path.join(OUTDIR, f"block_compression_{suffix}_n{N}.xml")


# ---------------------------------------------------------------------------
# 1. Write .geo
# ---------------------------------------------------------------------------
def write_geo(path, n, lq, hq, full=False):
    lc = lq / n
    model = "Full block" if full else "Quarter model (X1>=0, X2>=0)"
    geo = f"""// ============================================================
//  Nearly Incompressible Block Under Compression
//  {model}
//
//  Block: [{lq:.4f} x {lq:.4f} x {hq:.4f}]
//  (Full: L={lq*2 if not full else lq:.4f}, H={hq:.4f})
//  Divisions: n={n}  (element size lc={lc:.6f})
//
//  Physical markers:
//   {MARKER_SYM_X1}   sym_X1       (X1=0,  dir 0 = 0)
//   {MARKER_SYM_X2}   sym_X2       (X2=0,  dir 1 = 0)
//   {MARKER_BOTTOM}   bottom       (X3=0,  dir 2 = 0)
//   {MARKER_TOP}   top          (X3=H,  Neumann p0)
//   {MARKER_FREE_X1}   free_X1_pos  (X1=Lq, free)
//   {MARKER_FREE_X2}   free_X2_pos  (X2=Lq, free)
//   {MARKER_DOMAIN}  domain       (volume)
// ============================================================

lc = {lc:.6f};
n  = {n};

// Bottom face corners (Z=0)
Point(1) = {{ 0,   0,   0, lc }};   // origin
Point(2) = {{ {lq}, 0,   0, lc }};   // X1=Lq, X2=0
Point(3) = {{ {lq}, {lq}, 0, lc }};   // X1=Lq, X2=Lq
Point(4) = {{ 0,   {lq}, 0, lc }};   // X1=0,  X2=Lq

// Bottom face edges
Line(1) = {{1, 2}};   // X2=0 bottom   -> sym_X2 face
Line(2) = {{2, 3}};   // X1=Lq bottom  -> free_X1 face
Line(3) = {{3, 4}};   // X2=Lq bottom  -> free_X2 face
Line(4) = {{4, 1}};   // X1=0 bottom   -> sym_X1 face

Curve Loop(1)    = {{1, 2, 3, 4}};
Plane Surface(1) = {{1}};

Transfinite Curve{{1}} = n+1;
Transfinite Curve{{2}} = n+1;
Transfinite Curve{{3}} = n+1;
Transfinite Curve{{4}} = n+1;
Transfinite Surface{{1}};
Recombine Surface{{1}};

// Extrude in X3 (Z direction)
out[] = Extrude {{0, 0, {hq}}} {{
    Surface{{1}};
    Layers{{n}};
    Recombine;
}};

// out[0]=top(Z=H), out[1]=vol, out[2]=lat(l1/X2=0), out[3]=lat(l2/X1=Lq),
// out[4]=lat(l3/X2=Lq), out[5]=lat(l4/X1=0)

Physical Volume("domain",      {MARKER_DOMAIN}) = {{out[1]}};
Physical Surface("bottom",     {MARKER_BOTTOM}) = {{1}};
Physical Surface("top",        {MARKER_TOP})    = {{out[0]}};
Physical Surface("sym_X2",     {MARKER_SYM_X2}) = {{out[2]}};
Physical Surface("free_X1",    {MARKER_FREE_X1})= {{out[3]}};
Physical Surface("free_X2",    {MARKER_FREE_X2})= {{out[4]}};
Physical Surface("sym_X1",     {MARKER_SYM_X1}) = {{out[5]}};
"""
    with open(path, "w") as f:
        f.write(geo)
    print(f"[1/3] Wrote {path}")


# ---------------------------------------------------------------------------
# 2. Generate mesh via GMSH Python API
# ---------------------------------------------------------------------------
def generate_mesh(geo_path, msh_path, n, lq, hq, full=False):
    try:
        import gmsh
    except ImportError:
        print("ERROR: pip install gmsh")
        sys.exit(1)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm",    8)
    gmsh.option.setNumber("Mesh.Algorithm3D",  4)

    gmsh.model.add("block_compression")
    fac = gmsh.model.geo
    lc  = lq / n

    p1 = fac.addPoint(0.0, 0.0, 0.0, lc)
    p2 = fac.addPoint(lq,  0.0, 0.0, lc)
    p3 = fac.addPoint(lq,  lq,  0.0, lc)
    p4 = fac.addPoint(0.0, lq,  0.0, lc)

    l1 = fac.addLine(p1, p2)   # X2=0
    l2 = fac.addLine(p2, p3)   # X1=Lq
    l3 = fac.addLine(p3, p4)   # X2=Lq
    l4 = fac.addLine(p4, p1)   # X1=0

    cl  = fac.addCurveLoop([l1, l2, l3, l4])
    sur = fac.addPlaneSurface([cl])

    for line in [l1, l2, l3, l4]:
        fac.mesh.setTransfiniteCurve(line, n + 1)
    fac.mesh.setTransfiniteSurface(sur)
    fac.mesh.setRecombine(2, sur)
    fac.synchronize()

    out = fac.extrude([(2, sur)], 0, 0, hq, numElements=[n], recombine=True)
    fac.synchronize()

    top_tag  = out[0][1]
    vol_tag  = out[1][1]
    lat_sym2 = out[2][1]   # X2=0 (from l1)
    lat_fx1  = out[3][1]   # X1=Lq (from l2)
    lat_fx2  = out[4][1]   # X2=Lq (from l3)
    lat_sym1 = out[5][1]   # X1=0 (from l4)

    gmsh.model.addPhysicalGroup(3, [vol_tag],  MARKER_DOMAIN);  gmsh.model.setPhysicalName(3, MARKER_DOMAIN,  "domain")
    gmsh.model.addPhysicalGroup(2, [sur],      MARKER_BOTTOM);  gmsh.model.setPhysicalName(2, MARKER_BOTTOM,  "bottom")
    gmsh.model.addPhysicalGroup(2, [top_tag],  MARKER_TOP);     gmsh.model.setPhysicalName(2, MARKER_TOP,     "top")
    gmsh.model.addPhysicalGroup(2, [lat_sym2], MARKER_SYM_X2);  gmsh.model.setPhysicalName(2, MARKER_SYM_X2,  "sym_X2")
    gmsh.model.addPhysicalGroup(2, [lat_fx1],  MARKER_FREE_X1); gmsh.model.setPhysicalName(2, MARKER_FREE_X1, "free_X1")
    gmsh.model.addPhysicalGroup(2, [lat_fx2],  MARKER_FREE_X2); gmsh.model.setPhysicalName(2, MARKER_FREE_X2, "free_X2")
    gmsh.model.addPhysicalGroup(2, [lat_sym1], MARKER_SYM_X1);  gmsh.model.setPhysicalName(2, MARKER_SYM_X1,  "sym_X1")

    gmsh.model.mesh.generate(3)
    gmsh.write(msh_path)
    gmsh.finalize()
    print(f"[2/3] Wrote {msh_path}")


# ---------------------------------------------------------------------------
# 3. Parse .msh (format 2 or 4)
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
    node_lines = _block(content, "Nodes").splitlines()
    idx = 0
    nb  = int(node_lines[idx].split()[0]); idx += 1
    for _ in range(nb):
        n_in = int(node_lines[idx].split()[3]); idx += 1
        ids  = [int(node_lines[idx+k])-1 for k in range(n_in)]; idx += n_in
        for nid in ids:
            p = node_lines[idx].split(); idx += 1
            nodes[nid] = (float(p[0]), float(p[1]), float(p[2]))
    elem_lines = _block(content, "Elements").splitlines()
    idx = 0
    nb  = int(elem_lines[idx].split()[0]); idx += 1
    eid = 0
    for _ in range(nb):
        hdr    = elem_lines[idx].split(); idx += 1
        entity = int(hdr[1])
        etype  = int(hdr[2])
        n_el   = int(hdr[3])
        for _ in range(n_el):
            p    = elem_lines[idx].split(); idx += 1
            conn = [int(x)-1 for x in p[1:]]
            if   etype == 5: hexes.append((eid, conn, entity))
            elif etype == 3: quads.append((eid, conn, entity))
            eid += 1
    return nodes, hexes, quads


# ---------------------------------------------------------------------------
# 4. Write XML — matching cook_hex*_k1.xml format exactly
# ---------------------------------------------------------------------------
def write_xml(xml_path, nodes, hexes, quads,
              mu, kappa, p0, ninc, full=False):

    nodes_sorted = sorted(nodes.items())
    n_nodes = len(nodes_sorted)
    n_hexes = len(hexes)

    # Collect node sets per marker
    sym_x1_nodes = set()
    sym_x2_nodes = set()
    bottom_nodes = set()
    top_quads    = []

    for _, conn, marker in quads:
        if marker == MARKER_SYM_X1:  sym_x1_nodes.update(conn)
        if marker == MARKER_SYM_X2:  sym_x2_nodes.update(conn)
        if marker == MARKER_BOTTOM:  bottom_nodes.update(conn)
        if marker == MARKER_TOP:     top_quads.append(conn)

    lines = []

    # Header — same style as reference files
    lines.append("<!--")
    lines.append("      Nonlinear elasticity test problem")
    lines.append("      Nearly incompressible block under compression")
    lines.append("-->")
    lines.append('<?xml version="1.0"?>')
    lines.append("<problem>")

    # --- elasticity ---
    lines.append('<elasticity type="THREE_DIM">')
    lines.append('  <parameters> ')
    lines.append('    <material>NeoHookean</material>')
    lines.append(f'    <coefficients>{mu}, {kappa}, {kappa}</coefficients>')
    lines.append(f'    <ninc>{ninc}</ninc>')
    lines.append('  </parameters>')
    lines.append('')

    # Neumann: pressure on top face in -X3 direction (t2 = -p0)
    lines.append('  <neumann>')
    lines.append(f'    <node id="1" marker="{MARKER_TOP}" t0="0.0" t1="0.0" t2="-{p0}" />')
    lines.append('  </neumann>')
    lines.append('')

    # prescribed_displacement
    lines.append('  <prescribed_displacement>')

    if not full:
        # Symmetry X1=0: fix direction 0
        for nid in sorted(sym_x1_nodes):
            lines.append(f'    <node id="{nid}" direction="0" value="0.000000" />')
        # Symmetry X2=0: fix direction 1
        for nid in sorted(sym_x2_nodes):
            lines.append(f'    <node id="{nid}" direction="1" value="0.000000" />')

    # Bottom X3=0: fix direction 2
    for nid in sorted(bottom_nodes):
        lines.append(f'    <node id="{nid}" direction="2" value="0.000000" />')

    lines.append('  </prescribed_displacement>')
    lines.append('</elasticity>')
    lines.append('')

    # --- mesh ---
    lines.append('<mesh celltype="hexahedron" dim="3">')

    lines.append(f'  <nodes size="{n_nodes}">')
    for nid, (x, y, z) in nodes_sorted:
        lines.append(f'    <node id="{nid}" x="{x:.6f}" y="{y:.6f}" z="{z:.6f}" />')
    lines.append('  </nodes>')
    lines.append('')

    lines.append(f'  <elements size="{n_hexes}">')
    for i, (_, conn, _marker) in enumerate(hexes):
        verts = " ".join(f'v{j}="{v}"' for j, v in enumerate(conn))
        lines.append(f'    <element id="{i}" {verts} />')
    lines.append('  </elements>')
    lines.append('')

    lines.append('  <element_data type="fiber_isotropic">')
    lines.append('  </element_data>')
    lines.append('')

    # boundary: top face quads only (marker 40), matching reference style
    lines.append('  <boundary celltype="quadrilateral" dim="2">')
    for i, conn in enumerate(top_quads):
        verts = " ".join(f'v{j}="{v}"' for j, v in enumerate(conn))
        lines.append(f'    <element id="{i}" marker="{MARKER_TOP}" {verts} />')
    lines.append('  </boundary>')
    lines.append('</mesh>')
    lines.append('</problem>')

    with open(xml_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"[3/3] Wrote {xml_path}")
    print(f"      Nodes          : {n_nodes}")
    print(f"      Hex elements   : {n_hexes}")
    if not full:
        print(f"      Sym X1 nodes   : {len(sym_x1_nodes)}  (dir 0 fixed, marker {MARKER_SYM_X1})")
        print(f"      Sym X2 nodes   : {len(sym_x2_nodes)}  (dir 1 fixed, marker {MARKER_SYM_X2})")
    print(f"      Bottom nodes   : {len(bottom_nodes)}  (dir 2 fixed, marker {MARKER_BOTTOM})")
    print(f"      Top quads      : {len(top_quads)}  (Neumann p0={p0}, marker {MARKER_TOP})")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print(f"\nNearly Incompressible Block Under Compression")
    print(f"  Model          : {'Full block' if FULL else 'Quarter symmetry'}")
    print(f"  Divisions      : {N}x{N}x{N}  ({N**3} elements)")
    print(f"  Quarter size   : {Lq} x {Lq} x {Hq}  (full: {L} x {L} x {H})")
    print(f"  Material       : NeoHookean  mu={MU}  kappa={KAPPA}")
    print(f"  Pressure p0    : {P0} N/mm^2  ({NINC} increments)")
    print(f"  Output dir     : {os.path.abspath(OUTDIR)}\n")

    write_geo(GEO_FILE, N, Lq, Hq, full=FULL)
    generate_mesh(GEO_FILE, MSH_FILE, N, Lq, Hq, full=FULL)

    nodes, hexes, quads = parse_msh(MSH_FILE)
    write_xml(XML_FILE, nodes, hexes, quads,
              mu=MU, kappa=KAPPA, p0=P0, ninc=NINC, full=FULL)

    print("\nDone!")
    for fp in [GEO_FILE, MSH_FILE, XML_FILE]:
        size = os.path.getsize(fp) if os.path.exists(fp) else 0
        print(f"  {fp}  ({size} bytes)")
