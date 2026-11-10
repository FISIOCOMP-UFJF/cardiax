"""
Poisson Problem Mesh Generator
================================
Generates structured meshes for the Poisson problem with exact solution
on the unit square (2D) or unit cube (3D). No external dependencies.

Supported cell types:
  2D: quad (quadrilateral), tri (triangle)
  3D: hex  (hexahedron),    tet (tetrahedron)

XML format produced matches reference files exactly:
  unitsquare_quad*.xml  ->  celltype="quadrilateral"  dim="2"
  unitsquare_tri*.xml   ->  celltype="triangle"        dim="2"
  unitcube_tet*.xml     ->  celltype="tetrahedron"     dim="3"
  (hex)                 ->  celltype="hexahedron"      dim="3"

Boundary markers:
  2D: marker=7   (all 4 edges,   celltype="line"         dim="1")
  3D: marker=27  (all 6 faces,   celltype="triangle"     dim="2"  for tet)
                                  celltype="quadrilateral" dim="2"  for hex)

<poisson> block: single Dirichlet node id=0 with the boundary marker.

Usage:
  python poisson_mesh.py --type quad --n 5
  python poisson_mesh.py --type tri  --n 5
  python poisson_mesh.py --type hex  --n 4
  python poisson_mesh.py --type tet  --n 4
  python poisson_mesh.py --type quad --n 10 --out my_mesh.xml --out_dir ./results

  --type   : quad | tri | hex | tet  (required)
  --n      : divisions per edge (default 5)
  --out    : output filename (default: auto, e.g. unitsquare_quad_n5.xml)
  --out_dir: output directory (default: current dir)

Element counts:
  quad : n^2        (e.g. n=5 -> 25)
  tri  : 2*n^2      (e.g. n=5 -> 50)
  hex  : n^3        (e.g. n=4 -> 64)
  tet  : 6*n^3      (e.g. n=4 -> 384)

Examples:

  # 2D meshes
  python unit_cube_poisson.py --type quad --n 5    # 25 quads,  36 nodes
  python unit_cube_poisson.py --type tri  --n 5    # 50 tris,   36 nodes

  # 3D meshes
  python unit_cube_poisson.py --type hex  --n 4    # 64 hexes,  125 nodes
  python unit_cube_poisson.py --type tet  --n 4    # 384 tets,  125 nodes

  # Custom output
  python poisson_mesh.py --type quad --n 10 --out_dir ./meshes
"""

import argparse
import os

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Poisson problem mesh generator (no dependencies)")
parser.add_argument("--type",    required=True, choices=["quad","tri","hex","tet"],
                    help="Element type: quad | tri | hex | tet")
parser.add_argument("--n",       type=int, default=5,    help="Divisions per edge (default 5)")
parser.add_argument("--out",     type=str, default=None, help="Output filename (default: auto)")
parser.add_argument("--out_dir", type=str, default=".",  help="Output directory (default: .)")
args = parser.parse_args()

ETYPE  = args.type
N      = args.n
OUTDIR = args.out_dir
os.makedirs(OUTDIR, exist_ok=True)

# Auto filename matching reference naming
if args.out:
    XML_FILE = os.path.join(OUTDIR, args.out)
else:
    domain = "unitsquare" if ETYPE in ("quad","tri") else "unitcube"
    XML_FILE = os.path.join(OUTDIR, f"{domain}_{ETYPE}_n{N}.xml")

# Boundary markers matching reference files
MARKER_2D = 7
MARKER_3D = 27


# ===========================================================================
# 2D: QUAD mesh  [0,1]^2
# ===========================================================================
def make_quad_mesh(n):
    """
    Structured n x n quad mesh on [0,1]^2.
    Node ordering: row-major  id = j*(n+1) + i
    Element ordering: (BL, BR, TR, TL) for each cell
    Boundary: all 4 edges as line segments, marker=7
    """
    # Nodes
    nodes = {}
    for j in range(n+1):
        for i in range(n+1):
            nodes[j*(n+1)+i] = (i/n, j/n, 0.0)

    def idx(i,j): return j*(n+1)+i

    # Elements: quads BL->BR->TR->TL
    elements = []
    for j in range(n):
        for i in range(n):
            elements.append((idx(i,j), idx(i+1,j), idx(i+1,j+1), idx(i,j+1)))

    # Boundary edges (all 4 sides), marker=7
    boundary = []
    for i in range(n):   boundary.append((idx(i,0),   idx(i+1,0)))    # bottom
    for j in range(n):   boundary.append((idx(n,j),   idx(n,j+1)))    # right
    for i in range(n):   boundary.append((idx(n-i,n), idx(n-i-1,n)))  # top
    for j in range(n):   boundary.append((idx(0,n-j), idx(0,n-j-1)))  # left

    return (nodes, elements, boundary,
            "quadrilateral", 2, "line", 1, MARKER_2D)


# ===========================================================================
# 2D: TRI mesh  [0,1]^2
# ===========================================================================
def make_tri_mesh(n):
    """
    Structured n x n quad mesh, each quad split into 2 triangles -> 2*n^2 tris.
    Split: (BL,BR,TR) + (BL,TR,TL)  matches reference files.
    Boundary: all 4 edges as line segments, marker=7
    """
    nodes = {}
    for j in range(n+1):
        for i in range(n+1):
            nodes[j*(n+1)+i] = (i/n, j/n, 0.0)

    def idx(i,j): return j*(n+1)+i

    elements = []
    for j in range(n):
        for i in range(n):
            v0 = idx(i,   j)
            v1 = idx(i+1, j)
            v2 = idx(i+1, j+1)
            v3 = idx(i,   j+1)
            elements.append((v0, v1, v2))   # lower-right triangle
            elements.append((v0, v2, v3))   # upper-left  triangle

    boundary = []
    for i in range(n):   boundary.append((idx(i,0),   idx(i+1,0)))
    for j in range(n):   boundary.append((idx(n,j),   idx(n,j+1)))
    for i in range(n):   boundary.append((idx(n-i,n), idx(n-i-1,n)))
    for j in range(n):   boundary.append((idx(0,n-j), idx(0,n-j-1)))

    return (nodes, elements, boundary,
            "triangle", 2, "line", 1, MARKER_2D)


# ===========================================================================
# 3D: HEX mesh  [0,1]^3
# ===========================================================================
def make_hex_mesh(n):
    """
    Structured n x n x n hex mesh on [0,1]^3.
    Node id = k*(n+1)^2 + j*(n+1) + i
    Hex8 node ordering: bottom face (v0-v3) then top face (v4-v7)
    Boundary: all 6 faces as quads, marker=27
    """
    nodes = {}
    for k in range(n+1):
        for j in range(n+1):
            for i in range(n+1):
                nodes[k*(n+1)**2 + j*(n+1) + i] = (i/n, j/n, k/n)

    def idx(i,j,k): return k*(n+1)**2 + j*(n+1) + i

    elements = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                elements.append((
                    idx(i,   j,   k),   # v0
                    idx(i+1, j,   k),   # v1
                    idx(i+1, j+1, k),   # v2
                    idx(i,   j+1, k),   # v3
                    idx(i,   j,   k+1), # v4
                    idx(i+1, j,   k+1), # v5
                    idx(i+1, j+1, k+1), # v6
                    idx(i,   j+1, k+1), # v7
                ))

    boundary = []
    # Bottom k=0
    for j in range(n):
        for i in range(n):
            boundary.append((idx(i,j,0), idx(i+1,j,0), idx(i+1,j+1,0), idx(i,j+1,0)))
    # Top k=n
    for j in range(n):
        for i in range(n):
            boundary.append((idx(i,j,n), idx(i,j+1,n), idx(i+1,j+1,n), idx(i+1,j,n)))
    # Front j=0
    for k in range(n):
        for i in range(n):
            boundary.append((idx(i,0,k), idx(i+1,0,k), idx(i+1,0,k+1), idx(i,0,k+1)))
    # Back j=n
    for k in range(n):
        for i in range(n):
            boundary.append((idx(i,n,k), idx(i,n,k+1), idx(i+1,n,k+1), idx(i+1,n,k)))
    # Left i=0
    for k in range(n):
        for j in range(n):
            boundary.append((idx(0,j,k), idx(0,j,k+1), idx(0,j+1,k+1), idx(0,j+1,k)))
    # Right i=n
    for k in range(n):
        for j in range(n):
            boundary.append((idx(n,j,k), idx(n,j+1,k), idx(n,j+1,k+1), idx(n,j,k+1)))

    return (nodes, elements, boundary,
            "hexahedron", 3, "quadrilateral", 2, MARKER_3D)


# ===========================================================================
# 3D: TET mesh  [0,1]^3
# ===========================================================================
def make_tet_mesh(n):
    """
    Structured n x n x n hex mesh, each hex split into 6 tetrahedra.
    Uses the Freudenthal/Kuhn triangulation (consistent across shared faces):
      Given hex corners a=v0..g=v6 (same as hex ordering):
        tet0: (a, b, c, g)
        tet1: (a, c, d, g)
        tet2: (a, d, h, g)
        tet3: (a, e, f, g)
        tet4: (a, e, g, h)
        tet5: (a, b, f, g)
    Boundary: all 6 faces split into 2 triangles each, marker=27
    """
    nodes = {}
    for k in range(n+1):
        for j in range(n+1):
            for i in range(n+1):
                nodes[k*(n+1)**2 + j*(n+1) + i] = (i/n, j/n, k/n)

    def idx(i,j,k): return k*(n+1)**2 + j*(n+1) + i

    elements = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                a = idx(i,   j,   k)
                b = idx(i+1, j,   k)
                c = idx(i+1, j+1, k)
                d = idx(i,   j+1, k)
                e = idx(i,   j,   k+1)
                f = idx(i+1, j,   k+1)
                g = idx(i+1, j+1, k+1)
                h = idx(i,   j+1, k+1)
                elements.append((a, b, c, g))
                elements.append((a, c, d, g))
                elements.append((a, d, h, g))
                elements.append((a, e, f, g))
                elements.append((a, e, g, h))
                elements.append((a, b, f, g))

    # Boundary triangles (each face quad -> 2 tris), marker=27
    boundary = []

    def add_face_tris(face_quads):
        for (p0, p1, p2, p3) in face_quads:
            boundary.append((p0, p1, p2))
            boundary.append((p0, p2, p3))

    # Bottom k=0
    add_face_tris([(idx(i,j,0), idx(i+1,j,0), idx(i+1,j+1,0), idx(i,j+1,0))
                   for j in range(n) for i in range(n)])
    # Top k=n
    add_face_tris([(idx(i,j,n), idx(i,j+1,n), idx(i+1,j+1,n), idx(i+1,j,n))
                   for j in range(n) for i in range(n)])
    # Front j=0
    add_face_tris([(idx(i,0,k), idx(i+1,0,k), idx(i+1,0,k+1), idx(i,0,k+1))
                   for k in range(n) for i in range(n)])
    # Back j=n
    add_face_tris([(idx(i,n,k), idx(i,n,k+1), idx(i+1,n,k+1), idx(i+1,n,k))
                   for k in range(n) for i in range(n)])
    # Left i=0
    add_face_tris([(idx(0,j,k), idx(0,j,k+1), idx(0,j+1,k+1), idx(0,j+1,k))
                   for k in range(n) for j in range(n)])
    # Right i=n
    add_face_tris([(idx(n,j,k), idx(n,j+1,k), idx(n,j+1,k+1), idx(n,j,k+1))
                   for k in range(n) for j in range(n)])

    return (nodes, elements, boundary,
            "tetrahedron", 3, "triangle", 2, MARKER_3D)


# ===========================================================================
# XML Writer — matching reference format exactly
# ===========================================================================
def write_xml(path, nodes, elements, boundary,
              celltype, dim, bnd_celltype, bnd_dim, marker):

    nodes_sorted = sorted(nodes.items())
    lines = []
    lines.append("<!--")
    lines.append("      Problem with analytical solution for testing Poisson solver.")
    lines.append("-->")
    lines.append('<?xml version="1.0"?>')
    lines.append("")
    lines.append(f'<mesh celltype="{celltype}" dim="{dim}">')

    # nodes
    lines.append(f'  <nodes size="{len(nodes_sorted)}">')
    for nid, (x, y, z) in nodes_sorted:
        lines.append(f'    <node id="{nid}" x="{x:.6f}" y="{y:.6f}" z="{z:.6f}" />')
    lines.append('  </nodes>')
    lines.append('')

    # elements
    lines.append(f'  <elements size="{len(elements)}">')
    for i, conn in enumerate(elements):
        verts = " ".join(f'v{j}="{v}"' for j, v in enumerate(conn))
        lines.append(f'    <element id="{i}" {verts} />')
    lines.append('  </elements>')
    lines.append('')

    # boundary
    lines.append(f'  <boundary celltype="{bnd_celltype}" dim="{bnd_dim}">')
    for i, conn in enumerate(boundary):
        verts = " ".join(f'v{j}="{v}"' for j, v in enumerate(conn))
        lines.append(f'    <element id="{i}" marker="{marker}" {verts} />')
    lines.append('  </boundary>')
    lines.append('</mesh>')
    lines.append('')

    # poisson block
    lines.append('<poisson>')
    lines.append('  <dirichlet>')
    lines.append(f'    <node id="0" marker="{marker}" value="0.0" />')
    lines.append('  </dirichlet>')
    lines.append('</poisson>')

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ===========================================================================
# Main
# ===========================================================================
GENERATORS = {
    "quad": make_quad_mesh,
    "tri":  make_tri_mesh,
    "hex":  make_hex_mesh,
    "tet":  make_tet_mesh,
}

EXPECTED = {
    "quad": lambda n: f"{n}x{n} = {n**2} quads",
    "tri":  lambda n: f"{n}x{n}x2 = {2*n**2} triangles",
    "hex":  lambda n: f"{n}x{n}x{n} = {n**3} hexes",
    "tet":  lambda n: f"{n}x{n}x{n}x6 = {6*n**3} tets",
}

DOMAIN = {
    "quad": "2D unit square [0,1]^2",
    "tri":  "2D unit square [0,1]^2",
    "hex":  "3D unit cube  [0,1]^3",
    "tet":  "3D unit cube  [0,1]^3",
}

if __name__ == "__main__":
    print(f"\nPoisson mesh generator  (no external dependencies)")
    print(f"  Type   : {ETYPE}  —  {DOMAIN[ETYPE]}")
    print(f"  n      : {N}  ->  {EXPECTED[ETYPE](N)}")
    print(f"  Output : {XML_FILE}")

    nodes, elements, boundary, celltype, dim, bnd_ct, bnd_dim, marker = GENERATORS[ETYPE](N)

    write_xml(XML_FILE, nodes, elements, boundary,
              celltype, dim, bnd_ct, bnd_dim, marker)

    size = os.path.getsize(XML_FILE)
    print(f"\n  Nodes          : {len(nodes)}")
    print(f"  Elements       : {len(elements)}")
    print(f"  Boundary elems : {len(boundary)}  (marker={marker})")
    print(f"  File size      : {size:,} bytes")
    print(f"\nDone!")
