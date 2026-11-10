"""
Linear-elasticity SQUARE PLATE under uniform transverse load (bending).

A square plate of side A, thickness h, lying in the x-y plane, loaded by a
uniform pressure q on the top face (traction in -z). Edges are either clamped
or simply supported. Timoshenko plate theory gives the exact center
deflection:

    clamped (all four edges):          w_c = alpha_c * q A^4 / D_plate
    simply supported (all four edges): w_s = alpha_s * q A^4 / D_plate

with flexural rigidity  D_plate = E h^3 / (12 (1 - nu^2)),
and tabulated coefficients (Timoshenko & Woinowsky-Krieger):
    alpha_c ~ 0.00126   (clamped square)
    alpha_s ~ 0.00406   (simply supported square)

Setup (3D hex slab, THREE_DIM — matches the Elasticity solver schema):
  - top face z=h: uniform pressure q as Neumann traction t2 = -q (marker 10)
  - edges constrained via prescribed_displacement:
        clamped          -> u_x=u_y=u_z=0 on all four side faces
        simply_supported -> u_z=0 on the bottom edge lines of the four sides
  (A thin multi-layer hex slab approximates plate bending; use several
   elements through the thickness or a fine in-plane mesh for accuracy.)

Examples:
  python gen_plate_bending.py
  python gen_plate_bending.py --support clamped --na 20 --nz 2
  python gen_plate_bending.py --support simply --A 1 --h 0.01 --q 1.0
"""

import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description="Linear-elasticity square plate bending benchmark (3D hex slab).",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-A", "--side", type=float, default=1.0, help="plate side length")
parser.add_argument("--h", type=float, default=0.05, help="plate thickness")
parser.add_argument("--na", type=int, default=16, help="in-plane elements per side")
parser.add_argument("--nz", type=int, default=2,  help="elements through thickness")
parser.add_argument("--support", choices=["clamped", "simply"], default="clamped",
                    help="edge support condition")
parser.add_argument("--young",   type=float, default=1.0e4, help="Young's modulus E")
parser.add_argument("--poisson", type=float, default=0.3,   help="Poisson ratio nu")
parser.add_argument("--q", type=float, default=1.0, help="uniform transverse pressure")
parser.add_argument("-o", "--output", type=str, default="plate_bending.xml")
args = parser.parse_args()

A, h = args.side, args.h
nx = ny = args.na
nz = args.nz
support = args.support
E, nu = args.young, args.poisson
q = args.q
output_file = args.output

nnx, nny, nnz = nx+1, ny+1, nz+1
n_nodes = nnx*nny*nnz
n_elements = nx*ny*nz

def node_id(i, j, k):
    return i + j*nnx + k*nnx*nny

# ── Nodes ────────────────────────────────────────────────────────────────────────
nodes = []
for k in range(nnz):
    for j in range(nny):
        for i in range(nnx):
            nodes.append((node_id(i,j,k), i*A/nx, j*A/ny, k*h/nz))

# ── Elements ──────────────────────────────────────────────────────────────────────
elements = []
eid = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            v = [node_id(i,j,k), node_id(i+1,j,k), node_id(i+1,j+1,k), node_id(i,j+1,k),
                 node_id(i,j,k+1), node_id(i+1,j,k+1), node_id(i+1,j+1,k+1), node_id(i,j+1,k+1)]
            elements.append((eid, *v))
            eid += 1

# ── Top-face pressure quads (z=h, k=nz) marker 10, traction -z ───────────────────
top_faces = []
bid = 0
for j in range(ny):
    for i in range(nx):
        bv0 = node_id(i,   j,   nz)
        bv1 = node_id(i+1, j,   nz)
        bv2 = node_id(i+1, j+1, nz)
        bv3 = node_id(i,   j+1, nz)
        top_faces.append((bid, 10, bv0, bv1, bv2, bv3))
        bid += 1
top_nodes = [node_id(i, j, nz) for j in range(nny) for i in range(nnx)]

# ── Edge constraints ──────────────────────────────────────────────────────────────
def on_edge(i, j):
    return i in (0, nx) or j in (0, ny)

prescribed = []
seen = set()
if support == "clamped":
    # fully fix all nodes on the four side faces, through the thickness
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                if on_edge(i, j):
                    nid = node_id(i, j, k)
                    for d in (0, 1, 2):
                        prescribed.append((nid, d, 0.0))
else:  # simply supported: u_z = 0 along the four edges (all thickness layers),
       # plus minimal in-plane pins to remove rigid-body modes.
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                if on_edge(i, j):
                    nid = node_id(i, j, k)
                    prescribed.append((nid, 2, 0.0))   # u_z = 0 (support)
    # remove in-plane rigid modes: pin one corner in x,y and an adjacent in y
    c0 = node_id(0, 0, 0)
    c1 = node_id(nx, 0, 0)
    prescribed.append((c0, 0, 0.0))
    prescribed.append((c0, 1, 0.0))
    prescribed.append((c1, 1, 0.0))

# ── Analytical reference ──────────────────────────────────────────────────────────
Dp = E * h**3 / (12.0 * (1.0 - nu*nu))
alpha = 0.00126 if support == "clamped" else 0.00406
w_center = alpha * q * A**4 / Dp

# ── Write XML ──────────────────────────────────────────────────────────────────────
support_desc = ("all four edges fully clamped (u_x=u_y=u_z=0)"
                if support == "clamped"
                else "all four edges simply supported (u_z=0), plus 3 in-plane pins")

lines = []
lines.append("<!--")
lines.append("  ===================================================================")
lines.append("   LINEAR ELASTICITY BENCHMARK 6 — SQUARE PLATE BENDING")
lines.append("  ===================================================================")
lines.append("")
lines.append("   PROBLEM DEFINITION")
lines.append("     A square plate of side A and thickness h, lying in the x-y plane,")
lines.append("     loaded by a uniform transverse pressure q on its top face")
lines.append("     (traction tz = -q). Modelled as a thin 3D hex slab.")
lines.append(f"       A = {A}, h = {h}")
lines.append(f"       E = {E}, nu = {nu}")
lines.append(f"       q = {q}  (uniform pressure)")
lines.append(f"       support = {support}")
lines.append(f"       mesh = {nx} x {ny} x {nz} hex")
lines.append("")
lines.append("   BOUNDARY CONDITIONS")
lines.append(f"     edges : {support_desc}")
lines.append("     top face z=h : Neumann pressure  tz = -q  (marker 10)")
lines.append("")
lines.append("   EXACT / REFERENCE SOLUTION  (Timoshenko & Woinowsky-Krieger)")
lines.append("     Flexural rigidity  D = E h^3 / (12 (1-nu^2)) = %.6g" % Dp)
lines.append("     Center deflection  w_max = alpha * q A^4 / D")
lines.append("       clamped square        : alpha ~ 0.00126")
lines.append("       simply-supported square: alpha ~ 0.00406")
lines.append("     For this case (alpha = %.5f):  w_max = %.6g" % (alpha, w_center))
lines.append("     (Thin-plate theory; a 3D hex slab approximates it. Use a fine")
lines.append("      in-plane mesh and >=2 elements through h for a good match.")
lines.append("      Compare the computed u_z at the plate center against w_max.)")
lines.append("  -->")
lines.append('<?xml version="1.0"?>')
lines.append('<elasticity type="THREE_DIM">')
lines.append('  <parameters>')
lines.append(f'    <young>{E}</young>')
lines.append(f'    <poisson>{nu}</poisson>')
lines.append('  </parameters>')
lines.append('')
# Neumann defined once per marker; traction integrated over the marked faces.
lines.append('  <neumann>')
lines.append(f'    <node id="0" marker="10" t0="0.0" t1="0.0" t2="{-q:.8f}" />')
lines.append('  </neumann>')
lines.append('')
lines.append('  <prescribed_displacement>')
for (nid, d, val) in prescribed:
    lines.append(f'    <node id="{nid}" direction="{d}" value="{val:.6f}" />')
lines.append('  </prescribed_displacement>')
lines.append('</elasticity>')
lines.append('')
lines.append('<mesh celltype="hexahedron" dim="3">')
lines.append(f'  <nodes size="{n_nodes}">')
for (nid, x, y, z) in nodes:
    lines.append(f'    <node id="{nid}" x="{x:.6f}" y="{y:.6f}" z="{z:.6f}" />')
lines.append('  </nodes>')
lines.append('')
lines.append(f'  <elements size="{n_elements}">')
for e in elements:
    eid = e[0]; v = e[1:]
    lines.append(f'    <element id="{eid}" v0="{v[0]}" v1="{v[1]}" v2="{v[2]}" v3="{v[3]}" '
                 f'v4="{v[4]}" v5="{v[5]}" v6="{v[6]}" v7="{v[7]}" />')
lines.append('  </elements>')
lines.append('')
lines.append('  <element_data type="fiber_isotropic">')
lines.append('  </element_data>')
lines.append('')
lines.append('  <boundary celltype="quadrilateral" dim="2">')
for (bid, marker, bv0, bv1, bv2, bv3) in top_faces:
    lines.append(f'    <element id="{bid}" marker="{marker}" '
                 f'v0="{bv0}" v1="{bv1}" v2="{bv2}" v3="{bv3}" />')
lines.append('  </boundary>')
lines.append('</mesh>')

with open(output_file, "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"Generated: {output_file}")
print(f"  Nodes: {n_nodes}  Elements: {n_elements}")
print(f"  Support: {support}")
print(f"  Top pressure quads: {len(top_faces)}  Neumann nodes: {len(top_nodes)}")
print(f"  Prescribed-displacement entries: {len(prescribed)}")
print(f"  Reference center deflection w = {w_center:.6g}")
