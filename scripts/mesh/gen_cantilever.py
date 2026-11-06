"""
Linear-elasticity CANTILEVER BEAM under end shear (Timoshenko benchmark).

A beam clamped at x=0, length L along x, depth D along y, thin in z.
A shear traction in the -y direction is applied on the free-end face x=L
(Neumann, marker 10). The classical Timoshenko end-loaded cantilever has an
exact tip deflection

    delta_tip = P L^3 / (3 E I)  +  shear correction

with I = D^3 * t / 12 (t = thickness in z), P = total end shear force.
For a thin beam (L >> D) the bending term dominates and delta ~ P L^3/(3 E I).

Setup (3D hex, THREE_DIM — matches the Elasticity solver schema):
  - x = 0 face: fully clamped (u_x=u_y=u_z=0) via prescribed_displacement
  - x = L face: Neumann shear traction t1 = -(P / end_area) in y (marker 10)
  - z faces left free (thin slab approximates plane stress)

Examples:
  python gen_cantilever.py
  python gen_cantilever.py --nx 40 --ny 8 --nz 1 -L 10 -D 1 --load 1.0
  python gen_cantilever.py --young 1000 --poisson 0.3 --load 1.0
"""

import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description="Linear-elasticity Timoshenko cantilever beam (3D hex).",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-L", "--length", type=float, default=10.0, help="beam length (x)")
parser.add_argument("-D", "--depth",  type=float, default=1.0,  help="beam depth (y)")
parser.add_argument("-t", "--thick",  type=float, default=1.0,  help="thickness (z)")
parser.add_argument("--nx", type=int, default=40, help="elements along length")
parser.add_argument("--ny", type=int, default=8,  help="elements along depth")
parser.add_argument("--nz", type=int, default=1,  help="elements along thickness")
parser.add_argument("--young",   type=float, default=1000.0, help="Young's modulus E")
parser.add_argument("--poisson", type=float, default=0.3,    help="Poisson ratio nu")
parser.add_argument("--load", type=float, default=1.0,
                    help="total end shear force P (applied in -y as traction)")
parser.add_argument("-o", "--output", type=str, default="cantilever.xml")
args = parser.parse_args()

L, D, t = args.length, args.depth, args.thick
nx, ny, nz = args.nx, args.ny, args.nz
E, nu = args.young, args.poisson
P = args.load
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
            nodes.append((node_id(i,j,k), i*L/nx, j*D/ny, k*t/nz))

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

# ── Boundary quads on the free end x=L (marker 10) for the shear traction ─────────
end_faces = []
bid = 0
for k in range(nz):
    for j in range(ny):
        bv0 = node_id(nx, j,   k  )
        bv1 = node_id(nx, j+1, k  )
        bv2 = node_id(nx, j+1, k+1)
        bv3 = node_id(nx, j,   k+1)
        end_faces.append((bid, 10, bv0, bv1, bv2, bv3))
        bid += 1

# Shear traction magnitude: total force P spread over the end area (D*t), in -y.
end_area = D * t
ty = -P / end_area

# End-face nodes carry the Neumann traction (marker 10)
end_nodes = [node_id(nx, j, k) for k in range(nnz) for j in range(nny)]

# ── Prescribed displacement: clamp x=0 face fully ─────────────────────────────────
prescribed = []
for k in range(nnz):
    for j in range(nny):
        nid = node_id(0, j, k)
        prescribed.append((nid, 0, 0.0))
        prescribed.append((nid, 1, 0.0))
        prescribed.append((nid, 2, 0.0))

# ── Analytical reference (printed, not in the file) ───────────────────────────────
I = (D**3) * t / 12.0
delta_bending = P * L**3 / (3.0 * E * I)

# ── Write XML ──────────────────────────────────────────────────────────────────────
delta_shear = 0.0
# Timoshenko end-loaded cantilever, exact Euler-Bernoulli tip deflection with
# an optional shear correction. Shear part ~ 3 P L / (5 G A) for a rectangle
# (kappa=5/6), G = E/(2(1+nu)), A = D*t. Reported for reference.
G = E / (2.0 * (1.0 + nu))
A_sec = D * t
delta_shear = 1.2 * P * L / (G * A_sec)       # 6/5 * PL/(GA)
delta_total = delta_bending + delta_shear

lines = []
lines.append("<!--")
lines.append("  ===================================================================")
lines.append("   LINEAR ELASTICITY BENCHMARK 1 — CANTILEVER BEAM (end shear load)")
lines.append("  ===================================================================")
lines.append("")
lines.append("   PROBLEM DEFINITION")
lines.append("     A straight beam of length L (x), depth D (y), thickness t (z),")
lines.append("     clamped at x=0, loaded by a transverse shear force P at the free")
lines.append("     end x=L (applied as a uniform traction ty = -P/(D*t) in -y).")
lines.append(f"       L = {L}, D = {D}, t = {t}")
lines.append(f"       E = {E}, nu = {nu}")
lines.append(f"       P = {P}  (total end shear force)")
lines.append(f"       mesh = {nx} x {ny} x {nz} hex elements")
lines.append("")
lines.append("   BOUNDARY CONDITIONS")
lines.append("     x = 0 face : clamped, u_x = u_y = u_z = 0  (prescribed_displacement)")
lines.append("     x = L face : Neumann shear traction  ty = -P/(D*t)  (marker 10)")
lines.append("")
lines.append("   EXACT / REFERENCE SOLUTION  (Timoshenko)")
lines.append("     Second moment of area   I = D^3 t / 12 = %.6g" % I)
lines.append("     Bending tip deflection  d_b = P L^3 / (3 E I)      = %.6g" % delta_bending)
lines.append("     Shear correction        d_s = 6/5 * P L /(G A)     = %.6g" % delta_shear)
lines.append("     Total tip deflection    d   = d_b + d_s            = %.6g" % delta_total)
lines.append("     (For slender beams L>>D the bending term dominates; compare the")
lines.append("      computed free-end y-displacement against d above.)")
lines.append("  -->")
lines.append('<?xml version="1.0"?>')
lines.append('<elasticity type="THREE_DIM">')
lines.append('  <parameters>')
lines.append(f'    <young>{E}</young>')
lines.append(f'    <poisson>{nu}</poisson>')
lines.append('  </parameters>')
lines.append('')
# Neumann is defined ONCE per marker: the traction (force per unit area) is
# integrated over every boundary element carrying that marker. A single entry
# is enough (the solver keys neumann_map by marker); one <node> per face node
# would be redundant. This matches the working Cook's-membrane convention.
lines.append('  <neumann>')
lines.append(f'    <node id="0" marker="10" t0="0.0" t1="{ty:.8f}" t2="0.0" />')
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
for (bid, marker, bv0, bv1, bv2, bv3) in end_faces:
    lines.append(f'    <element id="{bid}" marker="{marker}" '
                 f'v0="{bv0}" v1="{bv1}" v2="{bv2}" v3="{bv3}" />')
lines.append('  </boundary>')
lines.append('</mesh>')

with open(output_file, "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"Generated: {output_file}")
print(f"  Nodes: {n_nodes}  Elements: {n_elements}")
print(f"  End (Neumann) nodes: {len(end_nodes)}  shear ty={ty:.6g}")
print(f"  Clamped prescribed-displacement entries: {len(prescribed)}")
print(f"  I = {I:.6g}")
print(f"  Reference tip deflection (bending only) = {delta_bending:.6g}")
