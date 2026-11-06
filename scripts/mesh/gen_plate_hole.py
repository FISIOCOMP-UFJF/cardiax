"""
Linear-elasticity PLATE WITH A CIRCULAR HOLE under uniaxial tension
(Kirsch stress-concentration benchmark).

A large plate with a central circular hole (radius R) under remote tension S
has the famous stress concentration factor of 3 at the hole edge
(sigma_theta = 3 S at the point perpendicular to the load). Kirsch's solution
gives the full stress field.

We model ONE QUARTER by symmetry (hole at the origin corner):
  domain  : quarter plate [0, W] x [0, W], minus the quarter disk r < R
  mesh    : structured "O-grid" — nr radial x nt angular cells mapping the
            quarter annulus between the hole (r=R) and the outer square.
  tension : applied in +x on the right edge x=W (Neumann, marker 10)
  symmetry:
     - left cut  (x=0, the y-axis):  u_x = 0
     - bottom cut (y=0, the x-axis): u_y = 0
  thin in z (one layer), u_z pinned on z=0 to stabilise out-of-plane.

The mesh is a single hex layer (plane-stress-like slab). The blend from the
circular hole boundary to the straight outer square is done by linearly
interpolating, for each angular ray, between the point on the hole and the
point on the outer square boundary.

Examples:
  python gen_plate_hole.py
  python gen_plate_hole.py --R 0.5 --W 5 --nr 10 --nt 12 --tension 1.0
"""

import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description="Linear-elasticity plate-with-hole (Kirsch) quarter model (3D hex).",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--R", type=float, default=0.5, help="hole radius")
parser.add_argument("--W", type=float, default=5.0, help="quarter-plate half-width")
parser.add_argument("--thick", type=float, default=0.1, help="thickness (z)")
parser.add_argument("--nr", type=int, default=10, help="radial elements (hole->edge)")
parser.add_argument("--nt", type=int, default=12, help="angular elements (0..90 deg)")
parser.add_argument("--young",   type=float, default=1000.0, help="Young's modulus E")
parser.add_argument("--poisson", type=float, default=0.3,    help="Poisson ratio nu")
parser.add_argument("--tension", type=float, default=1.0,
                    help="remote tension S applied on x=W edge (+x traction)")
parser.add_argument("-o", "--output", type=str, default="plate_hole.xml")
args = parser.parse_args()

R, W, t = args.R, args.W, args.thick
nr, nt = args.nr, args.nt
E, nu = args.young, args.poisson
S = args.tension
output_file = args.output

nnr, nnt, nnz = nr+1, nt+1, 2        # one hex layer in z
n_nodes = nnr * nnt * nnz
n_elements = nr * nt * 1

def node_id(ir, it, iz):
    return ir + it*nnr + iz*nnr*nnt

def outer_point(theta):
    """Point on the outer square boundary of the quarter [0,W]^2 along the ray
    at angle theta from +x axis. The ray hits either the right edge x=W or the
    top edge y=W, whichever is closer."""
    c, s = np.cos(theta), np.sin(theta)
    # distance to x=W plane along the ray: W/c ; to y=W: W/s
    cand = []
    if c > 1e-12: cand.append(W / c)
    if s > 1e-12: cand.append(W / s)
    d = min(cand)
    return d*c, d*s

# ── Nodes: for each angular ray, blend from hole (r=R) to outer square ──────────
nodes = []
coord = {}
for iz in range(nnz):
    z = iz * t
    for it in range(nnt):
        theta = (np.pi/2.0) * it / nt          # 0 .. 90 deg
        xin, yin = R*np.cos(theta), R*np.sin(theta)       # on the hole
        xout, yout = outer_point(theta)                   # on the square
        for ir in range(nnr):
            f = ir / nr
            x = xin + f*(xout - xin)
            y = yin + f*(yout - yin)
            nid = node_id(ir, it, iz)
            nodes.append((nid, x, y, z))
            coord[nid] = (x, y, z)

# ── Elements ──────────────────────────────────────────────────────────────────────
elements = []
eid = 0
for it in range(nt):
    for ir in range(nr):
        v0 = node_id(ir,   it,   0)
        v1 = node_id(ir+1, it,   0)
        v2 = node_id(ir+1, it+1, 0)
        v3 = node_id(ir,   it+1, 0)
        v4 = node_id(ir,   it,   1)
        v5 = node_id(ir+1, it,   1)
        v6 = node_id(ir+1, it+1, 1)
        v7 = node_id(ir,   it+1, 1)
        elements.append((eid, v0, v1, v2, v3, v4, v5, v6, v7))
        eid += 1

# ── Outer boundary faces on the x=W edge (marker 10) for the tension ──────────────
# The outer ring is ir=nr. A face there lies on the square boundary; we tag only
# those whose outward side is on x=W (theta <= 45 deg rays reach x=W).
# We apply tension as +x traction on ALL outer faces, but it is physically the
# remote tension; for the quarter model, tagging the x=W portion is standard.
tension_faces = []
bid = 0
for it in range(nt):
    # midpoint angle of this angular cell
    th_mid = (np.pi/2.0) * (it + 0.5) / nt
    xout, yout = outer_point(th_mid)
    on_right = abs(xout - W) < 1e-9     # this cell's outer edge is on x=W
    bv0 = node_id(nr, it,   0)
    bv1 = node_id(nr, it+1, 0)
    bv2 = node_id(nr, it+1, 1)
    bv3 = node_id(nr, it,   1)
    marker = 10 if on_right else 20     # 10 = loaded (x=W), 20 = top edge (free/other)
    tension_faces.append((bid, marker, bv0, bv1, bv2, bv3))
    bid += 1

# nodes on x=W outer edge get +x traction
tension_nodes = []
for it in range(nnt):
    nid = node_id(nr, it, 0)
    x, y, z = coord[nid]
    if abs(x - W) < 1e-9:
        tension_nodes.append(node_id(nr, it, 0))
        tension_nodes.append(node_id(nr, it, 1))

# ── Symmetry + stabilisation BCs ─────────────────────────────────────────────────
prescribed = []
TOL = 1e-9
for nid,(x,y,z) in coord.items():
    if abs(x) < TOL:                 # left cut (y-axis): u_x = 0
        prescribed.append((nid, 0, 0.0))
    if abs(y) < TOL:                 # bottom cut (x-axis): u_y = 0
        prescribed.append((nid, 1, 0.0))
    if abs(z) < TOL:                 # z=0 layer: u_z = 0 (thin-slab stabilisation)
        prescribed.append((nid, 2, 0.0))

# ── Analytical reference ──────────────────────────────────────────────────────────
# Kirsch: hoop stress at hole edge, theta=90deg (top of hole) = 3S (tension in x).
sigma_max = 3.0 * S

# ── Write XML ──────────────────────────────────────────────────────────────────────
lines = []
lines.append("<!--")
lines.append("  ===================================================================")
lines.append("   LINEAR ELASTICITY BENCHMARK 4 — PLATE WITH A HOLE (Kirsch)")
lines.append("  ===================================================================")
lines.append("")
lines.append("   PROBLEM DEFINITION")
lines.append("     A large plate with a central circular hole of radius R, under")
lines.append("     remote uniaxial tension S in the x-direction. By symmetry only")
lines.append("     one quarter is modelled (hole centred at the origin corner),")
lines.append("     meshed as a structured O-grid from the hole edge (r=R) out to")
lines.append("     the straight outer square boundary [0,W]x[0,W]. One hex layer")
lines.append("     in z (plane-stress-like slab).")
lines.append(f"       R = {R}, W = {W} (quarter half-width), thickness = {t}")
lines.append(f"       E = {E}, nu = {nu}")
lines.append(f"       S = {S}  (remote tension, +x)")
lines.append(f"       O-grid = {nr} radial x {nt} angular")
lines.append("")
lines.append("   BOUNDARY CONDITIONS")
lines.append("     x = 0 cut (y-axis) : symmetry, u_x = 0")
lines.append("     y = 0 cut (x-axis) : symmetry, u_y = 0")
lines.append("     x = W edge         : Neumann tension  tx = S  (marker 10)")
lines.append("     z = 0 layer        : u_z = 0  (thin-slab out-of-plane stabilisation)")
lines.append("     hole edge r=R      : traction-free")
lines.append("")
lines.append("   EXACT / REFERENCE SOLUTION  (Kirsch, 1898)")
lines.append("     For an INFINITE plate, the hoop stress on the hole boundary is")
lines.append("       sigma_theta(R, theta) = S (1 - 2 cos 2theta)")
lines.append("     giving the maximum at theta = 90deg (top of the hole, perpendicular")
lines.append("     to the load):")
lines.append("       STRESS CONCENTRATION FACTOR Kt = 3")
lines.append("       sigma_theta_max = 3 S = %.6g" % sigma_max)
lines.append("     and the minimum at theta = 0 (sides, along the load): -S.")
lines.append("     (W/R finite here, so expect Kt slightly below 3; it approaches 3")
lines.append("      as W/R grows. Compare sig_xx at the hole top node vs 3S.)")
lines.append("  -->")
lines.append('<?xml version="1.0"?>')
lines.append('<elasticity type="THREE_DIM">')
lines.append('  <parameters>')
lines.append(f'    <young>{E}</young>')
lines.append(f'    <poisson>{nu}</poisson>')
lines.append('  </parameters>')
lines.append('')
# Neumann defined once per marker; traction integrated over the marked faces
# (here the x=W outer edge, marker 10).
lines.append('  <neumann>')
lines.append(f'    <node id="0" marker="10" t0="{S:.8f}" t1="0.0" t2="0.0" />')
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
for (bid, marker, bv0, bv1, bv2, bv3) in tension_faces:
    lines.append(f'    <element id="{bid}" marker="{marker}" '
                 f'v0="{bv0}" v1="{bv1}" v2="{bv2}" v3="{bv3}" />')
lines.append('  </boundary>')
lines.append('</mesh>')

with open(output_file, "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"Generated: {output_file}")
print(f"  Nodes: {n_nodes}  Elements: {n_elements}")
print(f"  Tension (Neumann, x=W) nodes: {len(tension_nodes)}")
print(f"  Prescribed-displacement entries: {len(prescribed)}")
print(f"  Reference stress concentration factor = 3  ->  sigma_max = {sigma_max:.6g}")
