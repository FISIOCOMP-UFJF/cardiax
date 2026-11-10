"""
Linear-elasticity PATCH TEST (constant-stress test).

The most fundamental FE verification: impose a linear displacement field on
the boundary of a mesh whose interior nodes are deliberately distorted. A
correct element formulation must reproduce a CONSTANT stress/strain state
exactly (to machine precision), regardless of the interior distortion.

Setup (3D hex, THREE_DIM — matches the Elasticity solver schema):
  Unit cube [0,1]^3, 2x2x2 elements, with the 8 interior... (for a single
  interior node on a 2x2x2 grid) node(s) perturbed off the regular grid.
  We impose the linear field
        u = (a*x, b*y, c*z)         (here a uniaxial-like field)
  as prescribed displacement on ALL boundary nodes. The exact interior
  solution is the same linear field, giving constant strain
        eps_xx=a, eps_yy=b, eps_zz=c  (constant everywhere).

Verification: after solving, every element's stress must equal the constant
value D : eps, identical across elements. If interior nodes move to exactly
u(x) at their (distorted) positions, the patch test passes.

This generator writes the prescribed displacement on boundary nodes from the
linear field; you then check the computed interior displacement equals the
same field at those nodes.

Examples:
  python gen_patch_test.py
  python gen_patch_test.py --distort 0.15 --young 1.0 --poisson 0.3
"""

import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description="Linear-elasticity constant-stress patch test (3D hex).",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-L", "--length", type=float, default=1.0, help="cube side")
parser.add_argument("--distort", type=float, default=0.1,
                    help="interior-node perturbation (fraction of element size)")
parser.add_argument("--young", type=float, default=1.0, help="Young's modulus E")
parser.add_argument("--poisson", type=float, default=0.3, help="Poisson ratio nu")
# Linear field u = (a x, b y, c z)
parser.add_argument("--a", type=float, default=0.01, help="du_x/dx (strain xx)")
parser.add_argument("--b", type=float, default=0.0,  help="du_y/dy (strain yy)")
parser.add_argument("--c", type=float, default=0.0,  help="du_z/dz (strain zz)")
parser.add_argument("-o", "--output", type=str, default="patch_test.xml")
args = parser.parse_args()

L = args.length
nx = ny = nz = 2                      # 2x2x2 so there is one interior node
distort = args.distort
E, nu = args.young, args.poisson
a, b, c = args.a, args.b, args.c
output_file = args.output

nnx, nny, nnz = nx+1, ny+1, nz+1
n_nodes = nnx*nny*nnz
n_elements = nx*ny*nz

def node_id(i, j, k):
    return i + j*nnx + k*nnx*nny

def is_boundary(i, j, k):
    return (i in (0, nx)) or (j in (0, ny)) or (k in (0, nz))

# ── Nodes (interior nodes perturbed to make the test meaningful) ────────────────
rng = np.random.default_rng(42)
h = L / nx
nodes = []
coord = {}
for k in range(nnz):
    for j in range(nny):
        for i in range(nnx):
            x, y, z = i*L/nx, j*L/ny, k*L/nz
            if not is_boundary(i, j, k):
                # perturb interior node off the grid
                x += distort*h*(rng.random()-0.5)*2
                y += distort*h*(rng.random()-0.5)*2
                z += distort*h*(rng.random()-0.5)*2
            nid = node_id(i, j, k)
            nodes.append((nid, x, y, z))
            coord[nid] = (x, y, z)

# ── Elements ────────────────────────────────────────────────────────────────────
elements = []
eid = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            v = [node_id(i,j,k), node_id(i+1,j,k), node_id(i+1,j+1,k), node_id(i,j+1,k),
                 node_id(i,j,k+1), node_id(i+1,j,k+1), node_id(i+1,j+1,k+1), node_id(i,j+1,k+1)]
            elements.append((eid, *v))
            eid += 1

# ── Prescribed displacement on ALL boundary nodes: u = (a x, b y, c z) ──────────
prescribed = []
for k in range(nnz):
    for j in range(nny):
        for i in range(nnx):
            if is_boundary(i, j, k):
                nid = node_id(i, j, k)
                x, y, z = coord[nid]
                ux, uy, uz = a*x, b*y, c*z
                prescribed.append((nid, 0, ux))
                prescribed.append((nid, 1, uy))
                prescribed.append((nid, 2, uz))

# ── Write XML (Elasticity solver schema: young/poisson, prescribed_displacement) ─
# exact constant strain/stress for the imposed linear field u=(a x, b y, c z)
lam = E*nu/((1+nu)*(1-2*nu))
mu  = E/(2*(1+nu))
exx, eyy, ezz = a, b, c
tr = exx+eyy+ezz
sxx = lam*tr + 2*mu*exx
syy = lam*tr + 2*mu*eyy
szz = lam*tr + 2*mu*ezz

lines = []
lines.append("<!--")
lines.append("  ===================================================================")
lines.append("   LINEAR ELASTICITY BENCHMARK 5 — CONSTANT-STRESS PATCH TEST")
lines.append("  ===================================================================")
lines.append("")
lines.append("   PROBLEM DEFINITION")
lines.append("     The fundamental FE verification test. A unit cube is meshed with")
lines.append("     2x2x2 hexahedra and its single interior node is deliberately")
lines.append("     perturbed off the regular grid. A linear displacement field")
lines.append("         u = (a*x, b*y, c*z)")
lines.append("     is imposed on EVERY boundary node. A correct element must then")
lines.append("     reproduce a CONSTANT strain/stress state throughout, regardless")
lines.append("     of the interior distortion.")
lines.append(f"       E = {E}, nu = {nu}")
lines.append(f"       a = {a} (strain_xx), b = {b} (strain_yy), c = {c} (strain_zz)")
lines.append(f"       interior node perturbed by {distort}*h")
lines.append("")
lines.append("   BOUNDARY CONDITIONS")
lines.append("     all boundary nodes : prescribed u = (a*x, b*y, c*z)")
lines.append("")
lines.append("   EXACT SOLUTION  (constant everywhere)")
lines.append("     Lame parameters:  lambda = %.6g,  mu = %.6g" % (lam, mu))
lines.append("     Strain:  eps_xx=%.6g  eps_yy=%.6g  eps_zz=%.6g" % (exx,eyy,ezz))
lines.append("     Stress:  sig_xx=%.6g  sig_yy=%.6g  sig_zz=%.6g" % (sxx,syy,szz))
lines.append("     PASS criteria (to machine precision):")
lines.append("       - interior node displaces to u=(a*x,b*y,c*z) at its position")
lines.append("       - every element reports the constant stress above")
lines.append("  -->")
lines.append('<?xml version="1.0"?>')
lines.append('<elasticity type="THREE_DIM">')
lines.append('  <parameters>')
lines.append(f'    <young>{E}</young>')
lines.append(f'    <poisson>{nu}</poisson>')
lines.append('  </parameters>')
lines.append('')
lines.append('  <prescribed_displacement>')
for (nid, d, val) in prescribed:
    lines.append(f'    <node id="{nid}" direction="{d}" value="{val:.8f}" />')
lines.append('  </prescribed_displacement>')
lines.append('</elasticity>')
lines.append('')
lines.append('<mesh celltype="hexahedron" dim="3">')
lines.append(f'  <nodes size="{n_nodes}">')
for (nid, x, y, z) in nodes:
    lines.append(f'    <node id="{nid}" x="{x:.8f}" y="{y:.8f}" z="{z:.8f}" />')
lines.append('  </nodes>')
lines.append('')
lines.append(f'  <elements size="{n_elements}">')
for e in elements:
    eid = e[0]; v = e[1:]
    lines.append(f'    <element id="{eid}" v0="{v[0]}" v1="{v[1]}" v2="{v[2]}" v3="{v[3]}" '
                 f'v4="{v[4]}" v5="{v[5]}" v6="{v[6]}" v7="{v[7]}" />')
lines.append('  </elements>')
lines.append('')
lines.append('  <boundary celltype="quadrilateral" dim="2">')
lines.append('  </boundary>')
lines.append('</mesh>')

with open(output_file, "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"Generated: {output_file}")
print(f"  Nodes: {n_nodes}  Elements: {n_elements}")
print(f"  Boundary prescribed-displacement entries: {len(prescribed)}")
interior = [node_id(i,j,k) for k in range(nnz) for j in range(nny) for i in range(nnx)
            if not is_boundary(i,j,k)]
print(f"  Interior (free) nodes: {interior}")
for nid in interior:
    x,y,z = coord[nid]
    print(f"    node {nid} at ({x:.4f},{y:.4f},{z:.4f}) -> expected u=({a*x:.5f},{b*y:.5f},{c*z:.5f})")
