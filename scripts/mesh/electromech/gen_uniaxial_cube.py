"""
Generates an XML input file for uniaxial traction in the X direction
on a structured hexahedral cube mesh.

Cube: [0, L] x [0, L] x [0, L]
Mesh: nx x ny x nz hexahedral elements

Boundary conditions (uniaxial traction, direction X):
  - Face x=0 (left):  u_x = 0  (symmetry / fixed in X)
                       u_y = 0  (one node, to remove rigid body rotation)
                       u_z = 0  (one node, to remove rigid body rotation)
  - Face x=L (right): Neumann traction t0 applied in X direction (marker=10)

All nodes on x=0 face have u_x = 0.
One node on x=0 face also has u_y = 0 and u_z = 0 to prevent rigid body modes.
"""

import numpy as np

# ── Configuration ─────────────────────────────────────────────────────────────
L   = 1.0   # cube side length
nx  = 3     # elements in x
ny  = 3     # elements in y
nz  = 3     # elements in z
traction = 10.0   # traction magnitude in X direction (Neumann)

# Material (NeoHookean: mu, lambda, ?)
mu     = 80.1938
lam    = 40094.0
extra  = 40094.0
ninc   = 20

output_file = "uniaxial_cube.xml"
# ──────────────────────────────────────────────────────────────────────────────

# Number of nodes in each direction
nnx = nx + 1
nny = ny + 1
nnz = nz + 1
n_nodes    = nnx * nny * nnz
n_elements = nx * ny * nz

# Node index helper
def node_id(i, j, k):
    """i: x-index, j: y-index, k: z-index"""
    return i + j * nnx + k * nnx * nny

# ── Generate nodes ─────────────────────────────────────────────────────────────
nodes = []
for k in range(nnz):
    for j in range(nny):
        for i in range(nnx):
            x = i * L / nx
            y = j * L / ny
            z = k * L / nz
            nodes.append((node_id(i, j, k), x, y, z))

# ── Generate hexahedral elements ───────────────────────────────────────────────
# Node ordering matches the reference file:
# v0..v3 = bottom face (z=k), v4..v7 = top face (z=k+1)
# within each face: (i,j), (i+1,j), (i+1,j+1), (i,j+1)  [counter-clockwise from below]
elements = []
eid = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            v0 = node_id(i,   j,   k  )
            v1 = node_id(i+1, j,   k  )
            v2 = node_id(i+1, j+1, k  )
            v3 = node_id(i,   j+1, k  )
            v4 = node_id(i,   j,   k+1)
            v5 = node_id(i+1, j,   k+1)
            v6 = node_id(i+1, j+1, k+1)
            v7 = node_id(i,   j+1, k+1)
            elements.append((eid, v0, v1, v2, v3, v4, v5, v6, v7))
            eid += 1

# ── Boundary conditions ────────────────────────────────────────────────────────

# Face x=0: all nodes with i=0  →  fix u_x=0
left_nodes = [node_id(0, j, k) for k in range(nnz) for j in range(nny)]

# Face x=L: all nodes with i=nx  →  Neumann traction in X
right_nodes = [node_id(nx, j, k) for k in range(nnz) for j in range(nny)]

# Boundary quad faces on x=L face (Neumann marker=10)
# Each quad: (i=nx, j, k) → (i=nx, j+1, k) → (i=nx, j+1, k+1) → (i=nx, j, k+1)
boundary_faces = []
bid = 0
for k in range(nz):
    for j in range(ny):
        bv0 = node_id(nx, j,   k  )
        bv1 = node_id(nx, j+1, k  )
        bv2 = node_id(nx, j+1, k+1)
        bv3 = node_id(nx, j,   k+1)
        boundary_faces.append((bid, 10, bv0, bv1, bv2, bv3))
        bid += 1

# Prescribed displacements:
#   - All nodes on x=0: u_x = 0  (direction 0)
#   - One corner node on x=0: u_y = 0, u_z = 0  (remove rigid body modes)
anchor = node_id(0, 0, 0)  # corner at origin

prescribed = []
# Fix u_x = 0 on entire left face
for nid in left_nodes:
    prescribed.append((nid, 0, 0.0))  # direction 0 = x
# Fix u_y and u_z on one node to remove rigid body rotation
prescribed.append((anchor, 1, 0.0))  # direction 1 = y
prescribed.append((anchor, 2, 0.0))  # direction 2 = z

# ── Write XML ──────────────────────────────────────────────────────────────────
lines = []

lines.append("<!--")
lines.append("  Uniaxial traction test — cube mesh")
lines.append(f"  Cube: [0,{L}]^3  |  Mesh: {nx}x{ny}x{nz} hex elements")
lines.append("  Traction applied on face x=L in the X direction")
lines.append("-->")
lines.append('<?xml version="1.0"?>')
lines.append('<elasticity type="THREE_DIM">')
lines.append('  <parameters>')
lines.append(f'    <material>NeoHookean</material>')
lines.append(f'    <coefficients>{mu}, {lam}, {extra}</coefficients>')
lines.append(f'    <ninc>{ninc}</ninc>')
lines.append('  </parameters>')
lines.append('')

# Neumann: one entry per boundary node (traction in x = direction 0)
lines.append('  <neumann>')
for nid in right_nodes:
    lines.append(f'    <node id="{nid}" marker="10" t0="{traction:.6f}" t1="0.0" t2="0.0" />')
lines.append('  </neumann>')
lines.append('')

# Prescribed displacements
lines.append('  <prescribed_displacement>')
for (nid, direction, value) in prescribed:
    lines.append(f'    <node id="{nid}" direction="{direction}" value="{value:.6f}" />')
lines.append('  </prescribed_displacement>')
lines.append('</elasticity>')
lines.append('<!-- </problem> -->')
lines.append('')

# Mesh block
lines.append('<mesh celltype="hexahedron" dim="3">')
lines.append(f'  <nodes size="{n_nodes}">')
for (nid, x, y, z) in nodes:
    lines.append(f'    <node id="{nid}" x="{x:.6f}" y="{y:.6f}" z="{z:.6f}" />')
lines.append('  </nodes>')
lines.append('')
lines.append(f'  <elements size="{n_elements}">')
for (eid, v0, v1, v2, v3, v4, v5, v6, v7) in elements:
    lines.append(f'    <element id="{eid}" v0="{v0}" v1="{v1}" v2="{v2}" v3="{v3}" '
                 f'v4="{v4}" v5="{v5}" v6="{v6}" v7="{v7}" />')
lines.append('  </elements>')
lines.append('')
lines.append('  <element_data type="fiber_isotropic">')
lines.append('  </element_data>')
lines.append('')

# Boundary faces (Neumann surface)
lines.append('  <boundary celltype="quadrilateral" dim="2">')
for (bid, marker, bv0, bv1, bv2, bv3) in boundary_faces:
    lines.append(f'    <element id="{bid}" marker="{marker}" '
                 f'v0="{bv0}" v1="{bv1}" v2="{bv2}" v3="{bv3}" />')
lines.append('  </boundary>')
lines.append('</mesh>')

with open(output_file, "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"Generated: {output_file}")
print(f"  Nodes   : {n_nodes}")
print(f"  Elements: {n_elements}")
print(f"  Left face (x=0) nodes fixed in X: {len(left_nodes)}")
print(f"  Right face (x=L) nodes with traction: {len(right_nodes)}")
print(f"  Boundary quad faces (Neumann): {len(boundary_faces)}")
print(f"  Anchor node (rigid body removal): {anchor}")
