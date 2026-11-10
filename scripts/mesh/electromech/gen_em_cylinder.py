#!/usr/bin/env python3
"""
Thick-walled cylinder (hollow tube) with helical fibres, for testing the
coupled electromechanics solver (Electromechanic) before an LV geometry.

Geometry (micrometres, as in Cardiax):
    inner radius Ri, outer radius Ro, length L along z
    nr elements across the wall, nt around the circumference (multiple of 4),
    nz along the length. Hexahedra, same XML layout and local node order as
    coupled_electromech_Tomek.xml / make_em_slab.py.

Fibres (transversely isotropic), helix angle varying linearly across the wall:
    alpha = alpha_endo at r = Ri  ->  alpha_epi at r = Ro
    f = cos(alpha) e_theta + sin(alpha) e_z
    (0 deg = circumferential; defaults +60 endo, -60 epi, as in the LV)

Stimulus: a ring at the bottom end, z in [0, stim_height], whole cross-section.
The wave travels along z, and with the helical fibres the tube should shorten,
thicken its wall, narrow its cavity and twist.

Boundary conditions (statically determinate in-plane, so free to contract):
    bottom end z = 0                : u_z = 0          (symmetry plane / "base")
    bottom inner node at theta = 0  : u_y = 0
    bottom inner node at theta = pi : u_y = 0          (with the above: no
                                                        y translation, no rotation)
    bottom inner node at theta = pi/2: u_x = 0         (no x translation)
    top end: free.
Radial contraction and twist relative to the base satisfy these exactly, so
they remove the rigid-body motions without resisting the contraction.

No cavity pressure: the <pressure> entry is zero, as in the slab.

Usage:
    python3 make_em_cylinder.py                          # defaults below
    python3 make_em_cylinder.py --nt 64 --nz 16          # finer
    python3 make_em_cylinder.py --alpha-endo 0 --alpha-epi 0   # circumferential fibres
"""

import argparse
import math


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ri", type=float, default=2500.0, help="inner radius [um]")
    ap.add_argument("--ro", type=float, default=4000.0, help="outer radius [um]")
    ap.add_argument("--length", type=float, default=4000.0, help="length along z [um]")
    ap.add_argument("--nr", type=int, default=5, help="elements across the wall")
    ap.add_argument("--nt", type=int, default=48, help="elements around (multiple of 4)")
    ap.add_argument("--nz", type=int, default=12, help="elements along z")
    ap.add_argument("--alpha-endo", type=float, default=60.0, help="helix angle at r = Ri [deg]")
    ap.add_argument("--alpha-epi", type=float, default=-60.0, help="helix angle at r = Ro [deg]")
    ap.add_argument("--stim-height", type=float, default=400.0,
                    help="height of the stimulated ring at z = 0 [um]")
    ap.add_argument("--stim-value", type=float, default=1.0,
                    help="stimulus value (NP: positive depolarises; TT2/ToRORd: negative)")
    ap.add_argument("--stim-duration", type=float, default=2.0)
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()

    if a.nt % 4 != 0:
        raise SystemExit("--nt must be a multiple of 4 (supports at 0, pi/2, pi)")
    if not 0 < a.ri < a.ro:
        raise SystemExit("need 0 < ri < ro")

    nr, nt, nz = a.nr, a.nt, a.nz
    Ri, Ro, L = a.ri, a.ro, a.length
    out = a.out or f"em_cylinder_{nr}x{nt}x{nz}.xml"

    def nid(i, j, k):
        """i: radial 0..nr, j: circumferential 0..nt-1 (periodic), k: axial 0..nz"""
        return k * (nr + 1) * nt + (j % nt) * (nr + 1) + i

    def point(i, j, k):
        r = Ri + (Ro - Ri) * i / nr
        th = 2.0 * math.pi * j / nt
        return r * math.cos(th), r * math.sin(th), L * k / nz

    n_nodes = (nr + 1) * nt * (nz + 1)
    n_elems = nr * nt * nz

    # ---- boundary conditions --------------------------------------------
    bcs = set()
    for j in range(nt):
        for i in range(nr + 1):
            bcs.add((nid(i, j, 0), 2))
    bcs.add((nid(0, 0, 0), 1))             # theta = 0
    bcs.add((nid(0, nt // 2, 0), 1))       # theta = pi
    bcs.add((nid(0, nt // 4, 0), 0))       # theta = pi/2
    bcs = sorted(bcs, key=lambda b: (b[1], b[0]))

    # ---- elements and fibres ---------------------------------------------
    elements, fibres = [], []
    for k in range(nz):
        for j in range(nt):
            for i in range(nr):
                jn = j + 1
                # same local order as the slab: (r, theta, z) plays (x, y, z)
                v = [nid(i, j, k + 1), nid(i, j, k), nid(i, jn, k), nid(i, jn, k + 1),
                     nid(i + 1, j, k + 1), nid(i + 1, j, k), nid(i + 1, jn, k), nid(i + 1, jn, k + 1)]
                elements.append(v)

                t = (i + 0.5) / nr                                   # 0 endo -> 1 epi
                alpha = math.radians(a.alpha_endo + t * (a.alpha_epi - a.alpha_endo))
                th = 2.0 * math.pi * (j + 0.5) / nt
                fx = -math.cos(alpha) * math.sin(th)
                fy = math.cos(alpha) * math.cos(th)
                fz = math.sin(alpha)
                fibres.append((fx, fy, fz))

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
        for n, d in bcs:
            w(f'    <node id="{n}" direction="{d}" value="0.0" />\n')
        w('  </prescribed_displacement>\n')
        w('</elasticity>\n\n')

        # ---- mesh --------------------------------------------------------
        w('<mesh celltype="hexahedron" dim="3">\n')
        w(f'  <nodes size="{n_nodes}">\n')
        for k in range(nz + 1):
            for j in range(nt):
                for i in range(nr + 1):
                    x, y, z = point(i, j, k)
                    w(f'    <node id="{nid(i, j, k)}" x="{x:.6f}" y="{y:.6f}" z="{z:.6f}" />\n')
        w('  </nodes>\n')

        w(f'  <elements size="{n_elems}">\n')
        for e, v in enumerate(elements):
            attrs = " ".join(f'v{m}="{v[m]}"' for m in range(8))
            w(f'    <element id="{e}" {attrs}  />\n')
        w('  </elements>\n')

        w('  <element_data type="fiber_transversely_isotropic">\n')
        for e, (fx, fy, fz) in enumerate(fibres):
            w(f'    <element id="{e}">\n')
            w(f'        <fiber>{fx:.6f},{fy:.6f},{fz:.6f}</fiber>\n')
            w('    </element>\n')
        w('  </element_data>\n')
        w('</mesh>\n')

        # ---- electrophysiology -------------------------------------------
        eps = 1.0
        w('<electrophysiology>\n')
        w('  <stimuli number="1">\n')
        w(f'      <stim start="0.0" duration="{a.stim_duration:g}" value="{a.stim_value:g}"'
          f' x0="{-Ro - eps:.2f}" x1="{Ro + eps:.2f}"'
          f' y0="{-Ro - eps:.2f}" y1="{Ro + eps:.2f}"'
          f' z0="{-eps:.2f}" z1="{a.stim_height + eps:.2f}" />\n')
        w('  </stimuli>\n')
        w('</electrophysiology>\n')

    print(f"{out}: {nr}x{nt}x{nz} hexahedra, Ri = {Ri:g}, Ro = {Ro:g}, L = {L:g} um, "
          f"helix {a.alpha_endo:g} -> {a.alpha_epi:g} deg, "
          f"{n_nodes} nodes, {n_elems} elements, {len(bcs)} prescribed dofs")


if __name__ == "__main__":
    main()
