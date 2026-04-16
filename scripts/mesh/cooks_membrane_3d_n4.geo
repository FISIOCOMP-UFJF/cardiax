// ============================================================
//  Cook's Membrane — 3D GMSH geometry
//  n=4, h=1.0
//
//  Trapezoid corners at Z=0:
//    A=(0,0)   B=(48,44)   C=(48,60)   D=(0,44)
//  Extruded by h=1.0 in Z (1 layer -> plane strain).
//
//  Physical markers:
//   10   shear_load  (right face B-C, Neumann t1=6.25)
//   20   clamped     (left face  D-A, u0=u1=0)
//   100  domain      (volume)
// ============================================================

lc = 12.000000;
n  = 4;
h  = 1.0;

Point(1) = {  0,  0, 0, lc };   // A
Point(2) = { 48, 44, 0, lc };   // B
Point(3) = { 48, 60, 0, lc };   // C  (tip, reference point)
Point(4) = {  0, 44, 0, lc };   // D

Line(1) = {1, 2};   // bottom A->B
Line(2) = {2, 3};   // right  B->C  (shear)
Line(3) = {3, 4};   // top    C->D
Line(4) = {4, 1};   // left   D->A  (clamped)

Curve Loop(1)    = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve{1} = n+1;
Transfinite Curve{2} = n+1;
Transfinite Curve{3} = n+1;
Transfinite Curve{4} = n+1;
Transfinite Surface{1};
Recombine Surface{1};

// Extrude 1 layer in Z (plane strain / thin slab)
out[] = Extrude {0, 0, h} {
    Surface{1};
    Layers{1};
    Recombine;
};

// out[0]=back(Z=h), out[1]=vol, out[2]=lat(l1/bottom),
// out[3]=lat(l2/right/shear), out[4]=lat(l3/top), out[5]=lat(l4/left/clamped)

Physical Volume("domain",      100) = {out[1]};
Physical Surface("front",       50) = {1};
Physical Surface("back",        60) = {out[0]};
Physical Surface("bottom_edge", 30) = {out[2]};
Physical Surface("shear_load",  10) = {out[3]};
Physical Surface("top_edge",    40) = {out[4]};
Physical Surface("clamped",     20) = {out[5]};
