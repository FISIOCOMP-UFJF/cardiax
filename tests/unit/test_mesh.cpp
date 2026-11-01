/*
 * Unit tests of src/mesh: reading meshes in every supported format, writing results,
 * and the helpers built on the connectivity.
 *
 * What is covered and why:
 *  - Mesh::read_xml for each element type: node, element and boundary counts,
 *    markers, AHA numbers and fibres are all input to the physics, so a parsing
 *    mistake here silently changes the problem being solved.
 *  - set_tet_positive_volume: an inverted tetrahedron must be reordered, otherwise
 *    the Jacobian determinant is negative and every integral on it changes sign.
 *  - Mesh::read (.pts/.elem/.fib) and GmshIO (.msh): the other two input formats.
 *  - BoxMesh: the structured generator.
 *  - PurkinjeMesh: detection of bifurcations and Purkinje-muscle junctions (PMJs).
 *  - comm_matrix/greedy_coloring: the colouring must never give the same colour to two
 *    elements sharing a node, since that is what would make parallel assembly race.
 *  - WriterHDF5: fields go to the right dataset and step, the checkpoint round trip
 *    returns exactly what was saved, and the XDMF file points to the HDF5 file.
 *  - ReaderHDF5 and Mesh::read_hdf5: reading back what WriterHDF5 wrote.
 *  - WriterVTK: the XML file is produced with the right sizes.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <filesystem>
#include <set>

#include "mesh/mesh.hpp"
#include "mesh/box_mesh.hpp"
#include "mesh/purkinje_mesh.hpp"
#include "mesh/gmshio.hpp"
#include "mesh/graph_coloring.hpp"
#include "mesh/writer_hdf5.hpp"
#include "mesh/writer_vtk.hpp"
#include "mesh/reader_hdf5.hpp"
#include "util/util.hpp"

using namespace cardiax_test;

namespace
{

double tet_volume(const Mesh & m, int e)
{
  std::vector<int> p;
  m.get_element_pt_nums(e, p);
  arma::mat33 J;
  for (int i = 0; i < 3; i++)
    J.col(i) = m.get_point(p[i + 1]) - m.get_point(p[0]);
  return arma::det(J) / 6.0;
}

BoxSpec small_box()
{
  BoxSpec b;
  b.nx = 3; b.ny = 2; b.nz = 1;
  b.Lx = 3.0; b.Ly = 2.0; b.Lz = 1.0;
  return b;
}

} // namespace

// ---- read_xml ---------------------------------------------------------------

TEST_CASE("read_xml: tetrahedra with boundary triangles and orthotropic fibres")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");

  CHECK_EQ(m.get_n_dim(), 3);
  CHECK_EQ(m.get_nen(), 4);
  CHECK_EQ(m.get_n_points(), uint(4 * 3 * 2));
  CHECK_EQ(m.get_n_elements(), 3 * 2 * 1 * 6);
  // Two triangles per boundary quad: 2*(3*2 + 3*1 + 2*1) quads
  CHECK_EQ(m.get_n_boundary_elements(), 2 * 2 * (6 + 3 + 2));

  const arma::vec3 & p = m.get_point(box_node(b, 3, 2, 1));
  CHECK_NEAR(p(0), 3.0, 1e-14);
  CHECK_NEAR(p(1), 2.0, 1e-14);
  CHECK_NEAR(p(2), 1.0, 1e-14);

  const Element & e = m.get_element(0);
  CHECK_NEAR(e.get_fiber()(0), 1.0, 0.0);
  CHECK_NEAR(e.get_trans()(1), 1.0, 0.0);
  CHECK_NEAR(e.get_normal()(2), 1.0, 0.0);
}

TEST_CASE("read_xml: hexahedra with boundary quads")
{
  BoxSpec b = small_box();
  write_hex_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  CHECK_EQ(m.get_nen(), 8);
  CHECK_EQ(m.get_n_elements(), 6);
  CHECK_EQ(m.get_n_boundary_elements(), 2 * (6 + 3 + 2));
  std::vector<int> pn;
  m.get_boundary_element_pt_nums(0, pn);
  CHECK_EQ(pn.size(), std::size_t(4));
}

TEST_CASE("read_xml: 2D triangles and quadrilaterals with boundary lines")
{
  BoxSpec b;
  b.nx = 2; b.ny = 3;
  write_tri_square("t.xml", b);
  write_quad_square("q.xml", b);

  Mesh t, q;
  t.read_xml("t.xml");
  q.read_xml("q.xml");
  CHECK_EQ(t.get_n_dim(), 2);
  CHECK_EQ(t.get_nen(), 3);
  CHECK_EQ(t.get_n_elements(), 12);
  CHECK_EQ(q.get_nen(), 4);
  CHECK_EQ(q.get_n_elements(), 6);
  CHECK_EQ(q.get_n_boundary_elements(), 10);
}

TEST_CASE("read_xml: 1D line mesh")
{
  write_cable("c.xml", 5, 10.0);
  Mesh m;
  m.read_xml("c.xml");
  CHECK_EQ(m.get_n_dim(), 1);
  CHECK_EQ(m.get_nen(), 2);
  CHECK_EQ(m.get_n_elements(), 5);
  CHECK_NEAR(m.get_point(5)(0), 10.0, 1e-14);
}

TEST_CASE("read_xml: element markers and AHA numbers, default 0 when absent")
{
  write_text("mk.xml",
    "<mesh celltype=\"hexahedron\" dim=\"3\"><nodes size=\"8\">"
    "<node id=\"0\" x=\"0\" y=\"0\" z=\"0\"/><node id=\"1\" x=\"1\" y=\"0\" z=\"0\"/>"
    "<node id=\"2\" x=\"1\" y=\"1\" z=\"0\"/><node id=\"3\" x=\"0\" y=\"1\" z=\"0\"/>"
    "<node id=\"4\" x=\"0\" y=\"0\" z=\"1\"/><node id=\"5\" x=\"1\" y=\"0\" z=\"1\"/>"
    "<node id=\"6\" x=\"1\" y=\"1\" z=\"1\"/><node id=\"7\" x=\"0\" y=\"1\" z=\"1\"/>"
    "</nodes><elements size=\"2\">"
    "<element id=\"0\" marker=\"3\" aha=\"7\" v0=\"0\" v1=\"1\" v2=\"2\" v3=\"3\" v4=\"4\" v5=\"5\" v6=\"6\" v7=\"7\"/>"
    "<element id=\"1\" v0=\"0\" v1=\"1\" v2=\"2\" v3=\"3\" v4=\"4\" v5=\"5\" v6=\"6\" v7=\"7\"/>"
    "</elements>"
    "<boundary celltype=\"quadrilateral\" dim=\"2\">"
    "<element id=\"0\" marker=\"30\" v0=\"0\" v1=\"3\" v2=\"2\" v3=\"1\"/>"
    "</boundary></mesh>");
  Mesh m;
  m.read_xml("mk.xml");
  CHECK_EQ(m.get_element_index(0), 3);
  CHECK_EQ(m.get_element(0).get_aha_num(), 7);
  CHECK_EQ(m.get_element_index(1), 0);
  CHECK_EQ(m.get_element(1).get_aha_num(), 0);
  CHECK_EQ(m.get_boundary_element(0).get_index(), 30);
}

TEST_CASE("read_xml: transversely isotropic fibres are completed into an orthonormal basis")
{
  write_text("f.xml",
    "<mesh celltype=\"hexahedron\" dim=\"3\"><nodes size=\"8\">"
    "<node id=\"0\" x=\"0\" y=\"0\" z=\"0\"/><node id=\"1\" x=\"1\" y=\"0\" z=\"0\"/>"
    "<node id=\"2\" x=\"1\" y=\"1\" z=\"0\"/><node id=\"3\" x=\"0\" y=\"1\" z=\"0\"/>"
    "<node id=\"4\" x=\"0\" y=\"0\" z=\"1\"/><node id=\"5\" x=\"1\" y=\"0\" z=\"1\"/>"
    "<node id=\"6\" x=\"1\" y=\"1\" z=\"1\"/><node id=\"7\" x=\"0\" y=\"1\" z=\"1\"/>"
    "</nodes><elements size=\"1\">"
    "<element id=\"0\" v0=\"0\" v1=\"1\" v2=\"2\" v3=\"3\" v4=\"4\" v5=\"5\" v6=\"6\" v7=\"7\"/>"
    "</elements>"
    "<element_data type=\"fiber_transversely_isotropic\">"
    "<element id=\"0\"><fiber>3.0,4.0,0.0</fiber></element></element_data></mesh>");
  Mesh m;
  m.read_xml("f.xml");
  const Element & e = m.get_element(0);
  CHECK_NEAR(arma::norm(e.get_fiber()), 1.0, 1e-14);
  CHECK_NEAR(e.get_fiber()(0), 0.6, 1e-14);
  CHECK_NEAR(arma::norm(e.get_trans()), 1.0, 1e-14);
  CHECK_NEAR(arma::norm(e.get_normal()), 1.0, 1e-14);
  CHECK_NEAR(arma::dot(e.get_fiber(), e.get_trans()), 0.0, 1e-14);
  CHECK_NEAR(arma::dot(e.get_fiber(), e.get_normal()), 0.0, 1e-14);
}

TEST_CASE("read_xml: a file with no cell type throws")
{
  write_text("bad.xml", "<mesh dim=\"3\"></mesh>");
  Mesh m;
  CHECK_THROWS(m.read_xml("bad.xml"));
}

TEST_CASE("read_xml: inverted tetrahedra are reordered to positive volume")
{
  write_text("inv.xml",
    "<mesh celltype=\"tetrahedron\" dim=\"3\"><nodes size=\"4\">"
    "<node id=\"0\" x=\"0\" y=\"0\" z=\"0\"/><node id=\"1\" x=\"1\" y=\"0\" z=\"0\"/>"
    "<node id=\"2\" x=\"0\" y=\"1\" z=\"0\"/><node id=\"3\" x=\"0\" y=\"0\" z=\"1\"/>"
    "</nodes><elements size=\"2\">"
    "<element id=\"0\" v0=\"0\" v1=\"1\" v2=\"2\" v3=\"3\"/>"
    "<element id=\"1\" v0=\"0\" v1=\"2\" v2=\"1\" v3=\"3\"/>"
    "</elements></mesh>");
  Mesh m;
  m.read_xml("inv.xml");
  CHECK_NEAR(tet_volume(m, 0), 1.0 / 6.0, 1e-14);
  CHECK_NEAR(tet_volume(m, 1), 1.0 / 6.0, 1e-14);
}

TEST_CASE("generated tetrahedral box: positive volumes that add up to the box")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  double total = 0.0;
  for (int e = 0; e < m.get_n_elements(); e++)
  {
    const double v = tet_volume(m, e);
    CHECK(v > 0.0);
    total += v;
  }
  CHECK_NEAR(total, 3.0 * 2.0 * 1.0, 1e-12);
}

TEST_CASE_KNOWN_BUG("Mesh: the centroid of an element is the mean of its nodes",
                    "src/mesh/mesh.cpp, get_element_centroid: uses the first three nodes of "
                    "the MESH instead of those of the element, multiplies by 0.3 instead of "
                    "1/3, and leaves z uninitialised")
{
  BoxSpec b;
  b.nx = 2; b.ny = b.nz = 1;
  write_hex_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");

  const arma::vec3 c = m.get_element_centroid(1);
  CHECK_NEAR(c(0), 0.75, 1e-14);
  CHECK_NEAR(c(1), 0.5, 1e-14);
  CHECK_NEAR(c(2), 0.5, 1e-14);
}

TEST_CASE("Mesh: swap_nodes and update_point")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 1;
  write_hex_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");

  std::vector<int> before, after;
  m.get_element_pt_nums(0, before);
  m.swap_nodes(0, 0, 1);
  m.get_element_pt_nums(0, after);
  CHECK_EQ(after[0], before[1]);
  CHECK_EQ(after[1], before[0]);

  m.update_point(0, arma::vec3({9.0, 8.0, 7.0}));
  CHECK_NEAR(m.get_point(0)(1), 8.0, 0.0);
}

TEST_CASE_KNOWN_BUG("Mesh: the element mapping carries the region marker",
                    "src/fem/mapping.cpp: the Mapping constructor receives the marker but "
                    "does not store it, so get_index() returns uninitialised memory")
{
  write_text("mk.xml",
    "<mesh celltype=\"tetrahedron\" dim=\"3\"><nodes size=\"4\">"
    "<node id=\"0\" x=\"0\" y=\"0\" z=\"0\"/><node id=\"1\" x=\"1\" y=\"0\" z=\"0\"/>"
    "<node id=\"2\" x=\"0\" y=\"1\" z=\"0\"/><node id=\"3\" x=\"0\" y=\"0\" z=\"1\"/>"
    "</nodes><elements size=\"1\">"
    "<element id=\"0\" marker=\"12345\" v0=\"0\" v1=\"1\" v2=\"2\" v3=\"3\"/>"
    "</elements></mesh>");
  Mesh m;
  m.read_xml("mk.xml");
  CHECK_EQ(m.get_element_mapping(0).get_index(), 12345);
}

// ---- legacy formats ---------------------------------------------------------

TEST_CASE("Mesh::read: .pts/.elem/.fib format")
{
  write_text("leg.pts", "4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n");
  write_text("leg.elem", "1\nTt 0 1 2 3 7\n");
  write_text("leg.fib", "0\n1 0 0 0 1 0 0 0 1\n");
  Mesh m;
  m.read("leg");
  CHECK_EQ(m.get_n_points(), 4u);
  CHECK_EQ(m.get_n_elements(), 1);
  CHECK_EQ(m.get_n_dim(), 3);
  CHECK_EQ(m.get_nen(), 4);
  CHECK_EQ(m.get_element_index(0), 7);
  CHECK_NEAR(m.get_element(0).get_trans()(1), 1.0, 0.0);
}

TEST_CASE("Mesh::read: a missing file throws")
{
  Mesh m;
  CHECK_THROWS(m.read("does_not_exist"));
}

TEST_CASE("GmshIO: triangles and boundary lines, with 1-based indices converted")
{
  write_text("g.msh",
    "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
    "$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 1 1 0\n4 0 1 0\n$EndNodes\n"
    "$Elements\n4\n"
    "1 1 2 11 1 1 2\n"
    "2 1 2 12 2 2 3\n"
    "3 2 2 5 1 1 2 3\n"
    "4 2 2 6 1 1 3 4\n"
    "$EndElements\n");
  Mesh m;
  GmshIO(m).read("g.msh");
  CHECK_EQ(m.get_n_points(), 4u);
  CHECK_EQ(m.get_n_dim(), 2);
  CHECK_EQ(m.get_n_elements(), 2);
  CHECK_EQ(m.get_n_boundary_elements(), 2);
  std::vector<int> p;
  m.get_element_pt_nums(1, p);
  CHECK_EQ(p[0], 0);
  CHECK_EQ(p[2], 3);
  CHECK_EQ(m.get_element_index(0), 5);
  CHECK_EQ(m.get_boundary_element(1).get_index(), 12);
}

TEST_CASE("GmshIO: a tetrahedron in a 3D mesh")
{
  write_text("g3.msh",
    "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
    "$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n$EndNodes\n"
    "$Elements\n1\n1 4 2 9 1 1 2 3 4\n$EndElements\n");
  Mesh m;
  GmshIO(m).read("g3.msh");
  CHECK_EQ(m.get_n_dim(), 3);
  CHECK_EQ(m.get_n_elements(), 1);
  CHECK_EQ(m.get_nen(), 4);
}

// ---- BoxMesh ----------------------------------------------------------------

TEST_CASE("BoxMesh: counts and connectivity on a box starting at the origin")
{
  BoxMesh m(0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 2, 1, 1);
  CHECK_EQ(m.get_n_points(), uint(3 * 2 * 2));
  CHECK_EQ(m.get_n_elements(), 2);
  CHECK_EQ(m.get_nen(), 8);
  CHECK_EQ(m.get_n_dim(), 3);
  CHECK_NEAR(m.get_point(2)(0), 2.0, 1e-14);
}

TEST_CASE_KNOWN_BUG("BoxMesh: interior faces are not boundary elements",
                    "src/mesh/box_mesh.cpp: a face is taken as boundary when each of its "
                    "nodes touches some wall, even different walls, which also catches "
                    "interior faces (here 12 faces instead of 10)")
{
  BoxMesh m(0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 2, 1, 1);
  CHECK_EQ(m.get_n_boundary_elements(), 10);
}

TEST_CASE_KNOWN_BUG("BoxMesh: a box away from the origin keeps its position",
                    "src/mesh/box_mesh.cpp: p(0) = i*dx ignores xmin (same for y and z), "
                    "so every box is generated starting at the origin")
{
  BoxMesh m(1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1);
  CHECK_NEAR(m.get_point(0)(0), 1.0, 1e-14);
  CHECK_NEAR(m.get_point(0)(1), 1.0, 1e-14);
  CHECK_NEAR(m.get_point(0)(2), 1.0, 1e-14);
}

// ---- PurkinjeMesh -----------------------------------------------------------

TEST_CASE("PurkinjeMesh: finds the bifurcation and the terminal nodes")
{
  // root(0) -- 1 -- 2 and 1 -- 3: node 1 is a bifurcation, 2 and 3 are PMJs
  write_text("pk.xml",
    "<mesh celltype=\"line\" dim=\"1\"><nodes size=\"4\">"
    "<node id=\"0\" x=\"0\" y=\"0\" z=\"0\"/><node id=\"1\" x=\"1\" y=\"0\" z=\"0\"/>"
    "<node id=\"2\" x=\"2\" y=\"1\" z=\"0\"/><node id=\"3\" x=\"2\" y=\"-1\" z=\"0\"/>"
    "</nodes><elements size=\"3\">"
    "<element id=\"0\" v0=\"0\" v1=\"1\"/><element id=\"1\" v0=\"1\" v1=\"2\"/>"
    "<element id=\"2\" v0=\"1\" v1=\"3\"/></elements></mesh>");
  PurkinjeMesh p;
  p.read_xml("pk.xml");
  p.set_root(0);
  p.setup();

  CHECK_EQ(p.get_num_bifurcations(), std::size_t(1));
  REQUIRE(p.begin() != p.end());
  CHECK_EQ(p.begin()->first, 1);
  std::set<int> neighbours(p.begin()->second.begin(), p.begin()->second.end());
  CHECK(neighbours == std::set<int>({0, 2, 3}));

  std::set<int> pmj(p.get_pmj_nodes().begin(), p.get_pmj_nodes().end());
  CHECK_MSG(pmj == std::set<int>({2, 3}), "the root is not a PMJ");
}

// ---- Graph colouring --------------------------------------------------------

TEST_CASE("comm_matrix: marks exactly the pairs of elements that share a node")
{
  BoxSpec b;
  b.nx = 3; b.ny = 1; b.nz = 1;
  write_hex_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  const arma::umat L = comm_matrix(m);
  CHECK_EQ(L(0, 1), arma::uword(1));
  CHECK_EQ(L(1, 2), arma::uword(1));
  CHECK_EQ(L(0, 2), arma::uword(0));
  CHECK_EQ(L(1, 0), L(0, 1));
  CHECK_EQ(L(0, 0), arma::uword(0));
}

TEST_CASE("greedy_coloring: elements that share a node get different colours")
{
  BoxSpec b;
  b.nx = 3; b.ny = 2; b.nz = 1;
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  const arma::umat L = comm_matrix(m);
  const arma::umat C = greedy_coloring(L, 2);
  REQUIRE(C.n_elem == L.n_rows);
  for (arma::uword i = 0; i < L.n_rows; i++)
  {
    CHECK_MSG(C(i) != 0, "element " << i << " has no colour");
    for (arma::uword j = i + 1; j < L.n_rows; j++)
      if (L(i, j)) CHECK_MSG(C(i) != C(j), "elements " << i << " and " << j);
  }
}

// ---- WriterHDF5 -------------------------------------------------------------

TEST_CASE("WriterHDF5: geometry, topology and time are written when the file is opened")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  WriterHDF5 w(&m);
  w.open("res", 3, 0.5);

  const arma::mat coords = h5_read("res.h5", "/geometry/coordinates");
  CHECK_EQ(coords.n_rows, arma::uword(m.get_n_points()));
  CHECK_EQ(coords.n_cols, arma::uword(3));
  CHECK_NEAR(coords(box_node(b, 3, 2, 1), 0), 3.0, 1e-14);

  const arma::mat time = h5_read("res.h5", "/time");
  CHECK_EQ(time.n_rows, arma::uword(3));
  CHECK_NEAR(time(2, 0), 1.0, 1e-14);

  CHECK(file_exists("res.xmf"));
  CHECK(read_text("res.xmf").find("res.h5") != std::string::npos);
}

TEST_CASE("WriterHDF5: nodal, displacement and cell fields land in the right step")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  const int np = m.get_n_points(), ne = m.get_n_elements();

  WriterHDF5 w(&m);
  w.open("res", 3, 1.0);

  arma::vec vm(np), ta(np), cell(ne);
  arma::mat disp(3, np);
  for (int i = 0; i < np; i++) { vm(i) = i; ta(i) = 2.0 * i; }
  for (int e = 0; e < ne; e++) cell(e) = 3.0 * e;
  for (int i = 0; i < np; i++)
    for (int d = 0; d < 3; d++) disp(d, i) = 100.0 * i + d;

  w.write_vm_step(1, vm.memptr());
  w.write_point_field_step(2, ta.memptr(), "active_stress");
  w.write_cell_field_step(1, cell.memptr(), "stress");
  w.write_displ_step(2, disp.memptr());

  const arma::mat rvm = h5_read("res.h5", "/point_data/vm");
  CHECK_EQ(rvm.n_rows, arma::uword(3));
  CHECK_NEAR(rvm(1, 7), 7.0, 0.0);
  CHECK_NEAR(rvm(0, 7), 0.0, 0.0);

  const arma::mat rta = h5_read("res.h5", "/point_data/active_stress");
  CHECK_NEAR(rta(2, 5), 10.0, 0.0);

  const arma::mat rc = h5_read("res.h5", "/cell_data/stress");
  CHECK_NEAR(rc(1, 3), 9.0, 0.0);

  // [nsteps, np, 3] read as nsteps rows of np*3 values
  const arma::mat rd = h5_read("res.h5", "/point_data/displacements");
  CHECK_NEAR(rd(2, 3 * 4 + 1), 401.0, 0.0);
}

TEST_CASE("WriterHDF5: writing a field with no dataset warns and writes nothing")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  WriterHDF5 w(&m);
  w.open("res", 1, 1.0);
  arma::vec v(m.get_n_points(), arma::fill::ones);
  w.write_point_field_step(0, v.memptr(), "no_such_field");
  CHECK(!h5_has("res.h5", "/point_data/no_such_field"));
}

TEST_CASE("WriterHDF5: region markers and AHA numbers are saved per element")
{
  write_text("mk.xml",
    "<mesh celltype=\"tetrahedron\" dim=\"3\"><nodes size=\"4\">"
    "<node id=\"0\" x=\"0\" y=\"0\" z=\"0\"/><node id=\"1\" x=\"1\" y=\"0\" z=\"0\"/>"
    "<node id=\"2\" x=\"0\" y=\"1\" z=\"0\"/><node id=\"3\" x=\"0\" y=\"0\" z=\"1\"/>"
    "</nodes><elements size=\"1\">"
    "<element id=\"0\" marker=\"4\" aha=\"11\" v0=\"0\" v1=\"1\" v2=\"2\" v3=\"3\"/>"
    "</elements></mesh>");
  Mesh m;
  m.read_xml("mk.xml");
  WriterHDF5 w(&m);
  w.open("res", 2, 1.0);
  CHECK_NEAR(h5_read("res.h5", "/cell_data/fibrosis")(1, 0), 4.0, 0.0);
  CHECK_NEAR(h5_read("res.h5", "/cell_data/aha_marker")(1, 0), 11.0, 0.0);
}

TEST_CASE("WriterHDF5: the EP checkpoint round trip returns exactly what was saved")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  const int np = m.get_n_points(), nv = 5;

  arma::vec vm = arma::randu<arma::vec>(np), sv = arma::randu<arma::vec>(np * nv);
  WriterHDF5 w(&m);
  w.write_checkpoint(123, 4.5, vm.memptr(), sv.memptr(), nv);
  REQUIRE(file_exists("checkpoint_t_4.50ms.h5"));

  int step = -1, nn = -1, nvars = -1;
  double t = -1.0;
  w.read_checkpoint_metadata("checkpoint_t_4.50ms.h5", step, t, nn, nvars);
  CHECK_EQ(step, 123);
  CHECK_NEAR(t, 4.5, 0.0);
  CHECK_EQ(nn, np);
  CHECK_EQ(nvars, nv);

  arma::vec vm2(np), sv2(np * nv);
  w.read_checkpoint_data("checkpoint_t_4.50ms.h5", vm2.memptr(), sv2.memptr());
  CHECK_NEAR(arma::abs(vm - vm2).max(), 0.0, 0.0);
  CHECK_NEAR(arma::abs(sv - sv2).max(), 0.0, 0.0);
}

TEST_CASE("WriterHDF5: the mechanics checkpoint round trip returns exactly what was saved")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  const int ndofs = 3 * m.get_n_points();

  arma::vec x = arma::randu<arma::vec>(ndofs), f = arma::randu<arma::vec>(ndofs);
  WriterHDF5 w(&m);
  w.write_mech_checkpoint(7, 2.0, 3, 0.25, x.memptr(), f.memptr(), ndofs);

  // The name is the writer's business: find what it produced
  std::string file;
  for (const auto & e : std::filesystem::directory_iterator("."))
    if (e.path().filename().string().find("mech") != std::string::npos ||
        e.path().filename().string().find("checkpoint") != std::string::npos)
      if (e.path().extension() == ".h5" &&
          e.path().filename().string() != "checkpoint_t_4.50ms.h5")
        file = e.path().string();
  REQUIRE(!file.empty());

  int step = -1, inc = -1, n = -1;
  double t = -1.0, load = -1.0;
  w.read_mech_checkpoint_metadata(file, step, t, inc, load, n);
  CHECK_EQ(step, 7);
  CHECK_NEAR(t, 2.0, 0.0);
  CHECK_EQ(inc, 3);
  CHECK_NEAR(load, 0.25, 0.0);
  CHECK_EQ(n, ndofs);

  arma::vec x2(ndofs), f2(ndofs);
  w.read_mech_checkpoint_data(file, x2.memptr(), f2.memptr());
  CHECK_NEAR(arma::abs(x - x2).max(), 0.0, 0.0);
  CHECK_NEAR(arma::abs(f - f2).max(), 0.0, 0.0);
}

TEST_CASE("WriterHDF5: the eikonal LAT output has its own file and geometry")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  arma::vec lat(m.get_n_points());
  for (arma::uword i = 0; i < lat.n_elem; i++) lat(i) = 0.5 * i;

  WriterHDF5 w(&m);
  w.write_eikonal_lat("eik", lat.memptr());
  REQUIRE(file_exists("eik_lat.h5"));
  CHECK(file_exists("eik_lat.xmf"));

  bool found = false;
  for (const std::string & d : h5_datasets("eik_lat.h5"))
    if (d.find("/point_data/") == 0)
    {
      const arma::mat r = h5_read("eik_lat.h5", d);
      CHECK_NEAR(r(r.n_elem - 1), lat(lat.n_elem - 1), 1e-14);
      found = true;
    }
  CHECK(found);
}

// ---- ReaderHDF5 -------------------------------------------------------------

TEST_CASE("ReaderHDF5 and Mesh::read_hdf5: read the mesh written by WriterHDF5")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  WriterHDF5 w(&m);
  w.open("res", 2, 1.0);

  ReaderHDF5 r;
  REQUIRE(r.open("res.h5"));
  CHECK_EQ(r.get_n_points(), int(m.get_n_points()));
  CHECK_EQ(r.get_n_elements(), m.get_n_elements());
  CHECK_EQ(r.get_nen(), 4);
  CHECK_EQ(r.get_time().size(), std::size_t(2));

  Mesh again;
  again.read_hdf5("res", 3);
  CHECK_EQ(again.get_n_points(), m.get_n_points());
  CHECK_EQ(again.get_n_elements(), m.get_n_elements());
  const arma::vec3 & p = again.get_point(box_node(b, 3, 2, 1));
  CHECK_NEAR(p(0), 3.0, 1e-14);
  std::vector<int> p1, p2;
  m.get_element_pt_nums(5, p1);
  again.get_element_pt_nums(5, p2);
  CHECK(p1 == p2);
}

TEST_CASE_KNOWN_BUG("ReaderHDF5: reads back the fields written by WriterHDF5",
                    "src/mesh/reader_hdf5.cpp looks for fields in /vertex_field and /cell_field, "
                    "but WriterHDF5 writes them in /point_data and /cell_data: no field is found")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  WriterHDF5 w(&m);
  w.open("res", 2, 1.0);
  arma::vec vm(m.get_n_points());
  for (arma::uword i = 0; i < vm.n_elem; i++) vm(i) = 1.5 * i;
  w.write_vm_step(1, vm.memptr());

  ReaderHDF5 r;
  REQUIRE(r.open("res.h5"));
  CHECK(r.get_n_fields() > 0);
  CHECK(r.find_field("vm") >= 0);
  std::vector<double> out;
  REQUIRE(r.read_field_step("vm", 1, out));
  REQUIRE(out.size() == vm.n_elem);
  CHECK_NEAR(out[4], 6.0, 0.0);
}

TEST_CASE("ReaderHDF5: a missing file fails cleanly")
{
  ReaderHDF5 r;
  CHECK(!r.open("not_here.h5"));
  Mesh m;
  CHECK_THROWS(m.read_hdf5("not_here", 3));
}

// ---- WriterVTK --------------------------------------------------------------

TEST_CASE("WriterVTK: writes an XML file with the mesh sizes")
{
  BoxSpec b = small_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  WriterVTK w(&m);
  w.open("res_vtk", 1, 1.0);
  w.close();

  bool found = false;
  for (const auto & e : std::filesystem::directory_iterator("."))
    if (e.path().filename().string().find("res_vtk") == 0)
    {
      const std::string text = read_text(e.path().string());
      found = text.find("VTKFile") != std::string::npos &&
              text.find(std::to_string(m.get_n_points())) != std::string::npos;
    }
  CHECK(found);
}
