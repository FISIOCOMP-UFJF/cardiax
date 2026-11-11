/*
 * Unit tests of the mechanics assembly (src/pdes/nonlinear_elasticity, updated and
 * total Lagrangian), with no Newton iterations.
 *
 * What is covered and why:
 *  - Consistent tangent: the matrix from jacobian() must be the derivative of the
 *    residual from evaluate(), compared here against centred finite differences of the
 *    residual on a deformed configuration. This is the strongest single check of the
 *    assembly: it involves the geometric and material parts, the mean dilatation
 *    (hexahedra), the following pressure load, the Robin springs and the active stress.
 *    If the tangent is wrong, Newton loses its quadratic rate and the runs get slow or
 *    diverge, with no other visible symptom.
 *  - Equilibrium: the undeformed, unloaded state has zero residual; a rigid body motion
 *    does too.
 *  - Boundary conditions: prescribed dofs have unit rows and columns in the matrix and
 *    are not moved by update().
 *  - The volume integrals: tissue volume, and the cavity volume with the divergence
 *    theorem on a closed surface and with the lid plane on an open one. The circulation
 *    models depend entirely on these numbers.
 *  - Mesh unit detection, which scales all volumes to mL.
 *  - Fibre stretch per element and per node, input of the mechano-electric feedback.
 *  - The active scale per material and per region.
 *
 * Each case builds its own solver object. The Newton solver of UpdatedLagrangian lives in
 * a function-local static, so complete solves are tested in the integration tests
 * (test_app_nonlinearelas.cpp), one per process.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <regex>

#include "pdes/updated_lagrangian.hpp"
#include "pdes/total_lagrangian.hpp"

using namespace cardiax_test;

namespace
{

template <typename Base>
struct Open : public Base
{
  using Base::x;
  using Base::x0;
  using Base::msh;
  using Base::ldgof;
  using Base::material;
  using Base::lc;

  void load(const std::string & file)
  {
    this->filename = file;
    this->basename = file.substr(0, file.rfind(".xml"));
    this->init_mesh();
    this->config(file, file);
    this->init_matvecs();
    this->pre_solve();
    // Whole load applied
    this->lc.set_nincs(1);
    this->lc.update();
  }

  arma::vec residual()
  {
    petsc::Vector r;
    r.create(this->size());
    this->evaluate(r);
    arma::vec out(this->size());
    r.get_data(out.memptr());
    return out;
  }

  arma::mat tangent()
  {
    petsc::Matrix K;
    K.create(this->size(), this->size(), 120);
    this->jacobian(K);
    const int n = this->size();
    arma::mat D(n, n);
    std::vector<int> idx(n);
    for (int i = 0; i < n; i++) idx[i] = i;
    MatGetValues(K.mat(), n, idx.data(), n, idx.data(), D.memptr());
    return D.t();
  }

  //! Smooth, non-trivial deformation of the free dofs
  void deform(double a)
  {
    const int nd = this->msh.get_n_dim();
    for (uint i = 0; i < this->msh.get_n_points(); i++)
    {
      const arma::vec3 & p = this->x0[i];
      arma::vec3 u = {a * (p(0) * p(0) + 0.3 * p(1)),
                      a * (0.5 * p(0) * p(1) - 0.2 * p(2)),
                      a * (0.4 * p(1) * p(2) + 0.1 * p(0))};
      for (int d = 0; d < nd; d++)
        if (this->ldgof[i * nd + d]) this->x[i](d) = p(d) + u(d);
    }
  }

  void move_dof(int dof, double h)
  {
    const int nd = this->msh.get_n_dim();
    this->x[dof / nd](dof % nd) += h;
  }
};

//! Total Lagrangian only fills its per-element data at the start of solve()
struct OpenTL : public Open<TotalLagrangian>
{
  void prepare() { this->pre_computation(); }
};

std::string elasticity_block(const std::string & material, const std::string & coefficients,
                             const std::string & fixed, const std::string & extra = "")
{
  return "<elasticity type=\"THREE_DIM\"><parameters><material>" + material + "</material>"
         "<coefficients>" + coefficients + "</coefficients><ninc>1</ninc></parameters>"
         + extra + "<prescribed_displacement>" + fixed + "</prescribed_displacement></elasticity>\n";
}

// Nodes on the x = 0 face of a box
std::vector<int> face_x0(const BoxSpec & b)
{
  std::vector<int> n;
  for (int k = 0; k <= b.nz; k++)
    for (int j = 0; j <= b.ny; j++)
      n.push_back(box_node(b, 0, j, k));
  return n;
}

const char * kGuccione = "2.0, 8.0, 2.0, 4.0, 0.0, 200.0";
const char * kHolzapfel = "2.0, 3.0, 2.0, 1.5, 4.0, 5.0, 3.0, 2.0, 0.0, 200.0";

//! -dR/dx by centred differences on the free dofs against the assembled tangent.
//! Returns the largest entry gap, relative to the largest tangent entry.
template <typename Problem>
double tangent_gap(Problem & p, double h = 1e-5)
{
  const int n = p.size();
  const arma::mat K = p.tangent();

  double gap = 0.0;
  const double scale = arma::abs(K).max();
  for (int j = 0; j < n; j++)
  {
    if (!p.ldgof[j]) continue;
    p.move_dof(j, h);
    const arma::vec rp = p.residual();
    p.move_dof(j, -2.0 * h);
    const arma::vec rm = p.residual();
    p.move_dof(j, h);

    const arma::vec col = -(rp - rm) / (2.0 * h);
    for (int i = 0; i < n; i++)
      if (p.ldgof[i])
        gap = std::max(gap, std::fabs(col(i) - K(i, j)));
  }
  return gap / scale;
}

template <typename Mesher>
double ul_gap(const std::string & name, Mesher write, const std::string & material,
              const std::string & coefs, const std::string & extra, double active,
              double amplitude = 0.02)
{
  BoxSpec b;
  b.nx = 2; b.ny = 1; b.nz = 1;
  b.prefix_xml = elasticity_block(material, coefs, fix_nodes_xml(face_x0(b), {0, 1, 2}), extra);
  write(name, b);

  Open<UpdatedLagrangian> p;
  p.load(name);
  if (active != 0.0)
  {
    arma::vec ta(p.msh.get_n_points());
    ta.fill(active);
    p.set_active_stress(ta);
  }
  p.deform(amplitude);
  return tangent_gap(p);
}

} // namespace

// ---- Consistency of the tangent ---------------------------------------------

TEST_CASE("Updated Lagrangian: tangent of Guccione on tetrahedra")
{
  const double g = ul_gap("a.xml", write_tet_box, "Guccione", kGuccione, "", 0.0);
  CHECK_MSG(g < 1e-3, "relative gap " << g);
}

TEST_CASE("Updated Lagrangian: tangent of Guccione on hexahedra (mean dilatation)")
{
  const double g = ul_gap("b.xml", write_hex_box, "Guccione", kGuccione, "", 0.0);
  CHECK_MSG(g < 1e-3, "relative gap " << g);
}

TEST_CASE("Updated Lagrangian: tangent of Holzapfel-Ogden on tetrahedra")
{
  const double g = ul_gap("c.xml", write_tet_box, "HolzapfelOgden", kHolzapfel, "", 0.0);
  CHECK_MSG(g < 1e-3, "relative gap " << g);
}

TEST_CASE("Updated Lagrangian: tangent of Holzapfel-Ogden on hexahedra")
{
  const double g = ul_gap("d.xml", write_hex_box, "HolzapfelOgden", kHolzapfel, "", 0.0);
  CHECK_MSG(g < 1e-3, "relative gap " << g);
}

TEST_CASE("Updated Lagrangian: tangent with active tension")
{
  const double g = ul_gap("e.xml", write_tet_box, "Guccione", kGuccione, "", 5.0);
  CHECK_MSG(g < 1e-3, "relative gap " << g);
}

TEST_CASE("Updated Lagrangian: tangent with a following pressure load")
{
  const std::string pressure =
    "<pressure><node id=\"1\" marker=\"2\" value=\"3.0\"/>"
    "<node id=\"2\" marker=\"1\" value=\"0.0\"/><node id=\"3\" marker=\"3\" value=\"0.0\"/>"
    "<node id=\"4\" marker=\"4\" value=\"0.0\"/><node id=\"5\" marker=\"5\" value=\"0.0\"/>"
    "<node id=\"6\" marker=\"6\" value=\"0.0\"/></pressure>";
  const double g = ul_gap("f.xml", write_tet_box, "Guccione", kGuccione, pressure, 0.0);
  CHECK_MSG(g < 1e-3, "relative gap " << g);
}

TEST_CASE("Updated Lagrangian: tangent with a following pressure on a closed quadrilateral surface")
{
  // The pressure stiffness of Bonet & Wood is symmetrised by an integration by parts, whose
  // boundary term vanishes on a closed surface. On bilinear faces of an open surface the
  // form is an approximation, so here the load covers the whole boundary.
  const std::string pressure =
    "<pressure><node id=\"1\" marker=\"1\" value=\"3.0\"/><node id=\"2\" marker=\"2\" value=\"3.0\"/>"
    "<node id=\"3\" marker=\"3\" value=\"3.0\"/><node id=\"4\" marker=\"4\" value=\"3.0\"/>"
    "<node id=\"5\" marker=\"5\" value=\"3.0\"/><node id=\"6\" marker=\"6\" value=\"3.0\"/></pressure>";
  const double g = ul_gap("g.xml", write_hex_box, "Guccione", kGuccione, pressure, 0.0);
  CHECK_MSG(g < 1e-3, "relative gap " << g);
}

TEST_CASE("Total Lagrangian: tangent of Guccione on hexahedra")
{
  BoxSpec b;
  b.nx = 2; b.ny = 1; b.nz = 1;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, fix_nodes_xml(face_x0(b), {0, 1, 2}));
  write_hex_box("tl.xml", b);
  OpenTL p;
  p.load("tl.xml");
  p.prepare();
  p.deform(0.02);
  const double g = tangent_gap(p);
  CHECK_MSG(g < 1e-3, "relative gap " << g);
}

TEST_CASE("Total Lagrangian: tangent of Holzapfel-Ogden on tetrahedra")
{
  BoxSpec b;
  b.nx = 2; b.ny = 1; b.nz = 1;
  b.prefix_xml = elasticity_block("HolzapfelOgden", kHolzapfel, fix_nodes_xml(face_x0(b), {0, 1, 2}));
  write_tet_box("tl3.xml", b);
  OpenTL p;
  p.load("tl3.xml");
  p.prepare();
  p.deform(0.02);
  const double g = tangent_gap(p);
  CHECK_MSG(g < 1e-3, "relative gap " << g);
}

// allocate_Ta used to read active_stress(0) of an empty vector, so the compressible laws
// threw in config(); fixed in 5668b24
TEST_CASE("Mechanics: a compressible material can be configured from XML")
{
  BoxSpec b;
  b.nx = 2; b.ny = 1; b.nz = 1;
  b.prefix_xml = elasticity_block("CompressibleNeoHookean", "1.5, 2.5", fix_nodes_xml(face_x0(b), {0, 1, 2}));
  write_tet_box("tl2.xml", b);
  OpenTL p;
  p.load("tl2.xml");
  p.prepare();
  p.deform(0.02);
  CHECK(tangent_gap(p) < 1e-3);
}

// ---- Equilibrium ------------------------------------------------------------

TEST_CASE("Mechanics: the undeformed state with no load has no residual")
{
  BoxSpec b;
  b.nx = 2; b.ny = 2; b.nz = 2;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, fix_nodes_xml(face_x0(b), {0, 1, 2}));
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  CHECK_NEAR(arma::abs(p.residual()).max(), 0.0, 1e-5);
}

TEST_CASE("Mechanics: a rigid body rotation has no residual")
{
  BoxSpec b;
  b.nx = 2; b.ny = 2; b.nz = 1;
  // Nothing fixed: the motion is free
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  const double a = 0.4;
  const arma::mat33 R = {{std::cos(a), -std::sin(a), 0}, {std::sin(a), std::cos(a), 0}, {0, 0, 1}};
  for (uint i = 0; i < p.msh.get_n_points(); i++)
    p.x[i] = R * p.x0[i] + arma::vec3({1.0, 2.0, 3.0});
  CHECK_NEAR(arma::abs(p.residual()).max(), 0.0, 1e-4);
}

TEST_CASE("Mechanics: internal forces of a free body add up to zero")
{
  BoxSpec b;
  b.nx = 2; b.ny = 1; b.nz = 1;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_hex_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  p.deform(0.03);
  const arma::vec r = p.residual();
  arma::vec3 total = {0, 0, 0};
  for (uint i = 0; i < p.msh.get_n_points(); i++)
    for (int d = 0; d < 3; d++)
      total(d) += r(3 * i + d);
  CHECK_MSG(arma::abs(total).max() < 1e-8 * (1.0 + arma::abs(r).max()), total.t());
  CHECK(arma::abs(r).max() > 1e-3);   // a real deformation, or the sum proves nothing
}

// ---- Boundary conditions ----------------------------------------------------

TEST_CASE("Mechanics: prescribed dofs have unit diagonal in the matrix and do not move")
{
  BoxSpec b;
  b.nx = 2; b.ny = 1; b.nz = 1;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, fix_nodes_xml(face_x0(b), {0, 1, 2}));
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  p.deform(0.01);

  const int n = p.size();
  CHECK_EQ(n, 3 * int(p.msh.get_n_points()));
  const arma::mat K = p.tangent();
  const int fixed = 3 * box_node(b, 0, 1, 1) + 2;
  REQUIRE(!p.ldgof[fixed]);
  CHECK_NEAR(K(fixed, fixed), 1.0, 0.0);
  for (int i = 0; i < n; i++)
    if (i != fixed)
    {
      CHECK_NEAR(K(fixed, i), 0.0, 0.0);
      CHECK_NEAR(K(i, fixed), 0.0, 0.0);
    }

  const arma::vec3 before = p.x[box_node(b, 0, 1, 1)];
  petsc::Vector du;
  du.create(n);
  du.set(0.1);
  p.update(du, 1.0);
  CHECK_NEAR(arma::abs(p.x[box_node(b, 0, 1, 1)] - before).max(), 0.0, 0.0);
  const int freen = box_node(b, 2, 0, 0);
  CHECK(p.ldgof[3 * freen]);
}

TEST_CASE("Mechanics: update(u, s) moves the free dofs by s*u")
{
  BoxSpec b;
  b.nx = 1; b.ny = 1; b.nz = 1;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  const int n = p.size();
  petsc::Vector du;
  du.create(n);
  du = 0.0;
  du.set(3 * 1 + 2, 0.5);
  du.assemble();
  const double z0 = p.x[1](2);
  p.update(du, 0.2);
  CHECK_NEAR(p.x[1](2) - z0, 0.1, 1e-14);
}

// ---- Volumes ----------------------------------------------------------------

namespace
{

std::string relabel(const std::string & xml, const std::map<int, int> & markers)
{
  // Only the <boundary> block has elements with marker attributes in the generated meshes
  std::string out = xml;
  for (const auto & kv : markers)
    out = std::regex_replace(out, std::regex("marker=\"" + std::to_string(kv.first) + "\""),
                             "marker=\"M" + std::to_string(kv.second) + "\"");
  out = std::regex_replace(out, std::regex("marker=\"M"), "marker=\"");
  return out;
}

} // namespace

TEST_CASE("calc_volume: tissue volume, reference and after a deformation")
{
  BoxSpec b;
  b.nx = 2; b.ny = 2; b.nz = 2;
  b.Lx = 2.0; b.Ly = 3.0; b.Lz = 4.0;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  CHECK_NEAR(p.calc_volume(), 24.0, 1e-12);

  for (uint i = 0; i < p.msh.get_n_points(); i++) p.x[i](2) *= 1.5;
  CHECK_NEAR(p.calc_volume(), 36.0, 1e-12);
}

TEST_CASE("Cavity volume: a closed surface gives the volume it encloses, in both element types")
{
  for (bool hex : {false, true})
  {
    BoxSpec b;
    b.nx = b.ny = b.nz = 2;
    b.Lx = 2.0; b.Ly = 3.0; b.Lz = 4.0;
    b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
    const std::string name = hex ? "h.xml" : "t.xml";
    if (hex) write_hex_box(name, b); else write_tet_box(name, b);
    write_text(name, relabel(read_text(name), {{1, 30}, {2, 30}, {3, 30}, {4, 30}, {5, 30}, {6, 30}}));

    Open<UpdatedLagrangian> p;
    p.load(name);
    p.set_mesh_units("cm");
    CHECK_NEAR(p.volume_LV(), 24.0, 1e-10);
    CHECK_NEAR(p.volume_RV(), 0.0, 0.0);

    // Same after translating the whole body: the surface flux formula must not depend
    // on where the origin is
    for (uint i = 0; i < p.msh.get_n_points(); i++) p.x[i] += arma::vec3({5.0, -3.0, 7.0});
    CHECK_NEAR(p.volume_LV(), 24.0, 1e-9);
  }
}

TEST_CASE("Cavity volume: closing the base with its own marker adds up to the box")
{
  // calc_cavity_volume integrates x.n on the requested marker only and does not apply the
  // lid plane, which only report_boundary_closure() uses: a cavity must be closed by its
  // own marker. Here the two markers of an open cavity add up to the closed volume.
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.Lx = 2.0; b.Ly = 3.0; b.Lz = 4.0;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("m.xml", b);
  write_text("m.xml", relabel(read_text("m.xml"),
                              {{1, 30}, {2, 30}, {3, 30}, {4, 30}, {5, 30}, {6, 20}}));
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  p.set_mesh_units("cm");
  CHECK_NEAR(p.volume_LV() + p.volume_RV(), 24.0, 1e-9);
  CHECK(p.volume_RV() > 0.0);

  p.detect_cavity_lid_plane(MARKER_RV);
  p.report_boundary_closure();
  CHECK_NEAR(p.volume_LV() + p.volume_RV(), 24.0, 1e-9);
}

TEST_CASE("Cavity volume: markers other than the requested one do not count")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  p.set_mesh_units("cm");
  CHECK_NEAR(p.volume_LV(), 0.0, 0.0);
  CHECK_NEAR(p.volume_RV(), 0.0, 0.0);
}

TEST_CASE("Mesh units: scale factors, detection by size, and an unknown unit is ignored")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 1;
  b.Lx = b.Ly = b.Lz = 0.1;   // 10 cm in metres
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");

  p.detect_mesh_units();
  CHECK_EQ(p.mesh_units(), std::string("m"));
  CHECK_NEAR(p.volume_scale_to_mL(), 1e6, 0.0);

  p.set_mesh_units("cm");
  CHECK_NEAR(p.volume_scale_to_mL(), 1.0, 0.0);
  p.set_mesh_units("mm");
  CHECK_NEAR(p.volume_scale_to_mL(), 1e-3, 0.0);
  p.set_mesh_units("furlong");
  CHECK_EQ(p.mesh_units(), std::string("mm"));
  CHECK_NEAR(p.volume_myocardium_mL(), 0.001 * 1e-3, 1e-15);
}

TEST_CASE("Mesh units: the size of the geometry selects the unit")
{
  auto guess = [](double L) {
    BoxSpec b;
    b.nx = b.ny = b.nz = 1;
    b.Lx = b.Ly = b.Lz = L;
    b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
    write_tet_box("u.xml", b);
    Open<UpdatedLagrangian> p;
    p.load("u.xml");
    p.detect_mesh_units();
    return p.mesh_units();
  };
  CHECK_EQ(guess(0.1), std::string("m"));
  CHECK_EQ(guess(10.0), std::string("cm"));
  CHECK_EQ(guess(100.0), std::string("mm"));
}

TEST_CASE_KNOWN_BUG("Mesh units: a mesh in micrometres gives volumes in mL",
                    "CLAUDE.md and cardiac_cycle state lengths in um, but set_mesh_units and "
                    "detect_mesh_units only know m, cm and mm: a heart in um is taken as mm, "
                    "so its volume comes out 1e9 times too large")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 1;
  b.Lx = b.Ly = b.Lz = 1.0e5;     // 10 cm in um: a volume of 1000 mL
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("um.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("um.xml");
  p.detect_mesh_units();
  CHECK_NEAR(p.volume_myocardium_mL(), 1000.0, 1e-6);
}

// ---- Fibre stretch ----------------------------------------------------------

TEST_CASE("Fibre stretch: the undeformed state is 1, per element and per node")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");

  arma::vec le, ln;
  p.fiber_stretch_elements(le);
  p.fiber_stretch_nodes(ln);
  CHECK_EQ(le.n_elem, arma::uword(p.msh.get_n_elements()));
  CHECK_EQ(ln.n_elem, arma::uword(p.msh.get_n_points()));
  CHECK_NEAR(arma::abs(le - 1.0).max(), 0.0, 1e-14);
  CHECK_NEAR(arma::abs(ln - 1.0).max(), 0.0, 1e-14);
}

TEST_CASE("Fibre stretch: a uniform stretch along the fibres is measured after the assembly")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  for (uint i = 0; i < p.msh.get_n_points(); i++) p.x[i](0) = 1.2 * p.x0[i](0);

  p.residual();      // the deformation gradient is computed in the element loop
  arma::vec le, ln;
  p.fiber_stretch_elements(le);
  p.fiber_stretch_nodes(ln);
  CHECK_NEAR(arma::abs(le - 1.2).max(), 0.0, 1e-12);
  CHECK_NEAR(arma::abs(ln - 1.2).max(), 0.0, 1e-12);
}

TEST_CASE("Fibre stretch: a stretch across the fibres does not count")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  for (uint i = 0; i < p.msh.get_n_points(); i++) p.x[i](1) = 1.3 * p.x0[i](1);
  p.residual();
  arma::vec le;
  p.fiber_stretch_elements(le);
  CHECK_NEAR(arma::abs(le - 1.0).max(), 0.0, 1e-12);
}

// ---- Active tension per region ----------------------------------------------

namespace
{

std::string two_materials(const std::string & scale0, const std::string & scale1)
{
  return "<elasticity type=\"THREE_DIM\"><parameters num_materials=\"2\"><regions>"
         "<marker id=\"0\" material=\"0\"/><marker id=\"1\" material=\"1\"/></regions>"
         "<material id=\"0\" type=\"Guccione\"" + scale0 + "><coefficients>" + kGuccione +
         "</coefficients></material>"
         "<material id=\"1\" type=\"Guccione\"" + scale1 + "><coefficients>" + kGuccione +
         "</coefficients></material>"
         "<ninc>1</ninc></parameters></elasticity>\n";
}

// The two tets of a pair of elements, with markers 0 and 1
void write_two_regions(const std::string & path, const std::string & prefix)
{
  write_text(path, prefix +
    "<mesh celltype=\"tetrahedron\" dim=\"3\"><nodes size=\"5\">"
    "<node id=\"0\" x=\"0\" y=\"0\" z=\"0\"/><node id=\"1\" x=\"1\" y=\"0\" z=\"0\"/>"
    "<node id=\"2\" x=\"0\" y=\"1\" z=\"0\"/><node id=\"3\" x=\"0\" y=\"0\" z=\"1\"/>"
    "<node id=\"4\" x=\"-1\" y=\"0\" z=\"0\"/>"
    "</nodes><elements size=\"2\">"
    "<element id=\"0\" marker=\"0\" v0=\"0\" v1=\"1\" v2=\"2\" v3=\"3\"/>"
    "<element id=\"1\" marker=\"1\" v0=\"4\" v1=\"0\" v2=\"2\" v3=\"3\"/>"
    "</elements>"
    "<element_data type=\"fiber_orthotropic\">"
    "<element id=\"0\"><fiber>1,0,0</fiber><sheet>0,1,0</sheet><normal>0,0,1</normal></element>"
    "<element id=\"1\"><fiber>1,0,0</fiber><sheet>0,1,0</sheet><normal>0,0,1</normal></element>"
    "</element_data></mesh>");
}

} // namespace

TEST_CASE("Active scale: effective tension follows the ta_scale of each region's material")
{
  write_two_regions("r.xml", two_materials(" ta_scale=\"1.0\"", " ta_scale=\"0.25\""));
  Open<UpdatedLagrangian> p;
  p.load("r.xml");
  arma::vec ta(5);
  ta.fill(100.0);
  arma::vec eff;
  p.effective_active_tension(ta, eff);
  REQUIRE(eff.n_elem == 2);
  CHECK_NEAR(eff(0), 100.0, 1e-12);
  CHECK_NEAR(eff(1), 25.0, 1e-12);

  arma::vec scale;
  p.active_scale_field(scale);
  CHECK_NEAR(scale(1), 0.25, 0.0);
}

TEST_CASE("Active scale: a passive region contributes no active force")
{
  write_two_regions("a.xml", two_materials(" ta_scale=\"1.0\"", " ta_scale=\"0.0\""));
  write_two_regions("b.xml", two_materials(" ta_scale=\"0.0\"", " ta_scale=\"0.0\""));
  Open<UpdatedLagrangian> act, pas;
  act.load("a.xml");
  pas.load("b.xml");
  arma::vec ta(5);
  ta.fill(100.0);
  act.set_active_stress(ta);
  pas.set_active_stress(ta);
  CHECK(arma::abs(act.residual()).max() > 1.0);
  CHECK_NEAR(arma::abs(pas.residual()).max(), 0.0, 1e-4);
}

TEST_CASE("Materials block: invalid declarations are rejected")
{
  const std::string bad_mat =
    "<elasticity type=\"THREE_DIM\"><parameters num_materials=\"1\"><regions>"
    "<marker id=\"0\" material=\"5\"/></regions>"
    "<material id=\"0\" type=\"Guccione\"><coefficients>1,1,1,1,0,1</coefficients></material>"
    "<ninc>1</ninc></parameters></elasticity>\n";
  write_two_regions("bad.xml", bad_mat);
  Open<UpdatedLagrangian> p;
  CHECK_THROWS(p.load("bad.xml"));
}

// ---- Pressure bookkeeping ---------------------------------------------------

TEST_CASE("cavity_pressure: follows the load factor unless asked for the target value")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 1;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "",
                                  "<pressure><node id=\"1\" marker=\"30\" value=\"8.0\"/></pressure>");
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  CHECK_NEAR(p.cavity_pressure(MARKER_LV, true), 8.0, 0.0);
  CHECK_NEAR(p.cavity_pressure(MARKER_RV, true), 0.0, 0.0);

  p.lc.reset();
  p.lc.set_nincs(4);
  p.lc.update();
  CHECK_NEAR(p.cavity_pressure(MARKER_LV, true), 2.0, 1e-14);
  CHECK_NEAR(p.cavity_pressure(MARKER_LV, false), 8.0, 0.0);

  p.set_pressure_Ta(MARKER_LV, 3.0, MARKER_RV, 1.0, arma::vec(8, arma::fill::zeros));
  CHECK_NEAR(p.cavity_pressure(MARKER_LV, false), 3.0, 0.0);
  CHECK_NEAR(p.cavity_pressure(MARKER_RV, false), 1.0, 0.0);
}

TEST_CASE("PV history: file format, repeated open and explicit pressures")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 1;
  b.prefix_xml = elasticity_block("Guccione", kGuccione, "");
  write_tet_box("m.xml", b);
  Open<UpdatedLagrangian> p;
  p.load("m.xml");
  p.set_mesh_units("cm");

  p.open_pv_history("hist");
  CHECK(p.pv_history_is_open());
  p.record_pv_history(1, 10.0, 20.0);
  p.open_pv_history("hist");     // must not truncate
  p.record_pv_history(2, 11.0, 21.0);
  p.close_pv_history();
  CHECK(!p.pv_history_is_open());

  const arma::mat h = read_table("hist_pv_history.dat");
  REQUIRE(h.n_rows == 2);
  CHECK_EQ(h.n_cols, arma::uword(5));
  CHECK_NEAR(h(1, 0), 2.0, 0.0);
  CHECK_NEAR(h(1, 1), 11.0, 0.0);
  CHECK_NEAR(h(1, 3), 21.0, 0.0);
}
