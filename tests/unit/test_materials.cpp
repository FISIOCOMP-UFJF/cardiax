/*
 * Unit tests of src/pdes/materials and src/pdes/mech_utils: the constitutive laws.
 *
 * What is covered and why:
 *  - MaterialData kinematics: C, b, E, J and the isochoric parts.
 *  - The material factory, with a valid name, an unknown name and no parameters.
 *  - Every material: no energy and no stress in the undeformed state, and no stress
 *    under a rigid rotation (frame indifference). A stress at F = I would load the
 *    mechanics with a force that does not exist.
 *  - The analytic stresses and tangents that exist (Saint Venant-Kirchhoff, compressible
 *    neo-Hookean) against the finite-difference ones used in production.
 *  - Guccione and Holzapfel-Ogden: energy against the published formula, and fibre terms
 *    that only act in tension (Holzapfel-Ogden).
 *  - The finite-difference elasticity tensor: major and minor symmetries.
 *  - Active stress: the active energy gives the stress Ta f x f at F = I (I4f of the full
 *    C; the isochoric split is for the passive part only), and the
 *    Regazzoni-Quarteroni stabilisation term behaves as an elastic force of stiffness Ka.
 *  - active_scale: the region -> contractility map, including the legacy rule. A region
 *    silently left with zero active tension is hard to notice in the results.
 *  - push_forward of the stress tensor.
 *  - mech_utils: Kronecker delta, Levi-Civita, polar decomposition, deviatoric operators.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <memory>

#include "pdes/materials/hyperelastic_material.hpp"
#include "pdes/materials/incompressible_material.hpp"
#include "pdes/materials/guccione.hpp"
#include "pdes/materials/holzapfel_ogden.hpp"
#include "pdes/mech_utils.hpp"

using namespace cardiax_test;

namespace
{

Element fibre_element()
{
  const arma::vec3 f = {1, 0, 0}, s = {0, 1, 0}, n = {0, 0, 1};
  return Element(ELEM_TETRA, {0, 1, 2, 3}, f, s, n, 0, 0);
}

// A rotated fibre frame, so that the local-to-global maps are exercised
Element rotated_fibre_element()
{
  const double c = std::cos(0.7), sn = std::sin(0.7);
  const arma::vec3 f = {c, sn, 0}, s = {-sn, c, 0}, n = {0, 0, 1};
  return Element(ELEM_TETRA, {0, 1, 2, 3}, f, s, n, 0, 0);
}

arma::mat33 general_F()
{
  arma::mat33 F = {{1.10, 0.05, 0.02}, {-0.03, 0.95, 0.04}, {0.01, -0.02, 1.02}};
  return F;
}

arma::mat33 rotation()
{
  // 0.5 rad about z, then 0.3 rad about x
  const double a = 0.5, b = 0.3;
  arma::mat33 Rz = {{std::cos(a), -std::sin(a), 0}, {std::sin(a), std::cos(a), 0}, {0, 0, 1}};
  arma::mat33 Rx = {{1, 0, 0}, {0, std::cos(b), -std::sin(b)}, {0, std::sin(b), std::cos(b)}};
  return Rx * Rz;
}

struct NamedMaterial { std::string name; std::vector<double> prm; };

// Stiffness parameters in the range used in testes/ (Ta = 0)
const std::vector<NamedMaterial> kMaterials = {
  {"SaintVenantKirchoff",    {1.0, 2.0}},
  {"CompressibleNeoHookean", {1.0, 2.0}},
  {"Guccione",               {0.65, 6.62, 3.65, 2.65, 0.0, 200.0}},
  {"HolzapfelOgden",         {0.15, 0.3, 0.2, 0.41, 6.0, 11.8, 5.16, 11.3, 0.0, 333.0}},
};

std::unique_ptr<HyperelasticMaterial> make(const NamedMaterial & m)
{
  ElasticityType t = THREE_DIM;
  return std::unique_ptr<HyperelasticMaterial>(HyperelasticMaterial::create(m.name, t, m.prm));
}

std::unique_ptr<HyperelasticMaterial> make(const std::string & name,
                                           const std::vector<double> & prm)
{
  return make(NamedMaterial{name, prm});
}

} // namespace

// ---- MaterialData -----------------------------------------------------------

TEST_CASE("MaterialData: kinematic tensors of a general deformation")
{
  const arma::mat33 F = general_F();
  MaterialData md(fibre_element(), F);

  CHECK_NEAR(md.get_J(), arma::det(F), 1e-14);

  const arma::mat33 C = F.t() * F, b = F * F.t();
  CHECK_NEAR(arma::abs(md.right_cauchy_green() - C).max(), 0.0, 1e-14);
  CHECK_NEAR(arma::abs(md.left_cauchy_green() - b).max(), 0.0, 1e-14);
  CHECK_NEAR(arma::abs(md.lagrangian_strain() - 0.5 * (C - arma::eye(3, 3))).max(), 0.0, 1e-14);

  CHECK_NEAR(arma::det(md.isochoric_def_grad()), 1.0, 1e-13);
  CHECK_NEAR(arma::det(md.isochoric_rcg()), 1.0, 1e-13);
  CHECK_NEAR(arma::det(md.isochoric_lcg()), 1.0, 1e-13);
}

TEST_CASE("MaterialData: carries the element fibres and marker, and starts with no active terms")
{
  Element e(ELEM_TETRA, {0, 1, 2, 3}, arma::vec3({1, 0, 0}), arma::vec3({0, 1, 0}),
            arma::vec3({0, 0, 1}), 4, 0);
  MaterialData md(e, arma::eye(3, 3));
  CHECK_EQ(md.get_marker(), 4);
  CHECK_NEAR(md.fiber()(0), 1.0, 0.0);
  CHECK_NEAR(md.sheet()(1), 1.0, 0.0);
  CHECK_NEAR(md.normal()(2), 1.0, 0.0);
  CHECK_NEAR(md.get_active_stress(), 0.0, 0.0);
  CHECK_NEAR(md.get_active_stiffness(), 0.0, 0.0);
}

// ---- Factory ----------------------------------------------------------------

TEST_CASE("HyperelasticMaterial::create: names, parameters and errors")
{
  for (const NamedMaterial & m : kMaterials)
  {
    auto mat = make(m);
    REQUIRE(mat != nullptr);
    CHECK_NEAR(mat->get_param(0), m.prm[0], 0.0);
  }
  CHECK(make(kMaterials[2])->is_incompressible());
  CHECK(!make(kMaterials[0])->is_incompressible());
  CHECK(make("NeoHookean", {1.0, 2.0, 50.0})->is_incompressible());

  ElasticityType t = THREE_DIM;
  CHECK_THROWS(HyperelasticMaterial::create("NoSuchMaterial", t, std::vector<double>{1.0}));
  CHECK_THROWS(HyperelasticMaterial::create("Guccione", t, std::vector<double>{}));
}

TEST_CASE("HyperelasticMaterial::create: multi-material variant only exists for the cardiac laws")
{
  ElasticityType t = THREE_DIM;
  std::vector<std::vector<double>> prm = {{0.65, 6.62, 3.65, 2.65, 0.0, 0.0},
                                          {6.50, 6.62, 3.65, 2.65, 0.0, 0.0}};
  std::vector<int> map = {0, 1};
  std::unique_ptr<HyperelasticMaterial> g(HyperelasticMaterial::create("Guccione", t, prm, 2, map));
  REQUIRE(g != nullptr);
  CHECK_THROWS(HyperelasticMaterial::create("NeoHookean", t, prm, 2, map));

  // The second material is ten times stiffer: same strain, ten times the energy
  arma::mat33 F = {{1.1, 0, 0}, {0, 1.0 / std::sqrt(1.1), 0}, {0, 0, 1.0 / std::sqrt(1.1)}};
  Element e0(ELEM_TETRA, {0, 1, 2, 3}, arma::vec3({1, 0, 0}), arma::vec3({0, 1, 0}),
             arma::vec3({0, 0, 1}), 0, 0);
  Element e1(ELEM_TETRA, {0, 1, 2, 3}, arma::vec3({1, 0, 0}), arma::vec3({0, 1, 0}),
             arma::vec3({0, 0, 1}), 1, 0);
  MaterialData m0(e0, F), m1(e1, F);
  const arma::mat33 E = m0.lagrangian_strain();
  CHECK_REL(g->strain_energy(&m1, E), 10.0 * g->strain_energy(&m0, E), 1e-10);
}

TEST_CASE_KNOWN_BUG("Guccione with several materials: the volumetric penalty is applied",
                    "src/pdes/materials/guccione.{hpp,cpp}: the multi-material constructor stores "
                    "only 4 entries in parameters, but strain_energy reads K = parameters[5], "
                    "past the end of the vector: the penalty uses whatever is in that memory")
{
  ElasticityType t = THREE_DIM;
  std::vector<std::vector<double>> prm = {{0.65, 6.62, 3.65, 2.65, 0.0, 200.0}};
  std::vector<int> map = {0};
  std::unique_ptr<HyperelasticMaterial> multi(HyperelasticMaterial::create("Guccione", t, prm, 1, map));
  auto single = make("Guccione", prm[0]);

  arma::mat33 F = 1.05 * arma::mat33(arma::eye(3, 3));
  MaterialData md(fibre_element(), F);
  const arma::mat33 E = md.lagrangian_strain();
  CHECK_REL(multi->strain_energy(&md, E), single->strain_energy(&md, E), 1e-12);
}

TEST_CASE_KNOWN_BUG("Neo-Hookean (incompressible): no stress in the undeformed state",
                    "src/pdes/materials/neo_hookean.cpp: strain_energy mixes the isochoric trace "
                    "with the -2 ln(J) term of the compressible law, so S = -mu I at F = I. The "
                    "updated Lagrangian assembly uses this energy through calc_fd_stress")
{
  auto mat = make("NeoHookean", {1.0, 2.0, 50.0});
  MaterialData md(fibre_element(), arma::eye(3, 3));
  arma::mat S(3, 3);
  mat->calc_fd_stress(0, &md, S);
  CHECK_NEAR(arma::abs(S).max(), 0.0, 1e-5);
}

// ---- Properties of every material -------------------------------------------

TEST_CASE("Every material: no energy and no stress in the undeformed state")
{
  for (const NamedMaterial & m : kMaterials)
  {
    auto mat = make(m);
    MaterialData md(fibre_element(), arma::eye(3, 3));
    const arma::mat E(3, 3, arma::fill::zeros);
    CHECK_NEAR(mat->strain_energy(&md, E), 0.0, 1e-13);

    arma::mat S(3, 3);
    mat->calc_fd_stress(0, &md, S);
    CHECK_MSG(arma::abs(S).max() < 1e-5, m.name << " max |S| = " << arma::abs(S).max());
  }
}

TEST_CASE("Every material: a rigid rotation produces no strain energy and no stress")
{
  for (const NamedMaterial & m : kMaterials)
  {
    auto mat = make(m);
    MaterialData md(rotated_fibre_element(), rotation());
    const arma::mat33 E = md.lagrangian_strain();
    CHECK_MSG(arma::abs(E).max() < 1e-14, m.name);
    CHECK_NEAR(mat->strain_energy(&md, E), 0.0, 1e-12);

    arma::mat S(3, 3);
    mat->calc_fd_stress(0, &md, S);
    CHECK_MSG(arma::abs(S).max() < 1e-5, m.name << " max |S| = " << arma::abs(S).max());
  }
}

TEST_CASE("Every material: strain energy grows when the deformation grows")
{
  for (const NamedMaterial & m : kMaterials)
  {
    auto mat = make(m);
    const auto energy = [&](double s) {
      arma::mat33 F = arma::eye(3, 3);
      F(0, 0) = 1.0 + 0.1 * s;
      F(1, 1) = 1.0 - 0.02 * s;
      F(0, 1) = 0.05 * s;
      MaterialData md(fibre_element(), F);
      return mat->strain_energy(&md, md.lagrangian_strain());
    };
    CHECK_MSG(energy(1.0) > 0.0, m.name);
    CHECK_MSG(energy(2.0) > energy(1.0), m.name);
  }
}

TEST_CASE("Every material: the finite-difference elasticity tensor has both symmetries")
{
  for (const NamedMaterial & m : kMaterials)
  {
    auto mat = make(m);
    MaterialData md(rotated_fibre_element(), general_F());
    Tensor4 A;
    A.zero();
    mat->calc_fd_elastensor(&md, A);

    double scale = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
            scale = std::max(scale, std::fabs(A(i, j, k, l)));
    REQUIRE(scale > 0.0);

    double major = 0.0, minor_ = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
          {
            major  = std::max(major,  std::fabs(A(i, j, k, l) - A(k, l, i, j)));
            minor_ = std::max(minor_, std::fabs(A(i, j, k, l) - A(j, i, k, l)));
          }
    CHECK_MSG(major < 1e-3 * scale, m.name << " major " << major << " of " << scale);
    CHECK_MSG(minor_ < 1e-3 * scale, m.name << " minor " << minor_ << " of " << scale);
  }
}

// ---- Analytic against finite differences ------------------------------------

TEST_CASE("Saint Venant-Kirchhoff: energy, analytic PK2 and tangent")
{
  const double mu = 1.5, lam = 2.5;
  auto mat = make("SaintVenantKirchoff", {mu, lam});
  MaterialData md(fibre_element(), general_F());
  const arma::mat33 E = md.lagrangian_strain();
  const double trE = arma::trace(E);

  CHECK_NEAR(mat->strain_energy(&md, E),
             0.5 * lam * trE * trE + mu * arma::accu(E % E), 1e-14);

  arma::mat Sa(3, 3), Sfd(3, 3);
  mat->piola2_stress(&md, Sa);
  mat->calc_fd_stress(0, &md, Sfd);
  const arma::mat Sexact = lam * trE * arma::eye(3, 3) + 2.0 * mu * arma::mat(E);
  CHECK_NEAR(arma::abs(Sa - Sexact).max(), 0.0, 1e-14);
  CHECK_NEAR(arma::abs(Sfd - Sa).max(), 0.0, 1e-6);

  arma::mat Da(6, 6), Dfd(6, 6);
  mat->set_ndim(3);
  mat->mt_elastensor(&md, Da);
  mat->calc_fd_elastensor(&md, Dfd);
  CHECK_NEAR(Da(0, 0), lam + 2 * mu, 0.0);
  CHECK_NEAR(Da(3, 3), mu, 0.0);
  CHECK_MSG(arma::abs(Da - Dfd).max() < 1e-2, "max gap " << arma::abs(Da - Dfd).max());
}

TEST_CASE("Compressible neo-Hookean: the Cauchy stress is the push-forward of dW/dE")
{
  auto mat = make("CompressibleNeoHookean", {1.5, 2.5});
  MaterialData md(fibre_element(), general_F());

  arma::mat Sfd(3, 3), sigma(3, 3);
  mat->calc_fd_stress(0, &md, Sfd);
  mat->cauchy_stress(&md, sigma);
  arma::mat33 S33 = Sfd, sigma_pf;
  const arma::mat33 F = md.get_F();
  mat->push_forward(F, S33, sigma_pf);
  CHECK_MSG(arma::abs(sigma - arma::mat(sigma_pf)).max() < 1e-5,
            "max gap " << arma::abs(sigma - arma::mat(sigma_pf)).max());
}

TEST_CASE("Compressible neo-Hookean: PK2 and tangent, where b and C coincide")
{
  auto mat = make("CompressibleNeoHookean", {1.5, 2.5});
  arma::mat33 F = arma::eye(3, 3);
  F(0, 0) = 1.1; F(1, 1) = 0.95; F(2, 2) = 1.02;
  MaterialData md(fibre_element(), F);

  arma::mat Sa(3, 3), Sfd(3, 3);
  mat->piola2_stress(&md, Sa);
  mat->calc_fd_stress(0, &md, Sfd);
  CHECK_MSG(arma::abs(Sa - Sfd).max() < 1e-5, "max gap " << arma::abs(Sa - Sfd).max());

  arma::mat Da(6, 6, arma::fill::zeros), Dfd(6, 6, arma::fill::zeros);
  mat->mt_elastensor(&md, Da);
  mat->calc_fd_elastensor(&md, Dfd);
  CHECK_MSG(arma::abs(Da - Dfd).max() < 1e-2, "max gap " << arma::abs(Da - Dfd).max());
}

TEST_CASE_KNOWN_BUG("Compressible neo-Hookean: PK2 of a general deformation",
                    "src/pdes/materials/compressible_neo_hookean.cpp, piola2_stress: "
                    "C = md->left_cauchy_green() takes b for C. Total Lagrangian uses the "
                    "finite-difference stress, so this analytic form is not on that path")
{
  auto mat = make("CompressibleNeoHookean", {1.5, 2.5});
  MaterialData md(fibre_element(), general_F());
  arma::mat Sa(3, 3), Sfd(3, 3);
  mat->piola2_stress(&md, Sa);
  mat->calc_fd_stress(0, &md, Sfd);
  CHECK_NEAR(arma::abs(Sa - Sfd).max(), 0.0, 1e-5);
}

TEST_CASE("push_forward: sigma = F S F^T / J, and the identity for F = I")
{
  auto mat = make(kMaterials[0]);
  arma::mat33 S = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
  arma::mat33 sigma;
  mat->push_forward(arma::eye(3, 3), S, sigma);
  CHECK_NEAR(arma::abs(sigma - S).max(), 0.0, 1e-14);

  const arma::mat33 F = general_F();
  mat->push_forward(F, S, sigma);
  const arma::mat33 expected = F * S * F.t() / arma::det(F);
  CHECK_NEAR(arma::abs(sigma - expected).max(), 0.0, 1e-13);
}

// ---- Cardiac materials against their formulas -------------------------------

TEST_CASE("Guccione: energy of an isochoric fibre stretch matches the published law")
{
  const double C = 2.0, bf = 8.0, bt = 2.0, bfs = 4.0;
  auto mat = make("Guccione", {C, bf, bt, bfs, 0.0, 1000.0});

  const double lam = 1.1;
  arma::mat33 F = arma::eye(3, 3);
  F(0, 0) = lam;
  F(1, 1) = F(2, 2) = 1.0 / std::sqrt(lam);
  MaterialData md(fibre_element(), F);
  const arma::mat33 E = md.lagrangian_strain();

  const double Eff = E(0, 0), Ett = E(1, 1), Enn = E(2, 2);
  const double Q = bf * Eff * Eff + bt * (Ett * Ett + Enn * Enn);
  CHECK_NEAR(mat->strain_energy(&md, E), 0.5 * C * (std::exp(Q) - 1.0), 1e-10);
}

TEST_CASE("Guccione: the fibre-sheet shear term uses bfs, and the frame follows the fibres")
{
  const double C = 2.0, bf = 8.0, bt = 2.0, bfs = 4.0;
  auto mat = make("Guccione", {C, bf, bt, bfs, 0.0, 0.0});

  // Pure strain E_fs = g in the fibre frame, with the fibres rotated in the xy plane
  const Element e = rotated_fibre_element();
  arma::mat33 R;
  R.col(0) = e.get_fiber();
  R.col(1) = e.get_trans();
  R.col(2) = e.get_normal();

  const double g = 0.05;
  arma::mat33 El(arma::fill::zeros);
  El(0, 1) = El(1, 0) = g;
  const arma::mat33 E = R * El * R.t();

  MaterialData md(e, arma::eye(3, 3));
  const double Q = bfs * (g * g + g * g);
  CHECK_NEAR(mat->strain_energy(&md, E), 0.5 * C * (std::exp(Q) - 1.0), 1e-12);
}

TEST_CASE("Holzapfel-Ogden: isotropic part matches the published law")
{
  const double a = 2.0, b = 3.0;
  auto mat = make("HolzapfelOgden", {a, 0, 0, 0, b, 0, 0, 0, 0.0, 1000.0});
  arma::mat33 F = arma::eye(3, 3);
  F(0, 0) = 1.1;
  F(1, 1) = F(2, 2) = 1.0 / std::sqrt(1.1);
  MaterialData md(fibre_element(), F);
  const double I1 = arma::trace(F.t() * F);
  CHECK_NEAR(mat->strain_energy(&md, md.lagrangian_strain()),
             a / (2 * b) * (std::exp(b * (I1 - 3.0)) - 1.0), 1e-10);
}

TEST_CASE("Holzapfel-Ogden: fibre and sheet terms act in tension only")
{
  const double af = 5.0, bf = 4.0, as = 3.0, bs = 2.0;
  auto base = make("HolzapfelOgden", {0.1, 0, 0, 0, 1.0, 0, 0, 0, 0.0, 0.0});
  auto mat  = make("HolzapfelOgden", {0.1, af, as, 0, 1.0, bf, bs, 0, 0.0, 0.0});

  const auto extra = [&](double lam) {
    arma::mat33 F = arma::eye(3, 3);
    F(0, 0) = lam;
    F(1, 1) = 1.0 / lam;     // det F = 1 with an inverse sheet stretch
    MaterialData md(fibre_element(), F);
    const arma::mat33 E = md.lagrangian_strain();
    return mat->strain_energy(&md, E) - base->strain_energy(&md, E);
  };

  const double lam = 1.15;
  const double x_f = lam * lam - 1.0;
  CHECK_NEAR(extra(lam), af / (2 * bf) * (std::exp(bf * x_f * x_f) - 1.0), 1e-10);

  const double x_s = 1.0 / (lam * lam) - 1.0;
  (void)x_s;  // compressed sheet direction: contributes nothing above
  const double l2 = 1.0 / lam;
  const double x_s2 = lam * lam - 1.0;
  CHECK_NEAR(extra(l2), as / (2 * bs) * (std::exp(bs * x_s2 * x_s2) - 1.0), 1e-10);
}

TEST_CASE("Holzapfel-Ogden: b = 0 reduces to the linear form, not to a division by zero")
{
  auto mat = make("HolzapfelOgden", {2.0, 0, 0, 0, 0.0, 0, 0, 0, 0.0, 0.0});
  arma::mat33 F = arma::eye(3, 3);
  F(0, 0) = 1.1;
  F(1, 1) = F(2, 2) = 1.0 / std::sqrt(1.1);
  MaterialData md(fibre_element(), F);
  const double I1 = arma::trace(F.t() * F);
  CHECK_NEAR(mat->strain_energy(&md, md.lagrangian_strain()), 1.0 * (I1 - 3.0), 1e-12);
}

TEST_CASE("Volumetric penalty: both cardiac laws charge K/2 (J - 1)^2 for a pure dilation")
{
  const double K = 100.0, s = 1.02, J = s * s * s;
  for (const NamedMaterial & m : {NamedMaterial{"Guccione", {0.65, 6.62, 3.65, 2.65, 0.0, K}},
                                  NamedMaterial{"HolzapfelOgden", {0.15, 0, 0, 0, 6.0, 0, 0, 0, 0.0, K}}})
  {
    auto mat = make(m);
    arma::mat33 F = s * arma::eye(3, 3);
    MaterialData md(fibre_element(), F);
    const double total = mat->strain_energy(&md, md.lagrangian_strain());

    // Isochoric part of the pure dilation vanishes in Holzapfel-Ogden; Guccione is
    // written on the full strain, so only require the penalty to be included
    if (m.name == "HolzapfelOgden")
      CHECK_NEAR(total, 0.5 * K * (J - 1.0) * (J - 1.0), 1e-10);
    else
      CHECK(total >= 0.5 * K * (J - 1.0) * (J - 1.0));
  }
}

// ---- Active stress ----------------------------------------------------------

TEST_CASE("Active stress: at F = I it is Ta f x f")
{
  const double Ta = 7.0;
  for (const NamedMaterial & m :
       {NamedMaterial{"Guccione", {0.65, 6.62, 3.65, 2.65, 0.0, 0.0}},
        NamedMaterial{"HolzapfelOgden", {0.15, 0, 0, 0, 6.0, 0, 0, 0, 0.0, 0.0}}})
  {
    auto mat = make(m);
    const Element e = rotated_fibre_element();
    MaterialData md(e, arma::eye(3, 3));
    md.set_active_stress(Ta);

    arma::mat S(3, 3);
    mat->calc_fd_stress(0, &md, S);

    const arma::vec3 f = e.get_fiber();
    const arma::mat expected = Ta * (f * f.t());
    CHECK_MSG(arma::abs(S - expected).max() < 1e-5,
              m.name << " max gap " << arma::abs(S - expected).max());
  }
}

TEST_CASE("Active stress: zero tension contributes nothing")
{
  auto mat = make("Guccione", {0.65, 6.62, 3.65, 2.65, 0.0, 200.0});
  MaterialData md(fibre_element(), general_F());
  const auto * im = static_cast<IncompressibleMaterial *>(mat.get());
  md.set_active_stress(0.0);
  CHECK_NEAR(im->active_strain_energy(&md, md.lagrangian_strain()), 0.0, 0.0);
}

TEST_CASE("Active stabilisation: psi = Ka/2 (lambda - lambda_prev)^2, and off when Ka = 0")
{
  auto mat = make("Guccione", {0.65, 6.62, 3.65, 2.65, 0.0, 200.0});
  const auto * im = static_cast<IncompressibleMaterial *>(mat.get());

  arma::mat33 F = arma::eye(3, 3);
  F(0, 0) = 1.1;
  MaterialData md(fibre_element(), F);
  const arma::mat33 E = md.lagrangian_strain();

  md.set_active_stabilization(0.0, 1.0);
  CHECK_NEAR(im->active_strain_energy(&md, E), 0.0, 0.0);

  const double Ka = 50.0, lprev = 1.04;
  md.set_active_stabilization(Ka, lprev);
  CHECK_NEAR(im->active_strain_energy(&md, E), 0.5 * Ka * (1.1 - lprev) * (1.1 - lprev), 1e-12);

  // Stationary at lambda = lambda_prev
  md.set_active_stabilization(Ka, 1.1);
  CHECK_NEAR(im->active_strain_energy(&md, E), 0.0, 1e-20);
}

TEST_CASE("HyperelasticMaterial: active fields are held and the stabilisation switch follows them")
{
  auto mat = make(kMaterials[2]);
  mat->allocate_Ta(4);
  CHECK_EQ(mat->get_Ta().n_elem, arma::uword(4));

  mat->set_Ta(arma::vec({1.0, 2.0}));
  CHECK_NEAR(mat->get_Ta()(1), 2.0, 0.0);

  CHECK(!mat->has_active_stabilization());
  mat->set_active_stabilization(arma::vec({1.0}), arma::vec({1.0}));
  CHECK(mat->has_active_stabilization());
  mat->set_active_stabilization(arma::vec({1.0, 1.0}), arma::vec({1.0}));
  CHECK(!mat->has_active_stabilization());
  mat->set_active_stabilization(arma::vec(), arma::vec());
  CHECK(!mat->has_active_stabilization());
}

TEST_CASE("active_scale: legacy rule when nothing is declared, table when it is")
{
  auto mat = make(kMaterials[2]);
  CHECK(!mat->has_active_scale());
  CHECK_NEAR(mat->active_scale(0), 1.0, 0.0);
  CHECK_NEAR(mat->active_scale(1), 0.0, 0.0);

  mat->set_active_scale({0, 2, 1}, {1.0, 0.5, 0.0});
  CHECK(mat->has_active_scale());
  CHECK_NEAR(mat->active_scale(0), 1.0, 0.0);
  CHECK_NEAR(mat->active_scale(1), 0.0, 0.0);   // marker 1 -> material 2
  CHECK_NEAR(mat->active_scale(2), 0.5, 0.0);   // marker 2 -> material 1

  // Markers outside the table are passive
  CHECK_NEAR(mat->active_scale(3), 0.0, 0.0);
  CHECK_NEAR(mat->active_scale(-1), 0.0, 0.0);
}

// ---- mech_utils -------------------------------------------------------------

TEST_CASE("mech_utils: Kronecker delta and the symmetric perturbation")
{
  CHECK_EQ(delta_kronecker(1, 1), 1);
  CHECK_EQ(delta_kronecker(1, 2), 0);

  arma::mat33 d;
  delta_tensor(0, 1, 1.0e-3, d);
  CHECK_NEAR(d(0, 1), 0.5e-3, 0.0);
  CHECK_NEAR(d(1, 0), 0.5e-3, 0.0);
  CHECK_NEAR(d(0, 0), 0.0, 0.0);

  delta_tensor(2, 2, 1.0e-3, d);
  CHECK_NEAR(d(2, 2), 1.0e-3, 0.0);
}

TEST_CASE("mech_utils: Levi-Civita symbol")
{
  CHECK_EQ(levi(0, 1, 2), 1);
  CHECK_EQ(levi(1, 2, 0), 1);
  CHECK_EQ(levi(2, 0, 1), 1);
  CHECK_EQ(levi(2, 1, 0), -1);
  CHECK_EQ(levi(1, 0, 2), -1);
  CHECK_EQ(levi(0, 2, 1), -1);
  CHECK_EQ(levi(0, 0, 1), 0);
}

TEST_CASE("mech_utils: polar decomposition F = R U with R orthogonal and U symmetric")
{
  const arma::mat33 F = general_F();
  arma::mat33 R, U;
  polar_decomposition(F, R, U);
  CHECK_NEAR(arma::abs(arma::mat33(R * U) - F).max(), 0.0, 1e-13);
  CHECK_NEAR(arma::abs(arma::mat33(R.t() * R) - arma::mat33(arma::eye(3, 3))).max(), 0.0, 1e-13);
  CHECK_NEAR(arma::det(R), 1.0, 1e-13);
  CHECK_NEAR(arma::abs(U - U.t()).max(), 0.0, 1e-13);

  // A pure rotation has U = I
  polar_decomposition(rotation(), R, U);
  CHECK_NEAR(arma::abs(U - arma::mat33(arma::eye(3, 3))).max(), 0.0, 1e-13);
}

TEST_CASE_KNOWN_BUG("mech_utils: the spatial deviatoric operator removes the trace",
                    "src/pdes/mech_utils.cpp: contract(A,B) returns the entrywise product "
                    "matrix instead of the scalar A:B, so sp_dev_operator(T) is "
                    "T - diag(T)/3 instead of T - tr(T) I/3. Used by the analytic "
                    "HolzapfelOgden::deviatoric_stress")
{
  arma::mat33 T = {{1, 2, 3}, {2, 5, 4}, {3, 4, 9}};
  const arma::mat33 dev = sp_dev_operator(T);
  CHECK_NEAR(arma::trace(dev), 0.0, 1e-14);
  CHECK_NEAR(dev(0, 0), 1.0 - 5.0, 1e-14);
  CHECK_NEAR(arma::abs(sp_dev_operator(arma::eye(3, 3))).max(), 0.0, 1e-14);
}

TEST_CASE_KNOWN_BUG("mech_utils: the material deviatoric operator removes the part along C^-1",
                    "src/pdes/mech_utils.cpp: same cause as sp_dev_operator, contract(T,C) "
                    "is not a scalar, so the deviator of C^-1 does not vanish")
{
  arma::mat33 F = general_F();
  const arma::mat33 C = F.t() * F;
  const arma::mat33 Cinv = arma::inv(C);
  CHECK_NEAR(arma::abs(mt_dev_operator(Cinv, C)).max(), 0.0, 1e-12);
}
