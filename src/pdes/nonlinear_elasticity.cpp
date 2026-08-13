#include "nonlinear_elasticity.hpp"
#include "util/pugixml.hpp"
#include <map>
#include <algorithm>
#include <cstddef>

//#define DEBUG_PRESSURE
#define INCs_TA 0

NonlinearElasticity::NonlinearElasticity()
  : timer(), lc(), output_step(false), material(0),
    log(), writer(&msh), first_step(true)
{
  // define default parameters
  parameters.add("tol_energy",1.0e-8);
  parameters.add("tol_residual",1.0e-6);
  parameters.add("tol_displacement",1.0e-16);

  //parameters.add("tol_energy",1.0e-6);
  //parameters.add("tol_residual",1.0e-5);
  //parameters.add("tol_displacement",1.0e-16);

  parameters.add("energy_norm0",0.0);
  parameters.add("residual_norm0",0.0);

  // create log file using the format:
  // its1 its2 ... itsN volume
  log.create("log_elasticity.dat");
}

NonlinearElasticity::~NonlinearElasticity()
{
  if (material != NULL)
    delete material;

  for(uint i=0; i<vecF.size(); i++)
    delete vecF[i];
}

void NonlinearElasticity::apply_boundary(petsc::Matrix & Ks)
{
  int nd = msh.get_n_dim();
  int np = msh.get_n_points();
  int nod, idx, dir;
  double val, dlamb = lc.load_step();
  NodalData inf;

  // prescribe dirichlet boundary conditions
  if (fixed_nodes_map.size())
  {
    std::map<int,double> boundary_values;
    std::map<int,NodalData>::iterator it;
    for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
    {
      nod = it->first;
      inf = it->second;
      dir = inf.first;
      val = dlamb * inf.second;
      //idx = (dir * np) + nod;
      idx = (nod*nd) + dir;
      assert(dir+1 <= nd);
      boundary_values.insert( std::pair<int,double>(idx,val) );
    }
    FETools::apply_boundary_values(boundary_values, Ks);
  }
}

void NonlinearElasticity::assemble_traction()
{
  MixedFiniteElement * bfe = fespace.create_boundary_FE();
  if (bfe != NULL && neumann_map.size() > 0)
  {
    cout << "Computing traction boundary conditions" << endl;
    int nu  = bfe->get_ndofs_u();
    int nb = msh.get_n_boundary_elements();
    arma::mat belmat(nu,nu);
    arma::vec belvec(nu);
    vector<int> bdof;
    for(int i=0; i<nb; i++)
    {
      calc_neumann_elvec(i, bfe, belvec);
      fespace.get_boundary_element_dofs_u(i,bdof);

      for(int k=0; k<nu; k++)
        fext(bdof[k]) += belvec(k);
    }
  }
  delete bfe;
}

void NonlinearElasticity::assemble_pressure()
{
  MixedFiniteElement * bfe = fespace.create_boundary_FE();

  if (bfe != NULL && pressure_map.size() > 0)
  {
    cout << " Assembling pressure component of stiffness matrix" << endl;
    int nu  = bfe->get_ndofs_u();
    int nb = msh.get_n_boundary_elements();
    double xlamb = lc.load();
    arma::mat belmat(nu,nu);
    arma::vec belvec(nu);
    vector<int> bdof;

    for(int i=0; i<nb; i++)
    {
      fespace.get_boundary_element_dofs_u(i,bdof);
      calc_pforce_kpress(i,bfe,bdof,belvec,belmat);

      // assembles the pressure component of the stiffness matrix
      for(int j=0; j<nu; j++)
        for(int k=0; k<nu; k++)
          K.add(bdof[j], bdof[k], belmat(j,k));

      // assembles the nodal forces due to normal pressure
      for(int k=0; k<nu; k++)
      {

        if(ldgof[bdof[k]])
        {
          r.add(bdof[k], xlamb * belvec(k));
          //fext(bdof[k]) += belvec(k);
          tload(bdof[k]) += belvec(k);
        }
        else
          react.add(bdof[k], -xlamb * belvec(k));
      }
    }
  }

  K.assemble();
  r.assemble();
  apply_boundary(K);


  delete bfe;
}

void NonlinearElasticity::set_mesh_units(const std::string & units)
{
  if (units == "m" || units == "meter" || units == "meters")
  {
    vol_to_mL = 1.0e6;   // 1 m^3   = 1e6 mL
    unit_name = "m";
  }
  else if (units == "cm" || units == "centimeter" || units == "centimeters")
  {
    vol_to_mL = 1.0;     // 1 cm^3  = 1 mL
    unit_name = "cm";
  }
  else if (units == "mm" || units == "millimeter" || units == "millimeters")
  {
    vol_to_mL = 1.0e-3;  // 1 mm^3  = 1e-3 mL
    unit_name = "mm";
  }
  else
  {
    cout << " Warning: unknown mesh unit '" << units
         << "'; expected m, cm or mm. Keeping " << unit_name << "." << endl;
    return;
  }
  units_detected = true;
  cout << " Mesh length unit set to " << unit_name
       << " (volume scale to mL: " << vol_to_mL << ")" << endl;
}

void NonlinearElasticity::detect_mesh_units()
{
  // Infer the length unit from the size of the reference geometry. A human
  // heart spans roughly 10 cm, so the bounding-box diagonal is about
  //   0.1  in m,  10  in cm,  100  in mm
  // These ranges are far enough apart to separate reliably for any
  // physiological heart size (roughly 5 to 25 cm).
  if (x0.empty())
  {
    cout << " Warning: cannot detect mesh units, reference coordinates empty."
         << endl;
    return;
  }

  arma::vec3 lo = x0[0], hi = x0[0];
  for (size_t i = 1; i < x0.size(); i++)
    for (int d = 0; d < 3; d++)
    {
      lo(d) = std::min(lo(d), x0[i](d));
      hi(d) = std::max(hi(d), x0[i](d));
    }

  double diag = arma::norm(hi - lo, 2);

  const char * guess;
  if      (diag <  1.0)  guess = "m";
  else if (diag < 30.0)  guess = "cm";
  else                   guess = "mm";

  cout << " Mesh bounding-box diagonal = " << diag
       << " -> assuming length unit '" << guess << "'" << endl;

  set_mesh_units(guess);

  // Sanity check: report the implied heart size in cm.
  double diag_cm = diag * (vol_to_mL == 1.0e6 ? 100.0
                          : vol_to_mL == 1.0  ?   1.0 : 0.1);
  cout << "   implied geometry size: " << diag_cm << " cm";
  if (diag_cm < 3.0 || diag_cm > 40.0)
    cout << "  <-- WARNING: outside the expected range for a heart,"
         << " unit detection may be wrong; use set_mesh_units() to override";
  cout << endl;
}

double NonlinearElasticity::volume_scale_to_mL()
{
  if (!units_detected)
    detect_mesh_units();
  return vol_to_mL;
}

double NonlinearElasticity::total_volume_cavity(const int cavity_marker)
{
  MixedFiniteElement * bfe = fespace.create_boundary_FE();
  double total_endo_volume=0;

  if (bfe != NULL)
  {
    int nb = msh.get_n_boundary_elements();
    vector<int> bdof;
    for(int i=0; i<nb; i++)
    {
      fespace.get_boundary_element_dofs_u(i, bdof);
      total_endo_volume += calc_cavity_volume(i, bfe, bdof, cavity_marker);
    }
  }
  delete bfe;

  // convert from the mesh length unit cubed to mL
  return total_endo_volume * volume_scale_to_mL();
}

double NonlinearElasticity::volume_myocardium_mL()
{
  return calc_volume() * volume_scale_to_mL();
}

void NonlinearElasticity::fiber_stretch_elements(arma::vec & lam_e)
{
  const int nelem = msh.get_n_elements();
  const int nint  = get_num_integration_points();

  lam_e.set_size(nelem);

  // vecF is sized nelem*nint in init(). If it is empty the problem has not
  // been initialised yet; report the undeformed state rather than reading
  // past the end.
  if (vecF.size() < static_cast<size_t>(nelem * nint) || nint < 1)
  {
    lam_e.ones();
    return;
  }

  for (int e = 0; e < nelem; e++)
  {
    // f0 is the fibre direction in the REFERENCE configuration, which is
    // what I4f = f0.C.f0 is defined with. Do not push it forward here.
    const arma::vec3 & f0 = msh.get_element(e).get_fiber();

    double lam = 0.0;
    for (int q = 0; q < nint; q++)
    {
      // lambda_f = sqrt(f0.F^T.F.f0) = ||F f0||
      const arma::vec3 Ff0 = (*vecF[e * nint + q]) * f0;
      lam += arma::norm(Ff0, 2);
    }
    lam_e(e) = lam / static_cast<double>(nint);
  }
}

void NonlinearElasticity::fiber_stretch_nodes(arma::vec & lam_n)
{
  arma::vec lam_e;
  fiber_stretch_elements(lam_e);

  const int nelem   = msh.get_n_elements();
  const int npoints = msh.get_n_points();
  const int nen     = msh.get_nen();

  lam_n.zeros(npoints);
  arma::vec count(npoints, arma::fill::zeros);

  std::vector<int> pnums(nen);
  for (int e = 0; e < nelem; e++)
  {
    msh.get_element_pt_nums(e, pnums);
    for (int i = 0; i < nen; i++)
    {
      lam_n(pnums[i]) += lam_e(e);
      count(pnums[i]) += 1.0;
    }
  }

  // A node with no incident element cannot be assigned a stretch; leaving it
  // at 1 keeps the cell model at its neutral, undeformed behaviour instead of
  // feeding it a division by zero.
  for (int i = 0; i < npoints; i++)
    lam_n(i) = (count(i) > 0.0) ? lam_n(i) / count(i) : 1.0;
}

double NonlinearElasticity::volume_LV()
{
  return total_volume_cavity(MARKER_LV);
}

double NonlinearElasticity::volume_RV()
{
  return total_volume_cavity(MARKER_RV);
}

double NonlinearElasticity::cavity_pressure(const int cavity_marker,
                                            bool apply_load_factor)
{
  // Pressure prescribed for the given cavity. Returns 0 if no pressure BC
  // is defined for that marker.
  //
  // apply_load_factor = true  : scales by lc.load(), i.e. the pressure
  //     actually applied at the current increment. Correct for a monotonic
  //     ramp such as passive filling driven by NonlinearElasticity::run().
  // apply_load_factor = false : returns the prescribed target pressure as
  //     set by set_pressure_Ta(). Correct when an outer driver (e.g.
  //     CardiacElectromechanic) prescribes the pressure at every step.
  std::map<int,double>::iterator it = pressure_map.find(cavity_marker);
  if (it == pressure_map.end())
    return 0.0;

  return apply_load_factor ? it->second * lc.load() : it->second;
}

void NonlinearElasticity::open_pv_history(const string & basename)
{
  // Guard against truncating the file when solve() is called repeatedly,
  // as happens when an outer driver steps through the cardiac cycle.
  if (pv_file.is_open())
    return;

  string aux = basename + string("_pv_history.dat");
  pv_file.open(aux.c_str());

  if (pv_file.is_open())
  {
    pv_file << "# increment  P_LV  V_LV  P_RV  V_RV\n";
    pv_file << std::scientific << std::setprecision(8);
    pv_counter = 0;
  }
  else
    cout << "Warning: could not open PV history file " << aux << endl;
}

void NonlinearElasticity::record_pv_history(int increment)
{
  // Pressure taken from the load ramp: valid for a monotonic increment
  // loop such as passive filling.
  record_pv_history(increment,
                    cavity_pressure(MARKER_LV, true),
                    cavity_pressure(MARKER_RV, true));
}

void NonlinearElasticity::record_pv_history(int increment,
                                            double plv, double prv)
{
  if (!pv_file.is_open())
    return;

  double vlv = volume_LV();
  double vrv = volume_RV();

  pv_file << increment << " "
          << plv << " " << vlv << " "
          << prv << " " << vrv << "\n";
  pv_file.flush();

  pv_counter++;

  cout << "  [inc " << increment << "]"
       << " LV: P=" << plv << " V=" << vlv
       << " | RV: P=" << prv << " V=" << vrv << endl;
}

void NonlinearElasticity::close_pv_history()
{
  if (pv_file.is_open())
    pv_file.close();
}


void NonlinearElasticity::body_forces()
{
  if(bforce.size() != 0)
  {
    int ndim  = msh.get_n_dim();
    int ndof  = msh.get_nen() * ndim;
    int nelem = msh.get_n_elements();
    int nubf  = ndof/ndim;
    arma::vec elvec(ndof);
    std::vector<int> dnums;

    MixedFiniteElement * fe = fespace.createFE();
    Quadrature * qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    //Quadrature * qd = Quadrature::create(2, fe->get_type());
    cout << " Computing body forces (or forcing term)" << endl;

    for(int i=0; i<nelem; i++)
    {
      double detJxW;
      arma::vec shape;
      arma::mat dshape;
      arma::vec3 pt, qpt;

      Mapping em = fe->get_mapping(i);

      elvec.zeros();
      fespace.get_element_dofs_u(i,dnums);

      for(int q=0; q<qd->get_num_ipoints(); q++)
      {
        qpt = qd->get_point(q);
        fe->calc_shape_u(qpt, shape);
        fe->calc_deriv_shape_u(qpt, dshape);
        em.calc_jacobian(dshape, nubf);
        detJxW = qd->get_weight(q) * em.get_det_jacobian();
        pt = em.map_point(qpt, shape);

        // compute forcing vector
        //arma::vec f = forcing_term(pt, *material);
        arma::vec3 f;
        for(int i=0; i<3; i++) f(i) = bforce[i];

        for(int i=0; i<ndim; i++)
          for(int j=0; j<nubf; j++)
            elvec(j+i*nubf) += f(i) * shape(j+i*nubf) * detJxW;
      }

      // assemble the body force into the external load vector
      //   if dof is active, adds the forces due
      //   to gravity/body force to the global vector fext
      for(int k=0; k<ndof; k++)
        if(ldgof[dnums[k]])
          fext(dnums[k]) += elvec(k);
    }
    delete qd;
    delete fe;
  }
}

double NonlinearElasticity::calc_energy()
{
  const int ndofs = fespace.get_ndofs();
  double energy_norm = 0.0;
  arma::vec resid(ndofs);
  arma::vec displ(ndofs);
  r.get_data(resid.memptr());
  u.get_data(displ.memptr());

  for(int i=0; i<ndofs; i++)
    if(ldgof[i])
      energy_norm += resid(i) * displ(i);

  return abs(energy_norm);
}

double NonlinearElasticity::calc_volume(bool update)
{
  int ndof, nnode;
  double totv = 0.0;
  MixedFiniteElement * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(0, fe->get_type());
  //Quadrature * qd = Quadrature::create(2, fe->get_type());
  ndof  = fe->get_ndofs_u();
  nnode = fe->get_nnode();

  for(int i=0; i<msh.get_n_elements(); i++)
  {
    double detJxW, evol = 0.0;
    arma::mat dshape;
    arma::vec3 qpt;
    std::vector<arma::vec3> xe(ndof);
    get_element_x (i, xe);

    Mapping em = fe->get_mapping(i, xe);

    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      qpt = qd->get_point(q);
      fe->calc_deriv_shape_u(qpt,dshape);
      em.calc_jacobian(dshape, nnode);
      detJxW = qd->get_weight(q) * em.get_det_jacobian();
      evol = evol + detJxW;
    }

    totv = totv + evol;

    // initial element volume
    if(update) vol0[i] = evol;
  }

  delete qd;
  delete fe;

  return totv;
}

void NonlinearElasticity::calc_neumann_elvec(const int eindex, const MxFE * fe,
                                             arma::vec & elvec)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;
  int index;
  double detJxW;
  arma::vec shape;
  arma::mat dshape;
  arma::vec3 f;
  std::map<int,arma::vec3>::iterator it;

  SurfaceMapping sm = fe->get_boundary_mapping(eindex);
  index = sm.get_index();
  it = neumann_map.find(index);

  f.zeros();
  elvec.zeros();

  // get traction vector on surface with number=index
  if(neumann_map.find(index) != neumann_map.end())
  {
    Quadrature * qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    //Quadrature * qd = Quadrature::create(2, fe->get_type());

    f = it->second;

    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe->calc_shape_u(qd->get_point(q),shape);
      fe->calc_deriv_shape_u(qd->get_point(q),dshape);
      sm.calc_jacobian(dshape, nubf);

      detJxW = qd->get_weight(q) * sm.get_det_jacobian();

      arma::vec3 pt = sm.map_point(qd->get_point(q)[0]);

      for(int i=0; i<ndim; i++)
        for(int j=0; j<nubf; j++)
	        elvec(j+i*nubf) += f(i) * shape(j+i*nubf) * detJxW;
    }
    delete qd;
  }
}

void NonlinearElasticity::elem_pforce(const int elem_id, const MxFE * fe,
                                      const std::vector<int> & bdof,
                                      arma::vec & belvec,
                                      bool & is_spring)
{
  is_spring = false;
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int neln = ndof/ndim;
  int index;
  double detJxW, press;
  arma::vec shape;
  arma::mat dshape, dxis(3,2);

  // get current coordinates of the element
  std::vector<arma::vec3> xe(neln);
  get_boundary_element_x (elem_id, xe);

  // boundary element mapping
  SurfaceMapping sm = fe->get_boundary_mapping(elem_id, xe);

  std::map<int,double>::iterator it;
  index = sm.get_index();
  it = pressure_map.find(index);

  if(pressure_map.find(index) != pressure_map.end())
  {
    belvec.zeros();
    press = it->second;

    // create quadrature rule (need to check if order is ok)
    Quadrature * qd;
    qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    //qd = Quadrature::create(2, fe->get_type());

    // quadrature loop
    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe->calc_shape_u(qd->get_point(q),shape);
      fe->calc_deriv_shape_u(qd->get_point(q),dshape);
      sm.calc_jacobian(dshape, neln);
      // trick
      detJxW = qd->get_weight(q);

      // evaluates dx/dxi, dx/deta
      // for 2D set dx/deta = (0,0,-1)
      dxis.zeros();
      if(ndim==2) dxis(2,1) = -1.0;
      for(int id=0; id<ndim; id++)
        for(int jd=0; jd<ndim-1; jd++)
          for(int in=0; in<neln; in++)
          {
            //int ip = bdof[in];
	    int ip = bdof[in]/ndim;
            dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
          }

      //
      // computes pressure force (nodal forces - subroutine pforce)
      //
      arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));
      
      //std::cout << xnorm <<std::endl;

      for(int in=0; in<neln; in++)
      {
        double vl = press * shape(in) * detJxW;
        // assemble boundary force to elemental vector (in + id*neln)
        // belvec data layout
        // [ bloco_x, bloco_y, bloco_z]
        // [n1_x, n2_x, ..., nN_x, n1_y, ..., nN_y, n1_z, ..., nN_z]
        for(int id=0; id<ndim; id++)
          belvec(in + id*neln) += xnorm(id) * vl;
      }
    }
    // end of integration loop
    delete qd;
  }
  else
  {
    it = spring_map.find(index);
    if(it != spring_map.end())
    {
      belvec.zeros();
      press = it->second;
      is_spring = true;
      int stype = spring_type_map[index]; // 0 = normal, 1 = isotropic

      // create quadrature rule
      Quadrature * qd;
      qd = Quadrature::create(fe->get_order_u(), fe->get_type());

      // quadrature loop
      for(int q=0; q<qd->get_num_ipoints(); q++)
      {
        fe->calc_shape_u(qd->get_point(q),shape);
        fe->calc_deriv_shape_u(qd->get_point(q),dshape);
        sm.calc_jacobian(dshape, neln);
        detJxW = qd->get_weight(q);

        // Compute dx0/dxi, dx0/deta using REFERENCE coordinates
        dxis.zeros();
        if(ndim==2) dxis(2,1) = -1.0;
        for(int id=0; id<ndim; id++)
          for(int jd=0; jd<ndim-1; jd++)
            for(int in=0; in<neln; in++)
            {
              int ip = bdof[in]/ndim;
              dxis(id,jd) = dxis(id,jd) + x0[ip][id] * dshape(in,jd);
            }

        arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));
        double norm_xnorm = arma::norm(xnorm, 2);

        // Interpolate displacement u_h at quadrature point
        arma::vec3 uh;
        uh.zeros();
        for(int jn=0; jn<neln; jn++)
        {
          int jp = bdof[jn]/ndim;
          for(int jd=0; jd<ndim; jd++)
            uh(jd) += shape(jn) * (x[jp][jd] - x0[jp][jd]);
        }

        if(stype == 0)
        {
          // Normal spring (restoring): traction = -K_epi * (u . n_hat_0) * n_hat_0
          // belvec = -K_epi * (u_h . cross_0) * cross_0 / ||cross_0|| * N_i * w
          double un = arma::dot(uh, xnorm);
          for(int in=0; in<neln; in++)
            for(int id=0; id<ndim; id++)
              belvec(in + id*neln) += -press * un * xnorm(id)
                                      / norm_xnorm * shape(in) * detJxW;
        }
        else
        {
          // Isotropic spring (restoring): PN + ku = 0  =>  traction = -ku
          // belvec = -k * u_h * N_i * dGamma_0
          double dGamma0 = norm_xnorm * detJxW;
          for(int in=0; in<neln; in++)
            for(int id=0; id<ndim; id++)
              belvec(in + id*neln) += -press * uh(id) * shape(in) * dGamma0;
        }
      }
      // end of integration loop
      delete qd;
    }
    else
    {
      error("Boundary index for pressure load not found");  
    }
  }
}

void NonlinearElasticity::elem_kpress(const int elem_id, const MxFE * fe,
                                      const std::vector<int> & bdof,
                                      arma::mat & belmat)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int neln = ndof/ndim;
  int index;
  double detJxW, press;
  arma::vec shape;
  arma::mat dshape, dxis(3,2);

  // get current coordinates of the element
  std::vector<arma::vec3> xe(neln);
  get_boundary_element_x (elem_id, xe);

  // boundary element mapping
  SurfaceMapping sm = fe->get_boundary_mapping(elem_id, xe);

  std::map<int,double>::iterator it;
  index = sm.get_index();
  //std::cout << "\n" << index << " ";
  it = pressure_map.find(index);

  if(it != pressure_map.end())
  {
    belmat.zeros();
    press = it->second;

    // create quadrature rule (need to check if order is ok)
    Quadrature * qd;
    qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    //qd = Quadrature::create(2, fe->get_type());

    // quadrature loop
    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe->calc_shape_u(qd->get_point(q),shape);
      fe->calc_deriv_shape_u(qd->get_point(q),dshape);
      sm.calc_jacobian(dshape, neln);
      // trick
      detJxW = qd->get_weight(q);

      // evaluates dx/dxi, dx/deta
      // for 2D set dx/deta = (0,0,-1)
      dxis.zeros();
      if(ndim==2) dxis(2,1) = -1.0;
      for(int id=0; id<ndim; id++)
        for(int jd=0; jd<ndim-1; jd++)
          for(int in=0; in<neln; in++)
          {
            //int ip = bdof[in];
	          int ip = bdof[in]/ndim;
            dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
          }

      //
      // computes pressure component for stiffness matrix (subroutine kpress)
      //
      double sum, val, val2, val3;
      double apress = press * lc.load_step();
      val = -(apress * detJxW)/2.0;

      for(int in=0; in<neln; in++)
        for(int jn=0; jn<neln; jn++)
        {
          val2 = val * ( shape(in)*dshape(jn,0)
                       - shape(jn)*dshape(in,0) );
          val3 = 0.0;
          if(ndim==3) val3 = val * ( shape(in)*dshape(jn,1)
                                   - shape(jn)*dshape(in,1) );
          for(int id=0; id<ndim; id++)
            for(int jd=0; jd<ndim; jd++)
            {
              sum = 0.0;
              for(int kd=0; kd<3; kd++)
                sum += levi(id,jd,kd) * (val2*dxis(kd,1) - val3*dxis(kd,0));

              // write in the elemental pressure stiffness matrix (in + id*neln)
              // belmat data layout
              // [ bloco_x, bloco_y, bloco_z] , [..]
              // [n1_x, n2_x, ..., nN_x, n1_y, ..., nN_y, n1_z, ..., nN_z], [..]
              belmat(in + id*neln , jn + jd*neln) += sum;
            }
        }
    }
    
    // end of integration loop
    delete qd;
  }
  else
  {
    it = spring_map.find(index);
    
    if(it != spring_map.end())
    {
      belmat.zeros();
      press = it->second;
      int stype = spring_type_map[index]; // 0 = normal, 1 = isotropic

      // create quadrature rule
      Quadrature * qd;
      qd = Quadrature::create(fe->get_order_u(), fe->get_type());

      // quadrature loop
      for(int q=0; q<qd->get_num_ipoints(); q++)
      {
        fe->calc_shape_u(qd->get_point(q),shape);
        fe->calc_deriv_shape_u(qd->get_point(q),dshape);
        sm.calc_jacobian(dshape, neln);
        detJxW = qd->get_weight(q);

        // Compute dx0/dxi, dx0/deta using REFERENCE coordinates
        dxis.zeros();
        if(ndim==2) dxis(2,1) = -1.0;
        for(int id=0; id<ndim; id++)
          for(int jd=0; jd<ndim-1; jd++)
            for(int in=0; in<neln; in++)
            {
              int ip = bdof[in]/ndim;
              dxis(id,jd) = dxis(id,jd) + x0[ip][id] * dshape(in,jd);
            }

        arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));
        double norm_xnorm = arma::norm(xnorm, 2);

        if(stype == 0)
        {
          // Normal spring (restoring): K = +K_epi * (n̂₀ ⊗ n̂₀) * Ni * Nj * dΓ₀
          //   Adds stiffness in the normal direction
          arma::mat33 TensorProd = xnorm * xnorm.t();
          for(int in=0; in<neln; in++)
            for(int jn=0; jn<neln; jn++)
              for(int id=0; id<ndim; id++)
                for(int jd=0; jd<ndim; jd++)
                {
                  double sum = press * shape(in) * shape(jn) * TensorProd(id,jd)
                               / norm_xnorm * detJxW;
                  belmat(in + id*neln, jn + jd*neln) += sum;
                }
        }
        else
        {
          // Isotropic spring (restoring): K = +k * delta_{id,jd} * Ni * Nj * dΓ₀
          //   Adds stiffness in all directions
          double dGamma0 = norm_xnorm * detJxW;
          for(int in=0; in<neln; in++)
            for(int jn=0; jn<neln; jn++)
            {
              double mass = press * shape(in) * shape(jn) * dGamma0;
              for(int id=0; id<ndim; id++)
                belmat(in + id*neln, jn + id*neln) += mass;
            }
        }
      }
      // end of integration loop
      delete qd;
    }
    else
    {
      error("Boundary index for pressure load not found");
    }
  }
}

void NonlinearElasticity::calc_pforce_kpress(const int elem_id, const MxFE * fe,
                                             const std::vector<int> & bdof,
                                             arma::vec & belvec,
                                             arma::mat & belmat)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int neln = ndof/ndim;
  int index;
  double detJxW, press;
  arma::vec shape;
  arma::mat dshape, dxis(3,2);

  // third order alternating tensor
  arma::cube eps = arma::zeros(3,3,3);
  eps(0,1,2) = 1; eps(1,2,0) = 1; eps(2,0,1) = 1;
  eps(0,2,1) = -1; eps(1,0,2) = -1; eps(2,1,0) = -1;

  // get current coordinates of the element
  std::vector<arma::vec3> xe(neln);
  get_boundary_element_x (elem_id, xe);

  // boundary element mapping
  SurfaceMapping sm = fe->get_boundary_mapping(elem_id, xe);

  std::map<int,double>::iterator it;
  index = sm.get_index();
  it = pressure_map.find(index);

  if(pressure_map.find(index) != pressure_map.end())
  {
    belvec.zeros();
    belmat.zeros();
    press = it->second;

    // create quadrature rule (need to check if order is ok)
    Quadrature * qd;
    qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    //qd = Quadrature::create(2, fe->get_type());

    // quadrature loop
    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe->calc_shape_u(qd->get_point(q),shape);
      fe->calc_deriv_shape_u(qd->get_point(q),dshape);
      sm.calc_jacobian(dshape, neln);

      // Ajuste !
      detJxW = qd->get_weight(q);

      //
      // evaluates dx/dxi, dx/deta
      //
      dxis.zeros();

      // for 2D set dx/deta = (0,0,-1)
      if(ndim==2) dxis(2,1) = -1.0;

      for(int id=0; id<ndim; id++)
	    for(int jd=0; jd<ndim-1; jd++)
	      for(int in=0; in<neln; in++)
	      {
	        //int ip = bdof[in];
		int ip = bdof[in]/ndim;
	        dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
	      }

      //
      // computes pressure force (nodal forces)
      // subroutine pforce
      //
      arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));
      for(int in=0; in<neln; in++)
      {
	    double vl = press * shape(in) * detJxW;
	    for(int id=0; id<ndim; id++)
          belvec(in + id*neln) += xnorm(id) * vl;
      }

      //
      // computes the pressure component of the stiffness matrix
      // subroutine kpress
      //
      double sum, val, val2, val3;
      double apress = press * lc.load();

      val = -(apress * detJxW)/2.0;

      for(int in=0; in<neln; in++)
	    for(int jn=0; jn<neln; jn++)
	    {
	      val2 = val * (shape(in)*dshape(jn,0) - shape(jn)*dshape(in,0));
	      val3 = 0.0;

	      if(ndim==3)
	        val3 = val * (shape(in)*dshape(jn,1) - shape(jn)*dshape(in,1));

	      for(int id=0; id<ndim; id++)
	        for(int jd=0; jd<ndim; jd++)
	        {
	          sum = 0.0;
	          for(int kd=0; kd<3; kd++)
		        sum = sum + eps(id,jd,kd) * (val2 * dxis(kd,1) - val3 * dxis(kd,0));

              // write in the elemental pressure stiffness matrix
	          belmat(in + id*neln, jn + jd*neln) += sum;
	        }
	    }
    }
    // end of integration loop
    std::cout << "belmat: \n" << belmat <<std::endl;
    delete qd;
  }

}

double NonlinearElasticity::calc_cavity_volume(const int elem_id, const MxFE * fe,
					                                     const std::vector<int> & bdof,
					                                     const int cavity_marker)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int neln = ndof/ndim;
  int index;
  double detJxW;
  arma::vec shape;
  arma::mat dshape, dxis(3,2);

  // get current coordinates of the element
  std::vector<arma::vec3> xe(neln);
  get_boundary_element_x (elem_id, xe);

  // boundary element mapping
  SurfaceMapping sm = fe->get_boundary_mapping(elem_id, xe);

  std::map<int,double>::iterator it;
  index = sm.get_index();

  // Return 0 for any marker that does not belong to the requested cavity.
  // Other markers (e.g. epicardium, base, spring/Robin surfaces) contribute
  // nothing to the cavity volume.
  if(index != cavity_marker)
	  return 0;

  // Divergence theorem, V = (1/3) int x.n dA, evaluated on the current
  // configuration. This requires the marker to bound a CLOSED surface on
  // its own -- run report_boundary_closure() to verify: ||int n dA||/area
  // must be ~0. If a cavity marker is open at the base, this result becomes
  // dependent on the position of the origin and a lid formula is needed
  // instead.
  double endo_volume = 0;

    // create quadrature rule (need to check if order is ok)
    Quadrature * qd;
    qd = Quadrature::create(2, fe->get_type());

    // quadrature loop
    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe->calc_shape_u(qd->get_point(q),shape);
      fe->calc_deriv_shape_u(qd->get_point(q),dshape);
      sm.calc_jacobian(dshape, neln);
      detJxW = qd->get_weight(q);

      //
      // evaluates dx/dxi, dx/deta
      //
      dxis.zeros();

      // for 2D set dx/deta = (0,0,-1)
      if(ndim==2) dxis(2,1) = -1.0;

      for(int id=0; id<ndim; id++)
        for(int jd=0; jd<ndim-1; jd++)
          for(int in=0; in<neln; in++)
	        {
		        int ip = bdof[in]/ndim;
	          dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
	        }

      // area-weighted normal in the current configuration: n dA
      arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));

      // position at the quadrature point
      arma::vec3 xq;
      xq.zeros();
      for(int in=0; in<neln; in++)
        xq += shape(in) * xe[in];

      endo_volume += arma::dot(xq, xnorm) * detJxW / 3.0;
    }
    // end of integration loop
    delete qd;

  return endo_volume;
}

void NonlinearElasticity::report_active_regions()
{
  const int nelem = msh.get_n_elements();
  if (material == NULL || nelem == 0)
    return;

  // Count elements per region marker.
  std::map<int, int> count;
  int max_mesh_marker = -1;
  for (int i = 0; i < nelem; i++)
  {
    const int mk = msh.get_element_index(i);
    count[mk]++;
    max_mesh_marker = std::max(max_mesh_marker, mk);
  }

  // Every element marker must have an entry in the region map, otherwise the
  // material classes read map_mat[marker] past the end of the vector.
  if (!region_material_map.empty() &&
      max_mesh_marker >= (int)region_material_map.size())
  {
    std::ostringstream oss;
    oss << "Mesh has element region marker " << max_mesh_marker
        << " but <regions> only declares markers up to "
        << (int)region_material_map.size() - 1
        << ". Add the missing <marker> entries.";
    throw runtime_error(oss.str());
  }

  cout << "\n Active tension per element region:" << endl;
  cout << "   source: "
       << (material->has_active_scale()
             ? "ta_scale declared per material"
             : "legacy rule (only marker 0 contracts)")
       << endl;

  int n_active = 0;
  for (std::map<int, int>::const_iterator it = count.begin();
       it != count.end(); ++it)
  {
    const double s = material->active_scale(it->first);
    cout << "   marker " << it->first << ": " << it->second
         << " elements, Ta x " << s;
    if (s == 0.0)
      cout << "   (passive)";
    cout << endl;
    if (s != 0.0)
      n_active += it->second;
  }

  if (n_active == 0)
    cout << " *** WARNING: no element receives active tension. Check the "
         << "region markers of the mesh against <regions>/<material>."
         << endl;
  else
    cout << "   " << n_active << " of " << nelem
         << " elements are contractile." << endl;
  cout << endl;
}

void NonlinearElasticity::report_boundary_closure()
{
  // For a CLOSED surface, the integral of the outward normal vanishes:
  //     || int n dA || / area  ~  0
  // For an OPEN surface it equals the area vector of the opening, so the
  // ratio is O(1). This tells whether a given marker, ON ITS OWN, bounds a
  // closed region -- which is what the cavity volume computation requires.
  MixedFiniteElement * bfe = fespace.create_boundary_FE();
  if (bfe == NULL)
    return;

  int ndim = bfe->get_ndim();
  int neln = bfe->get_ndofs_u()/ndim;
  int nb = msh.get_n_boundary_elements();

  std::map<int,arma::vec3> nsum;      // int n dA
  std::map<int,double>     area;      // int dA
  std::map<int,double>     vol_div;   // (1/3) int x.n dA
  std::map<int,double>     vol_lid;   // int ((x-x0).e)(e.n) dA
  std::map<int,int>        count;

  arma::vec shape;
  arma::mat dshape, dxis(3,2);
  vector<int> bdof;

  for (int i = 0; i < nb; i++)
  {
    std::vector<arma::vec3> xe(neln);
    get_boundary_element_x(i, xe);
    SurfaceMapping sm = bfe->get_boundary_mapping(i, xe);
    int mk = sm.get_index();

    if (nsum.find(mk) == nsum.end())
    {
      nsum[mk].zeros();
      area[mk] = 0.0; vol_div[mk] = 0.0; vol_lid[mk] = 0.0; count[mk] = 0;
    }
    count[mk]++;

    fespace.get_boundary_element_dofs_u(i, bdof);

    Quadrature * qd = Quadrature::create(2, bfe->get_type());
    for (int q = 0; q < qd->get_num_ipoints(); q++)
    {
      bfe->calc_shape_u(qd->get_point(q), shape);
      bfe->calc_deriv_shape_u(qd->get_point(q), dshape);
      sm.calc_jacobian(dshape, neln);
      double w = qd->get_weight(q);

      dxis.zeros();
      if (ndim == 2) dxis(2,1) = -1.0;
      for (int id = 0; id < ndim; id++)
        for (int jd = 0; jd < ndim-1; jd++)
          for (int in = 0; in < neln; in++)
          {
            int ip = bdof[in]/ndim;
            dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
          }

      arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));
      arma::vec3 xq;  xq.zeros();
      for (int in = 0; in < neln; in++)
        xq += shape(in) * xe[in];

      nsum[mk]    += xnorm * w;
      area[mk]    += arma::norm(xnorm,2) * w;
      vol_div[mk] += arma::dot(xq, xnorm) * w / 3.0;
      vol_lid[mk] += (arma::dot(xq, lid_normal) - lid_offset)
                     * arma::dot(lid_normal, xnorm) * w;
    }
    delete qd;
  }
  delete bfe;

  cout << "\n--- Boundary closure diagnostic ---\n";
  cout << " marker  nfaces        area   ||int n dA||/area"
       << "     V_(1/3)x.n [mL]     V_lid [mL]\n";
  for (auto & kv : area)
  {
    int mk = kv.first;
    double ratio = (kv.second > 0.0)
                   ? arma::norm(nsum[mk],2)/kv.second : 0.0;
    cout << setw(7) << mk << setw(8) << count[mk]
         << setw(12) << kv.second
         << setw(20) << ratio
         << setw(20) << vol_div[mk]*1e6
         << setw(15) << vol_lid[mk]*1e6 << "\n";
  }
  cout << " ratio ~ 0 => marker is closed on its own (both volumes agree)\n";
  cout << " ratio ~ O(1) => marker is open; only V_lid is meaningful,\n"
       << "                 and only if the opening is planar.\n";
  cout << "-----------------------------------\n\n";
}

void NonlinearElasticity::set_cavity_lid_plane(const arma::vec3 & normal,
                                               double offset)
{
  lid_normal = arma::normalise(normal);
  lid_offset = offset;
  lid_plane_set = true;

  cout << " Cavity lid plane set: n = (" << lid_normal(0) << ", "
       << lid_normal(1) << ", " << lid_normal(2) << "), offset = "
       << lid_offset << endl;
}

void NonlinearElasticity::detect_cavity_lid_plane(int base_marker)
{
  // Estimate the valve/base plane from the boundary elements carrying
  // base_marker: the plane normal is the (normalised) sum of area-weighted
  // normals, and the offset is the area-weighted mean of e.x over it.
  MixedFiniteElement * bfe = fespace.create_boundary_FE();
  if (bfe == NULL)
    return;

  int ndim = bfe->get_ndim();
  int neln = bfe->get_ndofs_u()/ndim;
  int nb = msh.get_n_boundary_elements();

  arma::vec3 nsum;  nsum.zeros();
  arma::vec3 csum;  csum.zeros();
  double asum = 0.0;
  int nfound = 0;

  arma::vec shape;
  arma::mat dshape, dxis(3,2);
  vector<int> bdof;

  for (int i = 0; i < nb; i++)
  {
    std::vector<arma::vec3> xe(neln);
    get_boundary_element_x(i, xe);
    SurfaceMapping sm = bfe->get_boundary_mapping(i, xe);
    if (sm.get_index() != base_marker)
      continue;

    nfound++;
    fespace.get_boundary_element_dofs_u(i, bdof);

    Quadrature * qd = Quadrature::create(2, bfe->get_type());
    for (int q = 0; q < qd->get_num_ipoints(); q++)
    {
      bfe->calc_shape_u(qd->get_point(q), shape);
      bfe->calc_deriv_shape_u(qd->get_point(q), dshape);
      sm.calc_jacobian(dshape, neln);
      double w = qd->get_weight(q);

      dxis.zeros();
      if (ndim == 2) dxis(2,1) = -1.0;
      for (int id = 0; id < ndim; id++)
        for (int jd = 0; jd < ndim-1; jd++)
          for (int in = 0; in < neln; in++)
          {
            int ip = bdof[in]/ndim;
            dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
          }

      arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));
      arma::vec3 xq;  xq.zeros();
      for (int in = 0; in < neln; in++)
        xq += shape(in) * xe[in];

      double da = arma::norm(xnorm, 2) * w;
      nsum += xnorm * w;
      csum += xq * da;
      asum += da;
    }
    delete qd;
  }
  delete bfe;

  if (nfound == 0 || asum <= 0.0 || arma::norm(nsum,2) <= 0.0)
  {
    cout << " Note: no lid plane detected from marker " << base_marker
         << "; using n = (0,0,1), offset = 0." << endl;
    cout << "       This is harmless if each cavity marker already bounds a"
         << " closed surface\n"
         << "       (the formula is exact for any e and x0 in that case)."
         << " Call report_boundary_closure()\n"
         << "       to verify." << endl;
    lid_plane_set = true;   // avoid retrying every call
    return;
  }

  arma::vec3 e  = arma::normalise(nsum);
  arma::vec3 x0 = csum / asum;
  set_cavity_lid_plane(e, arma::dot(e, x0));
}

void NonlinearElasticity::config(const string & mshfile, const string & parfile)
{
  cout << "Setup of non-linear elasticity problem" << endl;
  int num_increments;
  std::string mtype, etype;
  std::string extension = file_extension(parfile);
  std::vector<double> matprop;

  int nelem = msh.get_n_elements();

  if(extension == "par")
  {
    cout << "Reading parameters from .par file" << endl;

    InputFile ifile(parfile);
    ifile.read("num_increments", num_increments);
    ifile.read("problem_type", etype);
    ifile.read("material_type", mtype);
    ifile.read_array("material_coef", matprop);
    ifile.read_array("body_force", bforce);
    ifile.read_section("neumann", neumann_map);
    ifile.read_section("dirichlet", dirichlet_map);
    ifile.read_section("loads", nodal_loads_map);
    ifile.read_section("pressure", pressure_map);
    ifile.read_section("prescribed_displacement", fixed_nodes_map);
    ifile.close();

    // setup load control
    lc.set_nincs(num_increments);

    // setup elasticity problem type
    ElasticityType elastype;
    if (std::strncmp(etype.c_str(),"PLANE_STRAIN",12)==0)
      elastype = PLANE_STRAIN;
    else if(std::strncmp(etype.c_str(),"PLANE_STRESS",12)==0)
      elastype = PLANE_STRESS;
    else if(std::strncmp(etype.c_str(),"THREE_DIM",9)==0)
      elastype = THREE_DIM;
    else
      throw runtime_error("Unknown elasticity type.");

    // setup material type
    
    material = HyperelasticMaterial::create(mtype, elastype, matprop);
    cout << "Hyperelastic material: " << mtype << endl;

    // setup mesh filename
    assert(file_exists(mshfile.c_str()));
    std::string::size_type idx = mshfile.rfind(".msh");
    filename = mshfile;
    basename = filename.substr(0,idx);
  }

  if(extension == "xml")
  {
    cout << "Reading parameters from XML file" << endl;

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(parfile.c_str());

    if (!result) {
      std::cout << "XML parsed with errors\n";
      std::cout << "Error description: " << result.description() << "\n";
      std::cout << "Error offset: " << result.offset;
      std::cout << " (error at [..." << (result + result.offset) << "]\n\n";
    }

    // setup elasticity problem type
    ElasticityType elastype;
    pugi::xml_attribute type = doc.child("elasticity").attribute("type");
    if (std::strncmp(type.as_string(), "PLANE_STRAIN", 12) == 0)
      elastype = PLANE_STRAIN;
    else if (std::strncmp(type.as_string(), "PLANE_STRESS", 12) == 0)
      elastype = PLANE_STRESS;
    else if (std::strncmp(type.as_string(), "THREE_DIM", 9) == 0)
      elastype = THREE_DIM;

    // number of increments

    pugi::xml_node params = doc.child("elasticity").child("parameters");
    pugi::xml_node reg = params.child("regions");

    if (params.attribute("num_materials"))
    {
      int num_materials = params.attribute("num_materials").as_int();
      if (num_materials <= 0)
        throw runtime_error("num_materials must be positive.");

      // ---- regions: marker id -> material id ---------------------------
      // The map is indexed by the marker id itself. The previous version
      // pushed the material ids in document order and silently relied on the
      // markers being listed as 0,1,2,...; a gap or a reordering shifted
      // every region onto the wrong material without any warning.
      int max_marker_id = -1;
      for (pugi::xml_node m = reg.child("marker"); m;
           m = m.next_sibling("marker"))
      {
        int id = m.attribute("id").as_int(-1);
        if (id < 0)
          throw runtime_error("<marker> without a valid non-negative id.");
        max_marker_id = std::max(max_marker_id, id);
      }

      std::vector<int> map_matAHA;
      std::vector<bool> marker_given;
      if (max_marker_id >= 0)
      {
        map_matAHA.assign(max_marker_id + 1, 0);
        marker_given.assign(max_marker_id + 1, false);
      }

      for (pugi::xml_node m = reg.child("marker"); m;
           m = m.next_sibling("marker"))
      {
        int id = m.attribute("id").as_int();
        int mt = m.attribute("material").as_int();

        if (mt < 0 || mt >= num_materials)
          throw runtime_error("<marker> refers to a material id outside "
                              "[0, num_materials).");

        if (marker_given[id])
          cout << " *** WARNING: marker id " << id << " listed more than "
               << "once in <regions>; the last entry wins." << endl;

        map_matAHA[id] = mt;
        marker_given[id] = true;
      }

      for (std::size_t i = 0; i < marker_given.size(); i++)
        if (!marker_given[i])
          cout << " *** WARNING: marker id " << i << " is not listed in "
               << "<regions>; it falls back to material 0." << endl;

        // load control
        int num_increments = params.child("ninc").text().as_int();
        lc.set_nincs(num_increments);

        // ---- materials ------------------------------------------------
        // Indexed by the material id, for the same reason as the markers
        // above. ta_scale is the multiplier applied to the active tension of
        // the regions made of this material: 1.0 fully contractile, 0.0
        // purely passive (valve plugs, scar), anything in between for
        // hypocontractile tissue. Absent means 1.0.
        std::vector<std::vector<double>> vec_matprop(num_materials);
        std::vector<double> mat_ta_scale(num_materials, 1.0);
        std::vector<bool> material_given(num_materials, false);
        bool any_ta_scale = false;

        for (pugi::xml_node m = params.child("material"); m;
             m = m.next_sibling("material"))
        {
          int id = m.attribute("id").as_int(-1);
          if (id < 0 || id >= num_materials)
            throw runtime_error("<material> id outside [0, num_materials).");

          if (material_given[id])
            cout << " *** WARNING: material id " << id << " defined more "
                 << "than once; the last definition wins." << endl;

          mtype = m.attribute("type").as_string();

          std::string strprop = m.child("coefficients").text().as_string();
          matprop.clear();
          parse_to_vector(strprop, matprop);

          double ta_scale = 1.0;
          if (m.attribute("ta_scale"))
          {
            ta_scale = m.attribute("ta_scale").as_double();
            any_ta_scale = true;
          }

          std::cout << "Material " << id;
          if (m.attribute("ta_scale"))
            std::cout << " (ta_scale = " << ta_scale << ")";
          std::cout << std::endl;

          std::vector<double>::iterator it;
          for (it = matprop.begin(); it != matprop.end(); ++it)
            cout << *it << " ";
          cout << endl;

          vec_matprop[id] = matprop;
          mat_ta_scale[id] = ta_scale;
          material_given[id] = true;
        }

        for (int i = 0; i < num_materials; i++)
          if (!material_given[i])
            throw runtime_error("num_materials declares a material that has "
                                "no <material> block.");

      // setup material
      std::cout << mtype << std::endl;

      // Kept for the bounds check in report_active_regions(): the material
      // classes index map_mat[md->get_marker()] with no bounds check, so a
      // marker beyond the last <marker> id would read past the end of the
      // vector and silently pick up garbage material parameters. The mesh is
      // not loaded yet at this point, so the check cannot happen here.
      region_material_map = map_matAHA;

      material = HyperelasticMaterial::create(mtype, elastype,  vec_matprop, num_materials, map_matAHA);
      cout << " Hyperelastic material: " << mtype << endl;
      cout << " Material properties: ";
      cout << " Multiple materials were defined. ";
      cout << endl;

      // ---- active tension per region ----------------------------------
      // Only take over the active tension when the input actually asked for
      // it. Without any ta_scale attribute the material keeps the legacy
      // rule (region 0 contracts, everything else is passive), so existing
      // input files reproduce exactly.
      if (any_ta_scale)
        material->set_active_scale(map_matAHA, mat_ta_scale);
    }
    else
    {
      if (params)
      {
        // load control
        int num_increments = params.child("ninc").text().as_int();
        lc.set_nincs(num_increments);

        // material type
        mtype = params.child("material").text().as_string();

        // material properties
        std::string strprop = params.child("coefficients").text().as_string();
        parse_to_vector(strprop, matprop);

        // setup material
        material = HyperelasticMaterial::create(mtype, elastype, matprop);
        cout << " Hyperelastic material: " << mtype << endl;
        cout << " Material properties: ";
        std::vector<double>::iterator it;
        for (it = matprop.begin(); it != matprop.end(); ++it)
          cout << *it << " ";
        cout << endl;

        int nelem = msh.get_n_elements();
      }
    }

    // body force (ex: gravity)
    pugi::xml_node bf = doc.child("elasticity").child("body_forces");
    if (bf) {
      std::string strbf = bf.text().as_string();
      parse_to_vector(strbf, bforce);
    }

    //
    // reading boundary conditions
    //

    // neumann
    pugi::xml_node nodes = doc.child("elasticity").child("neumann");
    for (pugi::xml_node node = nodes.child("node"); node;
         node = node.next_sibling("node")) {
      //int idx = node.attribute("id").as_int();
      int marker = node.attribute("marker").as_int();
      arma::vec3 t;
      t(0) = node.attribute("t0").as_double();
      t(1) = node.attribute("t1").as_double();
      t(2) = node.attribute("t2").as_double();
      neumann_map.insert(std::pair<int, arma::vec3>(marker, t));
    }

    // dirichlet
    pugi::xml_node dnodes = doc.child("elasticity").child("dirichlet");
    for (pugi::xml_node node = dnodes.child("node"); node;
         node = node.next_sibling("node")) {
      int marker = node.attribute("marker").as_int();
      int dir = node.attribute("direction").as_int();
      dirichlet_map.insert(std::pair<int, int>(marker, dir));
    }

    // dirichlet - fixed nodes
    pugi::xml_node fnodes = doc.child("elasticity").child("prescribed_displacement");
    for (pugi::xml_node node = fnodes.child("node"); node;
         node = node.next_sibling("node"))
    {
      int idx = node.attribute("id").as_int();
      int dir = node.attribute("direction").as_int();
      double val = node.attribute("value").as_double();
      NodalData nd(dir, val);
      fixed_nodes_map.insert(std::pair<int, NodalData>(idx, nd));
    }

    // pressure (normal following load) boundary condition
    pugi::xml_node pnodes = doc.child("elasticity").child("pressure");
    for (pugi::xml_node node = pnodes.child("node"); node;
         node = node.next_sibling("node")) {
      int marker = node.attribute("marker").as_int();
      double val = node.attribute("value").as_double();
      pressure_map.insert(std::pair<int, double>(marker, val));
    }

    // spring (Robin) boundary condition
    // type="normal"    : FSn0 = K_epi (u . n0) n0  (penalizes normal displacement only)
    // type="isotropic" : PN + ku = 0               (penalizes all displacement components)
    // default is "normal" for backward compatibility
    pugi::xml_node snodes = doc.child("elasticity").child("spring");
    for (pugi::xml_node node = snodes.child("node"); node;
         node = node.next_sibling("node")) {
      int marker = node.attribute("marker").as_int();
      double val = node.attribute("value").as_double();
      spring_map.insert(std::pair<int, double>(marker, val));
      int stype = 0; // default: normal spring
      std::string type_str = node.attribute("type").as_string("normal");
      if (type_str == "isotropic") stype = 1;
      spring_type_map.insert(std::pair<int, int>(marker, stype));
    }

    // nodal loads
    pugi::xml_node lnodes = doc.child("elasticity").child("loads");
    for (pugi::xml_node node = lnodes.child("node"); node;
         node = node.next_sibling("node")) {
      int idx = node.attribute("id").as_int();
      int dir = node.attribute("direction").as_int();
      double val = node.attribute("value").as_double();
      NodalData nd(dir, val);
      nodal_loads_map.insert(std::pair<int, NodalData>(idx, nd));
    }

    // setup mesh filename
    assert(file_exists(mshfile.c_str()));
    std::string::size_type idx = mshfile.rfind(".xml");
    filename = mshfile;
    basename = filename.substr(0,idx);
  }

  // print some info
  cout << " Number of nodal loads: " << nodal_loads_map.size() << endl;
  cout << " Number of prescribed displ.: " << fixed_nodes_map.size() << endl;
  cout << " Number of traction (Neumann) loads: " << neumann_map.size() << endl;
  cout << " Number of dirichlet boundary conds:" << dirichlet_map.size() << endl;
  cout << " Number of normal pressure loads: " << pressure_map.size() << endl;
  cout << " Number of spring boundary conds: " << spring_map.size();
  if (spring_map.size() > 0)
  {
    cout << " (";
    for (auto &st : spring_type_map)
      cout << "marker " << st.first << ": " << (st.second == 0 ? "normal" : "isotropic") << " ";
    cout << ")";
  }
  cout << endl;

  int npoints= msh.get_n_points();
  material->allocate_Ta(npoints);

}

void NonlinearElasticity::elem_resid (const int iel, const MxFE * fe,
                                      const Quadrature * qd, arma::vec & Re)
{
  error("NonlinearElasticity is an abstract class.");
}

void NonlinearElasticity::elem_stiff (const int iel, const MxFE * fe,
                                      const Quadrature * qd, arma::mat & Ke)
{
  error("NonlinearElasticity is an abstract class.");
}

void NonlinearElasticity::evaluate_forces(petsc::Vector & R)
{
  int ndofs = fespace.get_ndofs();
  double dlamb = lc.load_step();

  //tload = fext;

  for(int i=0; i<ndofs; i++)
    if(ldgof[i])
      R.add(i, dlamb * fext[i]);
}

int NonlinearElasticity::get_num_nz_prescribed()
{
  int cont=0;
  double val;
  NodalData inf;
  std::map<int,NodalData>::const_iterator it;

  for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
  {
    inf = it->second;
    val = inf.second;
    if (val != 0.0) cont++;
  }

  return cont;
}

int NonlinearElasticity::get_num_integration_points()
{
  int nint = 0;

  if (msh.get_nen() == 3 && msh.get_n_dim() == 2)
    nint = 3;
  else if (msh.get_nen() == 4 && msh.get_n_dim() == 2)
    nint = 4;
  else if (msh.get_nen() == 4 && msh.get_n_dim() == 3)
    nint = 1;
  else if (msh.get_nen() == 8 && msh.get_n_dim() == 3)
    nint = 8;

  return nint;
}

void NonlinearElasticity::get_element_x(const int e,
                                        std::vector<arma::vec3> & xe)
{
  uint i, k;
  std::vector<int> ptnums;
  msh.get_element_pt_nums(e, ptnums);

  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    xe[i] = x[k];
  }
}

void NonlinearElasticity::get_element_x0(const int e,
                                         std::vector<arma::vec3> & xe)
{
  uint i, k;
  std::vector<int> ptnums;
  msh.get_element_pt_nums(e, ptnums);

  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    xe[i] = x0[k];
  }
}

void NonlinearElasticity::get_boundary_element_x(const int e,
                                                 std::vector<arma::vec3> & xe)
{
  uint i, k;
  std::vector<int> ptnums;
  msh.get_boundary_element_pt_nums(e, ptnums);
  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    xe[i] = x[k];
  }
}

void NonlinearElasticity::get_boundary_element_x0(const int e, std::vector<arma::vec3> & xe)
{
  uint i, k;
  std::vector<int> ptnums;
  msh.get_boundary_element_pt_nums(e, ptnums);
  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    xe[i] = x0[k];
  }
}

void NonlinearElasticity::get_displacements(arma::vec & u)
{
  for(uint i=0; i<msh.get_n_points(); i++)
  {
    double dx = x[i][0] - x0[i][0];
    double dy = x[i][1] - x0[i][1];
    double dz = x[i][2] - x0[i][2];
    u(0 + i*3) = dx;
    u(1 + i*3) = dy;
    u(2 + i*3) = dz;
  }
}

void NonlinearElasticity::get_displacements(arma::mat & umat)
{
  const int np = msh.get_n_points();
  const int nd = msh.get_n_dim();
  umat.reshape(np,nd);
  for(int d=0; d<nd; d++)
    for(int i=0; i<np; i++)
      umat(i,d) = x[i][d];
}

void NonlinearElasticity::init()
{
    init_mesh();
    init_matvecs();
}

void NonlinearElasticity::init_mesh()
{
  // initialize mesh
  std::string extension = file_extension(filename);

  if(extension == "msh")
  {
    cout << "Reading GMSH mesh" << endl;
    GmshIO gmshreader(msh);
    gmshreader.read(filename);
  }
  else if(extension == "xml")
  {
    cout << "Reading XML mesh" << endl;
    msh.read_xml(filename);
  }

  cout << msh;
  fespace.set_mesh(&msh);

  // setup size of Voigt vec/mat
  nvoig = (msh.get_n_dim() == 2 ? 3 : 6);

  // initialize coordinates
  const std::vector<arma::vec3> & pts = msh.get_points();
  for(uint i=0; i<pts.size(); i++)
  {
    const arma::vec3 & p = pts[i];
    x.push_back(p);
    x0.push_back(p);
  }

  // std::cout << "Output step: " << output_step << std::endl;
  // // setup data writer
  // if(output_step)
  // {
  //   std::string output = filename.substr(0, filename.length() - 4) + "_output_nl";
  //   writer.open(output, lc.get_nincs() + 1 + INCs_TA, 1);
  // }
}

void NonlinearElasticity::setup_data_writer(int size)
{
  std::cout << "Output step size: " << size << std::endl;
  std::string output = filename.substr(0, filename.length() - 4) + "_output_nl";
  writer.open(output, size + 1, 1);
}

void NonlinearElasticity::init_matvecs()
{
  int nod, dir, idx;
  int N = fespace.get_ndofs();
  int n_elem = msh.get_n_elements();

  NodalData inf;
  std::map<int,NodalData>::const_iterator it;

  num_dofs = N;

  fext.resize(N);  fext.zeros();
  fext0.resize(N); fext0.zeros();
  udisp.resize(N); udisp.zeros();
  U.resize(N); U.zeros();

  tload.resize(N); tload.zeros();

  // initialize ldgof array
  //  --> ldgof[i] = true, if dof i is active
  //  --> ldgof[i] = false, if dof i is prescribed
  for(int i=0; i<N; i++) ldgof.push_back(true);

  for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
  {
    nod = it->first;
    inf = it->second;
    dir = inf.first;
    //idx = nod + (dir * msh.get_n_points());
    idx = (nod*msh.get_n_dim()) + dir;
    ldgof[idx] = false;
  }

  // initialize nodal loads vector
  for(it=nodal_loads_map.begin(); it!=nodal_loads_map.end(); ++it)
  {
    nod = it->first;
    inf = it->second;
    dir = inf.first;
    //idx = nod + (dir * msh.get_n_points());
    idx = (nod*msh.get_n_dim()) + dir;
    fext(idx) = inf.second;
  }

  // Allocate matrices and vectors
  cout << "N: " << N << endl;
  K.create(N,N,120); K = 0.0;
  u.create(N);       u = 0.0;
  r.create(N);       r = 0.0;
  react.create(N);   react = 0.0;

  // initialize stress array
  int nint = get_num_integration_points();
  stressdb.set_size(msh.get_n_elements(), nint, 9);
  straindb.set_size(msh.get_n_elements(), nint, 9);

  // initialize vector of F's
  for(int i=0; i<msh.get_n_elements(); i++)
    for(int j=0; j<nint; j++)
    {
      arma::mat33 * F = new arma::mat33();
      F->eye(3,3);
      vecF.push_back(F);
    }

  // resize vectors for 3 field formulation and ALG
  ep.resize(n_elem);
  eJ.resize(n_elem);
  eL.resize(n_elem);
  eL0.resize(n_elem);
  eps.resize(n_elem);
  eps0.resize(n_elem);
  vol0.resize(n_elem);

  if( material->is_incompressible() )
  {
    double penalt;
    penalt = (static_cast<IncompressibleMaterial*>(material))->get_kappa();

    eL0.fill(1.0);
    eL.zeros();
    eps.fill(penalt);
    eps0.fill(penalt);
  }

  // Reported here and not in config(): the electromechanics path calls
  // config() before init(), so the mesh is still empty at parsing time.
  // init_matvecs() runs after the mesh is loaded in both entry points.
  report_active_regions();
}

void NonlinearElasticity::init_resid_stiff()
{
  K = 0.0;
  react = 0.0;
  double xlamb = lc.load();
  for(int i=0; i<fespace.get_ndofs(); i++)
  {
    r.set(i, fext0(i) + xlamb * fext(i) );
  }
}

void NonlinearElasticity::line_search(double & eta0, double & eta,
                                      double & rtu0, double & rtu)
{
  double rtu1, alfa, q;

  rtu = calc_energy();
  if (rtu == 0.0) return;

  rtu1 = (rtu - rtu0 * (1 - eta))/(eta*eta);
  alfa = rtu0/rtu1;
  eta0 = eta;

  // Take different actions, depending on the value of alfa
  if (alfa < 0)
  {
    q = (alfa - sqrt(alfa * (alfa - 4)))/2.0;
    eta = alfa / q;
  }
  else if (alfa < 2)
  {
    eta = alfa/2.0;
  }
  else
  {
    eta = 1.0;
    rtu = 0.0;
  }
}

void NonlinearElasticity::output_vtk(const int cont, const int step)
{
  int np = msh.get_n_points();
  //int ne = msh.get_n_elements();
  std::stringstream ss;
  std::string name;
  std::string vtuname;

  //arma::vec t11(ne), t22(ne), t33(ne);
  //arma::vec t12(ne), t13(ne), t23(ne);
  //t11.zeros(); t22.zeros(); t33.zeros();
  //t12.zeros(); t13.zeros(); t23.zeros();

  ss << cont << "_" << step;
  name = this->basename + "_" + ss.str();
  vtuname = name + ".vtu";

  // Displacements
  //cout << "Aloca vetor" << endl;
  double *u_field = new double[np*3];
  //cout << "Preenche vetor" << endl;
  for(int i=0; i<np; i++)
  {
    double dx = x[i][0] - x0[i][0];
    double dy = x[i][1] - x0[i][1];
    double dz = x[i][2] - x0[i][2];
    u_field[i*3 + 0] = dx;
    u_field[i*3 + 1] = dy;
    u_field[i*3 + 2] = dz;
  }


  //if(output_step)
  //{
  //cout << "Salva hdf5" << endl;
    writer.write_displ_step(step, u_field);
  //}

  //cout << "Libera vetor" << endl;
  delete [] u_field;

  /*
  // Average stress output
  int nint = get_num_integration_points();
  for(int i=0; i<msh.get_n_elements(); i++)
  {
    for(int j=0; j<nint; j++)
    {
      t11(i) += stressdb(i,j,0);
      t12(i) += stressdb(i,j,1);
      t13(i) += stressdb(i,j,2);
      t22(i) += stressdb(i,j,4);
      t23(i) += stressdb(i,j,5);
      t33(i) += stressdb(i,j,8);
    }
    t11(i) /= nint;
    t12(i) /= nint;
    t13(i) /= nint;
    t22(i) /= nint;
    t23(i) /= nint;
    t33(i) /= nint;
  }

  vtkout.write_point_data(ux, "ux");
  vtkout.write_point_data(uy, "uy");
  vtkout.write_point_data(uz, "uz");
  vtkout.write_point_data(um, "umagnitude");
  vtkout.write_cell_data(ep,  "pressure");
  vtkout.write_cell_data(t11, "cauchy11");
  vtkout.write_cell_data(t22, "cauchy22");
  vtkout.write_cell_data(t33, "cauchy33");
  vtkout.write_cell_data(t12, "cauchy12");
  vtkout.write_cell_data(t13, "cauchy13");
  vtkout.write_cell_data(t23, "cauchy23");
  vtkout.write_def_mesh(x, vtuname);
  */
}

void NonlinearElasticity::output_vtk(const int step, const arma::vec & v, const arma::vec & displ)
{
  // LUCAS: talvez a fonte de ter dois arquivos de saida
  cout << "Writing Data\n";
  writer.write_vm_step(step, v.memptr());
  writer.write_displ_step(step, displ.memptr());
}

void NonlinearElasticity::output_vtk(const int cont,
                                     const int step,
                                     const std::string & name,
                                     const arma::vec & v,
                                     const ArrayMat33 & matfib)
{
  int np = msh.get_n_points();

  // First fiber directions
  arma::vec3 f, s, n;
  std::vector<arma::vec3> vec_f, vec_s, vec_n;

  for(uint i=0; i<matfib.size(); i++)
  {
    f = (matfib[i])->col(0);
    s = (matfib[i])->col(1);
    n = (matfib[i])->col(2);
    vec_f.push_back(f);
    vec_s.push_back(s);
    vec_n.push_back(n);
  }

  double *u_field = new double[np*3];
  for(int i=0; i<np; i++)
  {
    double dx = x[i][0] - x0[i][0];
    double dy = x[i][1] - x0[i][1];
    double dz = x[i][2] - x0[i][2];
    u_field[i*3 + 0] = dx;
    u_field[i*3 + 1] = dy;
    u_field[i*3 + 2] = dz;
  }

  if(output_step)
  {
    cout << "Writing Data\n";
    writer.write_vm_step(step, v.memptr());
    writer.write_displ_step(step, u_field);
  }

  delete [] u_field;
}

void NonlinearElasticity::storeLVvolumes(string basename)
{
  ofstream file;
  string aux = basename + string("_volumes.dat");
  file.open(aux.c_str());

  file << "# V_LV V_RV V_myocardium\n";
  file << volume_LV() << "\n";
  file << volume_RV() << "\n";
  file << calc_volume();

  file.close();
}

void NonlinearElasticity::storeStress(int step)
{
  //ofstream sfile, strainfile;
  //string aux = basename + string("_stress.dat");
  //sfile.open(aux.c_str());
  //aux = basename + string("_strain.dat");
  //strainfile.open(aux.c_str());
  arma::mat33 sig, E;
  int nint = get_num_integration_points();
  std::string vtuname;
  int ne = arma::size(stressdb, 0);

  arma::vec t11(ne), t22(ne), t33(ne);
  arma::vec t12(ne), t13(ne), t23(ne);
  arma::vec fiber_stress(ne), fiber_strain(ne), long_strain(ne), circ_strain(ne), rad_strain(ne);
  t11.zeros(); t22.zeros(); t33.zeros();
  t12.zeros(); t13.zeros(); t23.zeros();
  fiber_stress.zeros(); fiber_strain.zeros();
  long_strain.zeros(); circ_strain.zeros(); rad_strain.zeros();

  //cout << "Pontos de integração: " << nint << endl;

  vtuname = basename + "_stress.vtu";

  for (int id_el = 0; id_el < ne; id_el++)
  {
    for (int kk = 0; kk < 9; kk++) // stress component
    {
      double mediastress = 0;
      double mediaE = 0;
      for (int ii = 0; ii < nint; ii++) // integration point
      {
        mediastress += stressdb(id_el, ii, kk);
        mediaE += straindb(id_el, ii, kk);
      }
      mediastress = mediastress / (double) nint;
      //sfile.setf(ios::scientific);
      //sfile << mediastress << "\t";
      sig[kk] = mediastress;

      mediaE = mediaE / (double) nint;
      //strainfile.setf(ios::scientific);
      //strainfile << mediaE << "\t";
      E[kk] = mediaE;
    }

    t11(id_el) = sig(0,0);
    t12(id_el) = sig(0,1);
    t13(id_el) = sig(0,2);
    t22(id_el) = sig(1,1);
    t23(id_el) = sig(1,2);
    t33(id_el) = sig(2,2);

    arma::mat33 F, *auxF;
    F.zeros();
    for (int ii = 0; ii < nint; ii++) // integration point
      {
        auxF = vecF[id_el*nint + ii];
        F = F + *auxF;
      }
    F = F/nint;

    arma::vec3 fib0 = msh.get_element(id_el).get_fiber();
    arma::vec3 fib = msh.get_element(id_el).get_fiber();
    arma::vec3 f_l = msh.get_element(id_el).get_long();
    arma::vec3 f_c = msh.get_element(id_el).get_circ();
    arma::vec3 f_r = msh.get_element(id_el).get_rad();


    //std::cout << f_l << " " << f_c << " " << f_r << std::endl;

    fib = F*fib; //F*fib; //Alterei aqui para testar stress longitudinal
    fib = fib/arma::norm(fib,2);

    fiber_stress(id_el) = fib(0) * (fib(0)*sig(0,0) + fib(1)*sig(0,1) + fib(2)*sig(0,2)) +
                          fib(1) * (fib(0)*sig(0,1) + fib(1)*sig(1,1) + fib(2)*sig(1,2)) +
                          fib(2) * (fib(0)*sig(0,2) + fib(1)*sig(1,2) + fib(2)*sig(2,2));
    //sfile << fiber_stress(id_el) << "\t";
    //sfile << endl;

    fiber_strain(id_el) = fib0(0) * (fib0(0)*E(0,0) + fib0(1)*E(0,1) + fib0(2)*E(0,2)) +
                          fib0(1) * (fib0(0)*E(0,1) + fib0(1)*E(1,1) + fib0(2)*E(1,2)) +
                          fib0(2) * (fib0(0)*E(0,2) + fib0(1)*E(1,2) + fib0(2)*E(2,2));

    long_strain(id_el) =  f_l(0) * (f_l(0)*E(0,0) + f_l(1)*E(0,1) + f_l(2)*E(0,2)) +
                          f_l(1) * (f_l(0)*E(0,1) + f_l(1)*E(1,1) + f_l(2)*E(1,2)) +
                          f_l(2) * (f_l(0)*E(0,2) + f_l(1)*E(1,2) + f_l(2)*E(2,2));

    circ_strain(id_el) =  f_c(0) * (f_c(0)*E(0,0) + f_c(1)*E(0,1) + f_c(2)*E(0,2)) +
                          f_c(1) * (f_c(0)*E(0,1) + f_c(1)*E(1,1) + f_c(2)*E(1,2)) +
                          f_c(2) * (f_c(0)*E(0,2) + f_c(1)*E(1,2) + f_c(2)*E(2,2));

    rad_strain(id_el) =   f_r(0) * (f_r(0)*E(0,0) + f_r(1)*E(0,1) + f_r(2)*E(0,2)) +
                          f_r(1) * (f_r(0)*E(0,1) + f_r(1)*E(1,1) + f_r(2)*E(1,2)) +
                          f_r(2) * (f_r(0)*E(0,2) + f_r(1)*E(1,2) + f_r(2)*E(2,2));

    //std::cout << fiber_strain(id_el) << " " << long_strain(id_el) << " " << circ_strain(id_el) << " " << rad_strain(id_el) << std::endl;

    //strainfile << fiber_strain(id_el) << "\t";
    //strainfile << endl;
  }
  //sfile.close();
  //strainfile.close();

  writer.write_cell_field_step(step, fiber_stress.memptr(), string("stress"));
  writer.write_cell_field_step(step, fiber_strain.memptr(), string("strain"));
  writer.write_cell_field_step(step, long_strain.memptr(), string("long_strain"));
  writer.write_cell_field_step(step, circ_strain.memptr(), string("circ_strain"));
  writer.write_cell_field_step(step, rad_strain.memptr(), string("rad_strain"));
}

void NonlinearElasticity::store_point_field(int step, const arma::vec & v,
                                            const std::string & name)
{
  // The writer copies exactly n_points values out of the buffer, so a longer
  // vector is harmless but a shorter one would read past the end.
  const arma::uword np = msh.get_n_points();
  if(v.n_elem < np)
  {
    cout << " Warning: nodal field '" << name << "' has " << v.n_elem
         << " entries but the mesh has " << np
         << " nodes; not written." << endl;
    return;
  }

  writer.write_point_field_step(step, v.memptr(), name);
}


void NonlinearElasticity::prescribe_displacements()
{
  int nod, dir, ndim = msh.get_n_dim();
  double val;
  NodalData inf;
  std::map<int,double> boundary_values;
  std::map<int,NodalData>::iterator it;

  // loop over nodes with prescribed displacement
  for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
  {
    nod = it->first;
    inf = it->second;
    dir = inf.first;
    val = inf.second;
    assert(dir+1 <= ndim);
    x[nod][dir] = x0[nod][dir] + lc.load() * val;
  }
}

void NonlinearElasticity::reset()
{
  K.assemble();

  K = 0.0;
  u = 0.0;
  udisp.zeros();
  lc.reset();
}

void NonlinearElasticity::set_pressure_Ta(int mlv, double plv, int mrv, double prv, arma::vec ta)
{
  pressure_map[mlv] = plv;
  pressure_map[mrv] = prv;
  material->set_Ta(ta);
}

void NonlinearElasticity::set_active_stress(arma::vec ta)
{
  material->set_Ta(ta); 
}

void NonlinearElasticity::run(const string & mshfile, const string & parfile)
{
  // TODO: esta fixo para formato XML, mas pensando em remover codigo do .par

  // setup mesh filename
  assert(file_exists(mshfile.c_str()));
  std::string::size_type idx = mshfile.rfind(".xml");
  filename = mshfile;
  basename = filename.substr(0,idx);

  // init(); // -> init_mesh(); init_matvecs();

  init_mesh(); 

  msg("Setting material and elasticity type");
  msg("Setting boundary conditions");
  config(mshfile, parfile);
  set_output_step(true);
  
  std::cout << "Output step: " << output_step << std::endl;
  // setup data writer
  if(output_step)
  {
    std::string output = filename.substr(0, filename.length() - 4) + "_output_nl";
    writer.open(output, lc.get_nincs() + 1 + INCs_TA, 1);
  }

  init_matvecs();

  pre_solve();
  cout << "Initial LV cavity volume: " << volume_LV() << "\n";
  cout << "Initial RV cavity volume: " << volume_RV() << "\n";

  // standalone run: the load ramp is monotonic, so per-increment PV data
  // is physically meaningful
  set_pv_record_increments(true);

  solve();
  
  storeStress(lc.get_nincs());
  
  storeLVvolumes(this->basename);
  
  close_pv_history();
  
  //cout << "Final cavity volume: " << total_volume_cavity() << "\n";
  
  timer.summary();
}

void NonlinearElasticity::update_geometry(double eta)
{
  const int np = msh.get_n_points();
  const int nd = msh.get_n_dim();
  arma::mat umat;

  // get data from PETSc vector and copy to udisp
  u.get_data( udisp.memptr() );
  umat = arma::reshape(udisp, np, nd);
  U = U + udisp;

  // update geometry
  for(int d=0; d<nd; d++)
    for(int i=0; i<np; i++)
    {
      //int idx = i + (d*np);
      int idx = (i*nd) + d;
      // update if dof is free (not fixed)
      if( ldgof[idx] )
        x[i][d] = x[i][d] + eta * umat(i,d);
    }
}

void NonlinearElasticity::update_vectors(const ArrayMat33 & matfib0,
                                         ArrayMat33 & matfib)
{
  int ndim, ndofs, nnode, nelem;
  MixedFiniteElement * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(0, fe->get_type());

  ndim  = msh.get_n_dim();
  nelem = msh.get_n_elements();
  ndofs = fe->get_ndofs_u();
  nnode = fe->get_nnode();

  for(int iel=0; iel<nelem; iel++)
  {
    arma::mat dshape;
    arma::mat gradn(ndofs,ndim);
    arma::mat jacinv(ndim,ndim);
    arma::vec3 centroid;
    std::vector<arma::vec3> xe(ndofs);
    std::vector<arma::vec3> x0(ndofs);
    get_element_x(iel, xe);
    get_element_x0(iel, x0);

    Mapping em = fe->get_mapping(iel, x0);

    // compute centroid
    //   works for quad
    //   works for hex
    centroid.zeros();

    // TODO: check for each reference element
    // TODO: need to map c to the reference element ?

    // compute shape and shape derivatives at centroid
    fe->calc_deriv_shape_u(centroid, dshape);
    em.calc_jacobian(dshape, nnode);
    jacinv = em.get_inv_jacobian();
    gradn = (dshape * jacinv);

    // compute F
    arma::mat33 F(arma::fill::zeros);

    F(2,2) = (ndim==2) ? 1.0 : 0.0;

    // F = dx/dX
    for(int id=0; id<ndim; id++)
      for(int jd=0; jd<ndim; jd++)
	      for(int k=0; k<nnode; k++)
	        F(id,jd) += xe[k][id] * gradn(k,jd);

    // update vector directions
    arma::vec3 f,s,n;

    f = F * (matfib0[iel])->col(0);
    s = F * (matfib0[iel])->col(1);
    n = F * (matfib0[iel])->col(2);

    // normalize vectors
    f = f / arma::norm(f,2);
    s = s / arma::norm(s,2);
    n = n / arma::norm(n,2);

    (matfib[iel])->col(0) = f;
    (matfib[iel])->col(1) = s;
    (matfib[iel])->col(2) = n;
  }

  delete qd;
  delete fe;
}

// ---------- NONLINEAR PROBLEM ------------------------------------------------
bool NonlinearElasticity::converged(petsc::Vector & du)
{
  double rnorm=0, enorm=0, dnorm=0;
  double etol   = parameters["tol_energy"];
  double rtol   = parameters["tol_residual"];
  double dtol   = parameters["tol_displacement"];

  // tload - external load including pressure
  // fext - external load

  if(first_step)
  {
    double f0 = react.l2norm() + arma::norm(lc.load()*tload,2);
    u.copy_values(du.size(), du);
    parameters["energy_norm0"]   = fabs(calc_energy());
    parameters["residual_norm0"] = f0;

    rnorm = r.l2norm();
    rnorm = rnorm/f0;
    enorm = fabs(calc_energy());
    dnorm = du.l2norm();

    first_step = false;
  }
  else
  {
    rnorm = r.l2norm() / ( react.l2norm() + arma::norm(tload,2) ); // R/F (F=force+react)
    enorm = fabs(calc_energy());
    dnorm = du.l2norm();

    cout << scientific << setprecision(3);
    cout << " rnorm " << rnorm;
    cout << " enorm " << enorm;
    cout << " dunorm " << dnorm;
    cout << endl;
 }

  // convergence criterion
  bool check;
  check  = enorm < etol;
  check |= rnorm < rtol;
  check |= dnorm < dtol;
  return check;
}

void NonlinearElasticity::evaluate(petsc::Vector & resid)
{
  timer.enter("Residual");

  int n_dofs = msh.get_n_dim() * msh.get_nen();
  int n_elem = msh.get_n_elements();

  // init resid
  double xlamb = lc.load();
  for(int i=0; i<fespace.get_ndofs(); i++)
    r.set(i, fext0(i) + xlamb * fext(i));

  // set tload (external loads incl. pressure) to current fext (external loads)
  tload = fext;

  //
  // assemble residual vector r
  //
  arma::vec Re(n_dofs);
  std::vector<int> dnums;
  MxFE * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(0, fe->get_type());

  for(int i=0; i<n_elem; i++)
  {
    elem_resid (i, fe, qd, Re);
    fespace.get_element_dofs_u (i, dnums);

    for(int k=0; k<n_dofs; k++)
    {
      if (ldgof[dnums[k]])
      {
        r.add(dnums[k], -Re(k)); // add -R
      }
      else
        react.add(dnums[k], Re(k));

    }
  }
  delete qd;
  delete fe;

  //
  // pressure forces contribution
  //
  MixedFiniteElement * bfe = fespace.create_boundary_FE();
  if (bfe != NULL && (pressure_map.size() > 0 || spring_map.size() > 0))
  {
    int nu  = bfe->get_ndofs_u();
    int nb = msh.get_n_boundary_elements();
    double xlamb = lc.load();
    arma::vec belvec(nu);
    vector<int> bdof;

    for(int i=0; i<nb; i++)
    {
      bool is_spring = false;
      fespace.get_boundary_element_dofs_u(i,bdof);
      elem_pforce(i,bfe,bdof,belvec,is_spring);

      // Pressure is an EXTERNAL load and is ramped by the load factor.
      // The spring is a displacement-dependent RESTORING force: it must not
      // be ramped, otherwise the residual and the (unramped) spring
      // stiffness are inconsistent and Newton loses quadratic convergence.
      double scale = is_spring ? 1.0 : xlamb;

      // assembles the nodal forces due to normal pressure or spring
      for(int k=0; k<nu; k++)
      {
        if(ldgof[bdof[k]])
        {
          r.add(bdof[k], scale * belvec(k));
          // only external loads contribute to the reference load norm
          if(!is_spring)
            tload(bdof[k]) += belvec(k);
        }
        else
          react.add(bdof[k], -scale * belvec(k));
      }
    }
  }

  r.assemble();

  // copy from r to resid
  resid.copy_values(r.size(), r);
  resid.assemble();

  timer.leave();
}

void NonlinearElasticity::jacobian(petsc::Matrix & Kstiff)
{

  timer.enter("Stiffness");

  int n_dofs = msh.get_n_dim() * msh.get_nen();
  int n_elem = msh.get_n_elements();
  arma::mat Ke(n_dofs,n_dofs);
  std::vector<int> dnums;

  Kstiff = 0.0;

  MxFE * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(0, fe->get_type());

  for(int i=0; i<n_elem; i++)
  {

    elem_stiff (i, fe, qd, Ke);
    fespace.get_element_dofs_u (i, dnums);
#ifndef USE_BFGS
    // Fast assembling
    int * pidx;
    pidx = &dnums[0];
    Ke = Ke.t();
    Kstiff.add(n_dofs, n_dofs, pidx, pidx, Ke.memptr());
#endif

#ifdef USE_BFGS
    for(int j=0; j<n_dofs; j++)
    {
      for(int k=0; k<n_dofs; k++)
      {
        int I = dnums[j];
        int J = dnums[k];
        if(J >= I) Kstiff.add(I, J, Ke(j,k));
      }
    }
#endif
  }

  //
  // *** PRESSURE COMPONENT OF THE STIFFNESS MATRIX ***
  //
  MxFE * bfe = fespace.create_boundary_FE();
  if (bfe != NULL && (pressure_map.size() > 0 || spring_map.size()>0))
  {
    //cout << " Assembling pressure component of stiffness matrix" << endl;
    int nu = bfe->get_ndofs_u();
    int nb = msh.get_n_boundary_elements();
    arma::mat belmat(nu, nu);
    vector<int> bdof;

    for (int i = 0; i < nb; i++)
    {
      fespace.get_boundary_element_dofs_u(i, bdof);
      elem_kpress(i, bfe, bdof, belmat);
      // assembles the pressure component of the stiffness matrix
      for (int j = 0; j < nu; j++)
        for (int k = 0; k < nu; k++)
          Kstiff.add(bdof[j], bdof[k], belmat(j, k));
    }
  }

  delete qd;
  delete fe;

  Kstiff.assemble();
  apply_boundary(Kstiff);

  timer.leave();
}

void NonlinearElasticity::update(petsc::Vector & uu, double s)
{
  const int n_dim   = msh.get_n_dim();
  const int n_nodes = msh.get_n_points();
  arma::mat umat;

    arma::vec udispB;

    // get data from PETSc vector and copy to udisp
    uu.get_data( udisp.memptr() );
    udispB.resize(fespace.get_ndofs());
    udispB.zeros();
    for(int i=0; i<n_nodes; i++)
    {
        udispB(i) = udisp(n_dim*i);
        udispB(i + n_nodes) = udisp(n_dim*i + 1);
        udispB(i + 2*n_nodes) = udisp(n_dim*i + 2);
    }
    umat = arma::reshape(udispB, n_nodes, n_dim);

    // update total displacement vector
    U = U + udisp;

  // copy uu to internal u
  u.copy_values(uu.size(), uu);

  // update geometry
  for(int d=0; d<n_dim; d++)
    for(int i=0; i<n_nodes; i++)
    {
      //int idx = i + (d * n_nodes);
      int idx = (i*n_dim) + d;

      // if degree of freedom is free (not fixed/prescribed)
      // then update x = x + u
      if( ldgof[idx] )
        x[i][d] = x[i][d] + s * umat(i,d);
    }
}