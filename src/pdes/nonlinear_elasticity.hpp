#ifndef NONLINEAR_ELASTICITY_HPP
#define NONLINEAR_ELASTICITY_HPP

#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "fem/fem.h"
#include "util/util.hpp"
#include "util/load_control.hpp"
#include "util/parameters.hpp"
#include "linalg/linalg.h"
#include "nls/nonlinear_problem.hpp"
#include "nls/nonlinear_solver.hpp"
#include "nls/bfgs.hpp"
#include "nls/newton.hpp"
#include "nls/newton_ls.hpp"
#include "materials/material_data.hpp"
#include "materials/incompressible_material.hpp"
#include "mech_utils.hpp"

// #define USE_BFGS

//! Enum with bitflags for assemble functions
enum AssembleFlags
{
  ASSEMBLE_RESID = 0,
  ASSEMBLE_STIFF = 1,
  ASSEMBLE_ALL   = 2
};

//! Some typedefs to make code cleaner
typedef std::vector<arma::mat33*> ArrayMat33;
typedef std::map<int,double> BCPressure;
typedef std::map<int,arma::vec3> BCTraction;
typedef std::multimap<int,NodalData> BCNode;

//! Boundary markers identifying the ventricular cavities.
//! Any other marker (epicardium, base, spring/Robin surfaces, ...)
//! contributes zero to the cavity volume computation.
const int MARKER_LV = 30;   //!< Left ventricle endocardium
const int MARKER_RV = 20;   //!< Right ventricle endocardium
const int MARKER_BASE = 10; //!< Base / valve plane (lid closing the cavities)

/*!
 *  Abstract base class for nonlinear elasticity problems    
 *  implemented using either the Updated Lagrangian formulation 
 *  or the Total Lagrangian formulation.
 */
class NonlinearElasticity : public NonlinearProblem
{
public:

  //! Constructor
  NonlinearElasticity();

  //! Destructor
  virtual ~NonlinearElasticity();

  //! Apply boundary conditions to system (impose prescribed displacements)
  void apply_boundary(petsc::Matrix & Ks);

  //! Computes the mesh volume
  double calc_volume(bool update=false);

  //! Configure problem with parameters from file
  void config(const string & mshfile, const string & parfile);

  //! Initialize problem (init_mesh, init_coords, etc...)
  void init();

  //! Getters
  int get_num_nz_prescribed();
  int get_num_integration_points();
  Mesh & get_mesh() { return msh; }
  const std::vector<arma::mat33*> & get_vec_F() { return vecF; };

  //! Fibre stretch lambda_f = sqrt(I4f) = sqrt(f0.C.f0) = ||F f0||, one
  //! value per ELEMENT (averaged over its integration points). Uses the
  //! deformation gradient left in vecF by the last solve, so it is only
  //! meaningful after solve() has run. Resizes lam_e to n_elements.
  void fiber_stretch_elements(arma::vec & lam_e);

  //! Same quantity projected onto the NODES by averaging over the elements
  //! sharing each node, WEIGHTED by the reference element volume vol0 (the
  //! lumped L2 projection). The cell models live on the nodes (one ODE
  //! system per mesh point), so this is the form the electrophysiology needs.
  //! Falls back to unweighted averaging if vol0 has not been filled yet
  //! (calc_volume(true), called from pre_solve()). Resizes lam_n to n_points.
  void fiber_stretch_nodes(arma::vec & lam_n);
  petsc::Matrix & get_K() { return K; }
  void get_displacements(arma::mat & umat);
  void get_displacements(arma::vec & u);

  //! Output VTK data of that step
  void output_vtk(const int cont, const int step);

  //! New output vtk
  void output_vtk(const int step, const arma::vec & v, const arma::vec & displ);

  void setup_data_writer(int size);
  //! Output VTK data of that step and write scalar field v
  void output_vtk(const int cont, const int step,
                  const std::string & name, const arma::vec & v,
                  const std::vector<arma::mat33*> & vecfsn);
  void storeStress(int step);

  //! Write a node-centred scalar field (e.g. "vm", "active_stress") into the
  //! HDF5/XDMF output for the given output step.
  //! Tensao ativa EFETIVA por elemento: a media nodal de Ta sobre os nos do
  //! elemento, vezes o ta_scale do material daquele elemento.
  //!
  //! E a grandeza que a montagem realmente usa
  //! (updated_lagrangian.cpp: md->set_active_stress(Ta_ip * lc.load() *
  //! ta_scale_el)), e NAO a que ta.max() imprime -- aquela e o valor bruto do
  //! modelo celular, antes da escala por material. Sao por elemento porque
  //! ta_scale e constante por elemento: um no na fronteira entre marcadores
  //! nao tem um valor efetivo unico, e forcar um campo nodal ali inventaria
  //! um numero que a mecanica nunca viu.
  void effective_active_tension(const arma::vec & ta_nodal,
                                arma::vec & ta_elem) const;

  //! Escala de tensao ativa de cada elemento, para inspecao/saida.
  void active_scale_field(arma::vec & scale_elem) const;

  //! Grava um campo por ELEMENTO (cell-centred) no HDF5/XDMF.
  void store_cell_field(int step, const arma::vec & v,
                        const std::string & name);

  void store_point_field(int step, const arma::vec & v,
                         const std::string & name);
  void storeLVvolumes(string basename);


  //! Reset stiffness matrix and displacement vector
  void reset();

  //! Run simulation
  void run(const string & mshfile, const string & parfile);

  void set_pressure_Ta(int mlv, double plv, int mrv, double prv, arma::vec ta);

  //! Hand the mechanics the fields needed to stabilise a velocity-dependent
  //! active tension: the nodal active stiffness Ka and the nodal fibre
  //! stretch at the previous converged solve. Passing empty vectors (the
  //! default state) disables the stabilisation entirely.
  void set_active_stabilization(const arma::vec & ka,
                                const arma::vec & lam_prev);

  void set_active_stress(arma::vec ta); 

  //! Compute total cavity volume for a given boundary marker.
  //! Boundary elements with other markers contribute zero.
  double total_volume_cavity(const int cavity_marker = MARKER_LV);

  //! Convenience wrappers for left/right ventricular cavity volumes, in mL
  double volume_LV();
  double volume_RV();

  //! Myocardial (tissue) volume in mL
  double volume_myocardium_mL();

  //! Force the mesh length unit: "m", "cm" or "mm".
  void set_mesh_units(const std::string & units);

  //! Infer the mesh length unit from the size of the reference geometry.
  //! Called automatically on the first volume evaluation.
  void detect_mesh_units();

  //! Factor converting (mesh length unit)^3 to mL
  double volume_scale_to_mL();

  //! Name of the current mesh length unit
  const std::string & mesh_units() const { return unit_name; }

  //! Set the lid (valve) plane used to close the open endocardial surfaces.
  //! offset is e.x0 for any point x0 on the plane.
  void set_cavity_lid_plane(const arma::vec3 & normal, double offset);

  //! Estimate the lid plane from the boundary elements of base_marker.
  //! Called automatically on the first cavity-volume evaluation.
  void detect_cavity_lid_plane(int base_marker = MARKER_BASE);

  //! Print, per boundary marker, whether the surface is closed on its own
  //! and what each volume formula yields. Useful to validate cavity setup.
  void report_boundary_closure();

  //! Print, per element region marker, how many elements it holds and which
  //! active tension scale it receives. A region that silently gets zero
  //! active tension is very hard to spot in the results, so the mapping is
  //! reported once, from init_matvecs(), when the mesh is already loaded.
  void report_active_regions();

  //! Region marker -> material id, as declared in <regions>. Kept only to
  //! validate the mesh markers once the mesh is available.
  std::vector<int> region_material_map;

  //! Pressure prescribed for a given cavity.
  //! apply_load_factor=true scales by the current load factor (monotonic
  //! ramp, e.g. passive filling); false returns the prescribed target
  //! pressure as set by set_pressure_Ta() (externally driven stepping).
  double cavity_pressure(const int cavity_marker,
                         bool apply_load_factor = true);

  //! Open the pressure-volume history file and write its header.
  //! Safe to call repeatedly: does nothing if already open.
  void open_pv_history(const string & basename);

  //! Append one row (increment, P_LV, V_LV, P_RV, V_RV) to the PV history,
  //! taking the pressures from the current load ramp.
  void record_pv_history(int increment);

  //! Same, with pressures supplied explicitly by the caller.
  void record_pv_history(int increment, double plv, double prv);

  //! Close the pressure-volume history file
  void close_pv_history();

  //! Whether the pressure-volume history file is currently open
  bool pv_history_is_open() { return pv_file.is_open(); }

  //! Enable/disable per-increment PV recording inside solve().
  //! Off by default: it is meaningful for a monotonic load ramp (passive
  //! filling), but not when an outer driver re-ramps the load at every
  //! cardiac-cycle step, where only the converged end-of-step state matters.
  void set_pv_record_increments(bool o) { pv_record_increments = o; }

  //! Save data at each timestep
  void set_output_step(bool o) { output_step = o; }

  //! Update vector direction (ex: f = F f0)
  void update_vectors(const ArrayMat33 & matfib0, ArrayMat33 & matfib);

  //!
  //! Interface for nonlinear elasticity problems formulations (TL and UL)
  //!
  virtual void solve() = 0;
  virtual void pre_solve() = 0;
  virtual void assemble_stiff() = 0;
  virtual void assemble_const() = 0;
  virtual void assemble_active(const arma::vec & s,
                               std::vector<arma::mat33*> & vs,
                               std::vector<arma::mat33*> & vf) = 0;

  //! Assemble pressure matrix and vector
  void assemble_pressure();

  //! NonlinearProblem stuff
  void evaluate  (petsc::Vector & r);
  void jacobian  (petsc::Matrix & K);
  void update    (petsc::Vector & u, double s);
  bool converged (petsc::Vector & du);
  uint size      () { return num_dofs; }

  //! Finite element space for the solution
  //! PUBLIC for TLSNES
  VecH1FESpace fespace;

  //! Timer for sections
  TimerSection timer;

  //! Controls loading (incremental Newton)
  LoadControl lc;

  void gambiarra(bool b) {
    cout << "PRESS MAP " << pressure_map.size() << endl;
    pressure_map.clear();
    cout << "PRESS MAP " << pressure_map.size() << endl;
  }

protected:
  
  int num_dofs;                   //!< Number of dofs (per node?)
  int nvoig;                      //!< Size of vector/matrix in Voigt notation
  bool output_step;               //!< Output at each increment

  Mesh msh;                       //!< Computational mesh
  std::string filename;           //!< Input filename
  std::string basename;           //!< Substring of the filename with basename

  //! Boundary conditions
  std::map<int,arma::vec3> neumann_map;         //!< Traction boundary cond
  std::map<int,double> pressure_map;            //!< Pressure boundary cond
  std::map<int,double> spring_map;              //!< Spring boundary cond
  std::map<int,int>    spring_type_map;          //!< Spring type: 0=normal, 1=isotropic
  std::multimap<int,int> dirichlet_map;         //!< Dirichlet boundary cond
  std::multimap<int,NodalData> fixed_nodes_map; //!< Fixed nodes boundary cond
  std::multimap<int,NodalData> nodal_loads_map; //!< Applied nodal loads

  //! Vectors
  std::vector<bool> ldgof;        //!< Which dof is fixed and free
  std::vector<arma::vec3> x;      //!< Current coordinates
  std::vector<arma::vec3> x0;     //!< Reference coordinates
  std::vector<arma::vec3> disp;   //!< Displacements
  std::vector<arma::mat33*> vecF; //!< Deformation gradient tensor
  std::vector<double> bforce;     //!< Body force (gravity))

  //! Auxiliary vectors
  arma::cube stressdb;            //!< 3-D array to store stresses
  arma::cube straindb;            //!< 3-D array to store strain
  arma::vec udisp;                //!< Auxiliary displacement vector
  arma::vec fext;                 //!< Nodal loads
  arma::vec fext0;                //!< Old nodal loads
  arma::vec tload;                //!< Total external load including pressure load
  arma::umat lbnod;               //!< Boundary element connect (Bonet: lbnod)
  arma::vec U;

  //! Mixed 3 field + Augmented Lagrangian vectors
  bool use_alg;                   //!< Boolean to check if AL is on
  arma::vec ep;                   //!< Element pressure
  arma::vec eJ;                   //!< Average element jacobian
  arma::vec eL;                   //!< Current Lagrange multiplier Lambda
  arma::vec eL0;                  //!< Previous Lagrange multiplier Lambda
  arma::vec eps;                  //!< Penalty parameter
  arma::vec eps0;                 //!< Initial penalty parameter
  arma::vec vol0;                 //!< Initial element volume

  //! PETSc vectors and matrix for system Ku = r
  petsc::Matrix K;                //!< Tangent stiffness matrix
  petsc::Vector u;                //!< Displacement vector
  petsc::Vector r;                //!< Residual vector (Fint - Fext)
  petsc::Vector react;            //!< Reaction forces

  //! More
  HyperelasticMaterial * material;//!< Constitutive law
  Log log;                        //!< Logger for history of newton iterations
  std::ofstream pv_file;          //!< Pressure-volume history (per load increment)
  int pv_counter = 0;             //!< Number of rows written to pv_file
  bool pv_record_increments = false; //!< Record PV at every load increment

  arma::vec3 lid_normal = {0.0, 0.0, 1.0}; //!< Unit normal of the valve plane
  double lid_offset = 0.0;                 //!< e.x0 for x0 on the valve plane
  bool lid_plane_set = false;              //!< Whether the plane was set/detected

  double vol_to_mL = 1.0e6;       //!< (mesh unit)^3 -> mL; default assumes m
  std::string unit_name = "m";    //!< Current mesh length unit
  bool units_detected = false;    //!< Whether units were set/detected

  WriterHDF5 writer;              //! Data writer (VTK,HDF5)

  Parameters parameters;

  //double rnorm0;
  //double enorm0;
  //double rtol;
  //double etol;
  //double dtol;

  //! Is this the first step of the solution
  bool first_step;

  //!
  //! Methods
  //!

  //! Computes the traction forces and assemble into nodal loads vector
  void assemble_traction();

  //! Computes body forces (or forcing term for use with the MMS)
  void body_forces();

  //! Computes energy
  double calc_energy();

  //! Traction (Neumann) boundary condition
  void calc_neumann_elvec(const int eindex, const MxFE * fe,
                          arma::vec & elvec);

  //! Computes boundary elemental normal pressures
  void calc_pforce_kpress(const int elem_id, const MixedFiniteElement * fe,
                          const std::vector<int> & bdof, arma::vec & belvec,
                          arma::mat & belmat);
  //! Compute cavity volume contribution of a boundary element.
  //! Returns 0 if the element marker differs from cavity_marker.
  double calc_cavity_volume(const int elem_id, const MxFE * fe,
                            const std::vector<int> & bdof,
                            const int cavity_marker = MARKER_LV);

  //! Compute stiffness matrix at the element
  virtual void elem_resid (const int iel, const MxFE * fe,
                           const Quadrature * qd, arma::vec & Re);

  //! Compute residual at the element
  virtual void elem_stiff (const int iel, const MxFE * fe,
                           const Quadrature * qd, arma::mat & Ke);

  //! Boundary force vector. Sets is_spring=true when the element belongs to
  //! a spring (Robin) surface, so the caller can skip the load-factor ramp.
  void elem_pforce(const int elem_id, const MxFE * fe,
                   const std::vector<int> & bdof,
                   arma::vec & belvec,
                   bool & is_spring);

  void elem_kpress(const int elem_id, const MxFE * fe,
                   const std::vector<int> & bdof,
                   arma::mat & belmat);

  //! Evaluate nodal forces
  void evaluate_forces(petsc::Vector & R);

  //! Returns the current coordinates of element e
  void get_element_x (const int eidx, std::vector<arma::vec3> & x);

  //! Returns the reference coordinates of element e
  void get_element_x0(const int eidx, std::vector<arma::vec3> & x0);

  //! Returns the current coordinates of boundary element e
  void get_boundary_element_x (const int eidx, std::vector<arma::vec3> & x);

  //! Returns the reference coordinates of boundary element e
  void get_boundary_element_x0(const int eidx, std::vector<arma::vec3> & x0);

  //! Read mesh
  void init_mesh();

  //! Initialize vectors and matrices
  void init_matvecs();

  //! Clean stiffness and residual
  void init_resid_stiff();

  //! Perform quadratic line search
  void line_search(double & eta0, double & eta, double & rtu0, double & rtu);

  //! Prescribed nodal displacements (x = X + u_prescribed)
  void prescribe_displacements();

  //! Update the current coordinates: (x = x + u)
  void update_geometry(double eta);

};

#endif