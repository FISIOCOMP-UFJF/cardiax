#ifndef NONLINEAR_ELASTICITY_HPP
#define NONLINEAR_ELASTICITY_HPP

#include <map>
#include <iostream>
#include <sstream>
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

//#define USE_BFGS

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
  void storeLVvolumes(string basename);


  //! Reset stiffness matrix and displacement vector
  void reset();

  //! Run simulation
  void run();

  void set_pressure_Ta(int mlv, double plv, int mrv, double prv, const arma::vec &ta, const arma::vec &dta);

  //! Compute total cavity volume
  double total_volume_cavity();

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
  //! Compute cavity volume in a element
  double calc_cavity_volume(const int elem_id, const MxFE * fe,
                                                   const std::vector<int> & bdof);

  //! Compute stiffness matrix at the element
  virtual void elem_resid (const int iel, const MxFE * fe,
                           const Quadrature * qd, arma::vec & Re);

  //! Compute residual at the element
  virtual void elem_stiff (const int iel, const MxFE * fe,
                           const Quadrature * qd, arma::mat & Ke);

  void elem_pforce(const int elem_id, const MxFE * fe,
                   const std::vector<int> & bdof,
                   arma::vec & belvec);

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
