#ifndef CARDIAC_PROBLEM_HPP
#define CARDIAC_PROBLEM_HPP

#include "stimulus.hpp"
#include "linalg/linalg.h"
#include "odes/odes.h"
#include "mesh/writer_hdf5.hpp"
#include "util/parameters.hpp"

typedef std::vector<arma::mat33*> ArrayMat33;

/**
   Abstract base class for transient cardiac problems
   such as Bidomain or Monodomain models. 
*/

class CardiacProblem
{
public:

  //! Constructor
  CardiacProblem();

  //! Virtual destructor
  virtual ~CardiacProblem();

  //! Interface - perform one time step
  virtual void advance() = 0;

  //! Interface - initialization  
  virtual void init() = 0;
  
  //! Interface - solution
  virtual void solve() = 0;

  //! Configure initial conditions of the cells
  virtual void initial_conditions() = 0;

  //! Return an array with the cell types
  const arma::ivec & get_cell_types() const { return cells->get_cell_types(); }

  // -------------------------------------------------------------------
  // Coordenadas biventriculares consistentes (Cobiveco), por no.
  // Validas apenas se has_cobiveco() retornar verdadeiro.
  //   tv  transventricular  0 = VE, 1 = VD           (BINARIA)
  //   tm  transmural        0 = endocardio, 1 = epicardio
  //   rt  rotacional        [0,1)  CICLICA: 0 e 1 sao o mesmo angulo
  //   ab  apicobasal        0 = apice, 1 = base
  // -------------------------------------------------------------------

  //! Coordenada transmural (0 endo -> 1 epi)
  const arma::vec & get_tm() const { return cbv_tm; }

  //! Coordenada rotacional. CICLICA: nunca interpole nem faca media
  //! atravessando a costura em 0/1; use sin/cos se precisar.
  const arma::vec & get_rt() const { return cbv_rt; }

  //! Coordenada apicobasal (0 apice -> 1 base)
  const arma::vec & get_ab() const { return cbv_ab; }

  //! Coordenada transventricular, binaria (0 = VE, 1 = VD)
  const arma::ivec & get_tv() const { return cbv_tv; }

  //! Havia bloco <tm> no arquivo de malha?
  bool has_cobiveco() const { return cbv_loaded; }

  //! Le os blocos <tv> <tm> <rt> <ab> do arquivo de malha. Blocos ausentes
  //! sao ignorados; a malha continua valida sem nenhum deles.
  void read_cobiveco(const std::string & mshfile);

  //! Define o tipo celular de cada no (ENDO/MCELL/EPI) a partir de tm.
  //! Precisa ser chamado ANTES de initial_conditions(), porque Cells::init()
  //! ja usa o tipo para escolher as condicoes iniciais de cada celula.
  void set_cell_types_from_tm(double endo_mid = 0.3, double mid_epi = 0.7);

  //! Entrega a coordenada ab as celulas, ativando o gradiente apicobasal
  //! (no ToRORd-Land, GKs = GKs * 0.2^(2*ab - 1)). Sem o bloco <ab> na malha
  //! nao faz nada e o resultado fica identico ao de antes.
  void set_apicobasal_from_ab();

  //! Entrega as celulas o estiramento de fibra lambda_f = sqrt(I4f) por no,
  //! calculado pela mecanica, e opcionalmente sua taxa. E o acoplamento
  //! mecano-eletrico: no ToRORd-Land lambda entra em h(lambda) e em ca50
  //! (ativacao dependente do comprimento) e lambda_rate no modelo de
  //! distorcao (ZETAS/ZETAW). lam_rate deve estar na unidade de tempo NATIVA
  //! do modelo celular (1/ms no ToRORd); passe um vetor vazio para deixar a
  //! taxa em zero.
  void set_fiber_stretch(const arma::vec & lam,
                         const arma::vec & lam_rate = arma::vec());

  //! Return reference to the Cells object
  const Cells & get_cells() { return *cells; }

  Stimuli & get_stimuli() { return stimuli; }

  //! Return reference to the Mesh object
  const Mesh & get_mesh() { return *mesh; }

  //! Return pointer to the Mesh object
  const Mesh * get_mesh_ptr() { return mesh; }

  //! Reference to time parameters object
  const TimeParameters & get_time_parameters() { return tip; }

  //! Return the fiber direction of a given element
  const arma::vec3 & get_fiber(int elindex)
  { return mesh->get_element(elindex).get_fiber(); }

  //! Return the sheet/transverse direction of a given element
  const arma::vec3 & get_trans(int elindex)
  { return mesh->get_element(elindex).get_trans(); }

  //! Return the sheet-normal direction of a given element
  const arma::vec3 & get_normal(int elindex)
  { return mesh->get_element(elindex).get_normal(); }

  //! Change PDE time step
  void set_timestep (double s) { timestep = s; }

  //! Change total time
  void set_totaltime(double t) { totaltime = t; }

  //! Change print rate
  void set_printrate(double r) { printrate = r; }

  //! Configure cardiac problem with given data
  void setup(std::string & b, std::string & c, std::string & m,
             double dt, double T, double pr, double pa);

  //! Configure the type of each cell
  void setup_types(std::string & f);

  //! Need to update stiffness matrix?
  void update_matrix(bool r) { re_assembly_mats = r; }
  
  //! Write solution to output file
  void write_data(const arma::vec & u, const std::string & s, int * step);

  //! Write solution and displacement vector to output file
  void write_data(const arma::vec & u, const arma::vec & displ, const std::string & s, int step);

  //! Write Vm to a text file for post-processing
  void write_data_text(const arma::vec & vm, int * step);
  
  //! To measure times
  TimerSection timer;
 
protected:

  bool re_assembly_mats;
  double timestep;
  double totaltime;
  double printrate;
  double printrate_apd;

  string cell_name;
  string odesolver;
  string mesh_filename;
  string stimuli_filename;

  Stimuli stimuli;
  Parameters parameters;
  TimeParameters tip;
  Cells * cells;
  CellModel * cellmodel;

  //! Coordenadas Cobiveco por no (ver getters acima)
  arma::vec  cbv_tm, cbv_rt, cbv_ab;
  arma::ivec cbv_tv;
  bool cbv_loaded = false;   //!< havia bloco <tm>?
  bool cbv_has_ab = false;   //!< havia bloco <ab>?

  Mesh * mesh;
  WriterHDF5 * writer;

};


#endif /* CARDIAC_PROBLEM_HPP_ */