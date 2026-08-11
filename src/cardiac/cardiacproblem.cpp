#include "cardiacproblem.hpp"
#include "util/pugixml.hpp"
#include <vector>
#include <cmath>

CardiacProblem::CardiacProblem()
{
  // setup default parameter values
  parameters.add("surface_to_volume", 0.14);
  parameters.add("theta_method", 0.5);
  parameters.add("pcgtol", 1.0e-16);
  parameters.add("maxnz", 40);
}

CardiacProblem::~CardiacProblem()
{
  writer->close();
  delete mesh;
  delete writer;
}

void CardiacProblem::setup(std::string & b, std::string & c, std::string & m,
                           double dt, double T, double pr, double pa)
{
  cout << "Setup of cardiac problem" << endl;
  const double theta = parameters["theta_method"];

  timestep  = dt;
  totaltime = T;
  printrate = pr;
  printrate_apd = pa;
  
  mesh_filename = b;
  stimuli_filename = b;

  cell_name = c;
  odesolver = m;

  cout << "Parabolic solver: ";
  if (theta == 0.0)      cout << "Explicit Euler" << endl;
  else if (theta == 1.0) cout << "Implicit Euler" << endl;
  else if (theta == 0.5) cout << "Crank-Nicolson" << endl;
}

void CardiacProblem::setup_types(std::string & f)
{
  if( !file_exists(f) )
    cout << "Cells: all cells are of the same type" << endl;
  else
  {
    cout << "Cells: configuring cell types" << endl;
   
    int aux, size;
    int * vtypes;
    ifstream in(f.c_str());
    in >> size;
    vtypes = new int[size];

    for(int i=0; i<size; i++)
    {
      in >> aux;
      vtypes[i] = aux;
    }

    cells->set_cell_types(size, vtypes);   
    delete [] vtypes;
  }
}


// ===========================================================================
// Coordenadas biventriculares consistentes (Cobiveco)
//
// Cada campo e um bloco proprio dentro de <mesh>, no mesmo estilo do
// <eikonal>, com um valor por NO:
//
//   <tm size="2840">
//       <node id="0" value="0.412345" />
//       ...
//   </tm>
//
// Blocos separados (em vez de um unico bloco com varios atributos por no)
// porque os campos tem origens diferentes -- o LAT vem de um solver eikonal,
// as coordenadas vem do Cobiveco -- e sao regerados em momentos diferentes.
// ===========================================================================

namespace {

//! Le um bloco <name><node id=".." value=".."/></name> de dentro de <mesh>.
//! Retorna o numero de nos lidos, ou 0 se o bloco nao existe.
int read_nodal_field(const pugi::xml_document & doc, const char * name,
                     arma::vec & v, int ndofs)
{
  v.zeros(ndofs);

  pugi::xml_node blk = doc.child("mesh").child(name);
  if (!blk) return 0;

  int n = 0, out_of_range = 0;
  for (pugi::xml_node nd = blk.child("node"); nd; nd = nd.next_sibling("node"))
  {
    const int i = nd.attribute("id").as_int();
    if (i < 0 || i >= ndofs) { out_of_range++; continue; }
    v(i) = nd.attribute("value").as_double();
    n++;
  }

  cout << "   " << name << ": " << n << "/" << ndofs << " nos";
  if (n > 0) cout << ", min " << v.min() << " max " << v.max();
  if (out_of_range > 0) cout << " (" << out_of_range << " id(s) fora de faixa)";
  cout << endl;

  if (n > 0 && n != ndofs)
    cout << " *** AVISO: bloco <" << name << "> incompleto; os nos sem valor"
         << " ficaram em zero." << endl;

  return n;
}

} // namespace

void CardiacProblem::read_cobiveco(const std::string & mshfile)
{
  const int ndofs = mesh->get_n_points();

  pugi::xml_document doc;
  if (!doc.load_file(mshfile.c_str()))
  {
    cout << " Cobiveco: nao foi possivel abrir " << mshfile << endl;
    return;
  }

  cout << " -- Reading Cobiveco coordinates from mesh file --" << endl;

  arma::vec tv;
  const int n_tv = read_nodal_field(doc, "tv", tv,     ndofs);
  const int n_tm = read_nodal_field(doc, "tm", cbv_tm, ndofs);
  const int n_rt = read_nodal_field(doc, "rt", cbv_rt, ndofs);
  const int n_ab = read_nodal_field(doc, "ab", cbv_ab, ndofs);
  (void) n_rt;

  cbv_has_ab = (n_ab > 0);

  // tv e binaria; guarda-se como inteiro para deixar explicito que esse campo
  // nao deve ser interpolado.
  cbv_tv.zeros(ndofs);
  if (n_tv > 0)
  {
    for (int i = 0; i < ndofs; i++) cbv_tv(i) = (tv(i) > 0.5) ? 1 : 0;
    cout << "   tv: " << arma::accu(cbv_tv == 0) << " nos no VE, "
         << arma::accu(cbv_tv == 1) << " nos no VD" << endl;
  }

  cbv_loaded = (n_tm > 0);
  if (!cbv_loaded)
    cout << " Sem coordenadas Cobiveco na malha; tipo celular inalterado."
         << endl;
}

void CardiacProblem::set_cell_types_from_tm(double endo_mid, double mid_epi)
{
  if (!cbv_loaded) return;

  if (cells == nullptr)
  {
    cout << " *** set_cell_types_from_tm() chamado antes de init();"
         << " tipo celular inalterado." << endl;
    return;
  }

  const int n = static_cast<int>(cbv_tm.n_elem);
  std::vector<int> vtypes(n);

  int n_endo = 0, n_mid = 0, n_epi = 0;
  for (int i = 0; i < n; i++)
  {
    const double tm = cbv_tm(i);
    if      (tm <  endo_mid) { vtypes[i] = ENDO;  n_endo++; }
    else if (tm <= mid_epi)  { vtypes[i] = MCELL; n_mid++;  }
    else                     { vtypes[i] = EPI;   n_epi++;  }
  }

  cells->set_cell_types(n, vtypes.data());

  cout << " Tipo celular a partir de tm (endo < " << endo_mid
       << " <= mid <= " << mid_epi << " < epi): "
       << n_endo << " endo, " << n_mid << " mid, " << n_epi << " epi" << endl;

  // Se quase nao ha celulas mid, a malha nao resolve a espessura da parede e
  // a heterogeneidade transmural fica reduzida a um degrau endo/epi.
  if (20 * n_mid < n)
    cout << " *** AVISO: menos de 5% das celulas sao mid. A malha"
         << " provavelmente nao tem resolucao transmural suficiente para"
         << " heterogeneidade endo/mid/epi." << endl;
}

void CardiacProblem::set_apicobasal_from_ab()
{
  if (!cbv_has_ab) return;

  if (cells == nullptr)
  {
    cout << " *** set_apicobasal_from_ab() chamado antes de init();"
         << " gradiente apicobasal desligado." << endl;
    return;
  }

  cells->set_apicobasal(cbv_ab);

  // O fator e o mesmo que o modelo celular aplica; imprimi-lo aqui deixa
  // visivel no log a faixa efetiva, que e o que importa conferir.
  const double f_min = std::pow(0.2, 2.0 * cbv_ab.min() - 1.0);
  const double f_max = std::pow(0.2, 2.0 * cbv_ab.max() - 1.0);
  cout << " Gradiente apicobasal ligado: ab de " << cbv_ab.min() << " a "
       << cbv_ab.max() << " -> fator de GKs de " << f_min
       << " (mais apical) a " << f_max << " (mais basal)" << endl;

  if (cbv_ab.min() < 0.0 || cbv_ab.max() > 1.0)
    cout << " *** AVISO: ab fora de [0,1]; o fator de GKs extrapola a faixa"
         << " 5x-0.2x prevista." << endl;
}


void CardiacProblem::write_data(const arma::vec & u, const std::string & s, int * step)
{
  if (tip.time2print())
  {
    std::string aux = s.substr(s.length() - 2);

    if (aux == "ve")
      writer->write_ve_step(*step, u.memptr());
    else if (aux == "vm")
      writer->write_vm_step(*step, u.memptr());
    else if(aux=="dv")
    {
      arma::vec current = u( arma::span(u.size()/2, u.size()-1) );
      writer->write_ve_step(*step, current.memptr());
    }
    
    *step = *step + 1;
  }
}

void CardiacProblem::write_data(const arma::vec & u, const arma::vec & displ,
                                const std::string & s, int step)
{
    if (tip.time2print())
    {
        std::string aux = s.substr(s.length()-2);

        // potentials scalar field
        if(aux=="ve")
          writer->write_ve_step(step, u.memptr());
        else if(aux=="vm")
            writer->write_vm_step(step, u.memptr());

        // displacement vector field
        writer->write_displ_step(step, displ.memptr());

  }
}

void CardiacProblem::write_data_text(const arma::vec & vm, int * step)
{
  //if (tip.time2print())
  //{
    char buffer[256];
    sprintf(buffer,"output/vm_text/vm_%05d.txt", *step);
    string name(buffer);   
    vm.save(name, arma::raw_ascii);
    *step = *step + 1;
  //}
}
