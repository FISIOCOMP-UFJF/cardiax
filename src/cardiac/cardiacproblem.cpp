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
  delete warmup;
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

  cout << "   " << name << ": " << n << "/" << ndofs << " nodes";
  if (n > 0) cout << ", min " << v.min() << " max " << v.max();
  if (out_of_range > 0) cout << " (" << out_of_range << " id(s) out of range)";
  cout << endl;

  if (n > 0 && n != ndofs)
    cout << " *** WARNING: block <" << name << "> incomplete; nodes without a"
         << " value were left at zero." << endl;

  return n;
}

} // namespace

void CardiacProblem::read_cobiveco(const std::string & mshfile)
{
  const int ndofs = mesh->get_n_points();

  pugi::xml_document doc;
  if (!doc.load_file(mshfile.c_str()))
  {
    cout << " Cobiveco: could not open " << mshfile << endl;
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
    cout << "   tv: " << arma::accu(cbv_tv == 0) << " nodes in the LV, "
         << arma::accu(cbv_tv == 1) << " nodes in the RV" << endl;
  }

  cbv_loaded = (n_tm > 0);
  if (!cbv_loaded)
    cout << " No Cobiveco coordinates in the mesh; cell types left unchanged."
         << endl;
}

void CardiacProblem::set_cell_types_from_tm(double endo_mid, double mid_epi)
{
  if (!cbv_loaded) return;

  if (cells == nullptr)
  {
    cout << " *** set_cell_types_from_tm() called before init();"
         << " cell types left unchanged." << endl;
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

  cout << " Cell type from tm (endo < " << endo_mid
       << " <= mid <= " << mid_epi << " < epi): "
       << n_endo << " endo, " << n_mid << " mid, " << n_epi << " epi" << endl;

  // If there are almost no mid cells the mesh does not resolve the wall
  // thickness, and transmural heterogeneity collapses to an endo/epi step.
  if (20 * n_mid < n)
    cout << " *** WARNING: fewer than 5% of the cells are mid. The mesh"
         << " probably lacks the transmural resolution needed for"
         << " endo/mid/epi heterogeneity." << endl;
}

void CardiacProblem::set_apicobasal_from_ab()
{
  if (!cbv_has_ab) return;

  if (cells == nullptr)
  {
    cout << " *** set_apicobasal_from_ab() called before init();"
         << " apicobasal gradient disabled." << endl;
    return;
  }

  cells->set_apicobasal(cbv_ab);

  // The factor is the one the cell model applies; printing it here puts the
  // effective range in the log, which is what is worth checking.
  const double f_min = std::pow(0.2, 2.0 * cbv_ab.min() - 1.0);
  const double f_max = std::pow(0.2, 2.0 * cbv_ab.max() - 1.0);
  cout << " Apicobasal gradient on: ab from " << cbv_ab.min() << " to "
       << cbv_ab.max() << " -> GKs factor from " << f_min
       << " (most apical) to " << f_max << " (most basal)" << endl;

  if (cbv_ab.min() < 0.0 || cbv_ab.max() > 1.0)
    cout << " *** WARNING: ab outside [0,1]; the GKs factor extrapolates"
         << " beyond the intended 5x-0.2x range." << endl;
}


void CardiacProblem::set_fiber_stretch(const arma::vec & lam,
                                       const arma::vec & lam_rate)
{
  if (cells == nullptr)
  {
    cout << " *** set_fiber_stretch() chamado antes de init();"
         << " acoplamento mecano-eletrico desligado." << endl;
    return;
  }

  // Uma celula por no da malha (ver init()). Um vetor de outro tamanho
  // significa que a mecanica e a EP estao em malhas diferentes, e escrever
  // assim mesmo daria estiramento no no errado, silenciosamente.
  if (lam.n_elem != static_cast<arma::uword>(cells->size()))
  {
    cout << " *** set_fiber_stretch(): lambda tem " << lam.n_elem
         << " valores mas ha " << cells->size() << " celulas; ignorado."
         << endl;
    return;
  }

  cells->set_stretch(lam);

  if (lam_rate.n_elem == lam.n_elem)
    cells->set_stretch_rate(lam_rate);
  else
    cells->set_stretch_rate(arma::vec());
}


// ===========================================================================
//  Pre-condicionamento ("warm up") do modelo celular
//
//  As condicoes iniciais publicadas de um modelo ionico sao um repouso
//  aproximado, nao o ciclo limite no BCL da simulacao. Comecar dali faz com
//  que os primeiros batimentos do modelo acoplado estejam integrando o
//  transitorio lento das concentracoes -- e portanto que a mecanica e a
//  circulacao desses ciclos sejam transitorias tambem.
//
//  Aqui o modelo celular 0D e resolvido isoladamente por varios batimentos e
//  o estado resultante substitui as condicoes iniciais de cada no. O custo e
//  controlado agrupando as celulas por (tipo, faixa apicobasal): o ciclo
//  limite depende dos parametros da celula, nao da sua posicao.
// ===========================================================================

void CardiacProblem::warmup_cells(const arma::vec & lat)
{
  if (cells == nullptr || cellmodel == nullptr) return;

  // ---- primeira chamada: le as opcoes e roda o pre-condicionamento ----
  if (!warmup_tried)
  {
    warmup_tried = true;

    // Os defaults saem da configuracao corrente do problema, convertidos
    // para milissegundos: o warm up usa, por default, o mesmo passo de tempo
    // e o mesmo periodo da simulacao acoplada. -warmup_bcl e -warmup_dt
    // sobrescrevem qualquer um dos dois.
    const double ms_per_unit = cells->get_solver_time_unit_ms();

    // BCL do warm up: warmup_bcl quando definido (caminho acoplado, onde
    // totaltime cobre varios batimentos e portanto NAO serve como BCL),
    // senao totaltime (caminho legado, onde -t e um batimento).
    const double default_bcl =
      (warmup_bcl > 0.0 ? warmup_bcl : totaltime) * ms_per_unit;

    CellWarmup::Options opt =
      CellWarmup::options_from_command_line(default_bcl,
                                            timestep * ms_per_unit,
                                            cell_name,
                                            mesh_filename);

    if (!opt.enabled) return;

    warmup = new CellWarmup(cellmodel, opt);

    if (!warmup->run(cells->get_cell_types(),
                     cells->get_apicobasal(),
                     cells->size()))
    {
      // O warm up falhou (divergiu, parametros invalidos). Mantem as
      // condicoes iniciais originais em vez de aplicar um estado suspeito.
      delete warmup;
      warmup = nullptr;
    }
  }

  if (warmup == nullptr || !warmup->ready()) return;

  // ---- aplica o estado do ciclo limite em cada no --------------------
  const int    n           = cells->size();
  const double ms_per_unit = cells->get_solver_time_unit_ms();
  const bool   have_lat    = (lat.n_elem == static_cast<arma::uword>(n));

  if (warmup->uses_lat() && !have_lat)
    cout << " *** warm up: -warmup_lat 1 requested, but there is no per-node"
         << " LAT; every cell starts at the same phase." << endl;

  arma::vec y(cells->get_ode_size());

  for (int i = 0; i < n; i++)
  {
    // lat arrives in the solver time unit; CellWarmup works in ms.
    const double lat_ms = have_lat ? lat(i) * ms_per_unit : 0.0;
    warmup->state_for_cell(i, lat_ms, y.memptr());
    cells->set_system_state(i, y.memptr());
  }

  cout << " Initial conditions of " << n << " cell(s) replaced by the"
       << " limit-cycle state (phase " << warmup->get_phase_ms()
       << " ms";
  if (warmup->uses_lat() && have_lat) cout << ", shifted by each node's LAT";
  cout << ")." << endl;
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
