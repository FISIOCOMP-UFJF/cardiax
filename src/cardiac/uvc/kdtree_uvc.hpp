#ifndef KDTREE_UVC_HPP
#define KDTREE_UVC_HPP

#include <string>
#include <vector>

#include "reader_hdf5.hpp"

// =============================================================================
//  Interface de transferencia de campos entre malhas cardiacas por UVC
// -----------------------------------------------------------------------------
//  Mesmo nucleo do transfere_campo_uvc.cpp, so que partido em duas funcoes:
//
//    kdtree_uvc_build()   CARA, chamada UMA vez. Monta o embedding cilindrico
//                         E = [ab, tm, ab*cos(rt), ab*sin(rt)], as KDTrees por
//                         ventriculo e, para cada no do ALVO, a lista de
//                         vizinhos da FONTE com os pesos ja normalizados.
//
//    kdtree_uvc_transf()  BARATA, chamada em LOOP. So percorre as listas
//                         prontas: media convexa (continuo) ou voto ponderado
//                         (categorico). Nenhuma busca, nenhuma alocacao de
//                         arvore.
//
//  A separacao se sustenta porque os vizinhos e os pesos dependem SO da
//  geometria e das UVC -- nunca dos valores do campo. Transferir M campos (ou
//  M passos de tempo) custa uma unica busca k-NN, feita no build.
//
//  Uso tipico (o que o cardiax faria dentro do laco de tempo):
//
//      MalhaTransf mf, ma;
//      ReaderHDF5  leitor_f, leitor_a;
//      carregar_malha_uvc(leitor_f, "Paciente_1", "fonte", mf);
//      carregar_malha_uvc(leitor_a, "Paciente_2", "alvo ", ma);
//
//      OpcoesUVC o;                 // k=12, peso gauss, clamp do ab ligado
//      KdtreeUVC est;
//      kdtree_uvc_build(mf, ma, o, est);          // 1 vez
//
//      std::vector<double> campo_alvo;
//      for (int passo = 0; passo < n_passos; passo++)
//      {
//        leitor_f.read_field_step("vm", passo, campo_fonte);
//        kdtree_uvc_transf(est, campo_fonte, false, campo_alvo);   // em loop
//      }
//
//  Dependencias: biblioteca padrao do C++, a API C do HDF5 (via ReaderHDF5)
//  e a kdtree-cpp (cdalitz), copiada sem modificacoes.
// =============================================================================

// =============================================================================
//  Malha + UVC
// =============================================================================

//! Geometria e UVC de uma malha, um valor por no. rt SEMPRE em radianos.
//! (Nome distinto do MalhaUVC do busca_uvc.hpp para os dois headers poderem
//! conviver no mesmo arquivo.)
struct MalhaTransf
{
  std::string rotulo;              //!< "fonte" / "alvo"

  int n_points;
  int n_elements;
  int nen;                         //!< nos por elemento

  std::vector<double> xyz;         //!< 3*n_points
  std::vector<int>    tets;        //!< nen*n_elements
  std::vector<double> ab, tm, rt, tv;

  MalhaTransf() : rotulo(), n_points(0), n_elements(0), nen(0),
                  xyz(), tets(), ab(), tm(), rt(), tv() {}
};

// =============================================================================
//  Parametros da transferencia (os mesmos do transfere_campo_uvc.cpp)
// =============================================================================

struct OpcoesUVC
{
  int         k;                   //!< vizinhos por no do alvo (12)
  std::string peso;                //!< "gauss" ou "idw"
  double      pot_idw;             //!< expoente do IDW (2.0)
  double      w_ab, w_tm, w_rt;    //!< pesos dos eixos do embedding (1.0)

  bool        tem_tv_split;        //!< false: calcula da media min/max de tv
  double      tv_split;
  bool        sem_clamp_ab;        //!< true: nao limita o ab do alvo

  bool        silencioso;          //!< true: build nao imprime nada

  OpcoesUVC()
    : k(12), peso("gauss"), pot_idw(2.0), w_ab(1.0), w_tm(1.0), w_rt(1.0),
      tem_tv_split(false), tv_split(0.0), sem_clamp_ab(false),
      silencioso(false) {}
};

// =============================================================================
//  Estrutura montada uma vez e consultada em loop
// =============================================================================

//! Tudo o que kdtree_uvc_transf() precisa: quem sao os vizinhos de cada no do
//! alvo e com que peso. Depois do build nao sobra nenhuma KDTree viva -- por
//! isso a struct e copiavel e barata de carregar.
//!
//! As listas sao guardadas em linhas de largura fixa k_max:
//!   viz  [k_max*i + j] -> no da FONTE (valido so para j < n_viz[i])
//!   pesos[k_max*i + j] -> peso ja normalizado (a linha soma 1)
struct KdtreeUVC
{
  int n_fonte;                     //!< nos da malha fonte
  int n_alvo;                      //!< nos da malha alvo
  int k_max;                       //!< largura das linhas de viz/pesos

  std::vector<int>    viz;         //!< k_max*n_alvo
  std::vector<double> pesos;       //!< k_max*n_alvo
  std::vector<int>    n_viz;       //!< n_alvo

  //! Nos do alvo que nao receberam vizinho nenhum (UVC nao-finita, ou
  //! ventriculo sem fonte). Copiam o valor de outro no do alvo, o mais
  //! proximo entre os que receberam.
  std::vector<int> copia_destino;
  std::vector<int> copia_origem;

  double tv_split;                 //!< limiar VE/VD efetivamente usado
  int    n_sem_valor;              //!< nos do alvo que ficam NaN mesmo assim

  KdtreeUVC()
    : n_fonte(0), n_alvo(0), k_max(0), viz(), pesos(), n_viz(),
      copia_destino(), copia_origem(), tv_split(0.5), n_sem_valor(0) {}
};

// =============================================================================
//  Leitura
// =============================================================================

//! Le geometria + UVC de um par XDMF/HDF5 e deixa rt em radianos.
//! O leitor fica ABERTO de proposito: o campo a transferir costuma vir do
//! mesmo arquivo e e lido depois, passo a passo.
//! Devolve false (com mensagem em cerr) se faltar algum campo UVC.
bool carregar_malha_uvc(ReaderHDF5 & leitor, const std::string & arquivo,
                        const std::string & rotulo, MalhaTransf & m,
                        const std::string & n_ab = "ab",
                        const std::string & n_tm = "tm",
                        const std::string & n_rt = "rt",
                        const std::string & n_tv = "tv");

//! Le um campo escalar como esta no arquivo (por no ou por elemento).
//! nome: so o nome ("tecido") ou o caminho ("vertex_field/tecido").
//! por_celula e nome_real saem preenchidos; nome_real e o CAMINHO do
//! dataset no arquivo ("/vertex_field/tecido"), unico mesmo se o nome
//! existir nos dois grupos.
bool ler_campo_cru_uvc(ReaderHDF5 & r, const std::string & nome,
                       std::vector<double> & val, bool & por_celula,
                       std::string & nome_real);

// =============================================================================
//  Interface principal
// =============================================================================

//! CHAMADA UMA VEZ. Monta embedding, arvores e listas de vizinhos.
//! Nao olha nenhum campo: so xyz/ab/tm/rt/tv das duas malhas.
//! Devolve false se as malhas estiverem incoerentes ou vazias.
bool kdtree_uvc_build(const MalhaTransf & fonte, const MalhaTransf & alvo,
                      const OpcoesUVC & o, KdtreeUVC & est);

//! CHAMADA EM LOOP. campo_fonte tem um valor por no da FONTE; campo_alvo sai
//! com um valor por no do ALVO (NaN onde nao houve como decidir).
//!
//!   categorico = false -> media convexa dos vizinhos (sem overshoot)
//!   categorico = true  -> voto ponderado, preserva 0/1 e rotulos
//!
//! Vizinhos com valor nao-finito sao descartados e os pesos restantes
//! renormalizados -- e por isso que o build nao precisa conhecer o campo.
bool kdtree_uvc_transf(const KdtreeUVC & est,
                       const std::vector<double> & campo_fonte,
                       bool categorico,
                       std::vector<double> & campo_alvo);

// =============================================================================
//  Conversao PointData <-> CellData (mantem o tipo do campo)
// =============================================================================

//! Valores distintos e finitos, em ordem crescente
std::vector<double> classes_de_uvc(const std::vector<double> & v);

//! Campo por CELULA -> por NO.
//! continuo   -> media das celulas incidentes
//! categorico -> voto majoritario (binario: fracao > 0.5; multiclasse: moda)
std::vector<double> celula_para_no(const MalhaTransf & m,
                                   const std::vector<double> & vals_cell,
                                   bool categorico);

//! Campo por NO -> por CELULA.
//! continuo   -> media dos nos da celula
//! categorico -> binario: marca a celula quando a fracao de nos da classe
//!               positiva atinge frac; multiclasse: moda dos nos
std::vector<double> no_para_celula(const MalhaTransf & m,
                                   const std::vector<double> & vals_no,
                                   bool categorico, double frac = 0.5);

// =============================================================================
//  Saida
// =============================================================================

struct CampoSaida
{
  std::string nome;
  std::vector<double> val;
  bool por_celula;                 //!< false = PointData, true = CellData

  CampoSaida() : nome(), val(), por_celula(false) {}
};

//! Grava a malha e os campos num .vtu ASCII, escrito a mao (sem VTK).
//! Cada campo vai para PointData ou CellData conforme por_celula.
bool salvar_vtu_uvc(const std::string & caminho, const MalhaTransf & m,
                    const std::vector<CampoSaida> & campos);

#endif /* KDTREE_UVC_HPP */
