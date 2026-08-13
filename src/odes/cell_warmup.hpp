#ifndef CELL_WARMUP_HPP
#define CELL_WARMUP_HPP

#include <string>
#include <vector>
#include <armadillo>

// Forward declaration
class CellModel;

/**
   Pre-condicionamento ("warm up") do modelo celular 0D.

   O problema: as condicoes iniciais publicadas de um modelo ionico
   (ToRORd-Land, ORd, TT2, ...) sao um estado de repouso aproximado, nao o
   ciclo limite do modelo no BCL da simulacao. Comecar dali significa que os
   primeiros batimentos do modelo acoplado estao integrando o transitorio
   lento das concentracoes (Na+, Ca2+ do SR, CaMKt), que leva centenas de
   batimentos para assentar. Todo o resultado mecanico e circulatorio desses
   primeiros ciclos e transitorio tambem.

   A solucao classica: resolver o modelo 0D isolado por N batimentos com
   estimulo periodico, e usar o estado do ciclo limite como condicao inicial
   de cada no.

   CUSTO. Uma celula do ToRORd com dt = 0.01 ms e BCL = 800 ms sao 80 mil
   passos por batimento. Fazer isso em cada no de uma malha biventricular e
   inviavel. Mas o estado do ciclo limite so depende dos PARAMETROS da
   celula, nao da sua posicao: duas celulas com o mesmo tipo (endo/mid/epi) e
   a mesma coordenada apicobasal tem exatamente o mesmo ciclo limite. Entao
   as celulas sao agrupadas por (tipo, faixa de ab), o pre-condicionamento
   roda uma vez por grupo e o resultado e distribuido para os nos.

   Por DEFAULT ha um grupo por TIPO CELULAR apenas (-warmup_abbins 1): o
   gradiente apicobasal nao entra no agrupamento e a celula representativa
   usa ab = 0.5, que e o valor neutro do gradiente. Medido no ToRORd-Land a
   200 batimentos, colapsar o gradiente assim custa ~1% nas variaveis lentas
   (cansr, CaMKt, Ta) e ~0.03% no Na+ intracelular -- contra dezenas de por
   cento das condicoes iniciais publicadas, que e o erro que o warm up
   existe para eliminar. Alem de ser 8x mais barato, o resultado nao depende
   da malha, entao o cache serve para qualquer malha com o mesmo protocolo.
   Quem quiser resolver o gradiente usa -warmup_abbins 8.

   FASE. Por padrao o estado guardado e o do REPOUSO, isto e, o instante
   imediatamente anterior ao proximo estimulo (fase 0 do batimento). Com
   -warmup_phase e possivel escolher qualquer instante do ciclo limite,
   medido em ms a partir do estimulo; valores negativos contam de tras para
   frente ("100 ms antes do repouso" = -warmup_phase -100).

   Com -warmup_lat 1 a fase deixa de ser a mesma em todos os nos: cada no
   recebe a fase (phase - lat_i), de modo que o tecido ja comeca com a
   sequencia de ativacao do eikonal impressa nele, em vez de todas as
   celulas em repouso simultaneo. Nesse caso o ultimo batimento inteiro e
   amostrado a cada -warmup_sample ms e o estado de cada no sai por
   interpolacao linear entre duas amostras.

   CACHE. O resultado e gravado em texto (-warmup_cache) junto com uma
   assinatura de todos os parametros que o determinam. Numa proxima rodada
   com os mesmos parametros o arquivo e lido e o pre-condicionamento e
   pulado. Qualquer diferenca na assinatura invalida o cache e ele e
   recalculado -- nunca se aproveita um cache de outra configuracao.

   UNIDADES. Todos os tempos desta classe estao em MILISSEGUNDOS, e sao
   convertidos internamente para a unidade nativa do modelo celular
   (CellModel::native_time_unit_ms). Isso vale inclusive para o LAT passado
   em state_for_cell().
*/
class CellWarmup
{
public:

  //! Parametros do pre-condicionamento. Todos os tempos em ms.
  struct Options
  {
    bool   enabled        = false;  //!< -warmup 1 liga (default: desligado)
    int    beats          = 100;    //!< -warmup_beats
    double bcl_ms         = 800.0;  //!< -warmup_bcl
    double dt_ms          = 0.01;   //!< -warmup_dt
    double stim_amp       = -53.0;  //!< -warmup_stimamp
    double stim_dur_ms    = 2.0;    //!< -warmup_stimdur
    double phase_ms       = 0.0;    //!< -warmup_phase (negativo = do fim)
    bool   phase_from_lat = false;  //!< -warmup_lat 1
    double sample_ms      = 1.0;    //!< -warmup_sample (so com -warmup_lat 1)
    int    ab_bins        = 1;      //!< -warmup_abbins (0 ou 1 ignora ab)
    double lambda         = 1.0;    //!< -warmup_lambda (estiramento fixo)
    double tol            = 0.0;    //!< -warmup_tol (0 = roda todos os beats)
    std::string cache;              //!< -warmup_cache ("none" desliga)
    std::string model_name;         //!< nome do modelo, so para o cache
  };

  //! Le as opcoes da linha de comando. Os defaults de BCL e dt vem da
  //! configuracao corrente do problema, para que o warm up use por padrao a
  //! mesma discretizacao e o mesmo periodo da simulacao.
  static Options options_from_command_line(double default_bcl_ms,
                                           double default_dt_ms,
                                           const std::string & model_name,
                                           const std::string & basename);

  //! model precisa ja ter passado por CellModel::setup() -- e o setup que
  //! cria o ODESolver usado aqui.
  CellWarmup(CellModel * m, const Options & o);

  //! Agrupa as celulas e roda (ou le do cache) o pre-condicionamento.
  //! types e ab podem vir vazios: sem tipos, todas as celulas caem num unico
  //! grupo com o tipo corrente do modelo; sem ab, o gradiente apicobasal nao
  //! entra no agrupamento.
  //! Retorna false se nada foi feito (e ai as condicoes iniciais originais
  //! devem ser mantidas).
  bool run(const arma::ivec & types, const arma::vec & ab, int n_cells);

  //! O pre-condicionamento terminou e state_for_cell() pode ser chamado?
  bool ready() const { return is_ready; }

  //! Escreve em y (num_state_vars valores) o estado do ciclo limite da
  //! celula `cell`. lat_ms so e usado com -warmup_lat 1.
  void state_for_cell(int cell, double lat_ms, double * y) const;

  bool   uses_lat()   const { return opt.phase_from_lat; }
  double get_bcl_ms() const { return opt.bcl_ms; }
  double get_phase_ms() const { return phase_wrapped; }
  int    num_groups() const { return (int) grp_type.size(); }
  int    group_of(int cell) const;

private:

  CellModel * model;
  Options     opt;

  int    nvars = 0;
  double unit  = 1.0;      //!< ms por unidade de tempo nativa do modelo
  double phase_wrapped = 0.0;  //!< opt.phase_ms dobrado em [0, bcl)

  // ---- grupos --------------------------------------------------------
  std::vector<int>    grp_type;   //!< tipo celular de cada grupo
  std::vector<double> grp_ab;     //!< coordenada ab representativa do grupo
  std::vector<int>    cell_group; //!< grupo de cada celula (vazio = grupo 0)

  // ---- resultado -----------------------------------------------------
  //! Sem -warmup_lat: um unico estado por grupo, ja na fase pedida.
  std::vector<arma::vec> rest_state;

  //! Com -warmup_lat: o ultimo batimento inteiro, nvars x nsnap por grupo.
  std::vector<arma::mat> traj;
  arma::vec snap_t;               //!< instantes das amostras [ms], [0, bcl]

  bool is_ready = false;

  void   build_groups(const arma::ivec & types, const arma::vec & ab,
                      int n_cells);

  //! Integra `length_ms` a partir do estado y, com um estimulo de duracao
  //! opt.stim_dur_ms no inicio do intervalo. Se rec != nullptr, guarda uma
  //! amostra a cada opt.sample_ms (e os instantes em snap_t, uma unica vez).
  void   pace(double * y, double length_ms, arma::mat * rec, bool fill_times);

  //! Diferenca relativa maxima entre dois estados, para medir a convergencia
  //! de batimento a batimento.
  static double drift(const arma::vec & a, const arma::vec & b);

  std::string fingerprint() const;
  bool   load_cache();
  void   save_cache() const;
  static const char * type_name(int t);
};

#endif
