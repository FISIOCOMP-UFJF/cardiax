#ifndef CELL_COUPLING_HPP
#define CELL_COUPLING_HPP

#include <string>
#include <armadillo>
#include "odes/cells.hpp"

// ======================================================================
//  Tudo o que depende de QUAL modelo celular esta em uso
// ======================================================================
//
//  O cardiac_cycle nao deve saber que a tensao ativa do ToRORdLand e a
//  variavel de estado 49. Esse conhecimento mora aqui, em um lugar so.

//! Como ler Ta, Vm e o estimulo de um modelo celular
struct CellProfile
{
  //! Variavel de estado que guarda Ta. < 0: valor monitorado 0 (modelos
  //! fenomenologicos como Kerckhoffs, que publicam active_stress).
  int    ta_var  = -1;

  //! Ta do modelo celular -> unidade de pressao do modelo. ToRORdLand entrega
  //! kPa (fator = kPa_to_p); Kerckhoffs usa o T_ref do <pvloop>.
  double ta_unit = 1.0;

  //! Variavel de estado do potencial de membrana, para a saida por no.
  //! Kerckhoffs nao tem Vm: -1 desliga a saida em vez de gravar uma variavel
  //! sem relacao (o estado 0 dele e o comprimento do elemento contratil).
  int    vm_var  = -1;

  //! Estados do submodelo de distorcao de Land; -1 = nao se aplica
  int land_xs = -1, land_xw = -1, land_zetas = -1, land_zetaw = -1;

  //! Estimulo aplicado a cada no quando seu LAT e atingido. Amplitude na
  //! unidade de corrente do modelo, duracao na unidade do solver. Amplitude 0
  //! o desliga, o correto para modelos dirigidos pelo LAT (Kerckhoffs).
  double stim_amplitude = 0.0;
  double stim_duration  = 0.002;

  //! So o submodelo de Land le o estiramento de fibra
  bool   uses_stretch   = false;

  //! Monta o perfil do modelo `name`; le -stimamp/-stimdur e imprime o resumo
  static CellProfile for_model(const std::string & name, double kPa_to_p,
                               double T_ref, double ms_to_solver);
};

//! Calculos que so fazem sentido para o submodelo de Land (ToRORdLand):
//! rigidez ativa para a estabilizacao da tensao dependente da velocidade e o
//! diagnostico de Ta negativa.
class LandCoupling
{
public:
  void set_state_indices(int xs, int xw, int zetas, int zetaw)
  { xs_ = xs; xw_ = xw; zetas_ = zetas; zetaw_ = zetaw; }

  bool has_stiffness_states() const { return xs_ >= 0 && xw_ >= 0; }

  //! Ka nodal do submodelo de Land [kPa], dado lambda_f de cada no
  void active_stiffness(const Cells & cells, const arma::vec & lambda,
                        arma::vec & ka);

  //! Diz QUAL termo da tensao ativa de Land ficou negativo, e a que taxa de
  //! estiramento. Diagnostico apenas: nada e modificado.
  void report_negative_ta(const Cells & cells, const arma::vec & ta,
                          const arma::vec & lambda_rate, double rate_clip);

private:
  int xs_ = -1, xw_ = -1, zetas_ = -1, zetaw_ = -1;

  //! Buffers do diagnostico, membros para nao realocar quatro vetores nodais
  //! a cada chamada.
  arma::vec buf_xs_, buf_xw_, buf_zetas_, buf_zetaw_;
};

#endif