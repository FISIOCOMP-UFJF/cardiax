// =============================================================================
//  teste_kdtree_uvc.cpp
//
//  Teste em loop pequeno da interface de kdtree_uvc.hpp: monta a estrutura
//  UMA vez (kdtree_uvc_build) e transfere o campo VARIAS vezes
//  (kdtree_uvc_transf), cronometrando as duas partes separadamente.
//
//  E a prova do que a separacao vale: o build paga o k-NN uma vez e cada
//  passo do loop vira uma media ponderada sobre listas prontas.
//
//  Se o campo for uma serie temporal, cada iteracao le um passo diferente do
//  HDF5 (passo % n_steps). Se for estatico, o mesmo vetor e reaproveitado --
//  o custo medido continua sendo o da transferencia.
//
//  Uso:
//      ./teste_kdtree_uvc Paciente_1 Paciente_2 --campo lat --passos 20
//      ./teste_kdtree_uvc P1 P2 --campo fecido --tipo categorico --passos 5
//      ./teste_kdtree_uvc P1 P2 --campo lat --saida alvo_lat.vtu --nos 783,5510
// =============================================================================

#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "kdtree_uvc.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// =============================================================================
//  Linha de comando
// =============================================================================

struct Opcoes
{
  string fonte;
  string alvo;
  string campo;
  string tipo;
  string saida;
  string nos;           //!< ids do alvo a imprimir, separados por virgula
  int    passos;
  OpcoesUVC uvc;

  Opcoes() : fonte(), alvo(), campo(), tipo("continuo"), saida(), nos(),
             passos(10), uvc() {}
};

static void ajuda(const char * prog)
{
  cout
    << "uso: " << prog << " <fonte> <alvo> --campo NOME [opcoes]\n\n"
    << "  fonte/alvo      basename ou arquivo .xmf/.xdmf/.h5 com UVC nos nos\n\n"
    << "opcoes:\n"
    << "  --campo NOME    campo a transferir (obrigatorio): so o nome\n"
    << "                  (tecido) ou o caminho (vertex_field/tecido,\n"
    << "                  cell_field/tecido) quando o nome for ambiguo\n"
    << "  --tipo T        continuo | categorico   (padrao continuo)\n"
    << "  --passos N      quantas chamadas de kdtree_uvc_transf (padrao 10)\n"
    << "  --saida ARQ     grava o resultado da ultima chamada num .vtu\n"
    << "  --nos A,B,C     imprime o valor transferido nesses ids do alvo\n"
    << "  --k N           vizinhos por no alvo (padrao 12)\n"
    << "  --peso P        gauss | idw            (padrao gauss)\n"
    << "  --pot-idw V     expoente do IDW        (padrao 2)\n"
    << "  --w-ab V --w-tm V --w-rt V   pesos dos eixos do embedding\n"
    << "  --tv-split V    limiar VE/VD (padrao: media entre min e max de tv)\n"
    << "  --sem-clamp-ab  nao limita o ab do alvo a faixa da fonte\n"
    << endl;
}

static bool ler_opcoes(int argc, char ** argv, Opcoes & o)
{
  vector<string> posicionais;

  for (int i = 1; i < argc; i++)
  {
    const string a = argv[i];

    if (a == "-h" || a == "--help" || a == "--ajuda") { ajuda(argv[0]); std::exit(0); }
    else if (a == "--sem-clamp-ab") o.uvc.sem_clamp_ab = true;
    else if (a.size() > 1 && a[0] == '-')
    {
      if (++i >= argc) { cerr << "[ERRO] " << a << " exige um valor." << endl; return false; }
      const string v = argv[i];

      if      (a == "--campo")    o.campo = v;
      else if (a == "--tipo")     o.tipo = v;
      else if (a == "--saida")    o.saida = v;
      else if (a == "--nos")      o.nos = v;
      else if (a == "--passos")   o.passos = atoi(v.c_str());
      else if (a == "--k")        o.uvc.k = atoi(v.c_str());
      else if (a == "--peso")     o.uvc.peso = v;
      else if (a == "--pot-idw")  o.uvc.pot_idw = atof(v.c_str());
      else if (a == "--w-ab")     o.uvc.w_ab = atof(v.c_str());
      else if (a == "--w-tm")     o.uvc.w_tm = atof(v.c_str());
      else if (a == "--w-rt")     o.uvc.w_rt = atof(v.c_str());
      else if (a == "--tv-split") { o.uvc.tv_split = atof(v.c_str()); o.uvc.tem_tv_split = true; }
      else { cerr << "[ERRO] opcao desconhecida: " << a << endl; return false; }
    }
    else posicionais.push_back(a);
  }

  if (posicionais.size() != 2) { ajuda(argv[0]); return false; }

  o.fonte = posicionais[0];
  o.alvo  = posicionais[1];

  if (o.campo.empty())
  {
    cerr << "[ERRO] --campo e obrigatorio neste teste." << endl;
    return false;
  }
  if (o.tipo != "continuo" && o.tipo != "categorico")
  {
    cerr << "[ERRO] --tipo deve ser 'continuo' ou 'categorico'." << endl;
    return false;
  }
  if (o.passos < 1) o.passos = 1;

  return true;
}

//! "783,5510,7133" -> {783, 5510, 7133}
static vector<int> ler_ids(const string & s)
{
  vector<int> ids;
  size_t a = 0;
  while (a <= s.size())
  {
    const size_t b = s.find(',', a);
    const string t = s.substr(a, (b == string::npos) ? string::npos : b - a);
    if (!t.empty()) ids.push_back(atoi(t.c_str()));
    if (b == string::npos) break;
    a = b + 1;
  }
  return ids;
}

static double segundos(clock_t a, clock_t b)
{
  return (double) (b - a) / (double) CLOCKS_PER_SEC;
}

// =============================================================================
//  Programa principal
// =============================================================================

int main(int argc, char ** argv)
{
  Opcoes o;
  if (!ler_opcoes(argc, argv, o)) return 1;

  const bool categorico = (o.tipo == "categorico");
  cout << std::fixed << std::setprecision(6);
  cout << string(64, '=') << endl;

  // ---------------------------------------------------------- 0) leitura
  cout << "Lendo malhas..." << endl;

  ReaderHDF5 leitor_f, leitor_a;
  MalhaTransf mf, ma;

  if (!carregar_malha_uvc(leitor_f, o.fonte, "fonte", mf)) return 1;
  if (!carregar_malha_uvc(leitor_a, o.alvo,  "alvo ", ma)) return 1;

  // ----------------------------------------------- 1) campo (passo 0) ----
  vector<double> campo_fonte;
  bool por_celula = false;
  string nome_real;

  if (!ler_campo_cru_uvc(leitor_f, o.campo, campo_fonte, por_celula, nome_real))
  {
    cerr << "[ERRO] campo '" << o.campo << "' nao encontrado na fonte." << endl;
    cerr << "       campos disponiveis:";
    for (int i = 0; i < leitor_f.get_n_fields(); i++)
      cerr << " " << leitor_f.get_field(i).path.substr(1)
           << (leitor_f.get_field(i).cell_centered ? "(cell)" : "");
    cerr << endl;
    return 1;
  }

  // nome_real e o caminho do dataset ("/vertex_field/tecido"): chave unica
  // para as leituras do loop. nome_curto ("tecido") nomeia a saida.
  const int idx = leitor_f.find_field(nome_real);
  const int n_steps = (idx >= 0) ? leitor_f.get_field(idx).n_steps : 1;
  const string nome_curto = (idx >= 0) ? leitor_f.get_field(idx).name
                                       : nome_real;

  cout << "\nCampo '" << nome_real << "': "
       << (por_celula ? "CellData" : "PointData")
       << ", " << n_steps << " passo(s), tipo " << o.tipo << endl;

  // ------------------------------------------------- 2) BUILD (uma vez) --
  cout << "\n--- kdtree_uvc_build (1x) ---" << endl;

  KdtreeUVC est;
  const clock_t t0 = clock();
  if (!kdtree_uvc_build(mf, ma, o.uvc, est)) return 1;
  const clock_t t1 = clock();

  cout << "  tempo do build: " << segundos(t0, t1) << " s"
       << "   (k_max=" << est.k_max << ", "
       << est.copia_destino.size() << " copias, "
       << est.n_sem_valor << " sem valor)" << endl;

  // -------------------------------------------- 3) TRANSF (dentro do loop)
  cout << "\n--- kdtree_uvc_transf (" << o.passos << "x) ---" << endl;

  vector<double> campo_alvo;
  vector<double> bruto;
  double t_transf = 0.0;

  for (int passo = 0; passo < o.passos; passo++)
  {
    // serie temporal: cada iteracao usa um passo diferente do arquivo
    if (n_steps > 1)
    {
      const int s = passo % n_steps;
      if (!leitor_f.read_field_step(nome_real, s, bruto)) return 1;
      campo_fonte = por_celula ? celula_para_no(mf, bruto, categorico) : bruto;
    }
    else if (passo == 0 && por_celula)
    {
      campo_fonte = celula_para_no(mf, campo_fonte, categorico);
    }

    if ((int) campo_fonte.size() != mf.n_points)
    {
      cerr << "[ERRO] campo com " << campo_fonte.size()
           << " valores; a fonte tem " << mf.n_points << " nos." << endl;
      return 1;
    }    
    


    const clock_t a = clock();
    if (!kdtree_uvc_transf(est, campo_fonte, categorico, campo_alvo)) return 1;
    const clock_t b = clock();
    t_transf += segundos(a, b);


  
/*
      vector<CampoSaida> campos;

      CampoSaida c_ab; c_ab.nome = "ab"; c_ab.val = ma.ab; campos.push_back(c_ab);
      CampoSaida c_tm; c_tm.nome = "tm"; c_tm.val = ma.tm; campos.push_back(c_tm);
      CampoSaida c_rt; c_rt.nome = "rt"; c_rt.val = ma.rt; campos.push_back(c_rt);
      CampoSaida c_tv; c_tv.nome = "tv"; c_tv.val = ma.tv; campos.push_back(c_tv);

      // o campo sai no MESMO tipo em que entrou
      CampoSaida c_out;
      c_out.nome = nome_curto;
      if (por_celula)
      {
        c_out.por_celula = true;
        c_out.val = no_para_celula(ma, campo_alvo, categorico);
      }
      else
      {
        c_out.val = campo_alvo;
      }
      campos.push_back(c_out);

      std::string novasaida = "passo_" + std::to_string(passo) + ".vtu";

      if (!salvar_vtu_uvc(novasaida, ma, campos)) return 1;
      cout << "\nArquivo salvo: " << o.saida << endl;
      
      */









    // --- resumo desta iteracao ---
    double lo = 1e300, hi = -1e300;
    int n_nan = 0;
    for (int i = 0; i < ma.n_points; i++)
    {
      const double v = campo_alvo[(size_t) i];
      if (!(v == v)) { n_nan++; continue; }
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }

    cout << "  passo " << std::setw(3) << passo
         << "  saida [" << lo << ", " << hi << "]"
         << "  NaN=" << n_nan << endl;
  }

  cout << "\n  tempo total das " << o.passos << " transferencias: "
       << t_transf << " s" << endl;
  cout << "  media por chamada               : "
       << t_transf / (double) o.passos << " s" << endl;
  if (t_transf > 0.0)
    cout << "  build / media por chamada       : "
         << segundos(t0, t1) / (t_transf / (double) o.passos) << "x" << endl;

  // ----------------------------------------- 4) nos pedidos na linha ----
  if (!o.nos.empty())
  {
    const vector<int> ids = ler_ids(o.nos);
    cout << "\n--- nos do alvo pedidos ---" << endl;
    for (size_t t = 0; t < ids.size(); t++)
    {
      const int i = ids[t];
      if (i < 0 || i >= ma.n_points)
      {
        cout << "  no " << i << ": fora do intervalo [0, "
             << ma.n_points - 1 << "]" << endl;
        continue;
      }
      cout << "  no " << i
           << "  " << nome_curto << "=" << campo_alvo[(size_t) i]
           << "  ab=" << ma.ab[(size_t) i]
           << "  tm=" << ma.tm[(size_t) i]
           << "  tv=" << ma.tv[(size_t) i]
           << "  rt=" << ma.rt[(size_t) i]
           << "  (" << est.n_viz[(size_t) i] << " vizinhos)" << endl;
    }
  }

  // ------------------------------------------------------- 5) saida ----
  if (!o.saida.empty())
  {
    vector<CampoSaida> campos;

    CampoSaida c_ab; c_ab.nome = "ab"; c_ab.val = ma.ab; campos.push_back(c_ab);
    CampoSaida c_tm; c_tm.nome = "tm"; c_tm.val = ma.tm; campos.push_back(c_tm);
    CampoSaida c_rt; c_rt.nome = "rt"; c_rt.val = ma.rt; campos.push_back(c_rt);
    CampoSaida c_tv; c_tv.nome = "tv"; c_tv.val = ma.tv; campos.push_back(c_tv);

    // o campo sai no MESMO tipo em que entrou
    CampoSaida c_out;
    c_out.nome = nome_curto;
    if (por_celula)
    {
      c_out.por_celula = true;
      c_out.val = no_para_celula(ma, campo_alvo, categorico);
    }
    else
    {
      c_out.val = campo_alvo;
    }
    campos.push_back(c_out);

    if (!salvar_vtu_uvc(o.saida, ma, campos)) return 1;
    cout << "\nArquivo salvo: " << o.saida << endl;
  }

  leitor_f.close();
  leitor_a.close();

  cout << string(64, '=') << endl;
  cout << "Done" << endl;
  return 0;
}
