#include "kdtree_uvc.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "kdtree.hpp"      // cdalitz/kdtree-cpp, sem modificacoes

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// =============================================================================
//  Auxiliares locais (namespace anonimo: nao vazam para o resto do projeto,
//  nem colidem com os homonimos do busca_uvc.cpp)
// =============================================================================

namespace
{
  const double PI_     = 3.1415926535897932384626433832795;
  const double DOIS_PI = 2.0 * PI_;
  const double NAO_NUM = 0.0 / 0.0;

  //! true se o valor nao e NaN nem infinito (sem depender do C++11)
  bool eh_finito(double v)
  {
    return (v == v) && (v > -1e300) && (v < 1e300);
  }

  //! Embedding cilindrico de um no: E = [ab*w_ab, tm*w_tm,
  //!                                     ab*cos(rt)*w_rt, ab*sin(rt)*w_rt]
  //! O cos/sin resolve a costura +-pi sem copias-fantasma, e o raio ab -> 0
  //! faz o rt deixar de pesar no apice, onde ele degenera.
  void embutir(double ab, double tm, double rt, const OpcoesUVC & o,
               double e[4])
  {
    double cx = ab * std::cos(rt);
    double cy = ab * std::sin(rt);
    if (!eh_finito(cx) || !eh_finito(cy)) { cx = 0.0; cy = 0.0; }  // rt NaN

    e[0] = ab * o.w_ab;
    e[1] = tm * o.w_tm;
    e[2] = cx * o.w_rt;
    e[3] = cy * o.w_rt;
  }

  //! Vizinho mais proximo de q numa nuvem de pontos de dimensao dim.
  //! Devolve a posicao do ponto na lista, ou -1.
  int mais_proximo(Kdtree::KdTree & arvore, const double * q, int dim)
  {
    Kdtree::KdNodeVector viz;
    arvore.k_nearest_neighbors(Kdtree::CoordPoint(q, q + dim), 1, &viz);
    if (viz.empty()) return -1;
    return viz[0].index;
  }

}  // namespace anonimo

// =============================================================================
//  Leitura
// =============================================================================

bool ler_campo_cru_uvc(ReaderHDF5 & r, const string & nome,
                       vector<double> & val, bool & por_celula,
                       string & nome_real)
{
  // nome pode ser so o nome ("tecido") ou o caminho ("vertex_field/tecido")
  const int idx = r.find_field(nome);
  if (idx < 0) return false;

  const FieldInfo & f = r.get_field(idx);
  if (f.n_comp != 1)
  {
    cerr << "[ERRO] campo '" << f.name << "' tem " << f.n_comp
         << " componentes; esperado escalar." << endl;
    return false;
  }

  if (!r.read_field_step(f.path, 0, val)) return false;

  por_celula = f.cell_centered;
  nome_real  = f.path;          // caminho: chave unica para as proximas leituras
  return true;
}

// -----------------------------------------------------------------------------

bool carregar_malha_uvc(ReaderHDF5 & leitor, const string & arquivo,
                        const string & rotulo, MalhaTransf & m,
                        const string & n_ab, const string & n_tm,
                        const string & n_rt, const string & n_tv)
{
  if (!leitor.open(arquivo))
  {
    cerr << "[ERRO] falha ao abrir '" << arquivo << "'" << endl;
    return false;
  }

  m.rotulo     = rotulo;
  m.n_points   = leitor.get_n_points();
  m.n_elements = leitor.get_n_elements();
  m.nen        = leitor.get_nen();
  m.xyz        = leitor.get_coordinates();
  m.tets       = leitor.get_connectivity();

  const string nomes[4] = { n_ab, n_tm, n_rt, n_tv };
  vector<double> * destino[4] = { &m.ab, &m.tm, &m.rt, &m.tv };
  string faltando;

  for (int c = 0; c < 4; c++)
  {
    bool por_celula = false;
    string nome_real;
    vector<double> v;

    if (!ler_campo_cru_uvc(leitor, nomes[c], v, por_celula, nome_real))
    {
      faltando += (faltando.empty() ? "" : ", ") + nomes[c];
      continue;
    }

    if (por_celula)
    {
      cout << "  [AVISO] UVC '" << nomes[c] << "' esta por ELEMENTO; "
           << "convertida para os nos por media." << endl;
      v = celula_para_no(m, v, false);
    }
    *destino[c] = v;
  }

  if (!faltando.empty())
  {
    cerr << "[ERRO] campos UVC ausentes em '" << arquivo << "': "
         << faltando << endl;
    cerr << "       campos disponiveis:";
    for (int i = 0; i < leitor.get_n_fields(); i++)
      cerr << " " << leitor.get_field(i).name;
    cerr << endl;
    return false;
  }

  // ------------------------------------------------------- unidade do rt
  // O embedding usa cos(rt)/sin(rt): rt precisa estar em RADIANOS.
  double lo = 1e300, hi = -1e300;
  for (int i = 0; i < m.n_points; i++)
  {
    if (!eh_finito(m.rt[(size_t) i])) continue;
    if (m.rt[(size_t) i] < lo) lo = m.rt[(size_t) i];
    if (m.rt[(size_t) i] > hi) hi = m.rt[(size_t) i];
  }

  if (!(hi > 1.6 || lo < -1.6))
  {
    cout << "  [INFO] " << rotulo << ": rt parece estar em voltas ["
         << lo << ", " << hi << "]; convertendo para radianos (rt * 2*pi)."
         << endl;
    for (int i = 0; i < m.n_points; i++)
      if (eh_finito(m.rt[(size_t) i])) m.rt[(size_t) i] *= DOIS_PI;
  }

  cout << "  " << rotulo << " : " << m.n_points << " nos, "
       << m.n_elements << " celulas (" << m.nen << " nos cada)"
       << endl;

  return true;
}

// =============================================================================
//  Montagem -- CHAMADA UMA VEZ
// =============================================================================

bool kdtree_uvc_build(const MalhaTransf & fonte, const MalhaTransf & alvo,
                      const OpcoesUVC & o, KdtreeUVC & est)
{
  const int nf = fonte.n_points;
  const int na = alvo.n_points;

  if (nf <= 0 || na <= 0)
  {
    cerr << "[ERRO] kdtree_uvc_build: malha vazia (fonte " << nf
         << " nos, alvo " << na << " nos)." << endl;
    return false;
  }

  if ((int) fonte.ab.size() != nf || (int) fonte.tm.size() != nf ||
      (int) fonte.rt.size() != nf || (int) fonte.tv.size() != nf ||
      (int) alvo.ab.size()  != na || (int) alvo.tm.size()  != na ||
      (int) alvo.rt.size()  != na || (int) alvo.tv.size()  != na)
  {
    cerr << "[ERRO] kdtree_uvc_build: UVC com tamanho diferente do numero "
         << "de nos." << endl;
    return false;
  }

  if (o.k < 1)
  {
    cerr << "[ERRO] kdtree_uvc_build: k deve ser >= 1." << endl;
    return false;
  }
  if (o.peso != "gauss" && o.peso != "idw")
  {
    cerr << "[ERRO] kdtree_uvc_build: peso deve ser 'gauss' ou 'idw'." << endl;
    return false;
  }

  est.n_fonte = nf;
  est.n_alvo  = na;

  // ------------------------------------------------------------ 1) mascaras
  // Repare que o campo NAO entra aqui: a validade e so das UVC. Vizinhos com
  // valor nao-finito sao descartados depois, em kdtree_uvc_transf().
  vector<char> ok_s((size_t) nf, 0);
  vector<char> ok_t((size_t) na, 0);
  int n_ok_s = 0, n_ok_t = 0;

  for (int i = 0; i < nf; i++)
  {
    ok_s[(size_t) i] = (char) (eh_finito(fonte.ab[(size_t) i]) &&
                               eh_finito(fonte.tm[(size_t) i]) &&
                               eh_finito(fonte.tv[(size_t) i]));
    n_ok_s += ok_s[(size_t) i];
  }
  for (int i = 0; i < na; i++)
  {
    ok_t[(size_t) i] = (char) (eh_finito(alvo.ab[(size_t) i]) &&
                               eh_finito(alvo.tm[(size_t) i]) &&
                               eh_finito(alvo.tv[(size_t) i]));
    n_ok_t += ok_t[(size_t) i];
  }

  if (n_ok_s == 0)
  {
    cerr << "[ERRO] kdtree_uvc_build: nenhum no valido na fonte." << endl;
    return false;
  }

  if (!o.silencioso)
    cout << "  nos fonte validos: " << n_ok_s << "/" << nf
         << " | alvo validos: " << n_ok_t << "/" << na << endl;

  // ----------------------------------------------------------- 2) tv_split
  if (o.tem_tv_split) est.tv_split = o.tv_split;
  else
  {
    double lo = 1e300, hi = -1e300;
    for (int i = 0; i < nf; i++)
      if (eh_finito(fonte.tv[(size_t) i]))
      {
        if (fonte.tv[(size_t) i] < lo) lo = fonte.tv[(size_t) i];
        if (fonte.tv[(size_t) i] > hi) hi = fonte.tv[(size_t) i];
      }
    for (int i = 0; i < na; i++)
      if (eh_finito(alvo.tv[(size_t) i]))
      {
        if (alvo.tv[(size_t) i] < lo) lo = alvo.tv[(size_t) i];
        if (alvo.tv[(size_t) i] > hi) hi = alvo.tv[(size_t) i];
      }
    est.tv_split = (lo > hi) ? 0.5 : 0.5 * (lo + hi);
  }
  const double tv_split = est.tv_split;

  if (!o.silencioso)
    cout << "  TV_SPLIT = " << tv_split << "  (VE: tv<split, VD: tv>=split)"
         << endl;

  // ------------------------------------------------- 3) clamp do ab do alvo
  vector<double> ab_t = alvo.ab;
  if (!o.sem_clamp_ab)
  {
    double lo = 1e300, hi = -1e300;
    for (int i = 0; i < nf; i++)
      if (ok_s[(size_t) i])
      {
        if (fonte.ab[(size_t) i] < lo) lo = fonte.ab[(size_t) i];
        if (fonte.ab[(size_t) i] > hi) hi = fonte.ab[(size_t) i];
      }

    for (int i = 0; i < na; i++)
      if (eh_finito(ab_t[(size_t) i]))
      {
        if (ab_t[(size_t) i] < lo) ab_t[(size_t) i] = lo;
        if (ab_t[(size_t) i] > hi) ab_t[(size_t) i] = hi;
      }

    if (!o.silencioso)
      cout << "  ab do alvo clampado a [" << lo << ", " << hi << "]." << endl;
  }

  // ------------------------------------------------ 4) embedding cilindrico
  vector<double> E_s((size_t) 4 * nf, 0.0);
  vector<double> E_t((size_t) 4 * na, 0.0);

  for (int i = 0; i < nf; i++)
    embutir(fonte.ab[(size_t) i], fonte.tm[(size_t) i], fonte.rt[(size_t) i],
            o, &E_s[(size_t) 4 * i]);

  for (int i = 0; i < na; i++)
    embutir(ab_t[(size_t) i], alvo.tm[(size_t) i], alvo.rt[(size_t) i],
            o, &E_t[(size_t) 4 * i]);

  // --------------------------------- 5) quantos nos da fonte por ventriculo
  // Precisa vir antes de alocar as linhas: k_max e o maior k efetivo dos
  // dois ventriculos.
  int n_src[2] = { 0, 0 };
  for (int i = 0; i < nf; i++)
  {
    if (!ok_s[(size_t) i]) continue;
    n_src[(fonte.tv[(size_t) i] < tv_split) ? 0 : 1]++;
  }

  int kk[2];
  for (int v = 0; v < 2; v++)
    kk[v] = (o.k < n_src[v]) ? o.k : n_src[v];

  est.k_max = (kk[0] > kk[1]) ? kk[0] : kk[1];
  if (est.k_max < 1)
  {
    cerr << "[ERRO] kdtree_uvc_build: nenhum ventriculo com nos na fonte."
         << endl;
    return false;
  }

  est.viz.assign((size_t) est.k_max * na, -1);
  est.pesos.assign((size_t) est.k_max * na, 0.0);
  est.n_viz.assign((size_t) na, 0);

  // ------------------------------------------------- 6) k-NN por ventriculo
  for (int vent = 0; vent < 2; vent++)
  {
    const string nome = vent ? "VD" : "VE";

    // --- nos da fonte deste ventriculo (KdNode::index = posicao em id_s) ---
    Kdtree::KdNodeVector pts;
    vector<int> id_s;
    for (int i = 0; i < nf; i++)
    {
      if (!ok_s[(size_t) i]) continue;
      const bool ve = (fonte.tv[(size_t) i] < tv_split);
      if (ve == (vent == 1)) continue;
      const double * e = &E_s[(size_t) 4 * i];
      pts.push_back(Kdtree::KdNode(Kdtree::CoordPoint(e, e + 4), NULL,
                                   (int) id_s.size()));
      id_s.push_back(i);
    }

    // --- nos do alvo deste ventriculo ---
    vector<int> id_t;
    for (int i = 0; i < na; i++)
    {
      if (!ok_t[(size_t) i]) continue;
      const bool ve = (alvo.tv[(size_t) i] < tv_split);
      if (ve == (vent == 1)) continue;
      id_t.push_back(i);
    }

    if (id_s.empty() || id_t.empty())
    {
      if (!o.silencioso)
        cout << "  " << nome << ": sem nos validos, pulando." << endl;
      continue;
    }

    // id_s nao esta vazio: a Kdtree::KdTree aceita a lista
    Kdtree::KdTree arvore(&pts);
    Kdtree::KdNodeVector().swap(pts);   // a arvore guarda a propria copia

    const int k = kk[vent];
    Kdtree::CoordPoint   q(4);
    Kdtree::KdNodeVector viz;           // em ordem crescente de distancia
    vector<double> d((size_t) k, 0.0);
    vector<double> w((size_t) k, 0.0);

    for (size_t t = 0; t < id_t.size(); t++)
    {
      const int i = id_t[t];
      for (int c = 0; c < 4; c++) q[(size_t) c] = E_t[(size_t) 4 * i + c];

      arvore.k_nearest_neighbors(q, (size_t) k, &viz);
      const int nv = (int) viz.size();
      if (nv == 0) continue;

      // a kdtree-cpp nao devolve as distancias: recalcula a partir do ponto
      for (int j = 0; j < nv; j++)
      {
        const Kdtree::CoordPoint & p = viz[(size_t) j].point;
        double s2 = 0.0;
        for (int c = 0; c < 4; c++)
        {
          const double dc = q[(size_t) c] - p[(size_t) c];
          s2 += dc * dc;
        }
        d[(size_t) j] = std::sqrt(s2);
      }

      // --- pesos ---
      if (o.peso == "gauss")
      {
        double media = 0.0;
        for (int j = 0; j < nv; j++) media += d[(size_t) j];
        media /= (double) nv;
        const double h = (media > 1e-9) ? media : 1e-9;
        for (int j = 0; j < nv; j++)
        {
          const double z = d[(size_t) j] / h;
          w[(size_t) j] = std::exp(-z * z);
        }
      }
      else
      {
        for (int j = 0; j < nv; j++)
        {
          const double dd = (d[(size_t) j] > 1e-12) ? d[(size_t) j] : 1e-12;
          w[(size_t) j] = 1.0 / std::pow(dd, o.pot_idw);
        }
      }

      // coincidencia exata: so o primeiro vizinho conta
      if (d[0] < 1e-12)
      {
        for (int j = 0; j < nv; j++) w[(size_t) j] = 0.0;
        w[0] = 1.0;
      }

      double soma = 0.0;
      for (int j = 0; j < nv; j++) soma += w[(size_t) j];
      if (soma <= 0.0) { w[0] = 1.0; soma = 1.0; }

      // --- guarda a linha ja normalizada ---
      const size_t base = (size_t) est.k_max * i;
      for (int j = 0; j < nv; j++)
      {
        est.viz[base + j]   = id_s[(size_t) viz[(size_t) j].index];
        est.pesos[base + j] = w[(size_t) j] / soma;
      }
      est.n_viz[(size_t) i] = nv;
    }

    if (!o.silencioso)
      cout << "  " << nome << ": " << id_t.size() << " nos do alvo (k="
           << k << ")." << endl;
  }

  // --------------------------- 7) nos do alvo sem vizinho: copia de um par
  // Mesmo papel do "preenche NaN por vizinho UVC" do script original, so que
  // resolvido aqui: quem copia de quem depende so da geometria.
  est.copia_destino.clear();
  est.copia_origem.clear();
  est.n_sem_valor = 0;

  {
    Kdtree::KdNodeVector pts_e;       // com valor, no espaco do embedding
    Kdtree::KdNodeVector pts_x;       // os mesmos, no espaco fisico
    vector<int> com_valor;
    vector<int> sem_valor;

    for (int i = 0; i < na; i++)
    {
      if (est.n_viz[(size_t) i] > 0)
      {
        const double * e = &E_t[(size_t) 4 * i];
        const double * x = &alvo.xyz[(size_t) 3 * i];
        pts_e.push_back(Kdtree::KdNode(Kdtree::CoordPoint(e, e + 4), NULL,
                                       (int) com_valor.size()));
        pts_x.push_back(Kdtree::KdNode(Kdtree::CoordPoint(x, x + 3), NULL,
                                       (int) com_valor.size()));
        com_valor.push_back(i);
      }
      else sem_valor.push_back(i);
    }

    if (!sem_valor.empty() && !com_valor.empty())
    {
      Kdtree::KdTree arv_e(&pts_e);

      // A arvore fisica so e montada se algum no orfao tiver embedding
      // nao-finito (ab ou tm NaN): consultar a arvore do embedding com NaN
      // devolveria lixo.
      Kdtree::KdTree * arv_x = 0;

      for (size_t t = 0; t < sem_valor.size(); t++)
      {
        const int i = sem_valor[t];
        const double * e = &E_t[(size_t) 4 * i];

        bool e_ok = true;
        for (int c = 0; c < 4; c++) if (!eh_finito(e[c])) e_ok = false;

        int pos = -1;
        if (e_ok)
        {
          pos = mais_proximo(arv_e, e, 4);
        }
        else
        {
          if (arv_x == 0) arv_x = new Kdtree::KdTree(&pts_x);
          pos = mais_proximo(*arv_x, &alvo.xyz[(size_t) 3 * i], 3);
        }

        if (pos < 0) { est.n_sem_valor++; continue; }

        est.copia_destino.push_back(i);
        est.copia_origem.push_back(com_valor[(size_t) pos]);
      }

      delete arv_x;

      if (!o.silencioso)
        cout << "  " << est.copia_destino.size()
             << " nos do alvo preenchidos por vizinho UVC." << endl;
    }
    else est.n_sem_valor = (int) sem_valor.size();
  }

  return true;
}

// =============================================================================
//  Transferencia -- CHAMADA EM LOOP
// =============================================================================

bool kdtree_uvc_transf(const KdtreeUVC & est,
                       const vector<double> & campo_fonte,
                       bool categorico,
                       vector<double> & campo_alvo)
{
  if ((int) campo_fonte.size() != est.n_fonte)
  {
    cerr << "[ERRO] kdtree_uvc_transf: campo com " << campo_fonte.size()
         << " valores; a fonte tem " << est.n_fonte << " nos." << endl;
    return false;
  }
  if (est.k_max < 1 || (int) est.n_viz.size() != est.n_alvo)
  {
    cerr << "[ERRO] kdtree_uvc_transf: estrutura nao montada "
         << "(chame kdtree_uvc_build antes)." << endl;
    return false;
  }

  campo_alvo.assign((size_t) est.n_alvo, NAO_NUM);

  const int kmax = est.k_max;

  for (int i = 0; i < est.n_alvo; i++)
  {
    const int nv = est.n_viz[(size_t) i];
    if (nv <= 0) continue;

    const size_t base = (size_t) kmax * i;

    // ----------------------------------------------------------- continuo
    if (!categorico)
    {
      double v = 0.0, soma = 0.0;
      for (int j = 0; j < nv; j++)
      {
        const double val = campo_fonte[(size_t) est.viz[base + j]];
        if (!eh_finito(val)) continue;          // vizinho sem valor: descarta
        v    += est.pesos[base + j] * val;
        soma += est.pesos[base + j];
      }
      if (soma > 0.0) campo_alvo[(size_t) i] = v / soma;   // renormaliza
      continue;
    }

    // --------------------------------------------------------- categorico
    // Voto ponderado sobre os valores que aparecem entre os vizinhos. Como
    // nv e pequeno (k ~ 12), varrer as repeticoes sai mais barato do que
    // manter a lista global de classes -- e dispensa passa-la a funcao.
    // Empate vai para o rotulo MENOR, como no script Python.
    double melhor_w = -1.0;
    double melhor_v = NAO_NUM;

    for (int j = 0; j < nv; j++)
    {
      const double val = campo_fonte[(size_t) est.viz[base + j]];
      if (!eh_finito(val)) continue;

      bool repetido = false;
      for (int j2 = 0; j2 < j && !repetido; j2++)
        if (campo_fonte[(size_t) est.viz[base + j2]] == val) repetido = true;
      if (repetido) continue;

      double sw = 0.0;
      for (int j2 = 0; j2 < nv; j2++)
        if (campo_fonte[(size_t) est.viz[base + j2]] == val)
          sw += est.pesos[base + j2];

      if (sw > melhor_w || (sw == melhor_w && val < melhor_v))
      {
        melhor_w = sw;
        melhor_v = val;
      }
    }

    if (melhor_w >= 0.0) campo_alvo[(size_t) i] = melhor_v;
  }

  // ------------------------------------- copias (nos do alvo sem vizinho)
  for (size_t t = 0; t < est.copia_destino.size(); t++)
    campo_alvo[(size_t) est.copia_destino[t]] =
        campo_alvo[(size_t) est.copia_origem[t]];

  return true;
}

// =============================================================================
//  Conversao PointData <-> CellData
// =============================================================================

vector<double> classes_de_uvc(const vector<double> & v)
{
  vector<double> c;
  for (size_t i = 0; i < v.size(); i++)
    if (eh_finito(v[i])) c.push_back(v[i]);

  std::sort(c.begin(), c.end());
  c.erase(std::unique(c.begin(), c.end()), c.end());
  return c;
}

// -----------------------------------------------------------------------------

vector<double> celula_para_no(const MalhaTransf & m,
                              const vector<double> & vals_cell,
                              bool categorico)
{
  const int np  = m.n_points;
  const int ne  = m.n_elements;
  const int nen = m.nen;

  const vector<double> classes = classes_de_uvc(vals_cell);

  // ------------------------------------------------------------- continuo
  if (!categorico)
  {
    vector<double> soma((size_t) np, 0.0);
    vector<double> cont((size_t) np, 0.0);

    for (int e = 0; e < ne; e++)
      for (int j = 0; j < nen; j++)
      {
        const int no = m.tets[(size_t) nen * e + j];
        soma[(size_t) no] += vals_cell[(size_t) e];
        cont[(size_t) no] += 1.0;
      }

    for (int i = 0; i < np; i++)
      if (cont[(size_t) i] > 0.0) soma[(size_t) i] /= cont[(size_t) i];
    return soma;
  }

  // --------------------------------------------------- categorico binario
  if (classes.size() <= 2)
  {
    const double menor = classes.empty() ? 0.0 : classes.front();
    const double maior = (classes.size() == 2) ? classes.back() : menor;

    vector<double> soma((size_t) np, 0.0);
    vector<double> cont((size_t) np, 0.0);

    for (int e = 0; e < ne; e++)
    {
      const double b = (classes.size() == 2 && vals_cell[(size_t) e] > menor)
                       ? 1.0 : 0.0;
      for (int j = 0; j < nen; j++)
      {
        const int no = m.tets[(size_t) nen * e + j];
        soma[(size_t) no] += b;
        cont[(size_t) no] += 1.0;
      }
    }

    vector<double> out((size_t) np, menor);
    for (int i = 0; i < np; i++)
    {
      const double f = cont[(size_t) i] ? soma[(size_t) i] / cont[(size_t) i]
                                        : 0.0;
      // Mesmo criterio do np.round do script Python: empate exato (f = 0.5)
      // vai para a classe MENOR, porque o numpy arredonda 0.5 para o par.
      out[(size_t) i] = (f > 0.5) ? maior : menor;
    }
    return out;
  }

  // ----------------------------------------------- categorico multiclasse
  const size_t nc = classes.size();
  vector<double> contagem((size_t) np * nc, 0.0);

  for (int e = 0; e < ne; e++)
  {
    const size_t c = (size_t) (std::lower_bound(classes.begin(), classes.end(),
                                                vals_cell[(size_t) e])
                               - classes.begin());
    if (c >= nc) continue;
    for (int j = 0; j < nen; j++)
      contagem[(size_t) m.tets[(size_t) nen * e + j] * nc + c] += 1.0;
  }

  vector<double> out((size_t) np, classes[0]);
  for (int i = 0; i < np; i++)
  {
    size_t melhor = 0;
    for (size_t c = 1; c < nc; c++)
      if (contagem[(size_t) i * nc + c] > contagem[(size_t) i * nc + melhor])
        melhor = c;
    out[(size_t) i] = classes[melhor];
  }
  return out;
}

// -----------------------------------------------------------------------------

vector<double> no_para_celula(const MalhaTransf & m,
                              const vector<double> & vals_no,
                              bool categorico, double frac)
{
  const int ne  = m.n_elements;
  const int nen = m.nen;

  vector<double> out((size_t) ne, 0.0);

  // ------------------------------------------------------------- continuo
  if (!categorico)
  {
    for (int e = 0; e < ne; e++)
    {
      double soma = 0.0;
      int    n    = 0;
      for (int j = 0; j < nen; j++)
      {
        const double v = vals_no[(size_t) m.tets[(size_t) nen * e + j]];
        if (!eh_finito(v)) continue;
        soma += v;
        n++;
      }
      out[(size_t) e] = n ? soma / (double) n : NAO_NUM;
    }
    return out;
  }

  const vector<double> classes = classes_de_uvc(vals_no);

  // ----- binario com zero como classe negativa: fracao de nos positivos ---
  const bool binario = (classes.size() <= 2 && !classes.empty() &&
                        classes.front() == 0.0);

  if (binario)
  {
    const double cls_pos = (classes.size() == 2) ? classes.back() : 1.0;

    for (int e = 0; e < ne; e++)
    {
      int pos = 0;
      for (int j = 0; j < nen; j++)
        if (vals_no[(size_t) m.tets[(size_t) nen * e + j]] == cls_pos) pos++;
      const double f = (double) pos / (double) nen;
      out[(size_t) e] = (f >= frac) ? cls_pos : 0.0;
    }
    return out;
  }

  // -------------------------------------------- multiclasse: moda dos nos
  for (int e = 0; e < ne; e++)
  {
    double melhor   = vals_no[(size_t) m.tets[(size_t) nen * e]];
    int    melhor_n = 0;

    for (int j = 0; j < nen; j++)
    {
      const double v = vals_no[(size_t) m.tets[(size_t) nen * e + j]];
      int n = 0;
      for (int j2 = 0; j2 < nen; j2++)
        if (vals_no[(size_t) m.tets[(size_t) nen * e + j2]] == v) n++;
      if (n > melhor_n) { melhor_n = n; melhor = v; }
    }
    out[(size_t) e] = melhor;
  }
  return out;
}

// =============================================================================
//  Escrita da malha alvo em .vtu (ASCII, escrito a mao -- sem VTK)
// =============================================================================

bool salvar_vtu_uvc(const string & caminho, const MalhaTransf & m,
                    const vector<CampoSaida> & campos)
{
  std::ofstream f(caminho.c_str());
  if (!f)
  {
    cerr << "[ERRO] nao foi possivel gravar '" << caminho << "'" << endl;
    return false;
  }

  // tipo de celula do VTK a partir do numero de nos por elemento
  int vtk_tipo = 10;                       // VTK_TETRA
  if      (m.nen == 3) vtk_tipo = 5;       // VTK_TRIANGLE
  else if (m.nen == 8) vtk_tipo = 12;      // VTK_HEXAHEDRON
  else if (m.nen == 2) vtk_tipo = 3;       // VTK_LINE

  f << std::setprecision(10);
  f << "<?xml version=\"1.0\"?>\n";
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    << "byte_order=\"LittleEndian\">\n";
  f << "  <UnstructuredGrid>\n";
  f << "    <Piece NumberOfPoints=\"" << m.n_points
    << "\" NumberOfCells=\"" << m.n_elements << "\">\n";

  // ----------------------------------------------------------- pontos
  f << "      <Points>\n";
  f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
    << "format=\"ascii\">\n";
  for (int i = 0; i < m.n_points; i++)
    f << "          " << m.xyz[(size_t) 3 * i + 0] << " "
      << m.xyz[(size_t) 3 * i + 1] << " " << m.xyz[(size_t) 3 * i + 2] << "\n";
  f << "        </DataArray>\n      </Points>\n";

  // ----------------------------------------------------------- celulas
  f << "      <Cells>\n";
  f << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
    << "format=\"ascii\">\n";
  for (int e = 0; e < m.n_elements; e++)
  {
    f << "          ";
    for (int j = 0; j < m.nen; j++) f << m.tets[(size_t) m.nen * e + j] << " ";
    f << "\n";
  }
  f << "        </DataArray>\n";

  f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
    << "          ";
  for (int e = 0; e < m.n_elements; e++) f << (m.nen * (e + 1)) << " ";
  f << "\n        </DataArray>\n";

  f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
    << "          ";
  for (int e = 0; e < m.n_elements; e++) f << vtk_tipo << " ";
  f << "\n        </DataArray>\n      </Cells>\n";

  // ------------------------------------------------------------ campos
  f << "      <PointData>\n";
  for (size_t c = 0; c < campos.size(); c++)
  {
    if (campos[c].por_celula) continue;
    f << "        <DataArray type=\"Float64\" Name=\"" << campos[c].nome
      << "\" NumberOfComponents=\"1\" format=\"ascii\">\n          ";
    for (size_t i = 0; i < campos[c].val.size(); i++)
      f << campos[c].val[i] << " ";
    f << "\n        </DataArray>\n";
  }
  f << "      </PointData>\n";

  f << "      <CellData>\n";
  for (size_t c = 0; c < campos.size(); c++)
  {
    if (!campos[c].por_celula) continue;
    f << "        <DataArray type=\"Float64\" Name=\"" << campos[c].nome
      << "\" NumberOfComponents=\"1\" format=\"ascii\">\n          ";
    for (size_t i = 0; i < campos[c].val.size(); i++)
      f << campos[c].val[i] << " ";
    f << "\n        </DataArray>\n";
  }
  f << "      </CellData>\n";

  f << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
  f.close();
  return true;
}
