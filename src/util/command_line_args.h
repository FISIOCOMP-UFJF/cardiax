#ifndef COMMAND_LINE_ARGS_H_
#define COMMAND_LINE_ARGS_H_

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <map>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

/**
   Command line argument parser.

   Alem de ler os argumentos, REGISTRA cada leitura: qual flag, qual valor
   resolvido, qual era o default e se o valor veio da linha de comando ou do
   default. O registro e preenchido pelo proprio read(), entao nunca
   desatualiza -- uma flag que nao aparece no dump() nao e lida por ninguem.

   O objetivo nao e produzir um texto de ajuda (para isso existe o usage() do
   app), e sim um registro de PROCEDENCIA da rodada. O modo de falha caro nao
   e desconhecer uma flag, e nao perceber que ela ficou no default.
*/
class CommandLineArgs
{
  static int argc;
  static const char ** argv;

  static char invalid[20];

  //! Uma entrada do registro.
  struct Entry
  {
    std::string value;      //!< valor resolvido, como texto
    std::string def;        //!< default
    std::string help;       //!< descricao, quando fornecida
    std::string group;      //!< agrupamento para a impressao
    bool from_cli;          //!< veio da linha de comando ou do default?

    Entry() : from_cli(false) {}
  };

  static std::map<std::string, Entry> registry;

  //! Registra uma leitura. Se a mesma flag for lida mais de uma vez, a
  //! chamada que trouxer help/group preenche os campos ainda vazios.
  static void note(const char * arg, const std::string & value,
                   const std::string & def, bool from_cli,
                   const char * help, const char * group)
  {
    Entry & e = registry[arg];
    e.value    = value;
    e.def      = def;
    e.from_cli = from_cli;
    if (help  && *help)  e.help  = help;
    if (group && *group) e.group = group;
  }

  template <typename T>
  static std::string to_text(const T & v)
  {
    std::ostringstream os;
    os << v;
    return os.str();
  }

public:

  /** Initialize command line arguments parser */
  static void init(int argc_, const char ** argv_)
  {
    strcpy(invalid, "____invalid____");
    argc = argc_;
    argv = argv_;
  }

  /** Read argument without value (switch) */
  static bool hasSwitch(const char * arg)
  {
    const char* thisArg;
    int tmp_argc = argc;
    const char ** tmp_argv = argv;

    while((thisArg = *++tmp_argv) && --tmp_argc){ // check all arguments
      if (!strcmp(thisArg, arg) && (tmp_argc > 0))
        return true; // return if match and ++tmp_argv exist
    }

    return false;
  }

  /** Busca crua. Devolve nullptr se a flag NAO esta na linha de comando.

      A sentinela "____invalid____" nao distinguia ausencia de um valor
      literalmente igual a ela, e o marcador [cli]/[default] do dump()
      depende dessa distincao para nao mentir. */
  static const char * read_raw(const char * arg)
  {
    const char* thisArg;
    int tmp_argc = argc;
    const char** tmp_argv = argv;

    while((thisArg = *++tmp_argv) && (--tmp_argc)) {
      if (!strcmp(thisArg, arg) && (tmp_argc > 0))
        return *++tmp_argv;
    }

    return nullptr;
  }

  /** Read argument with double value */
  static double read(const char * arg, const double & def,
                     const char * help = "", const char * group = "")
  {
    const char * pc = read_raw(arg);
    const bool found = (pc != nullptr);
    const double v = found ? atof(pc) : def;
    note(arg, to_text(v), to_text(def), found, help, group);
    return v;
  }

  /** Read argument with int value */
  static int read(const char * arg, const int & def,
                  const char * help = "", const char * group = "")
  {
    const char * pc = read_raw(arg);
    const bool found = (pc != nullptr);
    const int v = found ? atoi(pc) : def;
    note(arg, to_text(v), to_text(def), found, help, group);
    return v;
  }

  /** Read argument with string value */
  static string read(const char * arg, const string & def,
                     const char * help = "", const char * group = "")
  {
    const char * pc = read_raw(arg);
    const bool found = (pc != nullptr);
    const string v = found ? string(pc) : def;
    note(arg, v, def, found, help, group);
    return v;
  }

  /** Read argument with boolean value */
  static bool read(const char * arg, const bool & def,
                   const char * help = "", const char * group = "")
  {
    const char * pc = read_raw(arg);
    const bool found = (pc != nullptr);
    const bool v = found ? (strcmp(pc, "false") != 0) : def;
    note(arg, v ? "true" : "false", def ? "true" : "false",
         found, help, group);
    return v;
  }

  /** Read argument with C-string value. Mantido para as chamadas que passam
      um literal como default. */
  static const char* read(const char* arg, const char * def,
                          const char * help = "", const char * group = "")
  {
    const char * pc = read_raw(arg);
    const bool found = (pc != nullptr);
    const char * v = found ? pc : def;
    note(arg, v ? v : "", def ? def : "", found, help, group);
    return v;
  }

  //! Havia -h / --help / -help na linha de comando?
  static bool wants_help()
  {
    return hasSwitch("-h") || hasSwitch("--help") || hasSwitch("-help");
  }

  /** Imprime tudo o que foi lido ate agora, agrupado, marcando o que veio da
      linha de comando e o que ficou no default.

      Chamar no FIM de config(), quando toda a configuracao ja passou por
      aqui: o registro so contem o que ja foi lido ate o momento da chamada. */
  static void dump(std::ostream & os)
  {
    if (registry.empty()) return;

    std::size_t w = 0, wv = 0;
    for (std::map<std::string,Entry>::const_iterator it = registry.begin();
         it != registry.end(); ++it)
    {
      w  = std::max(w,  it->first.size());
      wv = std::max(wv, it->second.value.size());
    }
    if (wv > 34) wv = 34;   // caminhos longos nao devem quebrar o alinhamento

    std::map<std::string, std::vector<std::string> > by_group;

    for (std::map<std::string,Entry>::const_iterator it = registry.begin();
         it != registry.end(); ++it)
    {
      const Entry & e = it->second;
      std::ostringstream line;
      line << "   " << std::left << std::setw((int)w + 2) << it->first
           << std::setw((int)wv + 2) << e.value
           << std::setw(11) << (e.from_cli ? "[cli]" : "[default]");
      if (!e.help.empty()) line << e.help;
      by_group[e.group.empty() ? "outros" : e.group].push_back(line.str());
    }

    os << "\n=== flags desta rodada ===" << std::endl;
    for (std::map<std::string,std::vector<std::string> >::const_iterator
           g = by_group.begin(); g != by_group.end(); ++g)
    {
      os << " [" << g->first << "]" << std::endl;
      for (std::size_t i = 0; i < g->second.size(); i++)
        os << g->second[i] << std::endl;
    }
    os << std::endl;
  }

};

#endif
