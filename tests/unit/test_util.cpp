/*
 * Unit tests of src/util: the helpers everything else relies on.
 *
 * What is covered and why:
 *  - CommandLineArgs: each overload of read() returns the CLI value when the flag is
 *    present and the default otherwise, and records the [cli]/[default] provenance
 *    that the "flags desta rodada" table prints. A wrong provenance would lie about
 *    which configuration a run used.
 *  - ExprScalarFunction: expressions of x, y, z and the constants pi and e, used for
 *    source terms and exact solutions in the Poisson solver.
 *  - LinearInterp/BilinearInterp: interpolation at nodes, between nodes, and the
 *    search routines locate() and hunt().
 *  - InputFile: the key/value and section parser of the legacy .par format.
 *  - LoadControl: the incremental load ramp of the mechanics. The load must reach
 *    exactly 1 after ninc increments, and restart from a checkpoint must keep the ramp.
 *  - Parameter/Parameters: typed parameter sets and the XML overrides.
 *  - TimeParameters: step counting and the print rate of the transient solvers. Also
 *    set_stop(), which guards against the EP clock stopping silently mid-run.
 *  - Timer/TimerSection/Log: elapsed time accumulates and logs reach the file.
 *  - String, path and small tensor helpers in util.hpp.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <map>
#include <sstream>
#include <unistd.h>

#include "util/util.hpp"
#include "util/command_line_args.h"
#include "util/expr_function.hpp"
#include "util/load_control.hpp"
#include "util/parameters.hpp"
#include "util/time_parameters.hpp"
#include "util/interp.h"
#include "util/input_file.h"
#include "util/log.hpp"
#include "util/timer.hpp"

using namespace cardiax_test;

// ---- CommandLineArgs --------------------------------------------------------

namespace
{
// Re-points CommandLineArgs to a fake command line and restores the real one, so
// these cases do not leak into later cases.
struct FakeCommandLine
{
  std::vector<const char *> argv;
  explicit FakeCommandLine(std::vector<const char *> a) : argv(std::move(a))
  {
    argv.insert(argv.begin(), "prog");
    argv.push_back(nullptr);
    CommandLineArgs::init(static_cast<int>(argv.size()) - 1, argv.data());
  }
  ~FakeCommandLine()
  {
    static const char * restore[] = {"prog", "-linsolver", "petsc", nullptr};
    CommandLineArgs::init(3, restore);
  }
};
}

TEST_CASE("CommandLineArgs: each type returns the CLI value when the flag is present")
{
  FakeCommandLine cl({"-d", "2.5", "-i", "7", "-s", "hello", "-b", "true", "-b2", "false"});
  CHECK_NEAR(CommandLineArgs::read("-d", 1.0), 2.5, 0.0);
  CHECK_EQ(CommandLineArgs::read("-i", 1), 7);
  CHECK_EQ(CommandLineArgs::read("-s", std::string("x")), std::string("hello"));
  CHECK_EQ(CommandLineArgs::read("-b", false), true);
  CHECK_EQ(CommandLineArgs::read("-b2", true), false);
  CHECK_EQ(std::string(CommandLineArgs::read("-s", "x")), std::string("hello"));
}

TEST_CASE("CommandLineArgs: each type returns the default when the flag is absent")
{
  FakeCommandLine cl({"-other", "1"});
  CHECK_NEAR(CommandLineArgs::read("-d", 1.5), 1.5, 0.0);
  CHECK_EQ(CommandLineArgs::read("-i", 3), 3);
  CHECK_EQ(CommandLineArgs::read("-s", std::string("def")), std::string("def"));
  CHECK_EQ(CommandLineArgs::read("-b", true), true);
  CHECK_EQ(std::string(CommandLineArgs::read("-s", "def")), std::string("def"));
}

TEST_CASE("CommandLineArgs: a flag at the end of the line, without a value, reads as absent")
{
  FakeCommandLine cl({"-a", "1", "-flag"});
  CHECK_EQ(CommandLineArgs::read("-flag", 42), 42);
  CHECK(CommandLineArgs::hasSwitch("-flag"));
  CHECK(!CommandLineArgs::hasSwitch("-nope"));
}

TEST_CASE("CommandLineArgs: wants_help recognises the three spellings")
{
  {
    FakeCommandLine cl({"-h", "x"});
    CHECK(CommandLineArgs::wants_help());
  }
  {
    FakeCommandLine cl({"--help", "x"});
    CHECK(CommandLineArgs::wants_help());
  }
  {
    FakeCommandLine cl({"-a", "1"});
    CHECK(!CommandLineArgs::wants_help());
  }
}

TEST_CASE("CommandLineArgs: dump() reports value and provenance of each flag read")
{
  FakeCommandLine cl({"-from_cli_flag", "9"});
  CommandLineArgs::read("-from_cli_flag", 1, "help text one", "grp_a");
  CommandLineArgs::read("-from_default_flag", 5, "help text two", "grp_a");

  std::ostringstream os;
  CommandLineArgs::dump(os);
  const std::string out = os.str();

  std::istringstream lines(out);
  std::string l, cli_line, def_line;
  while (std::getline(lines, l))
  {
    if (l.find("-from_cli_flag") != std::string::npos) cli_line = l;
    if (l.find("-from_default_flag") != std::string::npos) def_line = l;
  }
  CHECK_MSG(cli_line.find("[cli]") != std::string::npos, cli_line);
  CHECK_MSG(cli_line.find("9") != std::string::npos, cli_line);
  CHECK_MSG(cli_line.find("help text one") != std::string::npos, cli_line);
  CHECK_MSG(def_line.find("[default]") != std::string::npos, def_line);
  CHECK_MSG(def_line.find("5") != std::string::npos, def_line);
  CHECK(out.find("[grp_a]") != std::string::npos);
}

// ---- ExprScalarFunction -----------------------------------------------------

TEST_CASE("ExprScalarFunction: evaluates x, y, z and the constant pi")
{
  ExprScalarFunction f("2*x + y^2 - z");
  CHECK_NEAR(f(arma::vec3({1.0, 2.0, 3.0})), 3.0, 1e-14);

  ExprScalarFunction g("sin(pi*x)");
  CHECK_NEAR(g(arma::vec3({0.5, 0.0, 0.0})), 1.0, 1e-14);
}

TEST_CASE_KNOWN_BUG("ExprScalarFunction: the constant e is available, as the header says",
                    "src/util/expr_function.hpp documents the constant e, but "
                    "symbols.add_constants() only defines pi, epsilon and inf: 'e' throws")
{
  ExprScalarFunction h("e");
  CHECK_NEAR(h(arma::vec3({0.0, 0.0, 0.0})), std::exp(1.0), 1e-14);
}

TEST_CASE("ExprScalarFunction: repeated evaluations do not keep state")
{
  ExprScalarFunction f("x*x");
  CHECK_NEAR(f(arma::vec3({3.0, 0.0, 0.0})), 9.0, 1e-14);
  CHECK_NEAR(f(arma::vec3({-2.0, 0.0, 0.0})), 4.0, 1e-14);
}

TEST_CASE("ExprScalarFunction: an invalid expression throws instead of evaluating to garbage")
{
  CHECK_THROWS(ExprScalarFunction("x +* )"));
}

// ---- Interpolation ----------------------------------------------------------

TEST_CASE("LinearInterp: reproduces nodes and a linear function")
{
  arma::vec x = arma::linspace(0.0, 4.0, 9);
  arma::vec y = 3.0 * x + 1.0;
  LinearInterp li(x, y);

  for (arma::uword i = 0; i < x.n_elem; i++)
    CHECK_NEAR(li.interp(x(i)), y(i), 1e-13);

  CHECK_NEAR(li.interp(1.3), 3.0 * 1.3 + 1.0, 1e-13);
  CHECK_NEAR(li.interp(3.77), 3.0 * 3.77 + 1.0, 1e-13);
}

TEST_CASE("LinearInterp: nonlinear data is piecewise linear between the nodes")
{
  arma::vec x = {0.0, 1.0, 3.0};
  arma::vec y = {0.0, 10.0, 0.0};
  LinearInterp li(x, y);
  CHECK_NEAR(li.interp(0.5), 5.0, 1e-13);
  CHECK_NEAR(li.interp(2.0), 5.0, 1e-13);
}

TEST_CASE("BaseInterp: locate and hunt find the same interval")
{
  arma::vec x = arma::linspace(0.0, 1.0, 101);
  arma::vec y = x;
  LinearInterp a(x, y), b(x, y);
  for (double q : {0.005, 0.5, 0.123, 0.99, 0.1, 0.8})
    CHECK_EQ(a.locate(q), b.hunt(q));
}

TEST_CASE("BilinearInterp: reproduces a bilinear function")
{
  arma::vec x1 = arma::linspace(0.0, 1.0, 5);
  arma::vec x2 = arma::linspace(0.0, 2.0, 7);
  arma::mat z(5, 7);
  auto fn = [](double a, double b) { return 1.0 + 2.0 * a - b + 3.0 * a * b; };
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 7; j++)
      z(i, j) = fn(x1(i), x2(j));

  BilinearInterp bi(x1, x2, z);
  CHECK_NEAR(bi.interp(0.33, 1.27), fn(0.33, 1.27), 1e-12);
  CHECK_NEAR(bi.interp(0.8, 0.1), fn(0.8, 0.1), 1e-12);
}

// ---- InputFile --------------------------------------------------------------

TEST_CASE("InputFile: reads keys of each type and ignores comments")
{
  write_text("input.par",
             "# num_increments 999\n"
             "num_increments 12\n"
             "tolerance 1.5e-3\n"
             "problem_type THREE_DIM\n"
             "material_coef 3 1.0 2.0 3.0\n");
  InputFile f("input.par");
  int i = 0;
  double d = 0.0;
  std::string s;
  std::vector<double> arr;
  f.read("num_increments", i);
  f.read("tolerance", d);
  f.read("problem_type", s);
  f.read_array("material_coef", arr);
  f.close();

  CHECK_EQ(i, 12);
  CHECK_NEAR(d, 1.5e-3, 1e-18);
  CHECK_EQ(s, std::string("THREE_DIM"));
  REQUIRE(arr.size() == 3);
  CHECK_NEAR(arr[2], 3.0, 0.0);
}

TEST_CASE("InputFile: reads each kind of section")
{
  write_text("sections.par",
             "pressure\n11 2.5\n12 3.5\nend_pressure\n"
             "traction\n7 1.0 2.0 3.0\nend_traction\n"
             "prescribed\n4 0 0.1\n5 2 0.2\nend_prescribed\n");
  InputFile f("sections.par");

  std::map<int, double> pressure;
  f.read_section("pressure", pressure);
  CHECK_EQ(pressure.size(), std::size_t(2));
  CHECK_NEAR(pressure[12], 3.5, 0.0);

  std::map<int, arma::vec3> traction;
  f.read_section("traction", traction);
  REQUIRE(traction.count(7) == 1);
  CHECK_NEAR(traction[7](2), 3.0, 0.0);

  std::multimap<int, NodalData> presc;
  f.read_section("prescribed", presc);
  REQUIRE(presc.count(5) == 1);
  CHECK_EQ(presc.find(5)->second.first, 2);       // z
  CHECK_NEAR(presc.find(5)->second.second, 0.2, 0.0);
  f.close();
}

TEST_CASE_KNOWN_BUG("InputFile: a marker -> direction section keeps exactly its entries",
                    "src/util/input_file.cpp, read_section(multimap<int,int>): the end_ line "
                    "goes through the insertion too, adding a spurious entry 0 -> 0")
{
  write_text("dirichlet.par", "dirichlet\n1 2\n1 3\nend_dirichlet\n");
  InputFile f("dirichlet.par");
  std::multimap<int, int> dirichlet;
  f.read_section("dirichlet", dirichlet);
  f.close();
  CHECK_EQ(dirichlet.count(1), std::size_t(2));
  CHECK_EQ(dirichlet.size(), std::size_t(2));
}

// ---- LoadControl ------------------------------------------------------------

TEST_CASE("LoadControl: the ramp reaches exactly 1 after ninc increments")
{
  LoadControl lc(10);
  int n = 0;
  while (lc.has_load())
  {
    lc.update();
    n++;
    REQUIRE(n <= 10);
    CHECK_NEAR(lc.load(), n * 0.1, 1e-12);
  }
  CHECK_EQ(n, 10);
  CHECK_NEAR(lc.load(), 1.0, 1e-12);
  CHECK_EQ(lc.increment(), 10);
}

TEST_CASE("LoadControl: reset() restarts the ramp with the same increments")
{
  LoadControl lc;
  lc.set_nincs(4);
  lc.update();
  lc.update();
  lc.reset();
  CHECK_EQ(lc.increment(), 0);
  CHECK_NEAR(lc.load(), 0.0, 0.0);
  CHECK_NEAR(lc.load_step(), 0.25, 0.0);
  CHECK(lc.has_load());
}

TEST_CASE("LoadControl: the load never goes above 1")
{
  LoadControl lc(3);
  for (int i = 0; i < 6; i++) lc.update();
  CHECK_NEAR(lc.load(), 1.0, 0.0);
  CHECK(!lc.has_load());
}

TEST_CASE("LoadControl: restart from a checkpoint distributes what is left of the ramp")
{
  LoadControl lc(10);
  lc.restart_from_checkpoint(5, 0.5, 10);
  CHECK_EQ(lc.increment(), 5);
  CHECK_NEAR(lc.load(), 0.5, 0.0);
  CHECK_NEAR(lc.load_step(), 0.1, 1e-14);

  // Asking for more increments from the same point makes the steps smaller
  lc.restart_from_checkpoint(5, 0.5, 15);
  CHECK_NEAR(lc.load_step(), 0.05, 1e-14);

  // Nothing left: the load stays
  lc.restart_from_checkpoint(10, 1.0, 10);
  CHECK_NEAR(lc.load_step(), 0.0, 0.0);
  CHECK(!lc.has_load());
}

TEST_CASE("LoadControl: adapt() grows the step for easy increments and shrinks it for hard ones")
{
  LoadControl easy(10), hard(10), middle(10);
  easy.adapt(2);
  hard.adapt(12);
  middle.adapt(7);
  CHECK(easy.load_step() > 0.1);
  CHECK(hard.load_step() < 0.1);
  CHECK_NEAR(middle.load_step(), 0.1, 1e-14);

  LoadControl capped(3);
  for (int i = 0; i < 10; i++) capped.adapt(1);
  CHECK(capped.load_step() <= 0.4);
}

// ---- Parameters -------------------------------------------------------------

TEST_CASE("Parameters: typed add, read and overwrite")
{
  Parameters p("my");
  p.add("tol", 1e-5);
  p.add("maxits", 30);
  p.add("name", "cg");
  p.add("rebuild", false);

  CHECK_EQ(p.name(), std::string("my"));
  CHECK_NEAR(double(p["tol"]), 1e-5, 0.0);
  CHECK_EQ(int(p["maxits"]), 30);
  CHECK_EQ(std::string(p["name"]), std::string("cg"));
  CHECK_EQ(bool(p["rebuild"]), false);

  p["maxits"] = 50;
  CHECK_EQ(int(p["maxits"]), 50);
}

TEST_CASE("Parameters: has_key, remove, clear and key listings")
{
  Parameters p;
  p.add("a", 1);
  p.add("b", 2.0);
  CHECK(p.has_key("a"));
  CHECK(p.has_parameter("b"));
  CHECK(!p.has_key("c"));

  std::vector<std::string> keys;
  p.get_parameter_keys(keys);
  CHECK_EQ(keys.size(), std::size_t(2));

  p.remove("a");
  CHECK(!p.has_key("a"));
  p.clear();
  CHECK(!p.has_key("b"));
}

TEST_CASE("Parameters: nested sets and copies are independent")
{
  Parameters inner("solver");
  inner.add("tol", 1e-3);
  Parameters outer("all");
  outer.add(inner);

  CHECK(outer.has_parameter_set("solver"));
  CHECK_NEAR(double(outer("solver")["tol"]), 1e-3, 0.0);

  Parameters copy(outer);
  copy("solver")["tol"] = 1e-9;
  CHECK_NEAR(double(outer("solver")["tol"]), 1e-3, 0.0);
  CHECK_NEAR(double(copy("solver")["tol"]), 1e-9, 0.0);
}

TEST_CASE("Parameters: parse_xml overrides declared parameters and ignores the others")
{
  write_text("params.xml",
             "<electrophysiology><parameters>"
             "<sigma>0.25</sigma><maxnz>12</maxnz><flag>true</flag><kind>abc</kind>"
             "<undeclared>1</undeclared>"
             "</parameters></electrophysiology>");
  Parameters p;
  p.add("sigma", 1.0);
  p.add("maxnz", 3);
  p.add("flag", false);
  p.add("kind", "zzz");
  p.parse_xml("params.xml", "electrophysiology");

  CHECK_NEAR(double(p["sigma"]), 0.25, 0.0);
  CHECK_EQ(int(p["maxnz"]), 12);
  CHECK_EQ(bool(p["flag"]), true);
  CHECK_EQ(std::string(p["kind"]), std::string("abc"));
  CHECK(!p.has_key("undeclared"));
}

// ---- TimeParameters ---------------------------------------------------------

TEST_CASE("TimeParameters: number of steps, print rate and marching")
{
  TimeParameters t(0.5, 10.0, 2.0);
  CHECK_EQ(t.get_nsteps(), 20);
  CHECK_EQ(t.pr(), 4);
  CHECK_EQ(t.get_size(), 5);
  CHECK(!t.finished());
  CHECK(t.time2print());

  int n = 0;
  while (!t.finished()) { t.increase_time(); n++; REQUIRE(n <= 100); }
  CHECK_EQ(n, 20);
  CHECK_NEAR(t.time(), 10.0, 1e-12);
  t.reset();
  CHECK_EQ(t.it(), 0);
  CHECK_NEAR(t.time(), 0.0, 0.0);
}

TEST_CASE("TimeParameters: set_stop recomputes the step count and the end of the run")
{
  TimeParameters t(1.0, 5.0, 1.0);
  for (int i = 0; i < 5; i++) t.increase_time();
  CHECK(t.finished());

  t.set_stop(8.0);
  CHECK_EQ(t.get_nsteps(), 8);
  CHECK(!t.finished());
  CHECK_NEAR(t.get_stop(), 8.0, 0.0);
}

TEST_CASE("TimeParameters: restore_state continues from the saved step")
{
  TimeParameters t(0.1, 1.0, 0.1);
  t.restore_state(5, 0.5);
  CHECK_EQ(t.it(), 5);
  t.increase_time();
  CHECK_EQ(t.it(), 6);
  CHECK_NEAR(t.time(), 0.6, 1e-12);
}

// ---- Timers and log ---------------------------------------------------------

TEST_CASE("Timer: accumulates time and counts calls")
{
  Timer t;
  t.start();
  usleep(2000);
  t.stop();
  const double first = t.get_total_time();
  CHECK(first > 0.0);
  CHECK_EQ(t.get_num_calls(), 1);

  t.start();
  usleep(2000);
  t.stop();
  CHECK(t.get_total_time() > first);
  CHECK_EQ(t.get_num_calls(), 2);

  t.reset();
  CHECK_NEAR(t.get_total_time(), 0.0, 0.0);
}

TEST_CASE("TimerSection: sections accumulate independently and the summary does not crash")
{
  TimerSection ts;
  ts.enter("a");
  usleep(1000);
  ts.leave();
  ts.enter("b");
  usleep(1000);
  ts.leave();
  CHECK(ts.time("a") > 0.0);
  CHECK(ts.time("b") > 0.0);
  ts.summary();
}

TEST_CASE("Log: everything written reaches the file")
{
  {
    Log log("test.log");
    log << std::string("text ") << 2.5 << std::endl;
    log.write("value=%d\n", 42);
  }
  const std::string out = read_text("test.log");
  CHECK(out.find("text") != std::string::npos);
  CHECK(out.find("2.5") != std::string::npos);
  CHECK(out.find("value=42") != std::string::npos);
}

// ---- util.hpp ---------------------------------------------------------------

TEST_CASE("file helpers: file_exists and file_extension")
{
  write_text("some.file.xml", "x");
  CHECK(file_exists("some.file.xml"));
  CHECK(file_exists(std::string("some.file.xml")));
  CHECK(!file_exists("does_not_exist.xml"));
  CHECK_EQ(file_extension("a/b.c/mesh.xml"), std::string("xml"));
  CHECK_EQ(file_extension("noext"), std::string(""));
}

TEST_CASE("string helpers: split, trim and parsing")
{
  std::vector<std::string> parts;
  string_split("a b  c", " ", &parts);
  REQUIRE(parts.size() == 3);
  CHECK_EQ(parts[2], std::string("c"));

  std::string s = "  \t text here \t ";
  trim_white_spaces(s);
  CHECK_EQ(s, std::string("text here"));

  arma::vec3 v;
  parse_to_vec3("1.5,2.5, 3.5", v);
  CHECK_NEAR(v(0), 1.5, 0.0);
  CHECK_NEAR(v(1), 2.5, 0.0);
  CHECK_NEAR(v(2), 3.5, 0.0);

  std::vector<double> dv;
  parse_to_vector("1, 2.5,3", dv);
  REQUIRE(dv.size() == 3);
  CHECK_NEAR(dv[1], 2.5, 0.0);

  std::vector<int> iv;
  parse_to_vector("4 5 6", iv);
  REQUIRE(iv.size() == 3);
  CHECK_EQ(iv[2], 6);
}

TEST_CASE("print_error: throws carrying the reason")
{
  bool got = false;
  try { print_error("f.cpp", "task", "the reason"); }
  catch (const std::runtime_error & e)
  {
    got = std::string(e.what()).find("the reason") != std::string::npos;
  }
  CHECK(got);
}

TEST_CASE("read_xml_section: reads marker -> value pairs")
{
  write_text("sec.xml",
             "<elasticity><pressure>"
             "<node id=\"1\" marker=\"30\" value=\"5.5\"/>"
             "<node id=\"2\" marker=\"20\" value=\"1.5\"/>"
             "</pressure></elasticity>");
  std::map<int, double> m;
  read_xml_section("sec.xml", "elasticity", "pressure", m);
  CHECK_EQ(m.size(), std::size_t(2));
  CHECK_NEAR(m[30], 5.5, 0.0);
  CHECK_NEAR(m[20], 1.5, 0.0);
}

TEST_CASE("tensor helpers: contract, trace and Voigt vector")
{
  arma::mat A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  arma::mat I = arma::eye(3, 3);
  CHECK_NEAR(trace(A), 15.0, 0.0);
  CHECK_NEAR(contract(A, I), 15.0, 0.0);
  CHECK_NEAR(contract(A, A), arma::accu(A % A), 0.0);

  arma::mat S = {{1, 4, 6}, {4, 2, 5}, {6, 5, 3}};
  arma::mat v(6, 1);
  voigtvec(3, S, v);
  CHECK_NEAR(v(0), 1.0, 0.0);
  CHECK_NEAR(v(1), 2.0, 0.0);
  CHECK_NEAR(v(2), 3.0, 0.0);
  // Off-diagonal terms: all three components must come out, in some order
  CHECK_NEAR(v(3) + v(4) + v(5), 15.0, 0.0);
  CHECK_NEAR(v(3) * v(4) * v(5), 120.0, 0.0);
}
