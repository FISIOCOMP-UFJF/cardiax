//
// Created by rocha on 1/28/16.
//

#include <iostream>
#include "../parameter.hpp"
#include "../parameters.hpp"

using namespace std;

int main(int argc, char** argv)
{
  cout << "\nTesting Parameter class\n\n";

  // --------------------------------------------------------------------------
  //
  // TEST 1: create one Parameter of each type...int double etc..
  //
  // --------------------------------------------------------------------------

  IntParameter ip("param_inteiro");
  if(ip.is_set())
    cout << "Integer parameter is set." << endl;
  else
    cout << "Integer parameter is NOT set." << endl;

  ip = 10;

  if(ip.is_set())
    cout << "Integer parameter is set." << endl;
  else
    cout << "Integer parameter is NOT set." << endl;
  cout << "Type: " << ip.type_str() << endl;
  cout << "Value: " << ip.value_str() << endl;
  cout << "String: " << ip.str() << endl;
  cout << endl;


  DoubleParameter dp("param_real");
  if(dp.is_set())
    cout << "Double parameter is set." << endl;
  else
    cout << "Double parameter is NOT set." << endl;
  dp = 3.141516;
  if(dp.is_set())
    cout << "Double parameter is set." << endl;
  else
    cout << "Double parameter is NOT set." << endl;
  cout << "Type: " << dp.type_str() << endl;
  cout << "Value: " << dp.value_str() << endl;
  cout << "String: " << dp.str() << endl;
  double dvar = dp;
  cout << "Variable from cast: " << dvar << endl;
  cout << endl;


  BoolParameter bp("param_booleano");
  if(bp.is_set())
    cout << "Bool parameter is set." << endl;
  else
    cout << "Bool parameter is NOT set." << endl;
  bp = true;
  if(bp.is_set())
    cout << "Bool parameter is set." << endl;
  else
    cout << "Bool parameter is NOT set." << endl;
  cout << "Type: " << bp.type_str() << endl;
  cout << "Value: " << bp.value_str() << endl;
  cout << "String: " << bp.str() << endl;
  cout << endl;


  StringParameter sp("param_string");
  if(sp.is_set())
    cout << "String parameter is set." << endl;
  else
    cout << "String parameter is NOT set." << endl;
  sp = "some test string";
  if(sp.is_set())
    cout << "String parameter is set." << endl;
  else
    cout << "String parameter is NOT set." << endl;
  cout << "Type: " << sp.type_str() << endl;
  cout << "Value: " << sp.value_str() << endl;
  cout << "String: " << sp.str() << endl;
  cout << endl;

  cout << "All tests with Parameter class were done!" << endl;

  // --------------------------------------------------------------------------
  //
  // TEST 2: create one Parameters for testing
  //
  // --------------------------------------------------------------------------

  cout << "\nTesting Parameters class\n\n";

  cout << "Creating my_parameters" << endl;
  Parameters p("my_parameters");
  p.add("cg_relative_tolerance",  1e-15);
  p.add("cg_absolute_tolerance",  1e-15);
  p.add("cg_maxits",                 30);
  p.add("cg_monitor_convergence", true);

  cout << p.str() << endl;

  cout << "Test changing and retrieving" << endl;
  p["cg_maxits"] = 50;
  int its = p["cg_maxits"];
  p["cg_monitor_convergence"] = false;
  bool monitor = p["cg_monitor_convergence"];
  p["cg_absolute_tolerance"] = 12.345e-5;
  double dvalue = p["cg_absolute_tolerance"];
  cout << "Parameter 'cg_maxits'=" << its << endl;
  cout << "Parameter 'cg_monitor_convergence'=" << monitor << endl;
  cout << "Parameter 'cg_absolute_tolerance'=" << dvalue << endl;
  cout << endl;

  // Test nested parameter set
  Parameters pn("my_nested_parameters");
  pn.add("ep_teste", 123);
  pn.add("ep_rebuild", true);
  p.add(pn);
  cout << p.str() << endl;

  int teste = p("my_nested_parameters")["ep_teste"];
  cout << "Get nested value: " << teste << endl;

  cout << "Changing nested value to 987" << endl;
  p("my_nested_parameters")["ep_teste"] = 987;
  cout << "Get nested value again: " << (int) p("my_nested_parameters")["ep_teste"] << endl;

  // Test adding another parameter with the same key
  //p.add("cg_maxits",100);

  return 0;
}

