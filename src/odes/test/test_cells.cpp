/*

#include <iostream>
#include <cstdlib>
#include <fstream>
#include "datatype.h"
#include "cells.h"
#include "cellmodel.h"
#include "ode_solver.h"
#include "lr1.h"
#include "tt2.h"

using namespace std;

int main(int argc, char **argv)
{
  CellModel *p = new LuoRudy(); //TenTusscher();
  p->solver("ForwardEuler");
  p->setup(0.001, 10.0);

  Cells *c = new Cells(256,p);
  c->solve(); 

  delete p;
  delete c;
}

*/
