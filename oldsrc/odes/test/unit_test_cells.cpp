/*

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <typeinfo>
#include "../../util/datatype.hpp"
#include "cells.h"
#include "cellmodel.h"
#include "tt2.h"
#include "fhn.h"
#include "lr1.h"
#include "fe.h"
#include "rk4.h"

BOOST_AUTO_TEST_SUITE(test_suite_cells)

BOOST_AUTO_TEST_CASE(create_test_cells)
{
  cout << "\n*** Running Cells unit test " << endl;

  // Test 1
  // 10 cells * 19 sv
  CellModel *cell1 = new TenTusscher();
  cell1->solver("ForwardEuler");
  cell1->setup(0.001,60.0);
  Cells *c1 = new Cells(10,cell1);
  BOOST_CHECK_EQUAL(c1->get_size(),(int)190);
  BOOST_CHECK_EQUAL(c1->get_time(),(Real)0.0);

  // Test 2
  // 1024 cells * 2 sv
  CellModel *cell2 = new FitzHughNagumo();
  cell2->solver("ForwardEuler");
  cell2->setup(0.1,20.0);
  Cells *c2 = new Cells(1024,cell2);
  BOOST_CHECK_EQUAL(c2->get_size(),(int)2048);
  BOOST_CHECK_EQUAL(c2->get_time(),(Real)0.0);

  // Test 3
  // 100 cells * 8 sv
  CellModel *cell3 = new LuoRudy();
  cell3->solver("ForwardEuler");
  cell3->setup(0.01,600.0);
  Cells *c3 = new Cells(100,cell3);
  BOOST_CHECK_EQUAL(c3->get_size(),800);
  BOOST_CHECK_EQUAL(c3->get_time(),(Real)0.0);

}

BOOST_AUTO_TEST_SUITE_END()

*/
