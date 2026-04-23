#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "../../util/datatype.hpp"
#include "timestepper.h"

BOOST_AUTO_TEST_SUITE(test_suite_timestepper)

BOOST_AUTO_TEST_CASE(constructors_test)
{

  // Test 1
  TimeStepper ts0;
  BOOST_CHECK_EQUAL(ts0.number_of_time_steps(), (int)0);

  // Test 2 
  // timestep 0.1
  // final time 1.0
  // there should be 10 steps
  TimeStepper ts1(0.1,1.0);
  BOOST_CHECK_EQUAL(ts1.number_of_time_steps(), (int)10);  

  // Test 3
  // invalid 
  // nsteps should be 0
  TimeStepper ts2(0.1,-10.0);
  BOOST_CHECK_EQUAL(ts2.number_of_time_steps(), (int)0);

}

BOOST_AUTO_TEST_CASE(functions_test)
{
  
  // Test 1
  TimeStepper ts0(0.1,1.0);
  ts0.increase_time();
  BOOST_CHECK_EQUAL(ts0.time(), (Real)0.1);

}

BOOST_AUTO_TEST_SUITE_END()
