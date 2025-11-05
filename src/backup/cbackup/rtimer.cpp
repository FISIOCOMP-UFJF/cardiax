/**
 * @file   rtimer.cpp
 * @author ronanrmo
 * 
 * @brief  rtimer.cpp
 * 
 * 
 */

#include <sys/time.h>
#include <cstdlib>
#include "rtimer.h"

#ifdef GPU_OPENGL
#include <GL/gl.h>
#endif

RTimer::RTimer()
{
  total_time = 0.0;
  num_calls = 0;
}

void RTimer::start()
{
  num_calls++;

#ifdef GPU_OPENGL
  glFinish();
#endif

  gettimeofday(&start_time, NULL);
}

void RTimer::zero()
{
  total_time = 0.0;
  num_calls = 0;
}

void RTimer::stop()
{
#ifdef GPU_OPENGL
  glFinish();
#endif
  gettimeofday(&stop_time, NULL);
  total_time += stop_time.tv_sec-start_time.tv_sec + (stop_time.tv_usec - start_time.tv_usec)/1000000.0;
  start_time = stop_time;
}

double RTimer::get_total_time()
{
  return total_time;
}

int RTimer::get_num_calls()
{
  return num_calls;
}
