#include <iostream>
#include <cmath>
#include "load_control.hpp"

using namespace std;

LoadControl::LoadControl()
  : step(0), ninc(5), incs(0), xlamb(0.0), xlmax(1.0), dlamb(0.2)
{
  // do nothing
}

LoadControl::LoadControl(int ni)
  : step(0), ninc(ni), incs(0), xlamb(0.0), xlmax(1.0), dlamb(1.0/ni)
{
  // do nothing
}

void LoadControl::adapt(int nits)
{ 
  if(has_load())
  {
    // if num newton its smaller than 6 increase load
    // and decrease if it is larger than 9
    if(nits < 6)
    { 
      dlamb = dlamb * pow(10, 1./5.);
      cout << "Adapting load step: increasing load " << xlamb+dlamb << endl;
    }
    else if(nits >= 9)
    {
      dlamb = dlamb * pow(10, -1./5.);
      cout << "Adapting load step: decreasing load " << xlamb+dlamb << endl;
    }

    // keep dlamb limited
    if (dlamb > 0.4) dlamb = 0.4;
  }
}

bool LoadControl::has_load()
{
  return (xlamb <= xlmax) && (incs < ninc);
  
  // adaptive
  //return (xlamb < xlmax);
}

void LoadControl::reset()
{
  incs  = 0;
  xlamb = 0.0;
  xlmax = 1.0;
  dlamb = 1.0/ninc;
  //cout << "Ta factor: " << ta_factor << endl;
}


void LoadControl::set_nincs(int num)
{
  ninc  = num;
  xlmax = 1.0;
  dlamb = 1.0/ninc;
}

void LoadControl::update()
{
  // increment load
  incs += 1;
  xlamb = xlamb + dlamb;
  
  if(xlamb > 1.0) xlamb = 1.0;
  //if(xlamb > xlmax) xlamb = xlmax;
  
  std::cout << " Load increment " << incs << std::endl;
}

