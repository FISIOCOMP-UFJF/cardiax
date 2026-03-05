#include "interp.h"

int BaseInterp::locate(const double x)
{
  int ju,jm,jl;
  if (n < 2 || mm < 2 || mm > n) throw("locate size error");
  bool ascnd=(xx[n-1] >= xx[0]);
  jl=0;
  ju=n-1;
  while(ju-jl > 1){
    jm = (ju+jl) >> 1; // midpoint
    if((x >= xx[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  cor = fabs(jl-jsav) > dj ? 0 : 1;
  jsav = jl;
  return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
}

int BaseInterp::hunt(const double x)
{
  int jl=jsav, jm, ju, inc=1;
  if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
  bool ascnd = (xx[n-1] >= xx[0]);
  if (jl < 0 || jl > n-1){ // input not useful. go to bisection
    jl=0;
    ju=n-1;
  } else {
    if ((x >= xx[jl]) == ascnd) { // hunt up
      for(;;){
        ju = jl + inc;
        if (ju >= n-1) { ju = n-1; break; }    // off end of table
        else if ((x < xx[ju]) == ascnd) break; // found bracket
        else{
          jl = ju;
          inc += inc;
        }
      }
    } else { // hunt down
      ju = jl;
      for(;;) {
        jl = jl - inc;
        if (jl <= 0) { jl = 0; break; }
        else if ((x >= xx[jl]) == ascnd) break; // found bracket
        else {
          ju = jl;
          inc += inc;
        }
      }
    }
  }
  while(ju-jl > 1) { //hunt is done, so begin the final bisection
     jm = (ju+jl) >> 1;
     if ((x >= xx[jm]) == ascnd)
       jl=jm;
     else
       ju=jm;
  }
  // decide whether to use hunt or locate next time
  cor = fabs(jl-jsav) > dj ? 0 : 1;
  jsav = jl;
  return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
}
