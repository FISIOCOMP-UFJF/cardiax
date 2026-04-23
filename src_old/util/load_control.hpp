#ifndef LOADCONTROL_HPP
#define LOADCONTROL_HPP

#include <iostream>

class LoadControl
{
public:

  LoadControl();
  LoadControl(int ni);
  void adapt(int nits);
  bool has_load();
  int increment() const { return incs; }
  double load() const { return xlamb; }
  double load_step() const { return dlamb; }
  void reset();
  void set_nincs(int num);
  int get_nincs() { return ninc; }
  void update();

private:

  int step;     // number of steps
  int ninc;     // number of (initial) increments
  int incs;     // number of increments
  double xlamb; // current load
  double xlmax; // final load (default=1)
  double dlamb; // load step

};

#endif