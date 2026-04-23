#ifndef STIMULUS_H
#define STIMULUS_H

#include <set>
#include <map>
#include <vector>
#include <armadillo>
#include "mesh/mesh.hpp"

/** Stimulus file format:
      num_stimuli
      it val x0 x1 y0 y1 z0 z1
      it val x0 x1 y0 y1 z0 z1
      ...
      it val x0 x1 y0 y1 z0 z1

    Note: use it ordered in time
*/
struct Stimulus
{
  int it;
  double time;
  double start;
  double duration;
  double value;
  double x0, x1, y0, y1, z0, z1;
  
  // Constructor
  Stimulus(double st, double dur, double v, double xb, double xe,
           double yb, double ye, double zb, double ze)
      : it(0), start(st), duration(dur), value(v),
        x0(xb), x1(xe), y0(yb), y1(ye), z0(zb), z1(ze) {}
  
  // Destructor
  ~Stimulus(){}
  
};

/**
 *  Class to hold a set Stimulus
 *  File format for stimulus:
 *  start duration value x0 x1 y0 y1 z0 z1
 **/
class Stimuli
{
public:
  
  //! Default constructor
  Stimuli();

  //! Default destructor which calls clear()
  ~Stimuli();

  //! Add another stimulus to a specific node (e.g. an external applied current)
  void add(double start, double duration, double sval, double x0, double x1,
           double y0, double y1, double z0, double z1);

  //! Free Stimulus objects from vector of stimuli (svec) and then clear svec.
  void clear();

  //! Check and call apply in order to apply stimulus for Monodomain model.
  void check(const double time, const Mesh &msh, std::set<uint> & snodes,
             double * sval, bool * apply);

  //! Check if it has stimulus to apply and return a bool
  //! indicating whether we apply or not the stimulus,
  //! a std::set with the indexes of nodes to apply and the value
  //! of the stimulus.
  void check_bidomain(const double time, const Mesh &msh,
                      std::map<uint,uint> & tmap,
                      std::set<uint> & snodes, double * sval, bool * apply);

  //! Read data from .stim file
  void read(const std::string & filename);

  //! Read data from .xml file
  void read_xml(const std::string & filename);

private:

  //! Number of stimulus to be applied
  int num_stimuli;

  //! Number of stimulus already applied
  int num_applied;

  //! Vector of Stimulus pointers
  std::vector<Stimulus*> svec;

};

// ------------------- inline and template functions --------------

inline Stimuli::Stimuli() : num_stimuli(0), num_applied(0)
{
}

inline Stimuli::~Stimuli()
{
  clear();
}

#endif
