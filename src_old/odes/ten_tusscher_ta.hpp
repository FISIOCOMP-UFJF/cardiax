#ifndef TENTUSSCHERIITA_H
#define TENTUSSCHERIITA_H

#include <cmath>
#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

/** Modified ten Tusscher II model to generate an active tension using
    an ordinary differential equation proposed by Nash and Panfilov.

    References:
    - "A model for human ventricular tissue"
      K. H. W. J. ten Tusscher, D. Noble, P. J. Noble and A. V. Panfilov
      Am J Physiol Heart Circ Physiol. 286, H1573-H1589, 2004.

    - "Electromechanical model of excitable tissue to study cardiac arrhythmias"
      Progress in Biophysics & Molecular Biology, 85, (2004), 501-522

    Description of the currents:
    - IK1: inward rectifier potassium current
    - Ito: transient outward current
    - IKr: rapid time dependent potassium current
    - IKs: slow time dependent potassium current
    - ICaL: L-type calcium current
    - INaK: sodium potassium pump current
    - INa: fast sodium current
    - IbNa: sodium background current
    - IbCa: calcium background current
    - IpK: potassium pump current
    - IpCa: calcium pump current
    - Istim: external applied stimulus current    
*/

class TenTusscherTa : public CellModel
{
public:
  
    //! Default constructor
    TenTusscherTa();
  
    //! Set initial conditions for TT2 model
    virtual void init(double * values) const;

    //! Compute RHS equations for the TT2 model
    virtual void equation(const double t, const double * sv, double * values);

    //! Return the value of a variable
    virtual double variable(const int index) const { return -1; }
  

private:
  
    //! Auxiliary function to compute active tension.
    double active_stress_twitch(double v);
  
    double rateswitchNP(double v, double aux, double Vthreshold, double auxThreshold);
  
    double elapsedTime;
  
    bool marked;
  
};

#endif
