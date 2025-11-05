#include "mynash.hpp"
#include "cellmodel.hpp"
#include "ode_solver.hpp"

MyNashPanfilov::MyNashPanfilov() : CellModel(4)
{
  var_names.insert( std::pair<int, std::string>(0, "V") );
  var_names.insert( std::pair<int, std::string>(1, "r") );
  var_names.insert( std::pair<int, std::string>(2, "Ta") );
  var_names.insert( std::pair<int, std::string>(3, "aux2") );

  rlvars.insert(4);
}

void MyNashPanfilov::init(double * values) const
{
  if(values != 0)
  {
    values[0] = 0.0; // -80 ... 0.0
    values[1] = 0.0; // 0.1;
    values[2] = 0.0; // 0.0;
    values[3] = 0.0; // 0.0;
  }
}

void MyNashPanfilov::equation(const double time, const double *rY, double *rDY)
{
  // State variables
  const double V  = rY[0];
  const double r  = rY[1];
  const double Ta = rY[2];
  const double aux2= rY[3];
  
  // Parameters
  const double k   = 1.0;   //8.0;     
  const double a   = 0.075; //0.15
  const double eps = 0.002; //0.01; 
  const double mu1 = 0.2;   //0.12;   
  const double mu2 = 0.3;   
  const double kTa = 50.0;
  const double err = 5e-3;
  
  const double Vthreshold   = 0.05;
  const double auxThreshold = 0.1;
  
  double Istim = i_stim;
  
  // Calculations
  double d_dt_v	= -k*V*(V-a)*(V-1.0) - r*V + Istim;
  double d_dt_r	= (eps + (mu1 * r)/(mu2 + V)) * (-r - k*V*(V-a-1.0));  
  //double myCond   = rateswitchNP(V, aux2, Vthreshold, auxThreshold);
  double myCond = rateswitchNP(V, aux2, Vthreshold, auxThreshold, Ta/(kTa*0.06));
  double d_dt_Ta  = myCond * (kTa*V - Ta);  
  //double d_dt_aux = (V < Vthreshold) ? -aux2*0.1 : err;
  double d_dt_aux=  V<Vthreshold?-aux2*0.1:err;  
  
  rDY[0] = d_dt_v;
  rDY[1] = d_dt_r;
  rDY[2] = d_dt_Ta;
  rDY[3] = d_dt_aux;
}

double MyNashPanfilov::rateswitchNP(double V, double aux, double Vthreshold, double auxThreshold, double Ta)
{
    if(aux<auxThreshold && V>=0.01)
      return 0.0;
    
    if (V < Vthreshold){
      return /*0.1*//*0.01*/0.06/(Ta+1);
    }else{
      return /*0.01*//*0.02*/0.01;
    }
}

double E(double v)
{
  return 100.0*v - 80.0; 
}

