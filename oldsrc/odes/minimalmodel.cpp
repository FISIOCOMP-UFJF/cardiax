#include "minimalmodel.hpp"
#include <cmath>
#include <cassert>

MinimalModel::MinimalModel() : CellModel(6) {
    var_names.insert( std::pair<int, std::string>(0, "u") );
    var_names.insert( std::pair<int, std::string>(1, "v") );
    var_names.insert( std::pair<int, std::string>(2, "w") );
    var_names.insert( std::pair<int, std::string>(3, "s") );
    var_names.insert( std::pair<int, std::string>(4, "ca") );
    var_names.insert( std::pair<int, std::string>(5, "ta") );
    
    
    // variables for RL integration in Explicit Euler
    for(int i=1; i<=3; i++) rlvars.insert(i);

}



void MinimalModel::init(double *values) const {
    assert(values != nullptr);
    values[0] = -88.65; //u
    values[1] = 1.0; //v
    values[2] = 1.0; //w
    values[3] = 0.0; //s
    values[4] = 0.0; //ca
    values[5] = 0.0; //ta

}

void MinimalModel::equation(const double t, const double * sv, double * values){
    double u = sv[0];
    const double v = sv[1];
    const double w = sv[2];
    const double s = sv[3];
    const double ca = sv[4];
    const double ta = sv[5];

    const double dt = 0.01; 
    double stim_current = i_stim;

    //u = (u + 88.54) / 115.50;

    const double u_o = 0.000977;
    const double theta_v = 0.277344;
    const double theta_w = 0.222363; 
    const double tau_vplus = 1.451172;
    const double tau_si = 4.643555;
    const double tau_s1 = 1.567383;
    const double k_s = 0.224609;
    const double u_s = 1.380859;

    const double u_u = 1.016406;
    const double theta_vminus = 0.096387;
    const double theta_o = 0.002637;
    const double tau_v1minus = 15.722656;
    const double tau_v2minus = 1060.937500;
    const double tau_w1minus = 91.894531;
    const double tau_w2minus = 4.648438;
    const double k_wminus = 35.546875;
    const double u_wminus = 0.242188;
    const double tau_wplus = 158.496094;
    const double tau_fi = 0.044922;
    const double tau_o1 = 384.960938;
    const double tau_o2 = 1.005859;
    const double tau_so1 = 38.183594;

    const double tau_so2 = 1.175781;
    const double k_so = 3.969727;
    const double u_so = 0.823242;
    const double tau_s2 = 19.912109;
    
    const double tau_winf = 0.067383;
    const double w_ifstar = 0.820312;
    

    
    double H = (u - theta_v > 0) ? 1.0 : 0.0;
    double h_o = (u - theta_o > 0) ? 1.0 : 0.0;
    double h_w = (u - theta_w > 0) ? 1.0 : 0.0;
    double h_v_minus = (u - theta_vminus > 0) ? 1.0 : 0.0;

    double tau_o = (1.0 - h_o) * tau_o1 + h_o * tau_o2;
    double tau_so = tau_so1 + (tau_so2 - tau_so1) * (1.0 + tanh(k_so * (u - u_so))) * 0.5; 
    double tau_vminus = (1.0 - h_v_minus) * tau_v1minus + h_v_minus * tau_v2minus; 

    double J_fi = -v * H * (u - theta_v) * (u_u - u) / tau_fi;
    double J_so = (u - u_o) * (1.0 - h_w) / tau_o + h_w / tau_so;
    double J_si = -h_w * w * s / tau_si;

    values[0] = -(J_fi + J_so + J_si) + stim_current;
    //
    //values[0] = 115.50 * values[0] -88.54;

    double v_inf = (u < theta_vminus) ? 1.0 : 0.0;
    double tau_v_rl = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * H + tau_vminus * H);
    // values[1] = (v_inf -v) / tau_v;
    double v_inf_rl = (tau_vplus * v_inf * (1 - H)) / (tau_vplus - tau_vplus * H + tau_vminus * H);


    if(tau_v_rl > 1e-10) {
        values[1] = v_inf_rl - (v_inf_rl - v) * exp(-dt / tau_v_rl);
    } else {
        values[1] = dt * ((1.0 - H) * (v_inf - v) / tau_vminus - H * v / tau_vplus) + v;
    }

    double w_inf = (1.0 - h_o) * (1.0 - (u / tau_winf)) + h_o * w_ifstar;
    double tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0 + tanh(k_wminus * (u - u_wminus))) * 0.5;
    double tau_w_rl = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * h_w + tau_wminus * h_w);
    double w_inf_rl = (tau_wplus * w_inf * (1 - h_w)) / (tau_wplus - tau_wplus * h_w + tau_wminus * h_w);

    if(tau_w_rl > 1e-10) {
        values[2] = w_inf_rl - (w_inf_rl - w) * exp(-dt / tau_w_rl);
    } else {
        values[2] = dt * ((1.0 - h_w) * (w_inf - w) / tau_wminus - h_w * w / tau_wplus) + w;
    }

    double tau_s = (1.0 - h_w) * tau_s1 + h_w * tau_s2;
    double s_inf_rl = (1.0 + tanh(k_s * (u - u_s))) / 2;

    if(tau_s > 1e-10) {
        values[3] = s_inf_rl - (s_inf_rl - s) * exp(-dt / tau_s);
    } else {
        values[3] = dt *  (((1.0 + tanh(k_s * (u - u_s))) * 0.5 - s) / tau_s) + s;
    }

    // double tau_w = (158.496094 * 91.894531) / (158.496094 - H * (158.496094 - 91.894531));
    //JG

    double gamma0, gamma2, gamma3;
    double p1, p2;
    double a;
    
    gamma0 = 0.15332;
    gamma2 = 7.980469;
    p1 = 3.339844;
    p2 = 0.0;
    gamma3 = 8.099609;
    a = 0.355469;

    double fat_ca = pow(-J_fi, gamma0) / gamma2;
    double fina_ca = ((p1 + p2) / 2.0) - (p1 - p2) / 2.0 * tanh(-gamma3 * (J_si - a));
    double dca_dt = fat_ca - fina_ca*(ca);
    values[4] =  dca_dt ; // Euler

    double gamma7, gamma4, gamma5, gamma33;
    double p11, p22;
    double c1;
    double aa;

    gamma7 = 0.020996;
    gamma4 = 1.386719;
    gamma5 = 0.201172;
    p11 = 0.039453;
    p22 = 0.007715;
    gamma33 = 2.744141;
    aa = 0.054688;
    double fat_ta = gamma7*pow(ca, gamma4);

    if(pow(-J_fi, gamma5) > 0.001){
        c1 = 0;
    } else {
        c1 = 1;
    }

    double fina_ta = ((p11+p22/2.0) + (p11-p22)/2.0*tanh(-gamma33*(ca - aa)));

    double dta_dt = c1*fat_ta - fina_ta*(ta);
    values[5] = dta_dt ; //Euler

//outro teste
    // values[2] = (w_inf - w) / tau_w;

    // double s_inf = (1.0 + tanh(0.224609 * (u - 1.380859))) / 2.0;
    // values[3] = (s_inf - s) / 1.567383;

    // double fat_ca = pow(-J_fi, 0.15332) / 7.980469;
    // double fina_ca = ((p1 + 0.0) / 2.0) - ((p1 - 0.0) / 2.0) * tanh(-8.099609 * (J_si - 0.355469));
    // values[4] = fat_ca - fina_ca * ca;

    // double fat_ta = 0.020996 * pow(ca, 1.386719) ;
    // double fina_ta = ((p1 + p2) / 2.0) * ((p1 - p2) / 2.0) * tanh(-gamma3 * (ca - a));
    // values[5] = fat_ta - fina_ta * ta;


}