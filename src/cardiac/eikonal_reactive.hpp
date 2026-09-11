#ifndef EIKONAL_REACTIVE_HPP
#define EIKONAL_REACTIVE_HPP

#include "cardiacproblem.hpp" 
#include "pdes/eikonal.hpp" 
#include <vector>
#include <string>
#include <armadillo>
#include "fem/fem.h"


class EikonalReactive : public CardiacProblem {
private:
    Eikonal eikonal_solver; 
    arma::Col<double> lat;  

    // Variáveis que faltaram ser migradas do Eikonal antigo:
    uint ndofs;
    arma::Col<double> stim_values;
    std::set<unsigned int> stim_nodes; 
    double stim_val;
    bool stim_apply;
    bool stim_apply_nodes;
    CondTensorType condtype;
    H1FESpace fespace;

    void solve_odes();

public:
    EikonalReactive();
    ~EikonalReactive();

    void setup(std::string & b, std::string & c, std::string & m,
               double dt, double T, double pr, double pa);
    void set_conductivity(int cond);
    
    void init();
    void initial_conditions();
    void advance();
    
    void set_stimulus_value(int index, double val);
    void apply_lat_stimulus(double amplitude, double duration, double period);
    
    void read_eikonal_solution(const std::string &mshfile);

    void solve() override;

    const arma::vec& get_lat() const { return lat; }
};

#endif