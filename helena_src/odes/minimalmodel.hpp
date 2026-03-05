#ifndef MINIMAL_MODEL_HPP_
#define MINIMAL_MODEL_HPP_
#define INITIAL_V (-88.65)
#define NEQ 6

#include "cellmodel.hpp"
#include <cassert>
#include <unordered_map>

using namespace std;

class MinimalModel : public CellModel

{
public:
    MinimalModel();
    virtual void init(double *values) const override;
    virtual void equation(const double t, const double * sv, double * values);
    

// private:
//        double stim_current = 0.0;
// //     void RHS(const double *sv, double *rdY, double stim_current, double dt, int type_cell);
//        void equation(const double *sv, double *rdY, double stim_current, double dt, int type_cell) override;
};

#endif
