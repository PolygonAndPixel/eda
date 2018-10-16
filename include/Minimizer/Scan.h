#ifndef Scan_H_INCLUDED
#define Scan_H_INCLUDED

#include "helper/abbreviations.h"
#include "MinimizerResult.h"
#include "Minimizer.h"

class Scan : public Minimizer {
public:

    Scan(index_t n_points_per_dim, index_t max_points, index_t seed=1025, bool dump_points=false);

    virtual std::unique_ptr<Minimizer> clone() const override {
        return std::make_unique<Scan>(*this);
    }

    std::string get_name();
    
    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    void scan_space(index_t nDims);
    void scan_space_tmp(index_t nDims);

    /// Transform point from hypercube to physical space
    v_d to_physics(v_d cube, index_t nDims);

private:
    index_t n_points;
};

#endif