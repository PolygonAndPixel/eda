#ifndef Scan_H_INCLUDED
#define Scan_H_INCLUDED

#include "helper/abbreviations.h"
#include "MinimizerResult.h"
#include "Minimizer.h"

class Scan : public Minimizer {
public:

    Scan(int n_points_per_dim, int max_points, int seed=1025, bool dump_points=false);

    virtual std::unique_ptr<Minimizer> clone() const override {
        return std::make_unique<Scan>(*this);
    }

    std::string get_name();
    
    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    void scan_space(uint32_t nDims);

    /// Transform point from hypercube to physical space
    v_d to_physics(v_d cube, uint32_t nDims);

private:
    int n_points;
};

#endif