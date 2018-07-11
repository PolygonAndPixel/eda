#ifndef GradientDescent_H_INCLUDED
#define GradientDescent_H_INCLUDED

#ifdef OMP
#include <omp.h>
#endif

#include "helper/abbreviations.h"
#include "MinimizerResult.h"
#include "Minimizer.h"

class GradientDescent : public Minimizer {
public:

    GradientDescent(int max_iter, int min_iter, int max_points, 
                int n_gradients=1, int seed=1025, 
                double stepsize = 0.01, double conv_crit=0.0, 
                bool dump_points=false);

    virtual std::unique_ptr<Minimizer> clone() const override {
        return std::make_unique<GradientDescent>(*this);
    }

    std::string get_name();
    
    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    void descent(uint32_t nDims);
    
#ifdef OMP
    void descent_parallel(uint32_t nDims);
#endif

    v_d get_gradient(v_d & cube, uint32_t nDims, double llh);
    v_d get_gradient(v_d & cube, uint32_t nDims, double llh, 
        uint32_t & n_lh_calls);

    /// Transform point from hypercube to physical space
    v_d to_physics(v_d cube, uint32_t nDims);

private:
    double stepsize_;
    double conv_crit_;
    int min_iter_;
    int n_gradients_;
};

#endif
