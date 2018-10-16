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

    GradientDescent(index_t max_iter, index_t min_iter, index_t max_points, 
                index_t n_gradients=1, index_t seed=1025, 
                value_t stepsize = 0.01, value_t conv_crit=0.0, 
                bool dump_points=false);

    virtual std::unique_ptr<Minimizer> clone() const override {
        return std::make_unique<GradientDescent>(*this);
    }

    std::string get_name();
    
    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    void descent(index_t nDims);
    
#ifdef OMP
    void descent_parallel(index_t nDims);
#endif

    v_d get_gradient(v_d & cube, index_t nDims, value_t llh);
    v_d get_gradient(v_d & cube, index_t nDims, value_t llh, 
        index_t & n_lh_calls);

    /// Transform point from hypercube to physical space
    v_d to_physics(v_d cube, index_t nDims);

private:
    value_t stepsize_;
    value_t conv_crit_;
    index_t min_iter_;
    index_t n_gradients_;
};

#endif
