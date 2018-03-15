#ifndef MINIMIZER_H_INCLUDED
#define MINIMIZER_H_INCLUDED

#include <boost/random.hpp>
#include <boost/nondet_random.hpp>

#include "../helper/abbreviations.h"
#include "MinimizerResult.h"
#include "../likelihood/TestFunctions.h"

class Minimizer {
public:
    /// constructor for service factory
    Minimizer(double tolerance, uint32_t max_iter, uint32_t min_iter, 
              uint32_t max_points, uint32_t seed=1025, bool dump_points=false) {
        precision_criterion_ = tolerance;
        max_iter_ = max_iter;
        min_iter_ = min_iter;
        max_points_ = max_points;
        intgen.seed(seed);
        result.n_lh_calls = 0;
        result.lh_efficiency = 0;
        dump_points_ = dump_points;
    };

    uint32_t get_n_lh_calls() {return result.n_lh_calls;};
    double get_lh_efficiency() {return result.lh_efficiency;};
    void reset_calls(){result.n_lh_calls=0; result.lh_efficiency=0;};
    void set_bounds(v_d upper, v_d lower) {upper_bnds=upper; lower_bnds=lower;};
    
    void set_output(std::string path) {base_dir_ = path;};
    
    /// core method: minimizer a given function with given initial conditions
    virtual MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds, 
                             v_d upper_bounds) = 0;

    /** Function that evaluates the llh.
     *
     *  \param theta            Physical parameters of point that
     *                          shall be evaluated.
     *
     *  \return                 Log-likelihood
     * */
    double get_llh(v_d theta) {result.n_lh_calls++; 
        return test_func_->get_lh(theta);};
    /// Transform point from hypercube to physical space
    virtual v_d to_physics(v_d cube, uint32_t nDims) = 0;
    boost::mt19937 intgen;
    
    v_d upper_bnds, lower_bnds;
    v_d params_best_fit;
    double lh_bestFit_, lh_worstFit_;
    double precision_criterion_;
    uint32_t min_iter_, max_iter_, max_points_;
    
    TestFunctions *test_func_;
    MinimizerResult result;
    std::string base_dir_, file_name_;
    bool dump_points_;
    v_d results;
    
private:
    

};

#endif