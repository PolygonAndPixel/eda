#ifndef MINIMIZER_H_INCLUDED
#define MINIMIZER_H_INCLUDED

#include <boost/random.hpp>
#include <boost/nondet_random.hpp>

#include "MinimizerResult.h"
#include "../likelihood/TestFunctions.h"

// typedef std::vector<double> v_d;
// typedef std::vector<uint32_t>  v_i;
// typedef std::vector<v_d>  m_d;

class Minimizer {
public:
    /// constructor for service factory
    Minimizer(double tolerance, uint32_t max_iter, uint32_t min_iter, 
              uint32_t max_points, uint32_t seed=1025) {
        precision_criterion_ = tolerance;
        max_iter_ = max_iter;
        min_iter_ = min_iter;
        max_points_ = max_points;
        intgen.seed(seed);
        boost::uniform_real<> uni_dist(0,1);
        boost::variate_generator<boost::mt19937&, 
            boost::uniform_real<> > uf(intgen, uni_dist);
        n_lh_calls = 0;
        n_accepted = 0;
    };

    uint32_t get_n_lh_calls() {return n_lh_calls;};
    double get_lh_efficiency() {return n_accepted/n_lh_calls;};
    void reset_calls(){n_lh_calls=0; n_accepted=0;};
    void set_bounds(v_d upper, v_d lower) {upper_bnds=upper; lower_bnds=lower;};
    
    /// core method: minimizer a given function with given initial conditions
    virtual MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds, 
                             v_d upper_bounds) = 0;

    /// Get the likelihood of a current point
    virtual double get_llh(v_d theta) = 0;
    /// Transform point from hypercube to physical space
    virtual v_d to_physics(double *cube, uint32_t nDims) = 0;
    boost::mt19937 intgen;
    
    v_d upper_bnds, lower_bnds;
    v_d params_best_fit;
    double lh_bestFit_, lh_worstFit_;
    double precision_criterion_;
    uint32_t min_iter_, max_iter_, max_points_;
    
    uint32_t n_lh_calls, n_accepted;
    TestFunctions *test_func_;
    
private:

    

};

#endif
