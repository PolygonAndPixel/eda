#ifndef SampleSpace_H_INCLUDED
#define SampleSpace_H_INCLUDED

#include "MinimizerResult.h"
#include "Minimizer.h"

class SampleSpace : Minimizer {
public:
    /// constructor for service factory
    SampleSpace(uint32_t max_iter,  uint32_t max_points, uint32_t seed=1025);
    
//     virtual ~SampleSpace();

    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds, 
                             v_d upper_bounds);

    void sample_space(uint32_t nDims);

    /// Get the likelihood of a current point
    double get_llh(v_d theta);
    /// Transform point from hypercube to physical space
    v_d to_physics(double *cube, uint32_t nDims);
    
    void set_output(std::string path);

private:
    std::string base_dir_, file_name_;
    v_d results;
};

#endif
