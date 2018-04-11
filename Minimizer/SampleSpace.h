#ifndef SampleSpace_H_INCLUDED
#define SampleSpace_H_INCLUDED

#include "../helper/abbreviations.h"
#include "MinimizerResult.h"
#include "Minimizer.h"

class SampleSpace : public Minimizer {
public:

    SampleSpace(int max_iter,  int max_points, int seed=1025,
                bool dump_points=false);
    
    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds, 
                             v_d upper_bounds);

    void sample_space(uint32_t nDims);

    /// Transform point from hypercube to physical space
    v_d to_physics(v_d cube, uint32_t nDims);
    
private:

};

#endif
