/**
 * @brief
 * Interface of 
 * Accelerated Parameter Estimation with DALEχ
 * by Scott, Daniel F.
 * https://arxiv.org/abs/1705.02007
 * see https://github.com/danielsf/Dalex
 * 
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 */

#ifndef DalexMinimizer_H_INCLUDED
#define DalexMinimizer_H_INCLUDED

#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/nondet_random.hpp>
#include <fstream>
#include <limits>

#include "helper/abbreviations.h"
#include "MinimizerResult.h"
#include "Minimizer.h"

#include "dalex/include/dalex_driver.h"
#include "dalex/include/chisq.h"

class DalexMinimizer : public Minimizer {
public:
    
    DalexMinimizer(int max_iter, int max_points, int seed=1025, 
                bool dump_points=false);

    virtual std::unique_ptr<Minimizer> clone() const override {
        return std::make_unique<DalexMinimizer>(*this);
    }

    std::string get_name();

    void execute(uint32_t nDims);
    
    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    /// Transform point from hypercube to physical space
    v_d to_physics(v_d cube, uint32_t nDims);

private:
    int seed_;
};

#endif
