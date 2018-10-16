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
    
    DalexMinimizer(index_t max_iter, index_t max_points, index_t seed=1025, 
                bool dump_points=false);

    virtual std::unique_ptr<Minimizer> clone() const override {
        return std::make_unique<DalexMinimizer>(*this);
    }

    std::string get_name();

    void execute(index_t nDims);
    
    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    /// Transform point from hypercube to physical space
    v_d to_physics(v_d cube, index_t nDims);

private:
    index_t seed_;
};

#endif
