#ifndef MINIMIZERRESULT_H_INCLUDED
#define MINIMIZERRESULT_H_INCLUDED

#include <vector>
#include <string>
#include <boost/cstdint.hpp>

// A simple object to store the result of the minimization. 
struct MinimizerResult {
    double best_fit;
    std::vector<double> params_best_fit;
    uint32_t n_lh_calls;
    double lh_efficiency;
    std::string minimizer_name;
    std::string function_name;
};

#endif