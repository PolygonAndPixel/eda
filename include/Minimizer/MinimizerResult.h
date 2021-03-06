#ifndef MINIMIZERRESULT_H_INCLUDED
#define MINIMIZERRESULT_H_INCLUDED

#include <vector>
#include <string>

#include "helper/abbreviations.h"

// A simple object to store the result of the minimization.
struct MinimizerResult {
    double best_fit;
    v_d params_best_fit;
    uint32_t n_lh_calls;
    double lh_efficiency;
    std::string minimizer_name;
    std::string function_name;
    double run_time; // Not used yet. Could store the execution time in s.
    int runs;        // Not used yet. Could store the amount of runs for benchmarks.
};

#endif
