#ifndef NEUTRINOFUNCTION_H_INCLUDED
#define NEUTRINOFUNCTION_H_INCLUDED

#include <boost/cstdint.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <limits>
#include "helper/abbreviations.h"

class NeutrinoFunction {
public:
    NeutrinoFunction();

    NeutrinoFunction(std::string, uint32_t);

    /// function pointer.
    lh_pointer lh_p;
    std::string get_name();
    uint32_t get_ndims();
    /// Get the (log)likelihood
    double get_lh(v_d theta);

    // Called within get_lh
    EventHypothesis get_hypothesis_ptr(v_d &theta);
    VectorParticlePtr extract_hypothesis(const EventHypothesis &h);
    get_response_matrix();
    fit_statistics(); // <- Millipede Solver (Simplex?)

    // Called within get_response_matrix
    MillipedeAddOMSourcePairToMatrix(...);
    // Called within MillipedeAddOMSourcePairToMatrix
    get_probability_quantiles(...);

    // Called within get_hypothesis_ptr
    update_physics_variables(); // From I3SimpleParametrization

private:
    std::string name;
    uint32_t ndims_;

};

#endif
