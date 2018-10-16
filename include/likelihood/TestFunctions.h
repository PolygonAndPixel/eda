#ifndef TESTFUNCTIONS_H_INCLUDED
#define TESTFUNCTIONS_H_INCLUDED

#include <boost/cstdint.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <limits>
#include "helper/abbreviations.h"

class TestFunctions {
    typedef value_t(TestFunctions::*lh_pointer)(v_d & theta);
public:
    TestFunctions();
    /// Set the function and the number of dimensions for these functions.
    /// Last three parameters are the number of detectors for ICECUBE
    TestFunctions(std::string func_name, index_t ndims,
        index_t n_x=3, index_t n_y=3, index_t n_z=5);

    /// function pointer.
    lh_pointer lh_p;
    std::string get_name();
    index_t get_ndims();
    /// Generic likelihood function which calls the actual function which has
    /// been choosen during creation of this class.
    value_t get_lh(v_d & theta);
    value_t get_neg_lh(v_d & theta);
    value_t get_neg_llh(v_d & theta);
    /// The likelihood functions that take one vector of doubles and calculate
    /// the likelihood.
    value_t eggholder(v_d & theta);
    value_t paraboloid(v_d & theta);
    value_t townsend(v_d & theta);
    value_t rosenbrock(v_d & theta);
    value_t himmelblau(v_d & theta);
    value_t gauss_shell(v_d & theta);
    value_t icecube(v_d & theta);

    void set_func(std::string func_name, index_t ndims,
        index_t n_x=3, index_t n_y=3, index_t n_z=5);

private:
    std::string name;
    index_t ndims_;
    // Some variables for the gaussian shell
    value_t shell_width, r; // width and radius of the shells
};

#endif
