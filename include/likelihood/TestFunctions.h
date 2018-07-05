#ifndef TESTFUNCTIONS_H_INCLUDED
#define TESTFUNCTIONS_H_INCLUDED

#include <boost/cstdint.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <limits>
#include "helper/abbreviations.h"
#include "IceCubeToy/DOM.h"
#include "IceCubeToy/ESource.h"
#include "IceCubeToy/IceCube_helper.h"
#include "IceCubeToy/PulseMap.h"
#include "IceCubeToy/Track.h"

class TestFunctions {
    typedef double(TestFunctions::*lh_pointer)(v_d & theta);
public:
    TestFunctions();
    /// Set the function and the number of dimensions for these functions.
    /// Last three parameters are the number of detectors for ICECUBE
    TestFunctions(std::string func_name, uint32_t ndims,
        uint32_t n_x=3, uint32_t n_y=3, uint32_t n_z=5);

    /// function pointer.
    lh_pointer lh_p;
    std::string get_name();
    uint32_t get_ndims();
    /// Generic likelihood function which calls the actual function which has
    /// been choosen during creation of this class.
    double get_lh(v_d & theta);
    double get_neg_lh(v_d & theta);
    double get_neg_llh(v_d & theta);
    /// The likelihood functions that take one vector of doubles and calculate
    /// the likelihood.
    double eggholder(v_d & theta);
    double paraboloid(v_d & theta);
    double townsend(v_d & theta);
    double rosenbrock(v_d & theta);
    double himmelblau(v_d & theta);
    double gauss_shell(v_d & theta);
    double icecube(v_d & theta);

    void set_func(std::string func_name, uint32_t ndims,
        uint32_t n_x=3, uint32_t n_y=3, uint32_t n_z=5);

private:
    std::string name;
    uint32_t ndims_;
    // Some variables for the gaussian shell
    double shell_width, r; // width and radius of the shells
    // Some variables for IceCube
    PulseMap pulse_map;
    double length;
    double seg_length;
};

#endif
