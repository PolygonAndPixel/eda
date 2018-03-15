#ifndef TESTFUNCTIONS_H_INCLUDED
#define TESTFUNCTIONS_H_INCLUDED

#include <boost/cstdint.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <limits>
#include "../helper/abbreviations.h"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

typedef std::vector<double> v_d;

class TestFunctions {
    typedef double(TestFunctions::*lh_pointer)(v_d theta);
public:
    TestFunctions();
    /// Set the function and the number of dimensions for these functions.
    TestFunctions(std::string, uint32_t);
    
    /// function pointer.
    lh_pointer lh_p;
    std::string get_name();
    uint32_t get_ndims();
    /// Generic likelihood function which calls the actual function which has 
    /// been choosen during creation of this class.
    double get_lh(v_d theta);
    /// The likelihood functions that take one vector of doubles and calculate
    /// the likelihood.
    double eggholder(v_d theta);
    double townsend(v_d theta);
    double rosenbrock(v_d theta);
    double himmelblau(v_d theta);
    double gauss_shell(v_d theta);
    
    void set_func(std::string func_name, uint32_t ndims);

private:
    std::string name;
    uint32_t ndims_;
    // Some variables for the gaussian shell
    double shell_width, r; // width and radius of the shells
};



#endif