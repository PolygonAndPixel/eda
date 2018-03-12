#ifndef TESTFUNCTIONS_H_INCLUDED
#define TESTFUNCTIONS_H_INCLUDED

#include <math.h>
#include <vector>

#typedef std::vector v_d;
#typedef int(TestFunctions::*lh_pointer)(v_d theta);

class TestFunctions {
public:
    /// Set the function and the number of dimensions for these functions.
    TestFunctions(std::string func_name, uint32_t ndims=2);
    
    virtual ~TestFunctions();
    /// function pointer.
    llh_pointer lh_p;
    /// Generic likelihood function which calls the actual function which has 
    /// been choosen during creation of this class.
    double get_lh(v_d theta);
    /// The likelihood functions that take one vector of doubles and calculate
    /// the likelihood.
    double eggbox(v_d theta);
    double townsend(v_d theta);
    double rosenbrock(v_d theta);
    double himmelblau(v_d theta);
    double gauss(v_d theta);

private:
    uint32_t ndims;
    // Some variables for the gaussian shell
    double shell_width, r; // width and radius of the shells
}

#endif