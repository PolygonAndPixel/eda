/* Several different test functions for minimizing. See 
 * https://en.wikipedia.org/wiki/Test_functions_for_optimization
 * for a list of test functions and their definitions.
 * 
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */

#include "TestFunctions.h"

TestFunctions::TestFunctions(std::string func_name, uint32_t ndims) {
    
    switch(func_name) {
        case "egg":
            lh_p = &TestFunctions::eggholder;
            ndims = 2;
            break;
        case "town":
            lh_p = &TestFunctions::townsend;
            ndims = 2;
            break;
        case "rosenbrock":
            lh_p = &TestFunctions::rosenbrock;
            ndims = ndims;
            break;
        case "himmelblau":
            lh_p = &TestFunctions::himmelblau;
            ndims = 2;
            break;
        default:
            lh_p = &TestFunctions::gauss;
            ndims = ndims;
            // Values taken from
            // "MultiNest: an efficient and robust Bayesian inference tool 
            // for cosmology and particle physics"
            shell_width = 0.1;
            r = 2;
    }
}

/** Call the member function stored in lh_p. If you can use C++17, I 
 *  recommend using std::invoke instead.
 * 
 *  \param theta    Physical parameters of point that shall be evaluated.
 * 
 *  \return         Likelihood
 * */
double TestFunctions::get_lh(v_d theta) {
    
    return *lh_p(theta);
}

/** Calculate the likelihood by evaluating the eggholder function.
 *  If the dimensions are bounded to -512 and 512 then the minimum is at
 *  Minimum:        f(512, 404.2319) = -959.6407
 * 
 * 
 *  \param theta    Physical parameters of point that shall be evaluated.
 * 
 *  \return         Likelihood
 * */
double TestFunctions::eggholder(v_d theta) {
    double left = -(theta[1] + 47)*sin(sqrt(abs(theta[0]/2 + (theta[1]+47))));
    double right = theta[0] * sin(sqrt(abs(theta[0] - (theta[1]+47))));
    return left - right;
}

/** Calculate the likelihood by evaluating a modified townsend function.
 *  This function has some constraints that simply return infinity if not 
 *  matched.
 *  If the dimensions are bounded -2.25 <= theta[0] <= 2.5 and 
 *  -2.5 <= theta[1] <= 1.75 then the minimum is at
 *  Minimum:        f(2.0052938, 1.1944509) = -2,0239884
 * 
 *  \param theta    Physical parameters of point that shall be evaluated.
 * 
 *  \return         Likelihood
 * */
double TestFunctions::townsend(v_d theta) {
    double t = atan2(theta[0], theta[1]);
    double constraint = 2*cos(t) - 0.5*cos(2*t) 
        - 0.25*cos(3*t) - 0.125*cos(4*t);
    constraint *= constraint;
    constraint += (2*sin(t)) * (2*sin(t));
    if(theta[0]*theta[0] + theta[1]*theta[1] < constraint) {
        double left = cos((theta[0] - 0.1)*theta[1]);
        left *= left;
        return (-left - (theta[0]*sin(3*theta[0]+theta[1])));
    } else {
        return std::numeric_limits<double>::infinity();
    }
}

/** Calculate the likelihood by evaluating the rosenbrock function.
 *  
 *  Minimum:        f(1, 1, ..., 1) = 0
 * 
 *  \param theta    Physical parameters of point that shall be evaluated.
 * 
 *  \return         Likelihood
 * */
double TestFunctions::rosenbrock(v_d theta) {
    double lh = 0;
    for(uint32_t i = 0; i<theta.size()-1; ++i) {
        lh += (100 * (theta[i+1]-theta[i]*theta[i]) 
            * (theta[i+1]-theta[i]*theta[i]) + (theta[i]-1)*(theta[i]-1));
    }
    return lh;
}

/** Calculate the likelihood by evaluating Himmelblau's function.
 *  This has multiple minima.
 *  If the dimensions are bounded to -5 and 5 then the minima are at
 *  Minima:        f(3.0, 2.0)                  = 0.0
 *                 f(-2.805118, 3.131312)       = 0.0
 *                 f(-3.779310, -3.283186)      = 0.0
 *                 f(3.584428, -1.848126)       = 0.0
 * 
 *  \param theta    Physical parameters of point that shall be evaluated.
 * 
 *  \return         Likelihood
 * */
double TestFunctions::himmelblau(v_d theta) {
    double left = (theta[0]*theta[0] + theta[1] - 11);
    left *= left;
    double right = (theta[0] + theta[1]*theta[1] - 7);
    right *= right;
    return left+right;
}

/** Calculate the likelihood by evaluating teh gaussian shell function.
 *  See "MultiNest: an efficient and robust Bayesian inference tool for 
 *  cosmology and particle physics", Eq. 32 and 33.
 *  We always set the center of the shells at (-2.5, 0, 0, ..., 0)
 *  and (2.5, 0, 0, ..., 0).
 *  Dimensions should be bounded by -6 and 6 in all dimensions.
 * 
 *  \param theta    Physical parameters of point that shall be evaluated.
 * 
 *  \return         Likelihood
 * */
double TestFunctions::gauss(v_d theta) {
    double factor = 1/sqrt(2*M_PI*shell_width*shell_width);
    double left = (theta[0]-2.5)*(theta[0]-2.5);
    double right = (theta[0]+2.5)*(theta[0]+2.5);
    for(v_d::iterator theta_iter=theta.begin()+1; theta_iter<theta.end(); 
        theta_iter++) {
        
        left += (*theta_iter) * (*theta_iter);
        right += (*theta_iter) * (*theta_iter);
    }
    left = sqrt(left) - r;
    left *= left;
    left /= 2 * shell_width*shell_width;
    left = factor * exp(-left);
    
    right = sqrt(right) - r;
    right *= right;
    right /= 2 * shell_width*shell_width;
    right = factor * exp(right;
    
    return left + right;
}