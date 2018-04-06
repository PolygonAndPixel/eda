/* Several different test functions for minimizing. See
 * https://en.wikipedia.org/wiki/Test_functions_for_optimization
 * for a list of test functions and their definitions.
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */

#include "TestFunctions.h"

TestFunctions::TestFunctions() {


    lh_p = &TestFunctions::gauss_shell;
    ndims_ = 2;
    // Values taken from
    // "MultiNest: an efficient and robust Bayesian inference tool
    // for cosmology and particle physics"
    // radius made a bit smaller to avoid overlaps.
    shell_width = 0.1;
    r = 0.5;

}

TestFunctions::TestFunctions(
    std::string func_name,
    uint32_t ndims) {

    name = func_name;
    if(func_name == "eggholder") {
        lh_p = &TestFunctions::eggholder;
        ndims_ = 2;
    }
    else if(func_name == "townsend") {
        lh_p = &TestFunctions::townsend;
        ndims_ = 2;
    }
        else if(func_name == "rosenbrock") {
        lh_p = &TestFunctions::rosenbrock;
        ndims_ = ndims;
    }
    else if(func_name == "himmelblau") {
        lh_p = &TestFunctions::himmelblau;
        ndims_ = 2;
    }
    else {
        lh_p = &TestFunctions::gauss_shell;
        ndims_ = ndims;
        // Values taken from
        // "MultiNest: an efficient and robust Bayesian inference tool
        // for cosmology and particle physics"
        shell_width = 0.1;
        r = 0.5;
    }
}

/** Getter for the name of the used function.
 *
 *  \return         name
 * */
std::string TestFunctions::get_name() {
    return name;
}

/** Getter for the dimension.
 *
 *  \return         ndims_
 * */
uint32_t TestFunctions::get_ndims() {

    return ndims_;
}

/** Call the member function stored in lh_p. If you can use C++17, I
 *  recommend using std::invoke instead.
 *
 *  \param theta    Physical parameters of point that shall be evaluated.
 *
 *  \return         Likelihood
 * */
double TestFunctions::get_lh(
    v_d theta) {

    return CALL_MEMBER_FN(*this, lh_p)(theta);
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
double TestFunctions::eggholder(
    v_d theta) {

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
double TestFunctions::townsend(
    v_d theta) {

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
double TestFunctions::rosenbrock(
    v_d theta) {

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
double TestFunctions::himmelblau(
    v_d theta) {

    double left = (theta[0]*theta[0] + theta[1] - 11);
    left *= left;
    double right = (theta[0] + theta[1]*theta[1] - 7);
    right *= right;
    return left+right;
}

/** Calculate the negative log-likelihood by evaluating the
 *  gaussian shell function, where -log(0) is set to 0.
 *  See "MultiNest: an efficient and robust Bayesian inference tool for
 *  cosmology and particle physics", Eq. 32 and 33.
 *  We always set the center of the shells at (-2.5, 0, 0, ..., 0)
 *  and (2.5, 0, 0, ..., 0).
 *  Dimensions should be bounded by -6 and 6 in all dimensions.
 *
 *
 *  \param theta    Physical parameters of point that shall be evaluated.
 *
 *  \return         Likelihood
 * */
double TestFunctions::gauss_shell(
    v_d theta) {

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
    right = factor * exp(-right);

    if(isinf(log(left + right))) return 0;
    return log(left + right);
}

/** Change the used function and the number of dimensions.
 *
 *  \param func_name    The name of the desired function. Options are:
 *                      egg:        eggholder function
 *                      town:       townsend function
 *                      rosenbrock: rosenbrock function
 *                      himmelblau: himmelblau's function
 *                      else just use gaussian shell function
 *  \param ndims        The number of dimension.
 *
 * */
void TestFunctions::set_func(
    std::string func_name,
    uint32_t ndims) {

    name = func_name;
    if(func_name == "eggholder") {
        lh_p = &TestFunctions::eggholder;
        ndims_ = 2;
    }
    else if(func_name == "townsend") {
        lh_p = &TestFunctions::townsend;
        ndims_ = 2;
    }
        else if(func_name == "rosenbrock") {
        lh_p = &TestFunctions::rosenbrock;
        ndims_ = ndims;
    }
    else if(func_name == "himmelblau") {
        lh_p = &TestFunctions::himmelblau;
        ndims_ = 2;
    }
    else {
        lh_p = &TestFunctions::gauss_shell;
        ndims_ = ndims;
        // Values taken from
        // "MultiNest: an efficient and robust Bayesian inference tool
        // for cosmology and particle physics"
        shell_width = 0.1;
        r = 0.5;
    }
}
