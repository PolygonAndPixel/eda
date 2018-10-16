#include <stdint.h>
#include <opari2/pomp2_lib.h>
#include <opari2/pomp2_user_lib.h>
___POMP2_INCLUDE___
#line 1 "/home/mhierony/Documents/Uni/Masterthesis/git/github/eda/src/likelihood/TestFunctions.cpp"
/* Several different test functions for minimizing. See
 * https://en.wikipedia.org/wiki/Test_functions_for_optimization
 * for a list of test functions and their definitions.
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */
#include <iostream>
#include "likelihood/TestFunctions.h"

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
    uint32_t ndims,
    uint32_t n_x,
    uint32_t n_y,
    uint32_t n_z) {

    name = func_name;
    if(func_name == EGG) {
        lh_p = &TestFunctions::eggholder;
        ndims_ = 2;
    }
    else if(func_name == TOWN) {
        lh_p = &TestFunctions::townsend;
        ndims_ = 2;
    }
        else if(func_name == ROSEN) {
        lh_p = &TestFunctions::rosenbrock;
        ndims_ = ndims;
    }
    else if(func_name == HIMMEL) {
        lh_p = &TestFunctions::himmelblau;
        ndims_ = 2;
    }
    else if(func_name == PARABOLOID) {
        lh_p = &TestFunctions::paraboloid;
        ndims_ = ndims;
    }
    else if(func_name == ICECUBE) {
        // Generate detector
        uint32_t dist_xy = 50;
        uint32_t dist_z = 15;

        v_d dummy;
        v_d dummy2;
        std::mt19937 intgen(1025);
        std::normal_distribution<double> nf(0.2, 0.03);
        for(uint32_t i=0; i<n_x; i++) {
            for(uint32_t j=0; j<n_y; j++) {
                for(uint32_t k=0; k<n_z; k++) {
                    double x = (i-0.5*(n_x-1)) * dist_xy;
                    double y = (j-0.5*(n_y-1)) * dist_xy;
                    double z = (k-0.5*(n_z-1)) * dist_z;
                    DOM dom(x, y, z, 0.4+nf(intgen));
                    pulse_map.add_DOM(dom, dummy, dummy2);
                }
            }
        }
        // Generate track (single event)
        length = 5.0;
        seg_length = 5.0;
        Track track(35.0, -21.0, -5.0, 0.0, 0.6, 1.5, 5.0, 5.0);
        // Generate data
        pulse_map.fill(track);
        // track.fill_PulseMap(pulse_map);
        pulse_map.add_noise();
        lh_p = &TestFunctions::icecube;
        // x, y, z, t, theta, phi, energy
        ndims_ = 7;
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
 *  \return         negative Likelihood
 * */
double TestFunctions::get_lh(
    v_d & theta) {

    return CALL_MEMBER_FN(*this, lh_p)(theta);
}

/** Call the member function stored in lh_p. If you can use C++17, I
 *  recommend using std::invoke instead. 
 *  Returns the negative likelihood.
 *
 *  \param theta    Physical parameters of point that shall be evaluated.
 *
 *  \return         Negative likelihood
 * */
double TestFunctions::get_neg_lh(
    v_d & theta) {

    return -CALL_MEMBER_FN(*this, lh_p)(theta);
}

/** Call the member function stored in lh_p. If you can use C++17, I
 *  recommend using std::invoke instead. 
 *  Returns the negative log-likelihood which is used for real-world problems.
 *  However, most test functions here are usable with just the likelihood.
 *
 *  \param theta    Physical parameters of point that shall be evaluated.
 *
 *  \return         negative log-likelihood
 * */
double TestFunctions::get_neg_llh(
    v_d & theta) {

    return -log(CALL_MEMBER_FN(*this, lh_p)(theta));
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
    v_d & theta) {

    double left = -(theta[1] + 47)*sin(sqrt(abs(theta[0]/2 + (theta[1]+47))));
    double right = theta[0] * sin(sqrt(abs(theta[0] - (theta[1]+47))));

    return left - right;
}

/** Calculate the likelihood by evaluating a paraboloid.
 *  Minimum:        f(0, 0, ..., 0) = 0
 *
 *
 *  \param theta    Physical parameters of point that shall be evaluated.
 *
 *  \return         Likelihood
 * */
double TestFunctions::paraboloid(
    v_d & theta) {

    double sum = 0;
    for(auto v: theta) {
        sum += v*v;
    }
    return sum;
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
    v_d & theta) {

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
    v_d & theta) {

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
    v_d & theta) {

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
    v_d & theta) {

    double factor = 1.0/sqrt(2*M_PI*shell_width*shell_width);

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

    return -(left + right);
}

/** Calculate the negative log-likelihood by evaluating the
 *  likelihood of a neutrino event within an IceCube dummy detector with
 *  only one event modelled after Millipede.
 *
 *  \param theta    Physical parameters of point that shall be evaluated.
 *
 *  \return         Likelihood
 * */
double TestFunctions::icecube(
    v_d & theta) {

    double llh = 0.0;
    DOM dom;
    v_d charges;
    v_d times;
    Track test_track(theta[0], theta[1], theta[2], theta[3], theta[4],
        theta[5], theta[6], seg_length);
    while(pulse_map.get_next(dom, charges, times)) {
        double p_time_sum = 0.0;
        ESource source;
        // Sources in track is empty?
        while(test_track.get_next_source(source)) {
            double delta_r = dist(source.get_pos(), dom.get_pos());
            if(!charges.empty()) {
                for(uint32_t i=0; i<times.size(); i++) {
                    double delta_time = times[i] - source.get_time();
                    p_time_sum += hit_prob(delta_r, delta_time);
                }
            }
            llh -= hit_charge(delta_r);
        }
        if(!charges.empty()) {
            double charge_sum = 0.0;
            for(auto &i: charges) charge_sum += i;
            llh += charge_sum * log(p_time_sum + 0.3/80.0);
        }
    }
    std::cout << "llh before anything " << llh << " charge_sum: " << charge_sum << "\n";
    if(llh <= 0) return std::numeric_limits<double>::infinity();
    return llh;
}

/** Change the used function and the number of dimensions.
 *
 *  \param func_name    The name of the desired function. Options are:
 *                      egg:        eggholder function
 *                      town:       townsend function
 *                      rosenbrock: rosenbrock function
 *                      himmelblau: himmelblau's function
 *                      icecube:    IceCube toy model
 *                      else just use gaussian shell function
 *  \param ndims        The number of dimension.
 *
 * */
void TestFunctions::set_func(
    std::string func_name,
    uint32_t ndims,
    uint32_t n_x,
    uint32_t n_y,
    uint32_t n_z) {

    name = func_name;
    if(func_name == EGG) {
        lh_p = &TestFunctions::eggholder;
        ndims_ = 2;
    }
    else if(func_name == TOWN) {
        lh_p = &TestFunctions::townsend;
        ndims_ = 2;
    }
        else if(func_name == ROSEN) {
        lh_p = &TestFunctions::rosenbrock;
        ndims_ = ndims;
    }
    else if(func_name == HIMMEL) {
        lh_p = &TestFunctions::himmelblau;
        ndims_ = 2;
    }
    else if(func_name == PARABOLOID) {
        lh_p = &TestFunctions::paraboloid;
        ndims_ = ndims;
    }
    else if(func_name == ICECUBE) {
        // Generate detector
        uint32_t dist_xy = 50;
        uint32_t dist_z = 15;

        v_d dummy;
        v_d dummy2;
        std::mt19937 intgen(1025);
        std::normal_distribution<double> nf(0.2, 0.03);
        for(uint32_t i=0; i<n_x; i++) {
            for(uint32_t j=0; j<n_y; j++) {
                for(uint32_t k=0; k<n_z; k++) {
                    double x = (i-0.5*(n_x-1)) * dist_xy;
                    double y = (j-0.5*(n_y-1)) * dist_xy;
                    double z = (k-0.5*(n_z-1)) * dist_z;
                    DOM dom(x, y, z, 0.4+nf(intgen));
                    pulse_map.add_DOM(dom, dummy, dummy2);
                }
            }
        }
        // Generate track (single event)
        length = 5.0;
        seg_length = 5.0;
        // x, y, z, t, theta, phi, length, seg_length
        Track track(35.0, -21.0, -5.0, 0.0, 0.6, 1.5, length, seg_length);
        // Generate data
        pulse_map.fill(track);
        // track.fill_PulseMap(pulse_map);
        pulse_map.add_noise();
        lh_p = &TestFunctions::icecube;
        // x, y, z, t, theta, phi, energy
        ndims_ = 7;
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
