/**
 * @brief
 * It just samples randomly the space for plots. It does not minimize anything!
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 */
#include "Minimizer/GradientDescent.h"
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/nondet_random.hpp>
#include <fstream>

/** Constructor and destructor **/
GradientDescent::GradientDescent(
    int max_iter,
    int max_points,
    int seed,
    double step_size,
    double conv_crit,
    bool dump_points) : Minimizer(0, max_iter, 0, max_points, seed, dump_points)
{
    stepsize = step_size;
    this->conv_crit = conv_crit;
}

// GradientDescent::~GradientDescent() {}

/** Function that uses a simpel gradient descent.
 *  \param nDims             Dimensionallity of parameter space in terms of
 *                           free parameter for minimization
 * */
void GradientDescent::descent(
    uint32_t nDims){

    int iter = 0;
    int accepted = 0;
    bool converged = false;
    v_d cube(nDims);
    for(auto &v: cube) v = 0.5;
    v_d theta = to_physics(cube, nDims);
    double llh = get_llh(theta);
    v_d gradient = get_gradient(cube, nDims, llh);

    v_d cube_old(nDims);
    for(auto &v: cube_old) v = 0.5 + stepsize;
    v_d theta_old = to_physics(cube_old, nDims);
    double llh_old = get_llh(theta_old);
    v_d gradient_old = get_gradient(cube_old, nDims, llh_old);

    result.best_fit = llh;
    result.params_best_fit = theta;

    for(uint32_t iter = 0; iter < max_iter_; iter++) {
        cube_old = cube;
        llh_old = llh;
        gradient_old = gradient;
        // Get the next value
        for(uint32_t i = 0; i < nDims; i++) {
            cube[i] = cube_old[i] + stepsize*gradient_old[i];

        }
        theta = to_physics(cube, nDims);
        llh = get_llh(theta);
        result.n_lh_calls++;

        if(llh > result.best_fit) {
            result.best_fit = llh;
            result.params_best_fit = theta;
            accepted++;
        }
        converged = (fabs(llh - llh_old) < conv_crit);
        if(converged) break; 

        gradient = get_gradient(cube, nDims, llh);
    }
    result.lh_efficiency = (float) accepted / (float) result.n_lh_calls;
}

/** Calculate the gradient at the given point by evaluating the likelihood
 *  a tiny step in all other directions.
 * 
 * 
 * */
v_d GradientDescent::get_gradient(
    v_d cube, 
    uint32_t nDims, 
    double llh) {

    v_d gradient;
    for(uint32_t i = 0; i < nDims; i++) {
        v_d cube_tmp = cube;
        cube_tmp[i] += stepsize;
        v_d phys = to_physics(cube_tmp, nDims);
        double llh_tmp = get_llh(phys);
        result.n_lh_calls++;
        gradient.push_back( llh_tmp - llh );
    }
    return gradient;

}

/** Return the name of this class.
 *
 *  \return     Name of this class.
 */
std::string GradientDescent::get_name() {
    return ("GradientDescent");
}

/** Function to map from the unit hypercube to Theta in the physical space.
 *
 * \param cube              Hypercube coordinates of point that
 *                          shall be evaluated.
 * \param nDims             Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 * */
v_d GradientDescent::to_physics(
    v_d cube,
    uint32_t nDims) {

    v_d theta;
    
    for (int i=0; i<nDims; i++) {
        theta.push_back(this->lower_bnds[i]
        + (this->upper_bnds[i] - this->lower_bnds[i])
        * cube[i]);
    }

    return theta;
}

/** Required Minimize() function for every minimizer. Sets the bounds.
 *
 *  \param test_func        The function which shall be minimized
 *  \param lower_bounds     The lower bounds for each dimension
 *  \param upper_bounds     The upper bounds for each dimension
 *
 *  \return                 The result of the minimization
 * */
MinimizerResult
GradientDescent::Minimize(
    TestFunctions test_func,
    v_d lower_bounds,
    v_d upper_bounds ) {

    reset_calls();
    results.clear();
    upper_bnds = upper_bounds;
    lower_bnds = lower_bounds;
    test_func_ = &test_func;
    file_name_ = test_func_->get_name();
    descent(test_func_->get_ndims());
    result.minimizer_name = "GradientDescent";
    result.function_name = test_func_->get_name();
    return result;
}
