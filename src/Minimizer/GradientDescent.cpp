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
    int min_iter,
    int max_points,
    int n_gradients,
    int seed,
    double stepsize,
    double conv_crit,
    bool dump_points) : Minimizer(0, max_iter, 0, max_points, seed, dump_points)
{
    stepsize_ = stepsize;
    conv_crit_ = conv_crit;
    min_iter_ = min_iter;
    n_gradients_ = n_gradients;
}

// GradientDescent::~GradientDescent() {}

/** Function that uses a simple gradient descent.
 *  \param nDims             Dimensionallity of parameter space in terms of
 *                           free parameter for minimization
 * */
void GradientDescent::descent(
    uint32_t nDims){


    int accepted = 0;

    for(uint32_t k = 0; k < n_gradients_; k++) {
        int iter = 0;
        bool converged = false;
        v_d cube(nDims);
        for(auto &v: cube) v = uf(intgen);;

        v_d theta = to_physics(cube, nDims);
        double llh = get_llh(theta);
        v_d gradient = get_gradient(cube, nDims, llh);

        v_d cube_old(nDims);
        for(uint32_t i = 0; i < nDims; i++)
            cube_old[i] = cube[i] + stepsize_;
        v_d theta_old = to_physics(cube_old, nDims);
        double llh_old = get_llh(theta_old);
        v_d gradient_old = get_gradient(cube_old, nDims, llh_old);

        double tmp_best_fit = llh;
        v_d tmp_params_best_fit = theta;
        
        if(dump_points_) {
            results.insert(results.end(), theta.begin(), theta.end());
            results.push_back(llh);
        }
        accepted++;

        for(uint32_t iter = 0; iter < max_iter_; iter++) {
            cube_old = cube;
            llh_old = llh;
            gradient_old = gradient;
            // Get the next value
            for(uint32_t i = 0; i < nDims; i++) {
                cube[i] = cube_old[i] + stepsize_*gradient_old[i];

            }
            theta = to_physics(cube, nDims);
            llh = get_llh(theta);
            result.n_lh_calls++;

            if(llh > tmp_best_fit) {
                tmp_best_fit = llh;
                tmp_params_best_fit = theta;
                accepted++;
            }
            if(iter > min_iter_) {
                converged = (fabs(llh - llh_old) < conv_crit_);
            }

            if(dump_points_) {
                results.insert(results.end(), theta.begin(), theta.end());
                results.push_back(llh);

                if(iter%max_points_ == 0 || iter == max_iter_-1 || converged) {
                    uint32_t d = 1;
                    std::ofstream ofile((base_dir_+file_name_).c_str(),
                        std::ofstream::out  | std::ofstream::app);
                        
                    for(auto & p: results) {
                        ofile << p << "\t";
                        if(d%(nDims+1) == 0) ofile << std::endl;
                        d++;
                    }
                    ofile.close();
                    results.clear();
                }
            }
            if(converged) break; 

            gradient = get_gradient(cube, nDims, llh);
        }
        if(k == 0) {
            result.best_fit = tmp_best_fit;
            result.params_best_fit = tmp_params_best_fit;
        } else if(tmp_best_fit > result.best_fit) {
            result.best_fit = tmp_best_fit;
            result.params_best_fit = tmp_params_best_fit;
        }
    }
    result.lh_efficiency = (float) accepted / (float) result.n_lh_calls;
}

/** Function that uses a simple gradient descent but poorly parallelized.
 *  \param nDims             Dimensionallity of parameter space in terms of
 *                           free parameter for minimization
 * */
void GradientDescent::descent_parallel(
    uint32_t nDims){

    int accepted = 0;
    result.best_fit = std::numeric_limits<double>::lowest();
    uint32_t n_lh_calls = 0;
    #pragma omp parallel for reduction(+: n_lh_calls)
    for(uint32_t k = 0; k < n_gradients_; k++) {
        int iter = 0;
        bool converged = false;
        v_d cube(nDims);
        for(auto &v: cube) v = uf(intgen);;

        v_d theta = to_physics(cube, nDims);
        double llh = get_llh(theta);
        v_d gradient = get_gradient(cube, nDims, llh);

        v_d cube_old(nDims);
        for(uint32_t i = 0; i < nDims; i++)
            cube_old[i] = cube[i] + stepsize_;
        v_d theta_old = to_physics(cube_old, nDims);
        double llh_old = get_llh(theta_old);
        v_d gradient_old = get_gradient(cube_old, nDims, llh_old);

        double tmp_best_fit = llh;
        v_d tmp_params_best_fit = theta;
        
        if(dump_points_) {
            results.insert(results.end(), theta.begin(), theta.end());
            results.push_back(llh);
        }
        accepted++;

        for(uint32_t iter = 0; iter < max_iter_; iter++) {
            cube_old = cube;
            llh_old = llh;
            gradient_old = gradient;
            // Get the next value
            for(uint32_t i = 0; i < nDims; i++) {
                cube[i] = cube_old[i] + stepsize_*gradient_old[i];

            }
            theta = to_physics(cube, nDims);
            llh = get_llh(theta);
            n_lh_calls++;

            if(llh > tmp_best_fit) {
                tmp_best_fit = llh;
                tmp_params_best_fit = theta;
                accepted++;
            }
            if(iter > min_iter_) {
                converged = (fabs(llh - llh_old) < conv_crit_);
            }

            if(dump_points_) {
                results.insert(results.end(), theta.begin(), theta.end());
                results.push_back(llh);

                if(iter%max_points_ == 0 || iter == max_iter_-1 || converged) {
                    uint32_t d = 1;
                    std::ofstream ofile((base_dir_+file_name_).c_str(),
                        std::ofstream::out  | std::ofstream::app);
                        
                    for(auto & p: results) {
                        ofile << p << "\t";
                        if(d%(nDims+1) == 0) ofile << std::endl;
                        d++;
                    }
                    ofile.close();
                    results.clear();
                }
            }
            if(converged) break; 

            gradient = get_gradient(cube, nDims, llh, n_lh_calls);
        }
        #pragma omp critical 
        {
            if(tmp_best_fit > result.best_fit) {
                result.best_fit = tmp_best_fit;
                result.params_best_fit = tmp_params_best_fit;
            }
        }
    }
    result.n_lh_calls = n_lh_calls;
    result.lh_efficiency = (float) accepted / (float) result.n_lh_calls;
}

/** Calculate the gradient at the given point by evaluating the likelihood
 *  a tiny step in all other directions.
 *  
 *  \param cube     The hypercube coordinates of the point to evaluate
 *  \param nDims    Dimensionallity of parameter space in terms of
 *                  free parameter for minimization 
 *  \param llh      The likelihood at the point cube
 * 
 *  \return         A gradient vector
 * 
 * */
v_d GradientDescent::get_gradient(
    v_d & cube, 
    uint32_t nDims, 
    double llh) {

    v_d gradient;
    for(uint32_t i = 0; i < nDims; i++) {
        v_d cube_tmp = cube;
        cube_tmp[i] += stepsize_;
        v_d phys = to_physics(cube_tmp, nDims);
        double llh_tmp = get_llh(phys);
        result.n_lh_calls++;
        gradient.push_back( llh_tmp - llh );
    }
    return gradient;

}

/** Calculate the gradient at the given point by evaluating the likelihood
 *  a tiny step in all other directions.
 *  
 *  \param cube     The hypercube coordinates of the point to evaluate
 *  \param nDims    Dimensionallity of parameter space in terms of
 *                  free parameter for minimization 
 *  \param llh      The likelihood at the point cube
 * 
 *  \return         A gradient vector
 * 
 * */
v_d GradientDescent::get_gradient(
    v_d & cube, 
    uint32_t nDims, 
    double llh,
    uint32_t & n_lh_calls) {

    v_d gradient;
    for(uint32_t i = 0; i < nDims; i++) {
        v_d cube_tmp = cube;
        cube_tmp[i] += stepsize_;
        v_d phys = to_physics(cube_tmp, nDims);
        double llh_tmp = get_llh(phys);
        n_lh_calls++;
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
    descent_parallel(test_func_->get_ndims());
    result.minimizer_name = "GradientDescent";
    result.function_name = test_func_->get_name();
    return result;
}
