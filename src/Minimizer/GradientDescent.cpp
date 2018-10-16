/**
 * @brief
 * A simple gradient descent.
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
    index_t max_iter,
    index_t min_iter,
    index_t max_points,
    index_t n_gradients,
    index_t seed,
    value_t stepsize,
    value_t conv_crit,
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
    index_t nDims){

    index_t accepted = 0;

    for(index_t k = 0; k < n_gradients_; k++) {
        index_t iter = 0;
        bool converged = false;
        v_d cube(nDims);
        for(auto &v: cube) v = uf(intgen);;

        v_d theta = to_physics(cube, nDims);
        value_t llh = get_llh(theta);
        v_d gradient = get_gradient(cube, nDims, llh);

        v_d cube_old(nDims);
        for(index_t i = 0; i < nDims; i++)
            cube_old[i] = cube[i] + stepsize_;
        v_d theta_old = to_physics(cube_old, nDims);
        value_t llh_old = get_llh(theta_old);
        v_d gradient_old = get_gradient(cube_old, nDims, llh_old);

        value_t tmp_best_fit = llh;
        v_d tmp_params_best_fit = theta;
        
        if(dump_points_) {
            results.insert(results.end(), theta.begin(), theta.end());
            results.push_back(llh);
        }
        accepted++;

        for(index_t iter = 0; iter < max_iter_; iter++) {
            cube_old = cube;
            llh_old = llh;
            gradient_old = gradient;
            // Get the next value
            for(index_t i = 0; i < nDims; i++) {
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
                    index_t d = 1;
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

#ifdef OMP
/** Function that uses a simple gradient descent but poorly parallelized.
 *  Several gradients starting at different random locations are run in
 *  parallel.
 *  \param nDims             Dimensionallity of parameter space in terms of
 *                           free parameter for minimization
 * */
void GradientDescent::descent_parallel(
    index_t nDims){

    index_t accepted = 0;
    result.best_fit = std::numeric_limits<value_t>::lowest();
    index_t n_lh_calls = 0;
    #pragma omp parallel for reduction(+: n_lh_calls)
    for(index_t k = 0; k < n_gradients_; k++) {
        index_t iter = 0;
        bool converged = false;
        v_d cube(nDims);
        for(auto &v: cube) v = uf(intgen);;

        v_d theta = to_physics(cube, nDims);
        value_t llh = get_llh(theta);
        v_d gradient = get_gradient(cube, nDims, llh);

        v_d cube_old(nDims);
        for(index_t i = 0; i < nDims; i++)
            cube_old[i] = cube[i] + stepsize_;
        v_d theta_old = to_physics(cube_old, nDims);
        value_t llh_old = get_llh(theta_old);
        v_d gradient_old = get_gradient(cube_old, nDims, llh_old);

        value_t tmp_best_fit = llh;
        v_d tmp_params_best_fit = theta;
        
        if(dump_points_) {
            results.insert(results.end(), theta.begin(), theta.end());
            results.push_back(llh);
        }
        accepted++;

        for(index_t iter = 0; iter < max_iter_; iter++) {
            cube_old = cube;
            llh_old = llh;
            gradient_old = gradient;
            // Get the next value
            for(index_t i = 0; i < nDims; i++) {
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
                    index_t d = 1;
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
    index_t nDims, 
    value_t llh,
    index_t & n_lh_calls) {

    v_d gradient;
    for(index_t i = 0; i < nDims; i++) {
        v_d cube_tmp = cube;
        cube_tmp[i] += stepsize_;
        v_d phys = to_physics(cube_tmp, nDims);
        value_t llh_tmp = get_llh(phys);
        n_lh_calls++;
        gradient.push_back( llh_tmp - llh );
    }
    return gradient;

}

#endif

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
    index_t nDims, 
    value_t llh) {

    v_d gradient;
    for(index_t i = 0; i < nDims; i++) {
        v_d cube_tmp = cube;
        cube_tmp[i] += stepsize_;
        v_d phys = to_physics(cube_tmp, nDims);
        value_t llh_tmp = get_llh(phys);
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
    index_t nDims) {

    v_d theta;
    
    for (index_t i=0; i<nDims; i++) {
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

#ifdef OMP
    descent_parallel(test_func_->get_ndims());
#else 
    descent(test_func_->get_ndims());
#endif

    result.minimizer_name = "GradientDescent";
    result.function_name = test_func_->get_name();
    return result;
}
