/**
 * @brief 
 * It just samples in a grid from the space for plots. 
 * It does not minimize anything!
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 */
#include "Minimizer/Scan.h"
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/nondet_random.hpp>
#include <fstream>

/** Constructor and destructor **/
Scan::Scan(
    index_t n_points_per_dim,
    index_t max_points,
    index_t seed,
    bool dump_points) : Minimizer(0, 0, 0, max_points, seed, dump_points)
{
    n_points = n_points_per_dim;
}

/** Function that samples in a grid and dumps every max_points_ the points.
 *  \param nDims             Dimensionallity of parameter space in terms of
 *                           free parameter for minimization
 * */
void Scan::scan_space(
    index_t nDims){

    v_d cube(nDims, 0.0);

    if(dump_points_) {
        std::cout << "Saving files to " << base_dir_ << file_name_ << std::endl;
        std::ofstream ofile((base_dir_+file_name_).c_str(),
            std::ofstream::out  | std::ofstream::app);

        for(index_t j=0; j<nDims; j++) ofile << "Param" << j << "\t";
        // ofile << "X\tY\tZ\tT\tZenith\tAzimuth\tEnergy";
        ofile << std::endl;
        ofile.close();
    }
    value_t delta = 1.0/n_points;
    index_t i = 0;
    while(cube[nDims-1] <= 1.0) {
        i++;
        v_d theta = to_physics(cube, nDims);
        value_t current_result = get_llh(theta);
        if(dump_points_) {
            results.insert(results.end(), theta.begin(), theta.end());
            results.push_back(current_result);
        }

        if(i == 0 || result.best_fit > current_result) {
            result.best_fit = current_result;
            result.params_best_fit = theta;
        }
        cube[0] += delta;
        for(index_t j=0; j<nDims-1; j++) {
            if(cube[j] > 1.0) {
                cube[j+1] += delta;
                for(index_t k=j; (k>=0 && k>=j-1); k--) cube[k] = 0;
            } else {
                break;
            }
        }

        if(dump_points_ && (i%max_points_ == 0 || cube[nDims-1] >= 1.0)) {
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
            
            for(auto &c: cube) std::cout << c << ",\t";
            std::cout << std::endl;
        }
    }
    result.lh_efficiency = 1.0/i;
        
}


/** Function that samples in a grid and dumps every max_points_ the points.
 *  \param nDims             Dimensionallity of parameter space in terms of
 *                           free parameter for minimization
 * */
void Scan::scan_space_tmp(
    index_t nDims){

    v_d cube(nDims, 0.0);
//     value_t x,
//     value_t y,
//     value_t z,
//     value_t t,
//     value_t theta,
//     value_t phi,
//     value_t length,
//     value_t seg_length,
//     index_t seed) {
// // track(35.0, -21.0, -5.0, 0.0, 0.6, 1.5, 5.0, 5.0);
    if(dump_points_) {
        std::cout << "Saving files to " << base_dir_ << file_name_ << std::endl;
        std::ofstream ofile((base_dir_+file_name_).c_str(),
            std::ofstream::out  | std::ofstream::app);

        for(index_t j=0; j<nDims; j++) ofile << "Param" << j << "\t";
        // ofile << "X\tY\tZ\tT\tZenith\tAzimuth\tEnergy";
        ofile << std::endl;
        ofile.close();
    }
    index_t i = 0;
    value_t delta = 1.0/n_points;
    for(value_t a = 0; a <= 1.0; a += delta) {
        std::cout << a << "\n";
        for(value_t b = 0; b <= 1.0; b += delta) {
            for(value_t c = 0; c <= 1.0; c += delta) {
                i++;
                cube[0] = c;
                cube[1] = b;
                cube[3] = a;
                v_d theta = to_physics(cube, nDims);
                theta[2] = -5.0;
                theta[4] = 0.6;
                theta[5] = 1.5;
                theta[6] = 5.0;

                value_t current_result = get_llh(theta);
                if(dump_points_) {
                    results.insert(results.end(), theta.begin(), theta.end());
                    results.push_back(current_result);
                }

                if(i == 0 || result.best_fit > current_result) {
                    result.best_fit = current_result;
                    result.params_best_fit = theta;
                }
                if(dump_points_ && (i%max_points_ == 0 || a >= 1.0)) {
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
        }
    }
    result.lh_efficiency = 1.0/i;
        
}

/** Return the name of this class.
 *
 *  \return     Name of this class.
 */
std::string Scan::get_name() {
    return ("Scan (Scanning the space in a grid)");
}

/** Function to map from the unit hypercube to Theta in the physical space.
 *
 * \param cube              Hypercube coordinates of point that
 *                          shall be evaluated.
 * \param nDims             Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 * */
v_d Scan::to_physics(
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
Scan::Minimize(
    TestFunctions test_func,
    v_d lower_bounds,
    v_d upper_bounds ) {

    reset_calls();
    results.clear();
    upper_bnds = upper_bounds;
    lower_bnds = lower_bounds;
    test_func_ = &test_func;
    file_name_ = test_func_->get_name();
    scan_space_tmp(test_func_->get_ndims());
    result.minimizer_name = "Scan";
    result.function_name = test_func_->get_name();
    return result;
}
