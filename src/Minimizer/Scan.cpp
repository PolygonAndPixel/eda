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
    int n_points_per_dim,
    int max_points,
    int seed,
    bool dump_points) : Minimizer(0, 0, 0, max_points, seed, dump_points)
{
    n_points = n_points_per_dim;
}

/** Function that samples in a grid and dumps every max_points_ the points.
 *  \param nDims             Dimensionallity of parameter space in terms of
 *                           free parameter for minimization
 * */
void Scan::scan_space(
    uint32_t nDims){

    v_d cube(nDims, 0.0);

    if(dump_points_) {
        std::cout << "Saving files to " << base_dir_ << file_name_ << std::endl;
        std::ofstream ofile((base_dir_+file_name_).c_str(),
            std::ofstream::out  | std::ofstream::app);

        for(uint32_t j=0; j<nDims; j++) ofile << "Param" << j << "\t";
        // ofile << "X\tY\tZ\tT\tZenith\tAzimuth\tEnergy";
        ofile << std::endl;
        ofile.close();
    }
    double delta = 1.0/n_points;
    uint32_t i = 0;
    while(cube[nDims-1] <= 1.0) {
        i++;
        v_d theta = to_physics(cube, nDims);
        double current_result = get_llh(theta);
        if(dump_points_) {
            results.insert(results.end(), theta.begin(), theta.end());
            results.push_back(current_result);
        }

        if(i == 0 || result.best_fit > current_result) {
            result.best_fit = current_result;
            result.params_best_fit = theta;
        }
        cube[0] += delta;
        for(uint32_t j=0; j<nDims-1; j++) {
            if(cube[j] > 1.0) {
                cube[j+1] += delta;
                for(uint32_t k=j; (k>=0 && k>=j-1); k--) cube[k] = 0;
            } else {
                break;
            }
        }

        if(dump_points_ && (i%max_points_ == 0 || cube[nDims-1] >= 1.0)) {
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
            
            for(auto &c: cube) std::cout << c << ",\t";
            std::cout << std::endl;
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
    uint32_t nDims) {

    v_d theta;
    
    for (uint32_t i=0; i<nDims; i++) {
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
    scan_space(test_func_->get_ndims());
    result.minimizer_name = "Scan";
    result.function_name = test_func_->get_name();
    return result;
}
