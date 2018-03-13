/**
 * @brief Interface to SampleSpace for Gulliver.
 * It just samples randomly the space for plots. It does not minimize anything!
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 */
#include "SampleSpace.h"
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/nondet_random.hpp>
#include <fstream>


/** Constructor and destructor **/
SampleSpace::SampleSpace(
    uint32_t max_iter, 
    uint32_t max_points,
    uint32_t seed) : Minimizer(0, max_iter, 0, max_points, seed)
{
    
}

// SampleSpace::~SampleSpace() {}

/** Function that samples randomly and dumps every max_points_ the points.
 *  \param nDims             Dimensionallity of parameter space in terms of
 *                           free parameter for minimization
 * */
void SampleSpace::sample_space(
    uint32_t nDims){

    uint32_t c = 1;
    double *cube = new double[nDims];
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<boost::mt19937&,
        boost::uniform_real<> > uf(intgen, uni_dist);

    std::ofstream ofile(base_dir_.c_str(),
        std::ofstream::out  | std::ofstream::app);
    for(int j=0; j<nDims; j++) ofile << "Param" << j << "\t";
    ofile << std::endl;
    for(int i=0; i<max_iter_; i++, c++) {
        if(i%5000 == 0) {printf("\r (>^.^)> %d     ", i); fflush(stdout);}
        else if(i%7500 == 0) {printf("\r (>°.°)> %d", i); fflush(stdout);}
        for(int j=0; j<nDims; j++) cube[j] =uf();

        v_d theta = to_physics(cube, nDims);
        results.insert(results.end(), theta.begin(), theta.end());
        results.push_back(get_llh(theta));

        if(c%max_points_ == 0 || i == max_iter_-1) {
            int d = 1;
            for(v_d::iterator p=results.begin();
                p != results.end(); ++p) {

                ofile << *p << "\t";
                if(d%(nDims+1) == 0) ofile << std::endl;
                d++;
            }
            results.clear();
        }
    }
    ofile.close();
}

/** Function that evaluates the llh.
 *
 *  \param theta            Physical parameters of point that
 *                          shall be evaluated.
 * */
double SampleSpace::get_llh(
    v_d theta) {

    return test_func_->get_lh(theta);
}

/** Function to map from the unit hypercube to Theta in the physical space.
 *
 * \param cube              Hypercube coordinates of point that
 *                          shall be evaluated.
 * \param nDims             Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 * */
v_d SampleSpace::to_physics(
     double *cube,
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
 *  \return                 The result of the minimization
 * */
MinimizerResult
SampleSpace::Minimize(
    TestFunctions test_func,
    v_d lower_bounds,
    v_d upper_bounds ) {
    
    upper_bnds = upper_bounds;
    lower_bnds = lower_bounds;
    test_func_ = &test_func;
    sample_space(test_func_->get_ndims());

    MinimizerResult result;
    return result;
}
