#include "Minimizer/gpu/multinest_gpu.cuh"

/** 
 * @brief: Start the minimizing process.
 *
 * \param ins                   Use importance nested sampling
 * \param mmodal                Use multimodal mode
 * \param const_eff             Enforce a constant efficiency
 * \param n_live                Number of live points
 * \param precision_criterion   Stop if improvement is lower
 * \param enlargement           Enlarge ellipsoids by this factor
 * \param efficiency
 * \param n_dims
 * \param n_cdims               Number of dimensions to cluster
 * \param max_modes             Maximum number of modes
 * \params z_tol                Tolerance to cluster for
 * \param seed
 * \param p_wrap
 * \param log_zero              The worst likelihood
 * \param max_iter              Maximum number of iterations
 * \param func                  Each number stands for a different function.
                                CUDA pointer to functions are not very 
                                efficient.
 * \param lower_bnds            Lower bounds of the parameters
 * \param upper_bnds            Upper bounds of the parameters
 */
template<typename value_t>
__host__
void run(
    bool ins, 
    bool mmodal, 
    bool const_eff, 
    int n_live, 
    value_t precision_criterion, 
    value_t enlargement,  
    value_t efficiency,
    int n_dims, 
    int n_cdims, 
    int max_modes, 
    value_t z_tol, 
    int seed, 
    value_t * p_wrap, 
    value_t log_zero, 
    int max_iter,
    const int func,
    value_t * lower_bnds,
    value_t * upper_bnds) 
{

}

/** Function to map from the unit hypercube to Theta in the physical space.
 *  Not used here but still implemented because baseclass Minimizer needs it.
 *
 * \param lower_bnds        Lower bounds in every dimension
 * \param upper_bnds        Upper bounds in every dimension
 * \param theta             On out: Physical coordinates
 * \param cube              Hypercube coordinates of point that
 *                          shall be evaluated
 * \param n_dims            Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 * */
template<typename value_t>
__device__
void to_physics(
    value_t * lower_bnds,
    value_t * upper_bnds,
    value_t * theta,
    value_t * cube,
    int n_dims)
{
    for(int i=0; i<n_dims; i++)
        theta[i] = lower_bnds[i] + (upper_bnds[i] - lower_bnds[i]) * cube[i];
}

// Needed functions:
// nestrun, initrandoms, rmarinns
// nestsample, gen_initial_live, getrandom
// clusterednest, setlimits, wraparound, samp




template<typename value_t>
__device__
void calc_ellipsoid()
{

}