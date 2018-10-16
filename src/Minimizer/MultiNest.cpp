/**
 * @brief Interface to Multinest minimization algorithm for Gulliver. (Use carefully! )
 * Based on an interface from Martin Leuermann <leuermann@physik.rwth-aachen.de>
 *
 * References:
 * 'MultiNest: an efficient and robust Bayesian inference tool for cosmology and particle physics'
 * F. Feroz, M.P. Hobson, M. Bridges. Sep 2008. 14 pp.
 * Published in Mon.Not.Roy.Astron.Soc. 398 (2009) 1601-1614
 * DOI: 10.1111/j.1365-2966.2009.14548.x
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 */
#include "Minimizer/MultiNest.h"

/** Constructor and destructor **/
MultiNest::MultiNest(
    value_t tolerance,
    index_t max_iter,
    bool ins,
    bool mode_separation,
    bool const_eff,
    index_t n_live,
    value_t enlargement,
    index_t feedback_interval,
    index_t max_modes,
    bool feedback,
    index_t seed,
    bool dump_points) : Minimizer(tolerance, max_iter, 0, 0, seed, dump_points)
{
    // Various variables for the minimization routine in fortran.
    n_live_                 = n_live;
    // Do importance nested sampling?
    ins_                    = ins;
    mode_separation_        = mode_separation;
    const_eff_              = const_eff;
    enlargement_             = enlargement;
    feedback_interval_      = feedback_interval;
    max_modes_              = max_modes;
    feedback_               = feedback;
    logZero_                = -std::numeric_limits<value_t>::max();
    base_dir_               = "tmp"; // temporary files are written here

}

/** Return the name of this class.
 *
 *  \return     Name of this class.
 */
std::string MultiNest::get_name() {
    return ("MultiNest");
}

/** Function that shall be minimized. This is the actual function that is seen
 *  by the MultiNest Fortran code which itself refers to a Minimizer
 *  Object to evaluate the llh.
 *
 *  \param cube             Physical parameters of point that
 *                          shall be evaluated.
 *  \param n_dims           Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 *  \param n_pars           Parameters to consider (usually n_dims)
 *  \param llh              The new likelihood
 *  \param misc             Any other stuff user wants to hand over to fct
 *                          (structure of function parameters determined by
 *                          MultiNest requirements).
 * */
void MultiNest::fortran_get_llh(
    value_t *cube,
    index_t &n_dims,
    index_t &n_pars,
    value_t &llh,
    void *misc) {

    MultiNest *multiBase = static_cast<MultiNest*>(misc);
    multiBase->result.n_lh_calls++;
    v_d phys = multiBase->to_physics(v_d(cube, cube+n_dims), n_dims);
    llh = multiBase->test_func_->get_lh(phys);
}

/** Dumps information about the minimization at the end.
 *  \param n_samples
 *  \param n_live                    Number of live points per cluster
 *  \param n_dims                    Number of parameters
 *  \param phys_live
 *  \param posterior
 *  \param param_constr
 *  \param llh_best_fit
 *  \param log_z
 *  \param ins_log_z
 *  \param log_z_err
 *  \param misc                     Any other stuff user wants to hand over to fct
 *                                  (structure of function parameters determined by
 *                                  MultiNest requirements).
 * */
void MultiNest::c_dumper(
    index_t &n_samples,
    index_t &n_live,
    index_t &n_dims,
    value_t **phys_live,
    value_t **posterior,
    value_t **param_constr,
    value_t &llh_best_fit,
    value_t &log_z,
    value_t &ins_log_z,
    value_t &log_z_err,
    index_t &n_accepted,
    void *misc)
{
    MultiNest *multiBase = static_cast<MultiNest*>(misc);
    multiBase->result.params_best_fit.resize(n_dims);
    v_d cube(n_dims);
    for(index_t i=0; i<n_dims; ++i) {
        cube[i] = param_constr[0][2*n_dims + i];
    }
    multiBase->result.params_best_fit = multiBase->to_physics(cube, n_dims);

    value_t worst_fit = llh_best_fit;
    for(index_t k=0; k<n_live; k++) {
        if (worst_fit > phys_live[0][n_dims*n_live + k]) {
            worst_fit = phys_live[0][n_dims*n_live + k];
        }
    }
    worst_fit = -worst_fit;
    multiBase->result.best_fit  = -llh_best_fit;
    multiBase->result.lh_efficiency = (value_t) n_accepted / (value_t) multiBase->result.n_lh_calls;
}

/** Function to map from the unit hypercube to Theta in the physical space.
 *  Not used here but still implemented because baseclass Minimizer needs it.
 *
 * \param cube              Hypercube coordinates of point that
 *                          shall be evaluated.
 * \param n_dims             Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 *
 *  \return                 Physical coordinates
 * */
v_d MultiNest::to_physics(
    v_d cube,
    index_t n_dims) {

    v_d theta(n_dims);

    for (index_t i=0; i<n_dims; i++) {
        theta[i] = (this->lower_bnds[i]
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
MultiNest::Minimize(
    TestFunctions test_func,
    v_d lower_bounds,
    v_d upper_bounds) {

    reset_calls();
    results.clear();
    upper_bnds = upper_bounds;
    lower_bnds = lower_bounds;
    test_func_ = &test_func;
    file_name_ = test_func_->get_name() + "_";
    params_best_fit.resize(test_func_->get_ndims());
	seed_ = intgen()%30081; // Specified by MultiNest

    writefiles_ = dump_points_;
    index_t n_dims = test_func_->get_ndims();
    /// run Fortran routines for minimization
    // Hack to pass a string to fortran.
    char *b_dir = new char[base_dir_.size()+file_name_.size()+1];
    std::copy(base_dir_.begin(), base_dir_.end(), b_dir);
    std::copy(file_name_.begin(), file_name_.end(), b_dir+base_dir_.size());
    b_dir[base_dir_.size()+file_name_.size()] = '\0';

    index_t pWrap[n_dims];
    // We don't use peridodicity
    for(index_t i=0; i<n_dims; ++i) {
        pWrap[i] = 0;
    }
    value_t Ztol = -1E90;
    nested::run(ins_, mmodal_, const_eff_, n_live_, precision_criterion_, //5
        enlargement_,  n_dims, n_dims, // 8
        n_dims, // Number of parameters on which clustering should be done
        max_modes_, feedback_interval_, // 11
        Ztol, b_dir, seed_, // 14
        pWrap, feedback_, resume_, // 17
        writefiles_, init_mpi_, logZero_, max_iter_, // 21
        MultiNest::fortran_get_llh,
        MultiNest::c_dumper, static_cast<void*>(this));
    result.function_name = file_name_;
    result.minimizer_name = "MultiNest";
    return result;
}
