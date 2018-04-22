/**
 * @brief Interface to PolyChord minimization algorithm for Gulliver.
 * Based on Martin Leuermann's interface to MultiNest.
 *
 * References:
 * 'PolyChord: next-generation nested sampling'
 * W.J. Handley, M.P. Hobson, A.N. Lasenby. May 2015 15 pp.
 * https://arxiv.org/abs/1506.00171
 * DOI: 10.1093/mnras/stv1911
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 */
#include "Minimizer/PolyChord.h"

/** Constructor and destructor **/
PolyChord::PolyChord(
    double tolerance,
    int max_iter,
    int min_iter,
    int max_points,
    int n_prior,
    int n_grade,
    double *grade_frac,
    int n_live,
    int feedback,
    int max_dead,
    double boost_posterior,
    int num_repeats,
    bool posteriors,
    bool equals,
    bool cluster_posteriors,
    bool do_clustering,
    int seed,
    bool dump_points) : Minimizer(tolerance, max_iter, min_iter, max_points,
        seed, dump_points)
{
    nprior_                 = n_prior;
    nGrade_                 = (n_grade>0) ? n_grade : 1; // The number of grades
    // The fraction of time spent in each grade
    if(grade_frac) {
        grade_frac_         = new double[nGrade_];
        for(uint32_t i=0; i<nGrade_; ++i) grade_frac_[i] = grade_frac[i];
    } else {
        grade_frac_         = new double[nGrade_];
        for(uint32_t i=0; i<nGrade_; ++i) grade_frac_[i] = 1.0;
    }


    // Various variables for the minimization routine in fortran.
    n_live_                 = n_live;
    // The number of slow chords to draw. Default is 5*n_dims
    num_repeats_            = num_repeats;
    max_ndead_              = max_dead;     // -1 -> no maximum number
    // What factor should we bulk up the posterior points by (using
    // inter-chain points)
    // set to <=0 to use all of them
    boost_posterior_        = boost_posterior;
    // Calculate weighted posteriors
    posteriors_             = posteriors;
    // Calculate equally weighted posteriors. Reads and writes a file to disk
    equals_                 = equals;
    cluster_posteriors_     = cluster_posteriors;
    do_clustering_          = do_clustering;
    feedback_               = 0;
    write_resume_           = false;
    write_paramnames_       = false;
    read_resume_            = false;
    write_stats_            = false;
    write_live_             = dump_points;
    write_dead_             = false;
    write_prior_            = false;
    update_files_           = -1;
    nDerived_               = 0;
    base_dir_               = "tmp"; // temporary files are written here
}

/** Return the name of this class.
 *
 *  \return     Name of this class.
 */
std::string PolyChord::get_name() {
    return ("PolyChord");
}

/** Function to map from the unit hypercube to Theta in the physical space.
 *
 *  \param cube             Hypercube coordinates of point that
 *                          shall be evaluated.
 *  \param theta            Physical parameters
 *  \param n_dims            Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 *  \param misc             Any other stuff user wants to hand over to fct
 *                          (structure of function parameters determined by
 *                          PolyChord requirements).
 * */
void PolyChord::fortran_to_physics(
    double *cube,
    double *theta,
    int n_dims,
    void *misc) {

    PolyChord *pBase = static_cast<PolyChord*>(misc);

    for (int i=0; i<n_dims; i++) {
        theta[i] = (pBase->lower_bnds[i]
            + (pBase->upper_bnds[i] - pBase->lower_bnds[i])
            * cube[i]);
    }
}

/** Function that shall be minimized. This is the actual function that is seen
 *  by the PolyChord Fortran code which itself refers to a Minimizer
 *  Object to evaluate the llh.
 *
 *  \param theta            Physical parameters of point that
 *                          shall be evaluated.
 *  \param n_dims            Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 *  \param phi              Derived parameters (usually none)
 *  \param nDerived         Number of derived params (usually 0)
 *  \param misc             Any other stuff user wants to hand over to fct
 *                          (structure of function parameters determined by
 *                          PolyChord requirements).
 * */
double PolyChord::fortran_get_llh(
    double *theta,
    int n_dims,
    double *phi,
    int nDerived,
    void *misc) {

    PolyChord *pBase = static_cast<PolyChord*>(misc);
    pBase->result.n_lh_calls++;
    double llh = pBase->test_func_->get_lh(v_d(theta, theta+n_dims));
    // PolyChord maximizes, hence we take the negative result
    return -llh;
}

/** Dumps information about the minimization at the end.
 *  \param log_evidence             log(Z) (Bayesian evidence)
 *  \param error_log_evidence       Error on log(Z)
 *  \param ndead                    Number of dead live points used
 *  \param n_likelihood_calls       Number of likelihood calls
 *  \param live_params              The parameters of all the live points
 *  \param n_live                   Number of live points per cluster
 *  \param n_cluster                Number of cluster
 *  \param llh_bestFit              Likelihood of the live points
 *  \param nPar                     Number of parameters
 *  \param misc                     Any other stuff user wants to hand over to fct
 *                                  (structure of function parameters determined by
 *                                  PolyChord requirements).
 * */
void PolyChord::c_dumper(
    double log_evidence,
    double error_log_evidence,
    double ndead,
    double n_likelihood_calls,
    double n_accepted,
    double *live_params,
    int n_cluster,
    double llh_best_fit,
    double llh_worst_fit,
    int n_dims,
    void *misc)
{
    PolyChord *pBase = static_cast<PolyChord*>(misc);
    pBase->result.params_best_fit.resize(n_dims);
    for(uint32_t i=0; i<n_dims; i++) {
        pBase->result.params_best_fit[i] = live_params[i];
    }

    pBase->result.best_fit  = -llh_best_fit;
    pBase->result.lh_efficiency = n_accepted/n_likelihood_calls;
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
v_d PolyChord::to_physics(
    v_d cube,
    uint32_t n_dims) {

    v_d theta(n_dims);

    for (int i=0; i<n_dims; i++) {
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
PolyChord::Minimize(
    TestFunctions test_func,
    v_d lower_bounds,
    v_d upper_bounds) {

    reset_calls();
    results.clear();
    upper_bnds = upper_bounds;
    lower_bnds = lower_bounds;
    test_func_ = &test_func;
    file_name_ = test_func_->get_name();
    params_best_fit.resize(test_func_->get_ndims());

    int n_dims = test_func_->get_ndims();
    int num_repeats = (num_repeats_ < 0) ? 5*n_dims : num_repeats_;
    /// run Fortran routines for minimization
    // Hack to pass a string to fortran.
    char *b_dir = new char[base_dir_.size()+1];
    char *f_root = new char[file_name_.size()+1];
    std::copy(base_dir_.begin(), base_dir_.end(), b_dir);
    b_dir[base_dir_.size()] = '\0';
    std::copy(file_name_.begin(), file_name_.end(), f_root);
    f_root[file_name_.size()] = '\0';
    write_live_ = dump_points_;

    // The numbers are the indices of the last argument of each row.
    // It is easier to debug this way.
    interfaces::polychord_c_interface(n_live_, num_repeats, nprior_, // 4
        do_clustering_, feedback_, precision_criterion_, max_ndead_, // 8
        boost_posterior_, posteriors_, equals_, cluster_posteriors_, // 12
        write_resume_, write_paramnames_, read_resume_, write_stats_, write_live_, // 17
        write_dead_, write_prior_, update_files_, n_dims, nDerived_, // 22
        b_dir, f_root, nGrade_, grade_frac_, // 26
        &n_dims, PolyChord::fortran_get_llh, // 28
        PolyChord::fortran_to_physics, PolyChord::c_dumper,
        static_cast<void*>(this)); // 30
    result.function_name = file_name_;
    result.minimizer_name = "PolyChord";
    return result;
}
