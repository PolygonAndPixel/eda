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

 #include "PolyChord.h"

/** Constructor and destructor **/
PolyChord::PolyChord(

    ) : Minimizer(tolerance, max_iter, min_iter,
                  max_points, seed, dump_points)
{
    nprior_                 = -1;
    nGrade_                 = 1;    // The number of grades
    grade_frac_             = new double[1];
    grade_frac_[0]          = 1.0;    // The fraction of time spent in each grade
    // Various variables for the minimization routine in fortran.
    feedback_ = 0;
    write_resume_ = false;
    write_paramnames_ = false;
    read_resume_ = false;
    write_stats_ = false;
    write_live_ = false;
    write_dead_ = false;
    write_prior_ = false;
    update_files_ = -1;
    nDerived_ = 0;

    miniter_ = -1; // Needed for I3SingleServiceFactory. Not used by PolyChord.
    setbuf(stdout, NULL);
    base_dir_ = "/data/user/mhieronymus/polychord_dump/";
    file_root_ = "";

    // GetParameter("Tolerance", precision_criterion_);
    // GetParameter("Posteriors", posteriors_);
    // GetParameter("Equally", equals_);
    // GetParameter("ClusterPosteriors", cluster_posteriors_);
    // // What factor should we bulk up the posterior points by (using
    // // inter-chain points)
    // // set to <=0 to use all of them
    // GetParameter("BoostPosterior", boost_posterior_);
    // GetParameter("MaxDead", max_ndead_); // -1 -> no maximum number
    // GetParameter("NLive", nlive_);
    // GetParameter("NumRepeats", num_repeats_);
    // GetParameter("PeriodicParams", periodics_);
    // GetParameter("ModeSeparatedParams", modeSeparated_);
}

/** Function to map from the unit hypercube to Theta in the physical space.
 *
 * \param cube              Hypercube coordinates of point that
 *                          shall be evaluated.
 * \param ndims             Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 *
 *  \return                 Physical coordinates
 * */
double* PolyChord::fortran_to_physics(
    double* cube,
    uint32_t ndims) {

    double *theta = new double[ndims];

    for (int i=0; i<ndims; i++) {
        theta[i] = (this->lower_bnds[i]
            + (this->upper_bnds[i] - this->lower_bnds[i])
            * cube[i]);
    }
    return theta;
}

/** Function that shall be minimized. This is the actual function that is seen
 *  by the PolyChord Fortran code which itself refers to a Minimizer
 *  Object to evaluate the llh.
 *
 *  \param theta            Physical parameters of point that
 *                          shall be evaluated.
 *  \param ndims            Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 *  \param phi              Derived parameters (usually none)
 *  \param nDerived         Number of derived params (usually 0)
 *  \param misc             Any other stuff user wants to hand over to fct
 *                          (structure of function parameters determined by
 *                          PolyChord requirements).
 * */
double PolyChord::fortran_get_llh(
    double* theta,
    int ndims,
    double* phi,
    int nDerived,
    void *misc) {

    result.n_lh_calls++;
    return test_func_->get_lh(v_d(theta, theta+ndims));
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
 *                                  PolyChord requirements). For Gulliver, this is
 *                                  always an I3GulliverPolyChord object that contains
 *                                  all necessary informations about how to convert
 *                                  the standard PolyChord parameter ranges ([0,1])
 *                                  into physical ranges for the llh function
 *                                  (e.g. for millipede).
 * */
void PolyChord::c_dumper(
    double log_evidence,
    double error_log_evidence,
    double ndead,
    double n_likelihood_calls,
    double n_accepted,
    double *live_params,
    uint32_t n_cluster,
    double llh_best_fit,
    double llh_worst_fit,
    uint32_t ndims,
    void *misc)
{

    printf("Finished PoyChord with following paramters:\n"
           "log_evidence = %f\n"
           "error_log_evidence = %f\n"
           "ndead = %f\n"
           "n_likelihood_calls = %f\n"
           "n_cluster = %d\n"
           "best_llh = %f\n"
           "worst_llh = %f\n"
           "ndims = %d\n", log_evidence, error_log_evidence, ndead,
           n_likelihood_calls, n_cluster, llh_best_fit, llh_worst_fit, ndims);

    for(int i=0; i<nPar; i++) {
        printf("%f, ",live_params[i] );
    }
    printf("------\n");

    result.params_best_fit.clear();
    for(uint32_t i=0; i<ndims; i++) {
        result.params_best_fit[i] = live_params[i];
    }

    result.best_fit  = llh_best_fit;
    std::cout << "n_lh_calls according to PolyChord: " << n_likelihood_calls;
    std::cout << " and according to my code: " << result.n_lh_calls;
    std::cout << std::endl;
    result.n_lh_calls = n_likelihood_calls;
    lh_worstFit_ = llh_worst_fit;
    result.efficiency = n_accepted/n_likelihood_calls;
}

/** Function to map from the unit hypercube to Theta in the physical space.
 *  Not used here but still implemented because baseclass Minimizer needs it.
 *
 * \param cube              Hypercube coordinates of point that
 *                          shall be evaluated.
 * \param ndims             Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 *
 *  \return                 Physical coordinates
 * */
v_d PolyChord::to_physics(
    v_d cube,
    uint32_t ndims) {

    v_d theta(ndims);

    for (int i=0; i<ndims; i++) {
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
    if(dump_points_) {
        std::ofstream ofile((base_dir_+file_name_).c_str(),
            std::ofstream::out  | std::ofstream::app);

        for(int j=0; j<test_func_->get_ndims(); j++)
            ofile << "Param" << j << "\t";
        ofile << std::endl;
        ofile.close();
    }

    /// run Fortran routines for minimization
    // Hack to pass a string to fortran.
    char *b_dir = new char[base_dir_.size()+1];
    char *f_root = new char[file_name_.size()+1];
    std::copy(base_dir_.begin(), base_dir_.end(), b_dir);
    b_dir[base_dir_.size()] = '\0';
    std::copy(file_root_.begin(), file_root_.end(), f_root);
    f_root[file_root_.size()] = '\0';

    // The numbers are the indices of the last argument of each row.
    // It is easier to debug this way.
    interfaces::polychord_c_interface(nlive_, num_repeats_, nprior_, // 4
        do_clustering, feedback_, precision_criterion_, max_ndead_, // 8
        boost_posterior_, posteriors_, equals_, cluster_posteriors_, // 12
        write_resume_, write_paramnames_, read_resume_, write_stats_, write_live_, // 17
        write_dead_, write_prior_, update_files_, ndims, nDerived_, // 22
        b_dir, f_root, nGrade_, grade_frac_, // 26
        &grade_dims, I3GulliverPolyChord::fct, // 28
        PolyChord::to_physics, PolyChord::fortran_get_llh,
        static_cast<void*>(this)); // 30
    result.function_name = file_name_;
    result.minimizer_name = "PolyChord";
    return result;
}
