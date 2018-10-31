/**
 * @brief Interface to a GPU implementation of MultiNest.
 * 
 * References:
 * 'MultiNest: an efficient and robust Bayesian inference tool for cosmology and particle physics'
 * F. Feroz, M.P. Hobson, M. Bridges. Sep 2008. 14 pp.
 * Published in Mon.Not.Roy.Astron.Soc. 398 (2009) 1601-1614
 * DOI: 10.1111/j.1365-2966.2009.14548.x
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 */
#include "Minimizer/MultiNestGPU.cuh"

/** Constructor and destructor **/
MultiNestGPU::MultiNestGPU(
    double tolerance,
    int max_iter,
    bool ins,
    bool mode_separation,
    bool const_eff,
    int n_live,
    double enlargement,
    int feedback_interval,
    int max_modes,
    bool feedback,
    int seed,
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
    logZero_                = -std::numeric_limits<double>::max();
    base_dir_               = "tmp"; // temporary files are written here

}

/** Return the name of this class.
 *
 *  \return     Name of this class.
 */
std::string MultiNestGPU::get_name() {
    return ("MultiNestGPU");
}

/** Required Minimize() function for every minimizer. Sets the bounds.
 *  Copies data to GPU and calls the respective minimizer there.
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
    test_func_ = &test_func;
    file_name_ = test_func_->get_name() + "_";
    params_best_fit.resize(test_func_->get_ndims());
	seed_ = intgen()%30081; // Specified by MultiNest

    writefiles_ = dump_points_;
    int n_dims = test_func_->get_ndims();

    int pWrap[n_dims];
    // We don't use peridodicity
    for(uint32_t i=0; i<n_dims; ++i) {
        pWrap[i] = 0;
    }
    double Ztol = -1E90;
    // Copy data to GPU

    run(ins_, mmodal_, const_eff_, n_live_, precision_criterion_, //5
        enlargement_,  n_dims, n_dims, // 8
        n_dims, // Number of parameters on which clustering should be done
        max_modes_, feedback_interval_, // 11
        Ztol, b_dir, seed_, // 14
        pWrap, feedback_, resume_, // 17
        writefiles_, init_mpi_, logZero_, max_iter_,
        test_func_);

    // Get the result back
    result.function_name = file_name_;
    result.minimizer_name = "MultiNestGPU";
    return result;
}
