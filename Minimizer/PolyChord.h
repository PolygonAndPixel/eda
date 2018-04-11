#ifndef POLYCHORD_H_INCLUDED
#define POLYCHORD_H_INCLUDED

#include "MinimizerResult.h"
#include "Minimizer.h"
#include "../helper/abbreviations.h"

#include <random>
#include <lapacke.h>
#include <cblas.h>
#include <fstream>

class PolyChord : public Minimizer {
public:

    PolyChord(double tolerance, uint32_t max_iter=0, uint32_t min_iter=0,
        uint32_t max_points=0, int n_prior=-1, int n_grade=1, 
        double *grade_frac=nullptr, uint32_t n_live=250, uint32_t feedback=0, 
        int max_dead=-1, double boost_posterior=0.0, int num_repeats=-1,
        bool posteriors=false, bool equals=false, bool cluster_posteriors=false,
        uint32_t seed=1025, bool dump_points=false);

    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    void execute_polychord(uint32_t ndims);
    double* PolyChord::fortran_to_physics(double* cube, uint32_t ndims);
    double PolyChord::fortran_get_llh(double* theta, int ndims, double* phi,
        int nDerived, void *misc);
    /// save the results with this function
    void c_dumper(double log_evidence,
        double error_log_evidence, double ndead, double n_likelihood_calls,
        double n_accpeted, double *live_params, uint32_t n_cluster,
        double llh_best_fit, double llh_worst_fit, uint32_t ndims, void *misc);
    v_d PolyChord::to_physics(v_d cube, uint32_t ndims);

private:
    int n_live_, nprior_, max_ndead_, nGrade_;
    bool posteriors_, equals_, cluster_posteriors_;
    double boost_posterior_, *grade_frac_;

    // Various unimportant variables for the minimization routine in fortran.
    int feedback_, update_files_, nDerived_;
    bool write_resume_, write_paramnames_, read_resume_, write_stats_;
    bool write_live_, write_dead_, write_prior_;

};

#endif
