#ifndef POLYCHORD_H_INCLUDED
#define POLYCHORD_H_INCLUDED

#include "MinimizerResult.h"
#include "Minimizer.h"
#include "PolyChordInterface.h"
#include "helper/abbreviations.h"

#include <random>
#include <lapacke.h>
#include <cblas.h>
#include <fstream>

class PolyChord : public Minimizer {
public:

    PolyChord(value_t tolerance, index_t max_iter=0, index_t min_iter=0,
        index_t max_points=0, index_t n_prior=-1, index_t n_grade=1,
        value_t *grade_frac=nullptr, index_t n_live=250, index_t feedback=0,
        index_t max_dead=-1, value_t boost_posterior=0.0, index_t num_repeats=-1,
        bool posteriors=false, bool equals=false, bool cluster_posteriors=false,
        bool do_clustering=true, index_t seed=1025, bool dump_points=false);

    virtual std::unique_ptr<Minimizer> clone() const override {
        return std::make_unique<PolyChord>(*this);
    }

    std::string get_name();

    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    static void fortran_to_physics(value_t *cube, value_t *theta, index_t n_dims,
        void *misc);
    static value_t fortran_get_llh(value_t *theta, index_t n_dims, value_t *phi,
        index_t nDerived, void *misc);
    /// save the results with this function
    static void c_dumper(value_t log_evidence,
        value_t error_log_evidence, value_t ndead, value_t n_likelihood_calls,
        value_t n_accpeted, value_t *live_params, index_t n_cluster,
        value_t llh_best_fit, value_t llh_worst_fit, index_t n_dims, void *misc);
    v_d to_physics(v_d cube, index_t n_dims);

private:
    index_t n_live_, nprior_, max_ndead_, nGrade_, num_repeats_;
    bool posteriors_, equals_, cluster_posteriors_;
    value_t boost_posterior_, *grade_frac_;

    // Various unimportant variables for the minimization routine in fortran.
    index_t feedback_, update_files_, nDerived_;
    bool write_resume_, write_paramnames_, read_resume_, write_stats_;
    bool write_live_, write_dead_, write_prior_, do_clustering_;

};

#endif
