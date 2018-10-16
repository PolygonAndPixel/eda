#ifndef MULTINEST_H_INCLUDED
#define MULTINEST_H_INCLUDED

#include "MinimizerResult.h"
#include "Minimizer.h"
#include "MultiNestInterface.h"
#include "helper/abbreviations.h"

#include <fstream>

class MultiNest : public Minimizer {
public:

    MultiNest(value_t tolerance, index_t max_iter=0, bool ins=false,
        bool mode_separation=true, bool const_eff=false, index_t n_live=500,
        value_t enlargement=0.5, index_t feedback_interval=10, index_t max_modes=10,
        bool feedback=false, index_t seed=1025, bool dump_points=false);

    virtual std::unique_ptr<Minimizer> clone() const override {
        return std::make_unique<MultiNest>(*this);
    }

    std::string get_name();

    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    static void fortran_get_llh(value_t *cube, index_t &n_dims, index_t &n_pars,
        value_t &llh, void *misc);
    /// save the results with this function
    static void c_dumper(index_t &n_samples, index_t &n_live, index_t &n_dims,
        value_t **phys_live, value_t **posterior, value_t **param_constr,
        value_t &llh_best_fit, value_t &log_z, value_t &ins_log_z,
        value_t &log_z_err, index_t &n_accepted, void *misc);
    v_d to_physics(v_d cube, index_t n_dims);

private:
    index_t n_live_, max_modes_, feedback_, seed_;
    bool ins_, mode_separation_, const_eff_, mmodal_;
    value_t enlargement_;

    // Various unimportant variables for the minimization routine in fortran.
    index_t feedback_interval_, resume_, init_mpi_;
    value_t logZero_;
    index_t writefiles_;

};

#endif
