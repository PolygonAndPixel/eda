#ifndef MULTINEST_H_INCLUDED
#define MULTINEST_H_INCLUDED

#include "MinimizerResult.h"
#include "Minimizer.h"
#include "MultiNestInterface.h"
#include "helper/abbreviations.h"

#include <fstream>

class MultiNest : public Minimizer {
public:

    MultiNest(double tolerance, int max_iter=0, bool nis=false,
        bool mode_separation=true, bool const_eff=false, int n_live=250,
        double efficiency=1.0, int feedback_interval=100, int max_modes=100,
        bool feedback=false, uint32_t seed=1025, bool dump_points=false);

    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    static void fortran_get_llh(double *cube, int &n_dims, int &n_pars,
        double &llh, void *misc);
    /// save the results with this function
    static void c_dumper(int &n_samples, int &n_live, int &n_dims,
        double **phys_live, double **posterior, double **param_constr,
        double &llh_best_fit, double &log_z, double &ins_log_z,
        double &log_z_err, void *misc);
    v_d to_physics(v_d cube, uint32_t n_dims);

private:
    int n_live_, max_modes_, feedback_, seed_;
    bool ins_, mode_separation_, const_eff_, mmodal_;
    double efficiency_;

    // Various unimportant variables for the minimization routine in fortran.
    int feedback_interval_, resume_, init_mpi_;
    double logZero_;
    int writefiles_;

};

#endif
