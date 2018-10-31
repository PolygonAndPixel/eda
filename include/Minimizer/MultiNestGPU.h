#ifndef MULTINEST_H_INCLUDED
#define MULTINEST_H_INCLUDED

#include "MinimizerResult.h"
#include "Minimizer.h"
#include "MultiNestInterface.h"
#include "helper/abbreviations.h"

#include <fstream>

class MultiNestGPU : public Minimizer {
public:

    MultiNestGPU(double tolerance, int max_iter=0, bool ins=false,
        bool mode_separation=true, bool const_eff=false, int n_live=500,
        double enlargement=0.5, int feedback_interval=10, int max_modes=10,
        bool feedback=false, int seed=1025, bool dump_points=false);

    virtual std::unique_ptr<Minimizer> clone() const override {
        return std::make_unique<MultiNestGPU>(*this);
    }

    std::string get_name();

    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

private:
    int n_live_, max_modes_, feedback_, seed_;
    bool ins_, mode_separation_, const_eff_, mmodal_;
    double enlargement_;

    // Various unimportant variables for the minimization routine in fortran.
    int feedback_interval_, resume_, init_mpi_;
    double logZero_;
    int writefiles_;

};

#endif
