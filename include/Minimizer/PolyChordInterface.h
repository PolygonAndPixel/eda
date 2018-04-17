#ifdef __INTEL_COMPILER 			// if the PolyChord library was compiled with ifort
       #define POLYRUN poly_mp_polychord_c_interface_
#elif defined __GNUC__ 				// if the PolyChord library was compiled with gfortran
       #define POLYRUN __poly_MOD_polychord_c_interface
#else
       #error Do not know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) & edit example_eggbox_C++/eggbox.cc
#endif

#ifndef POLYCHORD_H
#define POLYCHORD_H

#ifdef __cplusplus

/***************************************** C++ Interface to PolyCHord *********/

#include <cstring>

namespace interfaces
{
	// module interfaces, function polychord_c_interface maps to
    // interfaces::polychord_c_interface

    extern "C" {
        void polychord_c_interface(
            int nlive, int num_repeats, int nprior, bool do_clustering,
            int feedback, double precision_criterion, int max_ndead,
            double boost_posterior, bool posteriors, bool equals,
            bool cluster_posteriors, bool write_resume,
            bool write_paramnames, bool read_resume, bool write_stats,
            bool write_live, bool write_dead, bool write_prior,
            int update_files, int nDims, int nDerived, char *base_dir,
            char *file_root, int nGrade, double *grade_frac, int *grade_dims,
            double (*Loglike)(double *theta, int nDims,
                    double *phi, int nDerived, void *),
            void (*prior)(double *cube, double *theta, int nDims, void *),
            void (*c_dumper)(double log_evidence, double error_log_evidence,
                double ndead, double n_likelihood_calls, double n_accepted,
                double *live_params,
                int n_cluster, double llh_best_fit, double llh_worst_fit,
                int nPar, void *),
            void *context);
    }
}
/******************************************************************************/

#else // ifdef __cplusplus

/***************************************** C Interface to PolyCHord ***********/

extern void polychord_c_interface(int, int, int,
    bool, int, double, int, double, bool, bool, bool, bool,
    bool, bool, bool, bool, bool, bool, int, int, int,
    char *, char *, int, double *, int *,
    double (*Loglike)(double *, int, double *, int, void *),
    void (*prior)(double *, double *, int, void *),
    void (*c_dumper)(double, double, double, double, double, double *, int,
        double, double, int, void *), void *context);

/******************************************************************************/

#endif // ifdef __cplusplus

#endif // POLYCHORD_H
