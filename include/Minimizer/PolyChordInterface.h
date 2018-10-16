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
            index_t nlive, index_t num_repeats, index_t nprior, bool do_clustering,
            index_t feedback, value_t precision_criterion, index_t max_ndead,
            value_t boost_posterior, bool posteriors, bool equals,
            bool cluster_posteriors, bool write_resume,
            bool write_paramnames, bool read_resume, bool write_stats,
            bool write_live, bool write_dead, bool write_prior,
            index_t update_files, index_t nDims, index_t nDerived, char *base_dir,
            char *file_root, index_t nGrade, value_t *grade_frac, index_t *grade_dims,
            value_t (*Loglike)(value_t *theta, index_t nDims,
                    value_t *phi, index_t nDerived, void *),
            void (*prior)(value_t *cube, value_t *theta, index_t nDims, void *),
            void (*c_dumper)(value_t log_evidence, value_t error_log_evidence,
                value_t ndead, value_t n_likelihood_calls, value_t n_accepted,
                value_t *live_params,
                index_t n_cluster, value_t llh_best_fit, value_t llh_worst_fit,
                index_t nPar, void *),
            void *context);
    }
}
/******************************************************************************/

#else // ifdef __cplusplus

/***************************************** C Interface to PolyCHord ***********/

extern void polychord_c_interface(index_t, index_t, index_t,
    bool, index_t, value_t, index_t, value_t, bool, bool, bool, bool,
    bool, bool, bool, bool, bool, bool, index_t, index_t, index_t,
    char *, char *, index_t, value_t *, index_t *,
    value_t (*Loglike)(value_t *, index_t, value_t *, index_t, void *),
    void (*prior)(value_t *, value_t *, index_t, void *),
    void (*c_dumper)(value_t, value_t, value_t, value_t, value_t, value_t *, index_t,
        value_t, value_t, index_t, void *), void *context);

/******************************************************************************/

#endif // ifdef __cplusplus

#endif // POLYCHORD_H
