#ifdef __INTEL_COMPILER 			// if the MultiNest library was compiled with ifort
       #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__ 				// if the MultiNest library was compiled with gfortran
       #define NESTRUN __nested_MOD_nestrun
#else
       #error Do not know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) & edit example_eggbox_C++/eggbox.cc
#endif

#ifndef MULTINEST_H
#define MULTINEST_H

#ifdef __cplusplus

/***************************************** C++ Interface to MultiNest **************************************************/

#include <cstring>

namespace nested
{

	// map the Fortran 90 entry points of libnest3.a to C++ functions

	// module nested, function nestRun maps to nested::run

	// the pass-by-reference nature of most of the Fortran is translated away
	// *apart* from the callbacks. The provided call back functions must still accept
	// references rather than values. There is also some confusion as to the type
	// of the first argument of LogLike.
	// Should it be a value_t * or an farray<value_t, 1> *? The former seems to
	// work and is simpler.

	// This structure is reverse engineered from looking
	// at gfortran stack traces. It is probably wrong

	template<typename type, index_t ndims> class farray_traits;

	template<> class farray_traits<value_t, 1> { public: static const index_t id = 537; };
	template<> class farray_traits<value_t, 2> { public: static const index_t id = 538; };
	template<> class farray_traits<index_t, 1> { public: static const index_t id = 265; };
	template<> class farray_traits<index_t, 2> { public: static const index_t id = 266; };

	// the extra data for f90 that defines how arrays are arranged.
	template<typename T, index_t ndim> class farray
	{
		public:
			farray(T *_data, index_t w, index_t h = 0) : data(_data), offset(0), type(farray_traits<T, ndim>::id),
			x_stride(1), x_lbound(1), x_ubound(w), y_stride(w), y_lbound(1), y_ubound(h) {};

			T *data;
			index_t offset;
			index_t type;
			index_t x_stride, x_lbound, x_ubound;
			index_t y_stride, y_lbound, y_ubound;
	};

	extern "C" {
		void NESTRUN(index_t &IS, index_t &mmodal, index_t &ceff, index_t &nlive, value_t &tol, value_t &efr, index_t &ndims,
			index_t &nPar, index_t &nClsPar, index_t &maxModes, index_t &updInt, value_t &Ztol, char *root, index_t &seed,
			index_t *pWrap, index_t &fb, index_t &resume, index_t &outfile, index_t &initMPI, value_t &logZero, index_t &maxiter,
			void (*Loglike)(value_t *Cube, index_t &n_dim, index_t &n_par, value_t &lnew, void *),
			void (*dumper)(index_t &, index_t &, index_t &, value_t **, value_t **, value_t **, value_t &, value_t &, value_t &, value_t &, index_t &, void *),
			void *context, index_t &root_len);
	}

	static void run(bool IS, bool mmodal, bool ceff, index_t nlive, value_t tol,
        value_t efr, index_t ndims, index_t nPar,
        index_t nClsPar,
        index_t maxModes, index_t updInt,
        value_t Ztol, const std::string & root, index_t seed,
        index_t *pWrap, bool fb, bool resume,
        bool outfile, bool initMPI, value_t logZero, index_t maxiter,
        void (*LogLike)(value_t *Cube, index_t &n_dim, index_t &n_par, value_t &lnew, void *),
		void (*dumper)(index_t &, index_t &, index_t &, value_t **, value_t **, value_t **, value_t &, value_t &, value_t &, value_t &, index_t &, void *), void *context)
	{
		char t_root[100];
		std::fill(t_root, t_root + 100, ' ');
		snprintf(t_root, 99, "%s", root.c_str());
		index_t root_len = strlen(t_root);
		t_root[strlen(t_root)] = ' ';

		index_t t_fb = fb;
		index_t t_resume = resume;
		index_t t_outfile = outfile;
		index_t t_initMPI = initMPI;
		index_t t_mmodal = mmodal;
		index_t t_IS = IS;
		index_t t_ceff = ceff;

		NESTRUN(t_IS, t_mmodal, t_ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, t_root, seed, pWrap, t_fb,
		t_resume, t_outfile, t_initMPI, logZero, maxiter, LogLike, dumper, context, root_len);
	}
}

/***********************************************************************************************************************/

#else // ifdef __cplusplus

/***************************************** C Interface to MultiNest **************************************************/

extern void NESTRUN(index_t *, index_t *, index_t *, index_t *, value_t *, value_t *, index_t *, index_t *, index_t *, index_t *, index_t *, value_t *,
char *, index_t *, index_t *, index_t *, index_t *, index_t *, index_t *, value_t *, index_t *, void (*Loglike)(value_t *, index_t *, index_t *,
value_t *, void *), void (*dumper)(index_t *, index_t *, index_t *, value_t **, value_t **, value_t **, value_t *,
value_t *, value_t *, value_t *, index_t *, void *), void *context);

void run(index_t IS, index_t mmodal, index_t ceff, index_t nlive, value_t tol, value_t efr, index_t ndims, index_t nPar, index_t nClsPar,
index_t maxModes, index_t updInt, value_t Ztol, char root[], index_t seed, index_t *pWrap, index_t fb, index_t resume, index_t outfile,
index_t initMPI, value_t logZero, index_t maxiter, void (*LogLike)(value_t *, index_t *, index_t *, value_t *, void *),
void (*dumper)(index_t *, index_t *, index_t *, value_t **, value_t **, value_t **, value_t *, value_t *, value_t *, value_t *, index_t *, void *),
void *context)
{
	index_t i;
	for (i = strlen(root); i < 100; i++) root[i] = ' ';

        NESTRUN(&IS, &mmodal, &ceff, &nlive, &tol, &efr, &ndims, &nPar, &nClsPar, &maxModes, &updInt, &Ztol,
        root, &seed, pWrap, &fb, &resume, &outfile, &initMPI, &logZero, &maxiter, LogLike, dumper, context);
}

/***********************************************************************************************************************/

#endif // ifdef __cplusplus

#endif // MULTINEST_H
