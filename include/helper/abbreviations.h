#ifndef ABBR_H_INCLUDED
#define ABBR_H_INCLUDED
// Some abbreviations

#include <math.h>
#include <vector>
#include <boost/cstdint.hpp>
#include <string>

typedef value_t value_t;
typedef index_t index_t;

#define SDIV(x,y)(((x)+(y)-1)/(y))
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
// Different methods
const std::string MAPSNAME      = "maps";
const std::string POLY          = "poly";
const std::string SAMPLE        = "sample";
const std::string MULTI         = "multi";
const std::string SCAN          = "scan";
const std::string GD            = "gradient descent";
const std::string DALEX         = "dalex";
// Different functions
const std::string HIMMEL        = "himmelblau";
const std::string TOWN          = "townsend";
const std::string ROSEN         = "rosenbrock";
const std::string EGG           = "eggholder";
const std::string GAUSS         = "gaussian_shell";
const std::string ICECUBE       = "icecube";
const std::string PARABOLOID    = "paraboloid";
const std::string ALL           = "all";

const value_t EULER_CONST = exp(1.0);

using v_d = std::vector<value_t>;
using v_i = std::vector<index_t>;
using m_d = std::vector<v_d>;



#define H2D (cudaMemcpyHostToDevice)
#define D2H (cudaMemcpyDeviceToHost)
#define H2H (cudaMemcpyHostToHost)
#define D2D (cudaMemcpyDeviceToDevice)
// For the CPU
struct splinetable {
	index_t ndim;
	index_t *order;     // Order the splines (ndim many entries)

	value_t **knots;    // Knots in each dimension
	long *nknots;       // Number of knots in each dimension

	value_t **extents;
	value_t *periods;

	value_t *coefficients;  // All coefficients for each spline
	long *naxes;            // Number of splines for each dimension
	unsigned long *strides; // stride to each dimension

	index_t naux;
	char ***aux;
};
// For the GPU
struct Splinetable {
	index_t *order;     // Order the splines (ndim many entries)

	value_t *knots;     // Knots in each dimension
	long *nknots;       // Number of knots in each dimension

	value_t *coefficients;  // All coefficients for each spline
	long *naxes;            // Number of splines for each dimension
	unsigned long *strides; // stride to each dimension
};

#endif
