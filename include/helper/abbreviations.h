#ifndef ABBR_H_INCLUDED
#define ABBR_H_INCLUDED
// Some abbreviations

#include <math.h>
#include <vector>
#include <boost/cstdint.hpp>
#include <string>

#define SDIV(x,y)(((x)+(y)-1)/(y))
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
// Different methods
const std::string MAPSNAME      = "maps";
const std::string POLY          = "poly";
const std::string SAMPLE        = "sample";
const std::string MULTI         = "multi";
// Different functions
const std::string HIMMEL        = "himmelblau";
const std::string TOWN          = "townsend";
const std::string ROSEN         = "rosenbrock";
const std::string EGG           = "eggholder";
const std::string GAUSS         = "gaussian_shell";
const std::string ICECUBE       = "icecube";
const std::string ALL           = "all";

const double EULER_CONST = exp(1.0);

using v_d = std::vector<double>;
using v_i = std::vector<uint32_t>;
using m_d = std::vector<v_d>;

#endif
