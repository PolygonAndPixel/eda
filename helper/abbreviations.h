#ifndef ABBR_H_INCLUDED
#define ABBR_H_INCLUDED
// Some abbreviations

#include <math.h>
#include <vector>
#include <boost/cstdint.hpp>

#define SDIV(x,y)(((x)+(y)-1)/(y))
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

const double EULER_CONST = exp(1.0);

using v_d = std::vector<double>;
using v_i = std::vector<uint32_t>;
using m_d = std::vector<v_d>;

#endif