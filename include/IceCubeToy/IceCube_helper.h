#ifndef ICECUBEHELPER_H_INCLUDED
#define ICECUBEHELPER_H_INCLUDED

#include <math.h>
#include <limits>
#include "helper/abbreviations.h"

extern double dist(const v_d &a, const v_d &b);
extern double hit_charge(double delta_r);
extern double hit_time(double delta_r, double delta_time);
extern double hit_prob(double delta_r, double delta_time);

#endif
