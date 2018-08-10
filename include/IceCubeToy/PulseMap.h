#ifndef PULSEMAP_H_INCLUDED
#define PULSEMAP_H_INCLUDED

#include <algorithm>
#include <random>
#include "helper/abbreviations.h"
#include "IceCubeToy/DOM.h"

class PulseMap {
public:
    PulseMap(int seed=1025);

    void add_DOM(DOM &dom, v_d &charge, v_d &time);
    bool get_next(DOM &dom, v_d &charge, v_d &time);
    void add_noise(double t_min = -5.0, double t_max = 10.0);

    std::mt19937 intgen;

private:
    std::vector<DOM> DOMS_;
    m_d charges_;
    m_d times_;
    uint32_t index_ = 0;
};

#endif
