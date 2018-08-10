#ifndef DOM_H_INCLUDED
#define DOM_H_INCLUDED

#include "helper/abbreviations.h"

class DOM {
public:
    DOM(){};
    DOM(double x, double y, double z, double noise_rate);

    v_d get_pos();
    double get_noise();

private:
    v_d pos_;
    double noise_rate_;
};

#endif
