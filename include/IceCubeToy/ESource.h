#ifndef ESource_H_INCLUDED
#define ESource_H_INCLUDED

#include "helper/abbreviations.h"

class ESource {
public:
    ESource(){};
    ESource(double x, double y, double z, double t, double theta,
            double phi, double en);

    v_d get_pos();
    double get_time();
    v_d get_direction();
    double get_energy();

private:
    v_d pos_;
    double time_;
    v_d direction_;
    double energy_;
};

#endif
