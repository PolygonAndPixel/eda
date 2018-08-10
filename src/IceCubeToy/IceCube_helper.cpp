#include "IceCubeToy/IceCube_helper.h"

double dist(
    const v_d &a,
    const v_d &b) {

    double distance = 0.0;
    for(uint32_t i=0; i<a.size(); i++) {
        distance += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(distance);
}

double hit_charge(
    double delta_r) {

    double scale = 2e3;
    double abs0 = 100;
    double geo = 1.0/(delta_r*delta_r);
    double absorption = exp(-delta_r/abs0);
    return scale*geo*absorption;
}

double hit_time(
    double delta_r,
    double delta_time) {

    double scat0 = 10;
    double a = delta_r/scat0;
    double x = delta_time;
    double gamma = tgamma(a);
    return ((pow(x, a-1) * exp(-x))/gamma);
}

double hit_prob(
    double delta_r,
    double delta_time) {

    double norm = hit_charge(delta_r);
    double gamma = hit_time(delta_r, delta_time - (delta_r/0.3));
    return norm*gamma;
}
