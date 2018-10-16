#ifndef TRACK_H_INCLUDED
#define TRACK_H_INCLUDED

#include <math.h>
#include <random>

#include "helper/abbreviations.h"
#include "IceCubeToy/DOM.h"
#include "IceCubeToy/ESource.h"
#include "IceCubeToy/IceCube_helper.h"

class Track {
public:
    Track(double x, double y, double z, double t, double theta, double phi,
        double length, double seg_length, int seed=1025);

    // void fill_PulseMap(PulseMap &void_map);

    ESource get_source(uint32_t i) {return sources_[i];};
    bool get_next_source(ESource &source);

    double minIo(double L) {return 0.222*L;}; // GeV

    std::mt19937 intgen;

private:
    double seg_length_;
    std::vector<ESource> sources_;
    double sum_length_;
    double nom_length_;
    uint32_t index_;
};

#endif
