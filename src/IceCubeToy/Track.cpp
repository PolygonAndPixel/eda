/* All files in this folder are based on Martin Leuermann's python script
 * for a toy model of the IceCube neutrino experiment.
 * Reference:
 * Martin Leuermann, May 2017
 * II. Physikalisches Institut RWTH Aachen
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */

#include "IceCubeToy/Track.h"

Track::Track(
    double x,
    double y,
    double z,
    double t,
    double theta,
    double phi,
    double length,
    double seg_length,
    int seed) {

    seg_length_ = seg_length;
    double current_length = 0;

    while(current_length < length) {
        double n_x = sin(phi)*sin(theta)*current_length + x;
        double n_y = cos(phi)*sin(theta)*current_length + y;
        double n_z = sin(theta)*current_length + z; // The DOMs look down...
        double n_time = current_length/0.3 + t;
        ESource source = ESource(n_x, n_y, n_z, n_time,
            theta, phi, minIo(seg_length_));
        sources.push_back(source);
        current_length += seg_length_;
    }
    sum_length = current_length;
    nom_length = length;
    intgen.seed(seed);
    index = 0;
}

/** Fill the pulse map with pulses from the energy sources along the track.
 *
 * \param void_map      The pulse map to be filled.
 * */
void Track::fill_PulseMap(
    PulseMap &void_map) {

    DOM dom;
    v_d charges;
    v_d times;
    std::uniform_real_distribution<double> uf2(0.0, 0.3);
    std::normal_distribution<double> nf(1.0, 0.2);
    for(ESource &source: sources) {
        while(void_map.get_next(dom, charges, times)) {
            double delta_r = dist(source.get_pos(), dom.get_pos());
            double abs_time = source.get_time() + delta_r/0.3;
            double mu = hit_charge(delta_r) * source.get_energy();
            std::poisson_distribution<uint32_t> pf(mu);
            uint32_t n_photons = pf(intgen);
            uint32_t count = 0;
            std::uniform_real_distribution<double> uf(0.0, 15.0*(delta_r/10.0));

            while(count < n_photons) {
                double t_time = uf(intgen);
                double t_y = uf2(intgen);
                if(hit_time(delta_r, t_time) > t_y) {
                    charges.push_back(nf(intgen));
                    times.push_back(t_time + t_y);
                    count++;
                }
            }
        }
    }
}

/** Get next Source.
 *
 * \param source       On out: the next ESource in the list.
 *
 * \return          False if no next element is available and resets counter.
 * */
bool Track::get_next_source(
    ESource &source) {

        if(index == sources.size()-1) {
            index = 0;
            return false;
        }
        index++;
        source = sources[index];
        return true;
}
