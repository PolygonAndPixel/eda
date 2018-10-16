/* All files in this folder are based on Martin Leuermann's python script
 * for a toy model of the IceCube neutrino experiment.
 * Reference:
 * Martin Leuermann, May 2017
 * II. Physikalisches Institut RWTH Aachen
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */

#include "IceCubeToy/PulseMap.h"
#include <iostream>

PulseMap::PulseMap(
    int seed) {

    intgen.seed(seed);
    index_ = 0;
}

/* Add a DOM with the charges it detected and the corresponding timings.
 *
 * \param dom       The DOM to be added.
 * \param charge    The charge of the DOM.
 * \param time      The time of the charge.
 * */
void PulseMap::add_DOM(
    DOM &dom,
    v_d &charge,
    v_d &time) {

    DOMS_.push_back(dom);
    charges_.push_back(charge);
    times_.push_back(time);
}

/** Get next DOM, charges and times.
 *
 * \param dom       On out: the next DOM in the list.
 * \param charge    On out: the next charges in the list.
 * \param time      On out: the next times for the charges in the list.
 *
 * \return          False if no next element is available and resets counter.
 * */
bool PulseMap::get_next(
    DOM & dom,
    v_d & charge,
    v_d & time) {

    if(index_ == DOMS_.size()) {
        index_ = 0;
        return false;
    }

    dom =  DOMS_[index_];
    charge =  charges_[index_];
    time =  times_[index_];
    index_++;
    return true;
}

/** Add noise to all DOMs in this class.
 *
 * \param t_min     Minimum time.
 * \param t_max     Maximum time.
 * */
void PulseMap::add_noise(
    double t_min,
    double t_max) {

    // DOM dom;
    // v_d dom_times;
    // v_d dom_charges;
    double tmp_min = std::numeric_limits<double>::lowest();
    double tmp_max = std::numeric_limits<double>::max();
    
    for(v_d t: times_) {
        if(t.empty()) continue;
        double tmp = *std::min_element(t.begin(), t.end());
        if(tmp < tmp_min) tmp_min = tmp;
        tmp = *std::max_element(t.begin(), t.end());
        if(tmp > tmp_max) tmp_max = tmp;
    }

    t_min = tmp_min - t_min;
    t_max = tmp_max + t_max;
    std::uniform_real_distribution<double> uf(t_min, t_max);

    for(uint32_t i = 0; i < DOMS_.size(); i++) {
        std::poisson_distribution<uint32_t> pf(DOMS_[i].get_noise());
        uint32_t n = pf(intgen);
        uint32_t cur = 0;
        while(cur < n) {
            double t = uf(intgen);
            charges_[i].push_back(1.0);
            times_[i].push_back(t);
            cur++;
        }
    }

    // while(get_next(dom, dom_times, dom_charges)) {
    //     std::poisson_distribution<uint32_t> pf(dom.get_noise());
    //     uint32_t n = pf(intgen);
    //     uint32_t cur = 0;
    //     while(cur < n) {
    //         double t = uf(intgen);
    //         dom_charges.push_back(1.0);
    //         dom_times.push_back(t);
    //         cur++;
    //     }
    // }
}

/** Fill the pulse map with pulses from the energy sources along a given track.
 * 
 * \param track     The track
 * */
void PulseMap::fill(
    Track & track) {

    ESource source;
    std::uniform_real_distribution<double> uf2(0.0, 0.3);
    std::normal_distribution<double> nf(1.0, 0.2);
    while(track.get_next_source(source)) {
        for(uint32_t i = 0; i < DOMS_.size(); i++) {
            double delta_r = dist(source.get_pos(), DOMS_[i].get_pos());
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
                    charges_[i].push_back(nf(intgen));
                    times_[i].push_back(t_time + t_y);
                    count++;
                }
            }
        }
    }
}

