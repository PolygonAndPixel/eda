/* All files in this folder are based on Martin Leuermann's python script
 * for a toy model of the IceCube neutrino experiment.
 * Reference:
 * Martin Leuermann, May 2017
 * II. Physikalisches Institut RWTH Aachen
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */

#include "IceCubeToy/PulseMap.h"

PulseMap::PulseMap(
    int seed) {

    intgen.seed(seed);
    index = 0;
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

    DOMS.push_back(dom);
    charges.push_back(charge);
    times.push_back(time);
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
    DOM &dom,
    v_d &charge,
    v_d &time) {

    if(index == DOMS.size()) {
        index = 0;
        return false;
    }
    dom = DOMS[index];
    charge = charges[index];
    time = times[index];
    index++;
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

    DOM dom;
    v_d dom_times;
    v_d dom_charges;
    std::uniform_real_distribution<double> uf(t_min, t_max);

    while(get_next(dom, dom_times, dom_charges)) {
        std::poisson_distribution<uint32_t> pf(dom.get_noise());
        uint32_t n = pf(intgen);
        uint32_t cur = 0;
        while(cur < n) {
            double t = uf(intgen);
            dom_charges.push_back(1.0);
            dom_times.push_back(t);
            cur++;
        }
    }
}
