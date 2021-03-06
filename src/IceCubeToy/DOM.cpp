/* All files in this folder are based on Martin Leuermann's python script
 * for a toy model of the IceCube neutrino experiment.
 * Reference:
 * Martin Leuermann, May 2017
 * II. Physikalisches Institut RWTH Aachen
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */

#include "IceCubeToy/DOM.h"

DOM::DOM(
    double x,
    double y,
    double z,
    double noise_rate) {

    pos_.push_back(x);
    pos_.push_back(y);
    pos_.push_back(z);
    noise_rate_ = noise_rate;
}

/** Getter for the position of the DOM.
 *
 *  \return         pos_
 * */
v_d DOM::get_pos() {
    return pos_;
}

/** Getter for the noise rate of the DOM.
 *
 *  \return         noise_rate
 * */
double DOM::get_noise() {
    return noise_rate_;
}
