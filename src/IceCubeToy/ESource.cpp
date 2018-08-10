/* All files in this folder are based on Martin Leuermann's python script
 * for a toy model of the IceCube neutrino experiment.
 * Reference:
 * Martin Leuermann, May 2017
 * II. Physikalisches Institut RWTH Aachen
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */

#include "IceCubeToy/ESource.h"

ESource::ESource(
    double x,
    double y,
    double z,
    double t,
    double theta,
    double phi,
    double en) {

    pos_.push_back(x);
    pos_.push_back(y);
    pos_.push_back(z);
    time_ = t;
    direction_.push_back(theta);
    direction_.push_back(phi);
    energy_ = en;
}

/** Getter for the position of the energy source.
 *
 *  \return         pos_
 * */
v_d ESource::get_pos() {
    return pos_;
}

/** Getter for the time of the energy source.
 *
 *  \return         time_
 * */
double ESource::get_time() {
    return time_;
}

/** Getter for the direction of the energy source.
 *
 *  \return         direction_ (theta, phi)
 * */
v_d ESource::get_direction() {
    return direction_;
}

/** Getter for the energy of the energy source.
 *
 *  \return         energy_
 * */
double ESource::get_energy() {
    return energy_;
}
