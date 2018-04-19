/* A test function which is used for neutrino reconstruction at IceCube.
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */

#include "likelihood/NeutrinoFunction.h"
#include <iostream>

NeutrinoFunction::NeutrinoFunction() {

}

NeutrinoFunction::NeutrinoFunction(
    std::string func_name,
    uint32_t ndims) {

    name = func_name;
}

/** Getter for the name of the used function.
 *
 *  \return         name
 * */
std::string NeutrinoFunction::get_name() {
    return name;
}

/** Getter for the dimension.
 *
 *  \return         ndims_
 * */
uint32_t NeutrinoFunction::get_ndims() {

    return ndims_;
}

/** Do some magic to get a value.
 *
 *  \param theta    Physical parameters of point that shall be evaluated.
 *
 *  \return         Likelihood
 * */
double NeutrinoFunction::get_lh(
    v_d theta) {

    EventHypothesis &h = *(get_hypothesis_ptr(theta));
    VectorParticlePtr sources = extract_hypothesis(h);
    get_response_matrix();
    fit_statistics();

    return 0;
}

/**
 *
 *  \return
 * */
EventHypothesis NeutrinoFunction::get_hypothesis_ptr(
    v_d &theta) {

    EventHypothesis h;
    update_physics_variables();
    return h;
}

/**
 *
 *  \return
 * */
VectorParticlePtr NeutrinoFunction::extract_hypothesis(
    const EventHypothesis &h) {

    VectorParticlePtr sources;
    return sources;
}

/**
 *
 *  \return
 * */
NeutrinoFunction::get_response_matrix() {

    MillipedeAddOMSourcePairToMatrix();
}

/**
 *
 *  \return
 * */
NeutrinoFunction::MillipedeAddOMSourcePairToMatrix() {

    get_probability_quantiles();
}

/**
 *
 *  \return
 * */
NeutrinoFunction::get_probability_quantiles() {
}


/**
 *
 *  \return
 * */
NeutrinoFunction::update_physics_variables() {

}

/**
 *
 *  \return
 * */
NeutrinoFunction::fit_statistics() {

}
