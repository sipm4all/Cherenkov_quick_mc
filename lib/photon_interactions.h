//  last update: 02/04/2025
//  author: Nicola Rubini - nicola.rubini@bo.infn.it
//  Scope: Cherenkov-related helper functions
//  TODO: -

#pragma once

#include "constants.h"
#include "probability.h"
#include "physics_vector.h"

namespace photons_interactions
{
    //  Basic step evolution
    std::array<double, 2> get_step_size(double absorption_coefficient, double scattering_coefficient);
    int get_interaction(double absorption_coefficient, double scattering_coefficient); // 0. Absorption
                                                                                       // 1. Scattering

    //  Interactions
    float scattering_interaction();
    float scattering_Rayleigh();
    float scattering_Mie();

};

std::array<double, 2> photons_interactions::get_step_size(double absorption_coefficient, double scattering_coefficient)
{
    double step_size = -std::log(probability::uniform_dist(probability::gen)) / (absorption_coefficient + scattering_coefficient);
    double time_step = step_size / _c0;
    return {step_size, time_step};
}
int photons_interactions::get_interaction(double absorption_coefficient, double scattering_coefficient)
{
    double absorption_fraction = (absorption_coefficient) / (absorption_coefficient + scattering_coefficient);
    double scattering_fraction = 1. - absorption_fraction;
    return probability::uniform_dist(probability::gen) <= scattering_fraction;
}
float photons_interactions::scattering_interaction()
{
    return scattering_Rayleigh();
}
float photons_interactions::scattering_Rayleigh()
{
    return 0.;
}