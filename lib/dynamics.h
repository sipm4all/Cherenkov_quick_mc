//  last update: 02/04/2025
//  author: Nicola Rubini - nicola.rubini@bo.infn.it
//  Scope: General dynamics and Cherenkov specifcs kinematic helper fuynctions
//  TODO: -

#pragma once

#include "constants.h"

namespace dynamics
{
    //  Kinematics
    //  --- Get Beta
    inline float get_relativistic_beta(float gamma) { return TMath::Sqrt(1 - gamma * gamma); }
    inline float get_relativistic_beta(float mass, float momentum) { return 1. / sqrt(1 + (mass * mass) / (momentum * momentum)); }
    //  --- Get Gamma
    inline float get_relativistic_gamma(float beta) { return 1. / (TMath::Sqrt(1 - beta * beta)); }
    inline float get_relativistic_gamma(float mass, float momentum) { return get_relativistic_gamma(get_relativistic_beta(mass, momentum)); }
}