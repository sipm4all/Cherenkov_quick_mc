// SPDX-License-Identifier: MIT
//
// cherenkov_mc/dynamics.h — relativistic kinematic helpers (beta / gamma).
//
// Header-only: all functions are `inline`, so they may be included in any
// number of translation units without ODR conflict.
//
#pragma once

#include <cmath>

#include <TMath.h>

#include "cherenkov_mc/constants.h"

namespace cherenkov_mc::dynamics
{
    //  Kinematics
    //  --- Get Beta
    inline float get_relativistic_beta(float gamma) { return TMath::Sqrt(1 - 1.f / (gamma * gamma)); }
    inline float get_relativistic_beta(float mass, float momentum) { return 1.f / std::sqrt(1 + (mass * mass) / (momentum * momentum)); }
    //  --- Get Gamma
    inline float get_relativistic_gamma(float beta) { return 1.f / (TMath::Sqrt(1 - beta * beta)); }
    inline float get_relativistic_gamma(float mass, float momentum) { return get_relativistic_gamma(get_relativistic_beta(mass, momentum)); }

} // namespace cherenkov_mc::dynamics
