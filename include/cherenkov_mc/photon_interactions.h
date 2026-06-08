// SPDX-License-Identifier: MIT
//
// cherenkov_mc/photon_interactions.h — optical photon interaction models:
// absorption / scattering step sampling, the parametrised transmission,
// scattering and absorption length curves, and the scattering / reflection
// kinematics.
//
// The transmission / scattering / absorption TF1s are exposed as Meyers
// singletons (lazy, one-per-process) instead of in-header globals; their
// parameter names are configured on first use.
//
#pragma once

#include <array>

#include <TF1.h>

#include "cherenkov_mc/physics_vector.h"

namespace cherenkov_mc::photon_interactions
{
    // Constants for interaction types
    constexpr int ABSORPTION = 0;
    constexpr int SCATTERING = 1;

    //  Parametrised optical curves (lazy, process-wide; parameter names set on
    //  first use). Replaces the former in-header global TF1* objects.
    TF1 &f_transmission();
    TF1 &f_scattering();
    TF1 &f_absorption();

    //  Basic step evolution
    std::array<double, 2> get_step_size(double absorption_coefficient, double scattering_coefficient);
    int get_interaction(double absorption_coefficient, double scattering_coefficient); // 0. Absorption
                                                                                       // 1. Scattering

    //  Interactions
    physics_vector scattering_interaction(physics_vector target_momentum);
    physics_vector scattering_Rayleigh(physics_vector target_momentum);
    physics_vector reflection(physics_vector target_momentum, physics_vector surface_normal, std::array<float, 2> refractive_index);

} // namespace cherenkov_mc::photon_interactions
