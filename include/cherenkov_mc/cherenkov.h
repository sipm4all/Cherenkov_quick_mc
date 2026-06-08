// SPDX-License-Identifier: MIT
//
// cherenkov_mc/cherenkov.h — Cherenkov emission helpers: emission angle,
// kinematic thresholds, and the (simplified) Frank-Tamm photon yield.
//
// Declarations live here; the non-trivial definitions live in
// src/cherenkov_mc/cherenkov.cxx. Trivial inline helpers stay in the header.
//
#pragma once

#include <TF1.h>
#include <TMath.h>

#include "cherenkov_mc/constants.h"
#include "cherenkov_mc/dynamics.h"

namespace cherenkov_mc::cherenkov
{
    //  Cherenkov
    //  --- Get Theta_c
    float get_ch_theta(float beta, float ref_index);

    //  --- Get thresholds
    inline float get_ch_thr_relativistic_beta(float ref_index) { return 1.f / ref_index; }
    inline float get_ch_thr_momentum(float ref_index, float mass) { return mass * get_ch_thr_relativistic_beta(ref_index) * dynamics::get_relativistic_gamma(get_ch_thr_relativistic_beta(ref_index)); }
    inline float get_ch_thr_mass(float ref_index, float momentum) { return momentum / (get_ch_thr_relativistic_beta(ref_index) * dynamics::get_relativistic_gamma(get_ch_thr_relativistic_beta(ref_index))); }

    //  --- Photons production
    //  Var: [0] beta, Par: [0] lambda (nm), [1] ref_index;
    //  Res: n_photons (1 / cm)
    double simp_frank_tamm(double *beta, double *parameters);

    //  Integral of the simplified Frank-Tamm yield between two wavelengths.
    inline double simp_frank_tamm_int(double lambda_min, double lambda_max, double beta, double ref_index)
    {
        return (1. / (lambda_min * 1.e-9) - 1. / (lambda_max * 1.e-9)) * 1.e-2 * 2 * TMath::Pi() * _a_EM * (1. - 1. / (beta * beta * ref_index * ref_index));
    }

    //  Process-wide TF1 wrapping simp_frank_tamm (lazy first-use construction,
    //  one instance per process — replaces the former in-header global).
    TF1 &f_simp_frank_tamm();

} // namespace cherenkov_mc::cherenkov
