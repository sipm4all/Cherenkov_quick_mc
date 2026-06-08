// SPDX-License-Identifier: MIT
#include "cherenkov_mc/cherenkov.h"

#include <stdexcept>

#include <TF1.h>
#include <TMath.h>

#include "cherenkov_mc/constants.h"

namespace cherenkov_mc::cherenkov
{
    float get_ch_theta(float beta, float ref_index)
    {
        // Ensure beta is within the valid range [0, 1]
        if (beta < 0.0f || beta > 1.0f)
            throw std::invalid_argument("[cherenkov::get_ch_theta] Beta must be in the range [0, 1].");

        float value = 1.0f / (beta * ref_index);
        // Ensure the input to TMath::ACos is within the valid range [-1, 1]
        if (value < -1.0f || value > 1.0f)
            throw std::domain_error("Invalid input for TMath::ACos: value out of range.");

        return TMath::ACos(value);
    }

    double simp_frank_tamm(double *beta, double *parameters)
    {
        // Ensure beta is within the valid range [0, 1]
        auto _beta = beta[0];
        if (_beta < 0.0f || _beta > 1.0f)
            throw std::invalid_argument("[cherenkov::simp_frank_tamm] Beta must be in the range [0, 1].");

        // Ensure ref_index is within the valid range (0, +inf)
        auto ref_index = parameters[1];
        if (ref_index <= 0.f)
            throw std::domain_error("[cherenkov::simp_frank_tamm] Invalid refractive index, must be greater than zero.");

        // Ensure lambda is within the valid range (0, +inf)
        // Convert lambda from nm to m
        auto lambda = parameters[0] * 1.e-9;
        if (lambda <= 0.f)
            throw std::domain_error("[cherenkov::simp_frank_tamm] Invalid photon wavelength, must be greater than zero.");

        // Frank-Tamm formula. The leading factor converts to units of 1/cm:
        //   2 * pi * alpha_EM * (1 / lambda^2) * (1 - 1 / (beta^2 n^2))
        auto result = 1.e-11 * 2 * TMath::Pi() * _a_EM * (1. / (lambda * lambda)) * (1. - 1. / (_beta * _beta * ref_index * ref_index));
        return result;
    }

    TF1 &f_simp_frank_tamm()
    {
        static TF1 instance("fsimp_frank_tamm", simp_frank_tamm, 0, 10, 2);
        return instance;
    }

} // namespace cherenkov_mc::cherenkov
