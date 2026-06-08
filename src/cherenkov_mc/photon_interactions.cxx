// SPDX-License-Identifier: MIT
#include "cherenkov_mc/photon_interactions.h"

#include <cmath>
#include <stdexcept>

#include <TF1.h>

#include <mist/hep/globals.h>

#include "cherenkov_mc/constants.h"

namespace cherenkov_mc::photon_interactions
{
    TF1 &f_transmission()
    {
        static TF1 instance = []
        {
            TF1 f("fTransmission", "[0]*TMath::Exp(-([2])/(TMath::Power(1000*x,4)))*TMath::Exp(-([1])/(TMath::Power(1000*x,8)))", 1e-9, 1e6);
            f.SetParName(0, "A");
            f.SetParName(1, "Bt");
            f.SetParName(2, "Ct");
            f.SetParLimits(0, 0., 1);
            return f;
        }();
        return instance;
    }

    TF1 &f_scattering()
    {
        static TF1 instance = []
        {
            TF1 f("fScattering", "1/([0]*TMath::Power(x,-4)+[1]*TMath::Power(x,-8))", 1e-9, 1e6);
            f.SetParName(0, "A");
            f.SetParName(1, "B");
            return f;
        }();
        return instance;
    }

    TF1 &f_absorption()
    {
        static TF1 instance = []
        {
            TF1 f("fAbsorption", "1/([0]+[3]/TMath::Power((x/1000),[4])+[1]*TMath::Exp(-[2]*(x)))", 1e-9, 1e6);
            f.SetParName(0, "A");
            f.SetParName(1, "B");
            f.SetParName(2, "C");
            f.SetParName(3, "D");
            f.SetParName(4, "E");
            return f;
        }();
        return instance;
    }

    std::array<double, 2> get_step_size(double absorption_coefficient, double scattering_coefficient)
    {
        double step_size = -std::log(mist::hep::rng().Uniform(0., 1.)) / (absorption_coefficient + scattering_coefficient);
        double time_step = step_size / _c0;
        return {step_size, time_step};
    }

    int get_interaction(double absorption_coefficient, double scattering_coefficient)
    {
        //  Check if the coefficients are valid
        if (absorption_coefficient < 0 || scattering_coefficient < 0)
            throw std::invalid_argument("[photon_interactions::get_interaction] Absorption and scattering coefficients must be non-negative.");

        // Calculate fraction of probability
        double absorption_fraction = (absorption_coefficient) / (absorption_coefficient + scattering_coefficient);
        double scattering_fraction = 1. - absorption_fraction;

        //  Return interaction type
        // 0. Absorption
        // 1. Scattering
        return mist::hep::rng().Uniform(0., 1.) <= scattering_fraction;
    }

    physics_vector scattering_interaction(physics_vector target_momentum)
    {
        return scattering_Rayleigh(target_momentum);
    }

    physics_vector scattering_Rayleigh(physics_vector target_momentum)
    {
        auto new_momentum = target_momentum;                                                  // Work on a copy
        auto energy = new_momentum.get_t();                                                   // Preserve the energy component
        new_momentum.set_t(0);                                                                // Rotate the spatial part only
        auto orth_base = new_momentum.get_spatial_orth_base();                                // Get orthogonal base
        double scattered_theta = std::acos(std::sqrt(mist::hep::rng().Uniform(0., 1.)));      // Rayleigh scattering: cos^2(theta) distribution
        double scattered_phi = 2.0 * M_PI * mist::hep::rng().Uniform(0., 1.);                 // Uniform azimuthal angle

        new_momentum.rotation_on_an_axis(scattered_theta, orth_base[1]);
        new_momentum.rotation_on_an_axis(scattered_phi, orth_base[0]);
        new_momentum.set_t(energy);                                                           // Restore the energy component

        return new_momentum;
    }

    physics_vector reflection(physics_vector target_momentum, physics_vector surface_normal, std::array<float, 2> refractive_index)
    {
        //  Copy momenta
        auto new_surface_normal = surface_normal.get_direction(); // Ensure the normal is a unit vector

        //  Calculate the reflected momentum
        auto new_momentum = target_momentum - new_surface_normal * 2 * (target_momentum * new_surface_normal); // Reflection formula

        //  Return result
        return new_momentum;
    }

} // namespace cherenkov_mc::photon_interactions
