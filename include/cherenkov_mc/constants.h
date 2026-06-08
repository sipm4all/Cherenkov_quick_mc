// SPDX-License-Identifier: MIT
//
// cherenkov_mc/constants.h — physical constants and particle masses.
//
// Header-only, ROOT-free: a flat set of inline-constexpr values living in the
// cherenkov_mc namespace. inline constexpr gives one definition across every
// translation unit (no ODR risk) while remaining usable in constant
// expressions.
//
#pragma once

namespace cherenkov_mc
{
    //  General physical constants
    inline constexpr float _c0 = 299792458.f, speed_of_light = 299792458.f;                      // [m / s]
    inline constexpr float _e = 1.60217663e-19f, elementary_charge = 1.60217663e-19f;            // [C]
    inline constexpr float _h = 6.62607015e-34f, planck_const = 6.62607015e-34f;                 // [m^2 kg / s]
    inline constexpr float _mu0 = 8.854e-12f, vacuum_mag_permittivity = 8.854e-12f;              // [F / m]
    inline constexpr float _a_EM = 1.f / 137.035999084f, fine_struct_const = 1.f / 137.035999084f; // [n.u.]

    //  Particle masses
    inline constexpr float kElectronMass = 0.0005110f; // [GeV / c^2]
    inline constexpr float kPhotonMass = 0.f;          // [GeV / c^2]
    inline constexpr float kMuonMass = 0.1057f;        // [GeV / c^2]
    inline constexpr float kPionMass = 0.1396f;        // [GeV / c^2]
    inline constexpr float kKaonMass = 0.4937f;        // [GeV / c^2]
    inline constexpr float kProtonMass = 0.9383f;      // [GeV / c^2]

} // namespace cherenkov_mc
