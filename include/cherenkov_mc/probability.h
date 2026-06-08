// SPDX-License-Identifier: MIT
//
// cherenkov_mc/probability.h — sampling helpers.
//
// Randomness is drawn from the process-wide mist::hep::rng() (a TRandom3),
// so the simulation shares one reproducible RNG stream with the rest of the
// mist-hep ecosystem instead of owning a private generator.
//
#pragma once

namespace cherenkov_mc::probability
{
    //  --- Poisson
    int poisson_sampling_Knuth(double lambda);    // low-lambda method
    int poisson_sampling_Atkinson(double lambda); // efficient across all lambda values

} // namespace cherenkov_mc::probability
