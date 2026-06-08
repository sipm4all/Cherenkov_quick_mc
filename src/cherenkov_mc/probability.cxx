// SPDX-License-Identifier: MIT
#include "cherenkov_mc/probability.h"

#include <cmath>

#include <TMath.h>

#include <mist/hep/globals.h>

namespace cherenkov_mc::probability
{
    int poisson_sampling_Knuth(double lambda)
    {
        int samples = 0;
        double target_probability = 1.0;
        double limit_probability = std::exp(-lambda);
        do
        {
            samples++;
            target_probability *= mist::hep::rng().Uniform(0., 1.);
        } while (target_probability > limit_probability);
        return samples - 1;
    }

    int poisson_sampling_Atkinson(double lambda)
    {
        if (lambda < 30)
            return poisson_sampling_Knuth(lambda); // Use Knuth for small lambda

        double c = 0.767 - 3.36 / lambda; // Empirical parameters
        double beta = M_PI / std::sqrt(3.0 * lambda);
        double alpha = beta * lambda;
        double k = std::log(c) - lambda - std::log(beta);

        while (true)
        {
            double u, x, n;
            u = mist::hep::rng().Uniform(0., 1.);
            x = (alpha - std::log((1.0 - u) / u)) / beta;
            n = std::floor(x + 0.5);

            if (n < 0)
                continue;

            double v = mist::hep::rng().Uniform(0., 1.);
            double y = alpha - beta * x;
            double lhs = y + std::log(v / (1.0 + std::exp(y)) * (1.0 + std::exp(y)));
            double rhs = k + n * std::log(lambda) - TMath::LnGamma(n + 1);

            if (lhs <= rhs)
                return n;
        }
    }

} // namespace cherenkov_mc::probability
