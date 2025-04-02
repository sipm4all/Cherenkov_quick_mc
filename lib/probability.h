//  last update: 02/04/2025
//  author: Nicola Rubini - nicola.rubini@bo.infn.it
//  Scope: Probability methods
//  TODO: -

#pragma once

namespace probability
{
    //  Random utility
    std::random_device rd;  // Seed for randomness
    std::mt19937 gen(rd()); // Mersenne Twister RNG
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    //  Probability
    //  --- Poisson
    int poisson_sampling_Knuth(double lambda);    // low lambda method
    int poisson_sampling_Atkinson(double lambda); // efficient accross all lambda values

};

int probability::poisson_sampling_Knuth(double lambda)
{
    int samples = 0;
    double target_probability = 1.0;
    double limit_probability = std::exp(-lambda);
    do
    {
        samples++;
        target_probability *= uniform_dist(probability::gen);
    } while (target_probability > limit_probability);
    return samples - 1;
}
int probability::poisson_sampling_Atkinson(double lambda)
{
    if (lambda < 30)
        return poisson_sampling_Knuth(lambda); // Use Knuth for small lambda

    double c = 0.767 - 3.36 / lambda; // Empirical parameters
    double beta = M_PI / sqrt(3.0 * lambda);
    double alpha = beta * lambda;
    double k = log(c) - lambda - log(beta);

    while (true)
    {
        double u, x, n;
        u = uniform_dist(gen);
        x = (alpha - log((1.0 - u) / u)) / beta;
        n = floor(x + 0.5);

        if (n < 0)
            continue;

        double v = gRandom->Uniform(0., 1.);
        double y = alpha - beta * x;
        double lhs = y + log(v / (1.0 + exp(y)) * (1.0 + exp(y)));
        double rhs = k + n * log(lambda) - TMath::LnGamma(n + 1);

        if (lhs <= rhs)
            return n;
    }
}