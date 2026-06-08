// tester_photon_interactions.cxx — optical interaction models.
#include "cherenkov_mc/photon_interactions.h"
#include "cherenkov_mc/physics_vector.h"

#include <cmath>
#include <cstdio>

#include <mist/hep/globals.h>

namespace pi = cherenkov_mc::photon_interactions;
using cherenkov_mc::physics_vector;

namespace {

int failures = 0;

void check(bool cond, const char *what)
{
    if (!cond) { std::printf("  FAIL: %s\n", what); ++failures; }
}

void check_close(double got, double want, double tol, const char *what)
{
    if (std::fabs(got - want) > tol)
    {
        std::printf("  FAIL: %s — got %.6g, want %.6g (tol %.3g)\n", what, got, want, tol);
        ++failures;
    }
}

} // namespace

int main()
{
    mist::hep::rng().SetSeed(2024);

    std::puts("[tester_photon_interactions] Rayleigh scatter preserves energy and |p|");
    {
        physics_vector p(2.5, 0.3, -0.4, 1.2); // (E, px, py, pz)
        const double p0 = p.get_r();
        for (int i = 0; i < 1000; ++i)
        {
            auto s = pi::scattering_Rayleigh(p);
            check_close(s.get_t(), 2.5, 1e-9, "energy (t) preserved");
            check_close(s.get_r(), p0, 1e-9, "|p| preserved");
        }
    }

    std::puts("[tester_photon_interactions] interaction selection fraction");
    {
        const double a = 1.0, s = 3.0; // P(scatter) = s/(a+s) = 0.75
        int scat = 0;
        const int N = 40000;
        for (int i = 0; i < N; ++i)
            if (pi::get_interaction(a, s) == pi::SCATTERING)
                ++scat;
        check_close(double(scat) / N, 0.75, 0.01, "P(scatter) = s/(a+s)");
    }

    std::puts("[tester_photon_interactions] step size is exponential with mean 1/(a+s)");
    {
        const double a = 1.0, s = 3.0; // mean free path = 1/(a+s) = 0.25
        double sum = 0;
        const int N = 40000;
        for (int i = 0; i < N; ++i)
        {
            auto step = pi::get_step_size(a, s);
            check(step[0] > 0, "step length positive");
            sum += step[0];
        }
        check_close(sum / N, 0.25, 0.01, "mean step = 1/(a+s)");
    }

    if (failures)
        std::printf("[tester_photon_interactions] %d FAILED\n", failures);
    else
        std::puts("[tester_photon_interactions] all passed");
    return failures ? 1 : 0;
}
