// tester_propagation.cxx — end-to-end sanity check of the propagator.
//
// Regression guard: before the world-boundary fix the propagator emitted
// ZERO photons (it stopped the electron before it entered any radiator). This
// test asserts a physically sensible non-zero yield.
#include "cherenkov_mc/fast_mc.h"
#include "cherenkov_mc/constants.h"

#include <cstdio>

#include <mist/hep/globals.h>

namespace {

int failures = 0;

void check(bool cond, const char *what)
{
    if (!cond)
    {
        std::printf("  FAIL: %s\n", what);
        ++failures;
    }
}

} // namespace

int main()
{
    using namespace cherenkov_mc;

    mist::hep::rng().SetSeed(12345); // deterministic

    fast_mc_geometry geom;
    auto *world = new box(physics_vector(0, 0, 0, 2), std::nullopt,
                          {{"side_x", 10.}, {"side_y", 10.}, {"side_z", 10.}});
    geom.add_world_container(world);

    ch_radiator r1(physics_vector(0, 0, 0, 1), 2, 5, 5, 1.0210, 1. / 6., 1. / 100.);
    r1.set_absorption_length(0.00753252, 10, 2e6, 1.49469e-05, 8.2902);
    r1.set_scattering_length(7.26886e+09, 0);
    geom.add_radiator(r1);

    ch_radiator r2(physics_vector(0, 0, 0, 3), 2, 5, 5, 1.0207, 1. / 6., 1. / 100.);
    r2.set_absorption_length(0.0136385, 10, 2e6, 1.29517e-05, 7.93124);
    r2.set_scattering_length(7.62127e+09, 0);
    geom.add_radiator(r2);

    const int n_events = 20;
    long total_photons = 0, total_absorbed = 0;
    for (int i = 0; i < n_events; ++i)
    {
        fast_mc_event event;
        event.set_geometry(geom);
        propagator prop(event);

        particle electron;
        electron.set_mass(kElectronMass);
        electron.set_momentum(0, 0, std::sqrt(11.5 * 11.5 - kElectronMass * kElectronMass));
        event.add_particle(0, electron);

        prop.run_event();

        for (auto &[id, p] : event.get_particles())
            if (p.get_mass() == 0)
            {
                ++total_photons;
                if (p.get_absorbed())
                    ++total_absorbed;
            }
    }

    const double mean_photons = double(total_photons) / n_events;
    const double absorbed_frac = total_photons ? double(total_absorbed) / total_photons : 0.0;
    std::printf("[tester_propagation] mean photons/event = %.1f, absorbed frac = %.3f\n",
                mean_photons, absorbed_frac);

    check(total_photons > 0, "propagator emits photons (regression: was 0)");
    check(mean_photons > 50.0, "11.5 GeV e- in aerogel yields > 50 photons/event");
    check(absorbed_frac >= 0.0 && absorbed_frac < 0.9, "absorbed fraction in a sane range");

    if (failures)
        std::printf("[tester_propagation] %d FAILED\n", failures);
    else
        std::puts("[tester_propagation] all passed");
    return failures ? 1 : 0;
}
