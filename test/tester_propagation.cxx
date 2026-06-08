// tester_propagation.cxx — end-to-end propagation: yield validation + parentage.
//
// Guards:
//  * Regression: before the world-boundary fix the propagator emitted ZERO
//    photons (it stopped the electron before it entered any radiator).
//  * Validation: the mean emitted yield must agree with the analytic
//    Frank-Tamm integral (sum over radiators of thickness * dN/dx).
//  * Parentage: every emitted photon is registered as a daughter of the
//    primary electron (regression for the dead event_id / relationships path).
#include "cherenkov_mc/cherenkov.h"
#include "cherenkov_mc/constants.h"
#include "cherenkov_mc/fast_mc.h"

#include <cstdio>

#include <mist/hep/globals.h>

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
    using namespace cherenkov_mc;

    mist::hep::rng().SetSeed(12345); // deterministic

    // Two aerogel radiators; a tight world (z in [0,4]) so photons exit
    // promptly and the test stays fast.
    const double thickness = 2.0;
    const double n1 = 1.0210, n2 = 1.0207;

    fast_mc_geometry geom;
    auto *world = new box(physics_vector(0, 0, 0, 2), std::nullopt,
                          {{"side_x", 10.}, {"side_y", 10.}, {"side_z", 4.}});
    geom.add_world_container(world);

    ch_radiator r1(physics_vector(0, 0, 0, 1), thickness, 5, 5, n1, 1. / 6., 1. / 100.);
    r1.set_absorption_length(0.00753252, 10, 2e6, 1.49469e-05, 8.2902);
    r1.set_scattering_length(7.26886e+09, 0);
    geom.add_radiator(r1);

    ch_radiator r2(physics_vector(0, 0, 0, 3), thickness, 5, 5, n2, 1. / 6., 1. / 100.);
    r2.set_absorption_length(0.0136385, 10, 2e6, 1.29517e-05, 7.93124);
    r2.set_scattering_length(7.62127e+09, 0);
    geom.add_radiator(r2);

    // Reference electron to read its beta (E = 11.5 GeV).
    particle ref;
    ref.set_mass(kElectronMass);
    ref.set_momentum(0, 0, std::sqrt(11.5 * 11.5 - kElectronMass * kElectronMass));
    const double beta = ref.get_beta();

    // Analytic expectation: thickness[cm] * dN/dx[1/cm] summed over radiators.
    const double expected =
        thickness * cherenkov::simp_frank_tamm_int(280, 850, beta, n1) +
        thickness * cherenkov::simp_frank_tamm_int(280, 850, beta, n2);

    const int n_events = 60;
    long total_photons = 0, total_absorbed = 0, total_daughters_of_primary = 0;
    for (int i = 0; i < n_events; ++i)
    {
        fast_mc_event event;
        event.set_geometry(geom);
        propagator prop(event);

        particle electron;
        electron.set_mass(kElectronMass);
        electron.set_momentum(0, 0, std::sqrt(11.5 * 11.5 - kElectronMass * kElectronMass));
        const int electron_id = event.add_particle(0, electron); // primary; id == 1

        prop.run_event();

        long photons_this_event = 0;
        for (auto &[id, p] : event.get_particles())
            if (p.get_mass() == 0)
            {
                ++total_photons;
                ++photons_this_event;
                if (p.get_absorbed())
                    ++total_absorbed;
            }

        // Parentage: photons are daughters of the primary electron.
        auto &rel = event.get_relationships();
        if (rel.count(electron_id))
            total_daughters_of_primary += rel[electron_id].size();
        check(rel.count(electron_id) && rel[electron_id].size() == size_t(photons_this_event),
              "all photons registered as daughters of the primary");
    }

    const double mean_photons = double(total_photons) / n_events;
    const double absorbed_frac = total_photons ? double(total_absorbed) / total_photons : 0.0;
    std::printf("[tester_propagation] mean photons/event = %.1f (analytic %.1f), absorbed frac = %.3f\n",
                mean_photons, expected, absorbed_frac);

    check(total_photons > 0, "propagator emits photons (regression: was 0)");
    check(total_daughters_of_primary == total_photons, "parentage count matches photon count");
    // Statistical tolerance: ~sqrt(expected/n_events) stderr on the mean; 5% is comfortable.
    check_close(mean_photons, expected, 0.05 * expected, "mean yield matches Frank-Tamm");
    check(absorbed_frac >= 0.0 && absorbed_frac < 0.9, "absorbed fraction in a sane range");

    if (failures)
        std::printf("[tester_propagation] %d FAILED\n", failures);
    else
        std::puts("[tester_propagation] all passed");
    return failures ? 1 : 0;
}
