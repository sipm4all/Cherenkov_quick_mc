// tester_kinematics.cxx — particle kinematics and dynamics helpers.
#include "cherenkov_mc/constants.h"
#include "cherenkov_mc/dynamics.h"
#include "cherenkov_mc/particle.h"

#include <cmath>
#include <cstdio>

namespace dyn = cherenkov_mc::dynamics;
using cherenkov_mc::particle;

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
    std::puts("[tester_kinematics] dynamics beta/gamma round-trip");
    {
        const float beta = 0.8f;
        const float gamma = dyn::get_relativistic_gamma(beta);
        check_close(gamma, 1.0 / std::sqrt(1 - beta * beta), 1e-4, "gamma(beta)");
        check_close(dyn::get_relativistic_beta(gamma), beta, 1e-4, "beta(gamma) round-trip");
    }
    {
        const float m = cherenkov_mc::kMuonMass, p = 1.0f;
        check_close(dyn::get_relativistic_beta(m, p), 1.0 / std::sqrt(1 + (m * m) / (p * p)), 1e-4, "beta(mass,momentum)");
    }

    std::puts("[tester_kinematics] particle energy/momentum consistency");
    {
        particle e;
        e.set_mass(cherenkov_mc::kElectronMass);
        e.set_momentum(0, 0, 5.0);
        check_close(e.get_E(), std::sqrt(5.0 * 5.0 + cherenkov_mc::kElectronMass * cherenkov_mc::kElectronMass), 1e-9, "E = sqrt(p^2+m^2)");
        check_close(e.get_beta(), e.get_p() / e.get_E(), 1e-12, "beta = p/E");
    }

    std::puts("[tester_kinematics] energy_shift conserves mass, reduces |p|");
    {
        particle pi_;
        pi_.set_mass(cherenkov_mc::kPionMass);
        pi_.set_momentum(0, 0, 3.0);
        const double m0 = pi_.get_mass();
        const double p0 = pi_.get_p();
        pi_.energy_shift(-0.5); // lose 0.5 GeV
        check(pi_.get_p() < p0, "|p| decreased after energy loss");
        check_close(pi_.get_mass(), m0, 1e-4, "mass conserved across energy_shift");
    }

    std::puts("[tester_kinematics] straight-line step displacement");
    {
        particle ph;
        ph.set_mass(0);
        ph.set_momentum(0, 0, 1.0); // beta = 1 along +z
        ph.set_position(0, 0, 0, 0);
        const double dt = 1e-12;
        ph.step(dt);
        check_close(ph.get_z(), 1.0 * cherenkov_mc::_c0 * 100 * dt, 1e-9, "z = beta*c*100*dt");
        check_close(ph.get_t(), dt, 1e-18, "t advanced by dt");
    }

    if (failures)
        std::printf("[tester_kinematics] %d FAILED\n", failures);
    else
        std::puts("[tester_kinematics] all passed");
    return failures ? 1 : 0;
}
