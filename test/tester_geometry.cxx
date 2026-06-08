// tester_geometry.cxx — exercises box / cylinder containment and exit geometry.
#include "cherenkov_mc/geometry.h"
#include "cherenkov_mc/particle.h"
#include "cherenkov_mc/physics_vector.h"

#include <cmath>
#include <cstdio>
#include <vector>

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

    std::puts("[tester_geometry] box is_inside");
    {
        box b(physics_vector(0, 0, 0, 0), std::nullopt, {{"side_x", 4.}, {"side_y", 4.}, {"side_z", 4.}});
        physics_vector inside(0, 0, 0, 0), outside(0, 3, 0, 0);
        check(b.is_inside(inside), "origin inside box");
        check(!b.is_inside(outside), "(3,0,0) outside 4-cube");
    }

    std::puts("[tester_geometry] box exit distance and normal");
    {
        box b(physics_vector(0, 0, 0, 0), std::nullopt, {{"side_x", 4.}, {"side_y", 4.}, {"side_z", 4.}});
        particle pho;
        pho.set_mass(0);
        pho.set_momentum(0.3, 0.0, 1.0); // unit dir = (0.287,0,0.957)
        pho.set_position(0, 0, 0, 0);
        auto step = b.get_max_step_size(pho);
        check(step.has_value(), "box step has value");
        if (step)
            check_close((*step)[0], 2.0 / 0.95783, 2e-3, "box exit distance ~2.088 cm");
        auto n = b.get_surface_norm(pho);
        check_close(n.get_z(), 1.0, 1e-9, "box exit normal +z");
    }

    std::puts("[tester_geometry] cylinder exit (caps) distance and normal");
    {
        std::vector<physics_vector> axis = {physics_vector(0, 0, 0, 1)};
        cylinder c(physics_vector(0, 0, 0, 0), axis, {{"radius", 2.0}, {"side_x", 10.0}});
        particle pho;
        pho.set_mass(0);
        pho.set_momentum(0.3, 0.0, 1.0);
        pho.set_position(0, 0, 0, 0);
        auto step = c.get_max_step_size(pho);
        check(step.has_value(), "cylinder step has value");
        if (step)
            check_close((*step)[0], 5.0 / 0.95783, 5e-3, "cylinder exits top cap ~5.22 cm");
        auto n = c.get_surface_norm(pho);
        check_close(n.get_z(), 1.0, 1e-6, "cylinder exit normal +z (cap)");
    }

    if (failures)
        std::printf("[tester_geometry] %d FAILED\n", failures);
    else
        std::puts("[tester_geometry] all passed");
    return failures ? 1 : 0;
}
