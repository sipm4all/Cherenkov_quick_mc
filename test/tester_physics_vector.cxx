// tester_physics_vector.cxx — exercises the four-vector algebra.
#include "cherenkov_mc/physics_vector.h"

#include <cmath>
#include <cstdio>

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
    using cherenkov_mc::physics_vector;

    std::puts("[tester_physics_vector] spatial_cross_product");
    {
        physics_vector x(0, 1, 0, 0), y(0, 0, 1, 0);
        auto z = x.spatial_cross_product(y); // x_hat × y_hat = z_hat
        check_close(z.get_x(), 0.0, 1e-12, "cross x");
        check_close(z.get_y(), 0.0, 1e-12, "cross y");
        check_close(z.get_z(), 1.0, 1e-12, "cross z");
    }

    std::puts("[tester_physics_vector] metric dot product");
    {
        physics_vector p(2, 0, 0, 1); // (E=2, pz=1) -> p·p = -E^2 + pz^2 = -3
        check_close(p * p, -3.0, 1e-12, "minkowski self-dot");
    }

    std::puts("[tester_physics_vector] rotation preserves spatial norm");
    {
        physics_vector v(0, 1, 0.5, -0.3);
        const double r0 = v.get_r();
        v.rotation_on_an_axis(0.7, physics_vector(0, 0, 0, 1)); // rotate about z
        check_close(v.get_r(), r0, 1e-9, "norm preserved under rotation");
        check_close(v.get_z(), -0.3, 1e-9, "z component fixed under z-rotation");
    }

    std::puts("[tester_physics_vector] orthonormal base");
    {
        physics_vector v(0, 1, 2, 3);
        auto base = v.get_spatial_orth_base();
        check_close(base[0] * base[1], 0.0, 1e-9, "e0 . e1 = 0");
        check_close(base[0] * base[2], 0.0, 1e-9, "e0 . e2 = 0");
        check_close(base[1] * base[2], 0.0, 1e-9, "e1 . e2 = 0");
        check_close(base[0].get_r(), 1.0, 1e-9, "|e0| = 1");
        check_close(base[1].get_r(), 1.0, 1e-9, "|e1| = 1");
        check_close(base[2].get_r(), 1.0, 1e-9, "|e2| = 1");
    }

    if (failures)
        std::printf("[tester_physics_vector] %d FAILED\n", failures);
    else
        std::puts("[tester_physics_vector] all passed");
    return failures ? 1 : 0;
}
