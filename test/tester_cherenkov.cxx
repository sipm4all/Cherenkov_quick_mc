// tester_cherenkov.cxx — Cherenkov emission physics (angle, thresholds, yield).
//
// Validation strategy: each quantity is checked against an independent
// closed-form expression, not against the library's own helper.
#include "cherenkov_mc/cherenkov.h"

#include <cmath>
#include <cstdio>

namespace ck = cherenkov_mc::cherenkov;

namespace {

int failures = 0;

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
    const double n = 1.0210;
    const double beta = 0.99999;

    std::puts("[tester_cherenkov] emission angle == acos(1/(beta n))");
    check_close(ck::get_ch_theta(beta, n), std::acos(1.0 / (beta * n)), 1e-6, "theta_c");

    std::puts("[tester_cherenkov] thresholds");
    // Threshold beta is 1/n; at threshold the emission angle collapses to 0.
    check_close(ck::get_ch_thr_relativistic_beta(n), 1.0 / n, 1e-6, "beta_thr = 1/n");
    check_close(ck::get_ch_theta(1.0 / n + 1e-6, n), 0.0, 5e-3, "theta -> 0 at threshold");

    std::puts("[tester_cherenkov] Frank-Tamm yield integral vs trapezoidal sum");
    // simp_frank_tamm_int is the analytic integral of dN/dlambda over
    // [lambda_min, lambda_max]; compare to a fine trapezoidal sum of the
    // per-wavelength point function (independent numeric integration).
    const double lo = 280, hi = 850;
    const int N = 200000;
    const double d = (hi - lo) / N;
    double sum = 0;
    double pars0[2] = {lo, n};
    double prev = ([&]{ double x = beta; return ck::simp_frank_tamm(&x, pars0); })();
    for (int i = 1; i <= N; ++i)
    {
        double lambda = lo + i * d;
        double pars[2] = {lambda, n};
        double x = beta;
        double cur = ck::simp_frank_tamm(&x, pars);
        sum += 0.5 * (prev + cur) * d; // trapezoid, per nm
        prev = cur;
    }
    // simp_frank_tamm returns 1/cm per nm (1e-11 prefactor); integrate over nm.
    const double analytic = ck::simp_frank_tamm_int(lo, hi, beta, n);
    check_close(sum, analytic, analytic * 1e-3, "trapezoid == analytic integral");

    if (failures)
        std::printf("[tester_cherenkov] %d FAILED\n", failures);
    else
        std::puts("[tester_cherenkov] all passed");
    return failures ? 1 : 0;
}
