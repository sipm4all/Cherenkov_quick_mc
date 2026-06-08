// SPDX-License-Identifier: MIT
#include "cherenkov_mc/particle.h"

#include <cmath>

#include "cherenkov_mc/constants.h"

namespace cherenkov_mc
{
    double particle::get_p() const
    {
        return std::sqrt(get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz());
    }

    void particle::set_momentum(double px, double py, double pz)
    {
        set_px(px);
        set_py(py);
        set_pz(pz);
        set_E(std::sqrt(px * px + py * py + pz * pz + mass * mass));
    }
    void particle::set_momentum(double E, double px, double py, double pz)
    {
        momentum = {E, px, py, pz};
        mass = std::sqrt(E * E - (px * px + py * py + pz * pz));
    }
    void particle::set_momentum(const std::array<double, 3> &direction, double magnitude)
    {
        double norm = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
        if (norm > 0)
        {
            set_px(magnitude * direction[0] / norm);
            set_py(magnitude * direction[1] / norm);
            set_pz(magnitude * direction[2] / norm);
            set_E(std::sqrt(get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz() + mass * mass));
        }
        else
        {
            set_px(0);
            set_py(0);
            set_pz(0);
            set_E(mass);
        }
    }
    void particle::energy_shift(double energy_shift)
    {
        auto momentum = get_p();
        if ((get_E() + energy_shift) <= mass)
        {
            set_momentum(0., 0., 0.);
            return;
        }
        auto new_momentum = std::sqrt((get_E() + energy_shift) * (get_E() + energy_shift) - mass * mass);
        set_momentum(get_px() * (new_momentum) / (momentum), get_py() * (new_momentum) / (momentum), get_pz() * (new_momentum) / (momentum));
    }
    void particle::step(double time_step)
    {
        // Update the position of the particle based on its momentum and time step
        set_t(get_t() + time_step);
        set_x(get_x() + get_betax() * _c0 * 100 * time_step);
        set_y(get_y() + get_betay() * _c0 * 100 * time_step);
        set_z(get_z() + get_betaz() * _c0 * 100 * time_step);
    }

} // namespace cherenkov_mc
