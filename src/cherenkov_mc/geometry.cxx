// SPDX-License-Identifier: MIT
#include "cherenkov_mc/geometry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "cherenkov_mc/constants.h"

namespace cherenkov_mc
{
    // === geometry ===

    double geometry::get_parameter(const std::string &key) const
    {
        if (parameters.find(key) != parameters.end())
            return parameters.at(key);
        throw std::invalid_argument("Parameter key not found.");
    }

    // === box ===

    box::box(physics_vector p, std::optional<std::vector<physics_vector>> ax, std::map<std::string, double> params)
        : geometry(p, ax, params)
    {
        // Set default parameters if not provided
        if (parameters.find("side_x") == parameters.end())
            parameters["side_x"] = 1.0;
        if (parameters.find("side_y") == parameters.end())
            parameters["side_y"] = 1.0;
        if (parameters.find("side_z") == parameters.end())
            parameters["side_z"] = 1.0;
    }

    bool box::is_inside(physics_vector &point) const
    {
        double side_x = get_parameter("side_x");
        double side_y = get_parameter("side_y");
        double side_z = get_parameter("side_z");

        return (std::fabs(point.get_x() - position.get_x()) <= side_x * 0.5 &&
                std::fabs(point.get_y() - position.get_y()) <= side_y * 0.5 &&
                std::fabs(point.get_z() - position.get_z()) <= side_z * 0.5);
    }

    std::optional<std::array<double, 2>> box::get_max_step_size(particle target_particle) const
    {
        // Ray-box (slab method) intersection: find the time/distance at which the
        // particle, travelling along its current velocity, exits the box.
        const double half[3] = {get_parameter("side_x") * 0.5,
                                get_parameter("side_y") * 0.5,
                                get_parameter("side_z") * 0.5};
        const double centre[3] = {position.get_x(), position.get_y(), position.get_z()};
        const double origin[3] = {target_particle.get_x(), target_particle.get_y(), target_particle.get_z()};
        // Velocity in cm / s (matches particle::step convention)
        const double velocity[3] = {target_particle.get_betax() * _c0 * 100,
                                    target_particle.get_betay() * _c0 * 100,
                                    target_particle.get_betaz() * _c0 * 100};

        double t_enter = -std::numeric_limits<double>::infinity();
        double t_exit = std::numeric_limits<double>::infinity();
        for (int i = 0; i < 3; ++i)
        {
            const double lo = centre[i] - half[i];
            const double hi = centre[i] + half[i];
            if (std::fabs(velocity[i]) < 1.e-30)
            {
                // Travelling parallel to this pair of planes: must already be between them
                if (origin[i] < lo || origin[i] > hi)
                    return std::nullopt;
                continue;
            }
            double t1 = (lo - origin[i]) / velocity[i];
            double t2 = (hi - origin[i]) / velocity[i];
            if (t1 > t2)
                std::swap(t1, t2);
            t_enter = std::max(t_enter, t1);
            t_exit = std::min(t_exit, t2);
        }

        // No forward intersection
        if (t_exit < std::max(t_enter, 0.0))
            return std::nullopt;

        const double speed = target_particle.get_beta() * _c0 * 100;
        return std::array<double, 2>{t_exit * speed, t_exit};
    }

    physics_vector box::get_surface_norm(particle target_particle) const
    {
        // Outward normal of the face through which the particle exits the box.
        auto step = get_max_step_size(target_particle);
        if (!step)
            return physics_vector(0, 0, 0, 0);

        const double t_exit = (*step)[1];
        const double half[3] = {get_parameter("side_x") * 0.5,
                                get_parameter("side_y") * 0.5,
                                get_parameter("side_z") * 0.5};
        const double centre[3] = {position.get_x(), position.get_y(), position.get_z()};
        const double origin[3] = {target_particle.get_x(), target_particle.get_y(), target_particle.get_z()};
        const double velocity[3] = {target_particle.get_betax() * _c0 * 100,
                                    target_particle.get_betay() * _c0 * 100,
                                    target_particle.get_betaz() * _c0 * 100};

        // The exit face is the axis whose surface is reached exactly at t_exit
        int axis = 0;
        double best = std::numeric_limits<double>::infinity();
        double sign = 1.0;
        for (int i = 0; i < 3; ++i)
        {
            const double exit_coord = origin[i] + velocity[i] * t_exit;
            const double d = std::fabs(std::fabs(exit_coord - centre[i]) - half[i]);
            if (d < best)
            {
                best = d;
                axis = i;
                sign = (exit_coord - centre[i] >= 0) ? 1.0 : -1.0;
            }
        }

        physics_vector normal(0, 0, 0, 0);
        if (axis == 0)
            normal.set_x(sign);
        else if (axis == 1)
            normal.set_y(sign);
        else
            normal.set_z(sign);
        return normal;
    }

    std::string box::get_type() const { return "Box"; }

    // === cylinder ===

    bool cylinder::is_inside(physics_vector &point) const
    {
        double radius = get_parameter("radius");
        double side_x = get_parameter("side_x");

        physics_vector d = point - position;

        if (axes)
        {
            double projection_side_x = d * axes.value()[0];

            if (std::fabs(projection_side_x) > side_x * 0.5)
                return false;

            double px = d.get_x() - projection_side_x * axes.value()[0].get_x();
            double py = d.get_y() - projection_side_x * axes.value()[0].get_y();
            double pz = d.get_z() - projection_side_x * axes.value()[0].get_z();

            double dist2 = px * px + py * py + pz * pz;

            return dist2 <= radius * radius;
        }

        return false;
    }

    std::optional<std::array<double, 2>> cylinder::get_max_step_size(particle target_particle) const
    {
        // Geometry is solved with a unit direction so the parameter is directly a
        // path length (cm); the time is then length / speed.
        auto cylinder_axis = (axes.value()[0]).get_direction();
        cylinder_axis.set_t(0);
        auto cylinder_center = position;
        cylinder_center.set_t(0);
        auto direction = target_particle.get_momentum().get_direction();
        direction.set_t(0);
        auto particle_position = target_particle.get_position();
        particle_position.set_t(0);
        const double radius = get_parameter("radius");
        const double half_length = get_parameter("side_x") * 0.5; // axial half-extent (see is_inside)

        // Split position offset and direction into axial / radial parts
        auto offset = particle_position - cylinder_center;
        auto offset_parallel = cylinder_axis * (offset * cylinder_axis);
        auto direction_parallel = cylinder_axis * (direction * cylinder_axis);
        auto offset_perp = offset - offset_parallel;
        auto direction_perp = direction - direction_parallel;

        const double eps = 1.e-12;
        double exit = std::numeric_limits<double>::infinity();

        // --- Barrel surface: |offset_perp + s * direction_perp|^2 = radius^2 ---
        const double A = direction_perp * direction_perp;
        const double B = 2.0 * (direction_perp * offset_perp);
        const double C = offset_perp * offset_perp - radius * radius;
        if (A > eps)
        {
            const double disc = B * B - 4 * A * C;
            if (disc >= 0)
            {
                const double sqrt_disc = std::sqrt(disc);
                for (double s : {(-B + sqrt_disc) / (2 * A), (-B - sqrt_disc) / (2 * A)})
                {
                    if (s <= eps)
                        continue;
                    const double axial = (offset * cylinder_axis) + s * (direction * cylinder_axis);
                    if (std::fabs(axial) <= half_length && s < exit)
                        exit = s;
                }
            }
        }

        // --- End caps: planes at axial = +/- half_length ---
        const double dir_axial = direction * cylinder_axis;
        if (std::fabs(dir_axial) > eps)
        {
            const double pos_axial = offset * cylinder_axis;
            for (double cap : {half_length, -half_length})
            {
                const double s = (cap - pos_axial) / dir_axial;
                if (s <= eps || s >= exit)
                    continue;
                auto hit = offset + direction * s;
                auto hit_parallel = cylinder_axis * (hit * cylinder_axis);
                auto hit_perp = hit - hit_parallel;
                if (hit_perp * hit_perp <= radius * radius)
                    exit = s;
            }
        }

        if (!std::isfinite(exit))
            return std::nullopt;

        const double speed = target_particle.get_beta() * _c0 * 100;
        return std::array<double, 2>{exit, speed > 0 ? exit / speed : 0.0};
    }

    physics_vector cylinder::get_surface_norm(particle target_particle) const
    {
        auto step = get_max_step_size(target_particle);
        if (!step)
            return physics_vector(0, 0, 0, 0);

        auto cylinder_axis = (axes.value()[0]).get_direction();
        cylinder_axis.set_t(0);
        auto cylinder_center = position;
        cylinder_center.set_t(0);
        auto direction = target_particle.get_momentum().get_direction();
        direction.set_t(0);
        auto particle_position = target_particle.get_position();
        particle_position.set_t(0);
        const double half_length = get_parameter("side_x") * 0.5;

        const double s = (*step)[0];
        auto offset = (particle_position - cylinder_center) + direction * s;
        const double axial = offset * cylinder_axis;

        // On an end cap if the exit point sits at the axial extent, else on the barrel
        if (std::fabs(std::fabs(axial) - half_length) < 1.e-6)
            return cylinder_axis * (axial >= 0 ? 1.0 : -1.0);

        auto radial = offset - cylinder_axis * axial;
        return radial.get_direction();
    }

    std::string cylinder::get_type() const { return "Cylinder"; }

} // namespace cherenkov_mc
