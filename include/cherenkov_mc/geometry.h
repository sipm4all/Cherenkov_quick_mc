// SPDX-License-Identifier: MIT
//
// cherenkov_mc/geometry.h — abstract volume interface plus box and cylinder
// implementations used as world containers / radiator envelopes.
//
// ROOT-free. Definitions live in src/cherenkov_mc/geometry.cxx.
//
// TODO: introduce rotation and translation of the geometry, with matching
//       is_inside / step-size handling.
//
#pragma once

#include <array>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "cherenkov_mc/particle.h"
#include "cherenkov_mc/physics_vector.h"

namespace cherenkov_mc
{
    // === geometry class ===

    class geometry
    {
    protected:
        physics_vector position;
        std::optional<std::vector<physics_vector>> axes;
        std::map<std::string, double> parameters;

    public:
        geometry(physics_vector pos = physics_vector(0, 0, 0, 0),
                 std::optional<std::vector<physics_vector>> ax = std::nullopt,
                 std::map<std::string, double> params = {})
            : position(pos), axes(ax), parameters(params) {}

        // Getters
        physics_vector get_position() const { return position; }
        std::optional<std::vector<physics_vector>> get_axes() const { return axes; }
        double get_parameter(const std::string &key) const;
        virtual std::string get_type() const = 0;

        // Setters
        void set_position(const physics_vector &pos) { position = pos; }
        void set_axes(const std::optional<std::vector<physics_vector>> &ax) { axes = ax; }
        void set_parameter(const std::string &key, double value) { parameters[key] = value; }

        virtual bool is_inside(physics_vector &point) const = 0;
        virtual std::optional<std::array<double, 2>> get_max_step_size(particle target_particle) const = 0;
        virtual physics_vector get_surface_norm(particle target_particle) const = 0;

        virtual ~geometry() = default;
    };

    // === box class ===

    class box : public geometry
    {
    public:
        box(physics_vector p = physics_vector(0, 0, 0, 0),
            std::optional<std::vector<physics_vector>> ax = std::nullopt,
            std::map<std::string, double> params = {});

        bool is_inside(physics_vector &point) const override;
        std::optional<std::array<double, 2>> get_max_step_size(particle target_particle) const override;
        physics_vector get_surface_norm(particle target_particle) const override;
        std::string get_type() const override;
    };

    // === cylinder class ===

    class cylinder : public geometry
    {
    public:
        cylinder(physics_vector c = physics_vector(0, 0, 0, 0),
                 std::optional<std::vector<physics_vector>> ax = std::nullopt,
                 std::map<std::string, double> params = {})
            : geometry(c, ax, params) {}

        bool is_inside(physics_vector &point) const override;
        std::optional<std::array<double, 2>> get_max_step_size(particle target_particle) const override;
        physics_vector get_surface_norm(particle target_particle) const override;
        std::string get_type() const override;

    private:
        physics_vector main_axis; // Main axis of the cylinder, check that is provided or set a default one
    };

} // namespace cherenkov_mc
