// SPDX-License-Identifier: MIT
#include "cherenkov_mc/physics_vector.h"

#include <cmath>
#include <string>

#include <mist/logger/logger.h>

namespace cherenkov_mc
{
    //  Constructors
    physics_vector::physics_vector()
    {
        _vector.fill(0.0);
        _metric = {-1, 1, 1, 1};
    }
    physics_vector::physics_vector(double t, double x, double y, double z)
    {
        _vector = {t, x, y, z};
        _metric = {-1, 1, 1, 1};
    }

    //  Overloaded operators
    physics_vector &physics_vector::operator*=(double scale)
    {
        for (int i = 0; i < 4; ++i)
            _vector[i] *= scale;
        return *this;
    }
    physics_vector &physics_vector::operator+=(const physics_vector &other)
    {
        for (int i = 0; i < 4; ++i)
            _vector[i] += other._vector[i];
        return *this;
    }
    physics_vector &physics_vector::operator-=(const physics_vector &other)
    {
        for (int i = 0; i < 4; ++i)
            _vector[i] -= other._vector[i];
        return *this;
    }
    physics_vector physics_vector::operator+(const physics_vector &other)
    {
        physics_vector result;
        for (int i = 0; i < 4; ++i)
            result._vector[i] = _vector[i] + other._vector[i];
        return result;
    }
    physics_vector physics_vector::operator-(const physics_vector &other)
    {
        physics_vector result;
        for (int i = 0; i < 4; ++i)
            result._vector[i] = _vector[i] - other._vector[i];
        return result;
    }
    double physics_vector::operator*(const physics_vector &other)
    {
        double result = 0.;
        for (int i = 0; i < 4; ++i)
            result += _vector[i] * _metric[i] * other._vector[i];
        return result;
    }
    physics_vector physics_vector::operator*(double scale)
    {
        physics_vector result;
        for (int i = 0; i < 4; ++i)
            result._vector[i] = _vector[i] * scale;
        return result;
    }

    //  Getters
    double physics_vector::get_r() const
    {
        return std::sqrt(get_x() * get_x() + get_y() * get_y() + get_z() * get_z());
    }
    double physics_vector::get_theta() const
    {
        return std::acos(get_z() / get_r());
    }
    double physics_vector::get_phi() const
    {
        return std::atan2(get_y(), get_x());
    }
    double physics_vector::get_norm() const
    {
        double result = 0.;
        for (int i = 0; i < 4; ++i)
            result += _vector[i] * _metric[i] * _vector[i];
        return std::sqrt(result);
    }
    physics_vector physics_vector::get_direction() const
    {
        auto new_vector = *this;
        new_vector.set_t(0);
        return new_vector * (1. / get_r());
    }

    //  Setters
    void physics_vector::set_t_x_y_z(double t, double x, double y, double z)
    {
        set_t(t);
        set_x(x);
        set_y(y);
        set_z(z);
    }
    void physics_vector::set_x_y_z(double x, double y, double z)
    {
        set_x(x);
        set_y(y);
        set_z(z);
    }
    void physics_vector::set_r_theta_phi(double r, double theta, double phi)
    {
        set_x(r * std::sin(theta) * std::cos(phi));
        set_y(r * std::sin(theta) * std::sin(phi));
        set_z(r * std::cos(theta));
    }

    //  Vector combination methods
    physics_vector physics_vector::spatial_cross_product(const physics_vector target_2)
    {
        return physics_vector(0.,
                              get_y() * target_2.get_z() - get_z() * target_2.get_y(),
                              get_z() * target_2.get_x() - get_x() * target_2.get_z(),
                              get_x() * target_2.get_y() - get_y() * target_2.get_x());
    }

    //  Vector shift methods
    std::array<physics_vector, 3> physics_vector::get_spatial_orth_base()
    {
        //  Define result
        std::array<physics_vector, 3> result;

        //  First base versor is the vector itself
        result[0] = get_direction();
        result[0].set_t(0);

        //  Define a second perpendicular versor
        if ((std::fabs(result[0].get_x()) < 1.e-9) && (std::fabs(result[0].get_y()) < 1.e-9))
            result[1] = physics_vector(0., 1., 0., 0.);
        else
            result[1] = physics_vector(0., -result[0].get_y(), result[0].get_x(), 0.);

        //  Use cross product to find last versor
        result[2] = result[1].spatial_cross_product(result[0]);

        //  Impose normalisation
        for (auto &current_vector : result)
            current_vector = current_vector.get_direction();

        // result
        return result;
    }
    void physics_vector::rotation_on_an_axis(double angle, physics_vector wrt_the_vector)
    {
        // Based on Rodrigues' rotation formula
        // https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        physics_vector target_vector = (*this);
        physics_vector result;
        // Assure the direction vector is unitary
        physics_vector rotation_axis = wrt_the_vector.get_direction();
        // Calculate
        result += target_vector * std::cos(angle);
        result += target_vector.spatial_cross_product(rotation_axis) * std::sin(angle);
        result += rotation_axis * (rotation_axis * target_vector) * (1 - std::cos(angle));
        *this = result;
    }

    //  Utility
    void physics_vector::print() const
    {
        mist::logger::debug("[physics_vector] t: " + std::to_string(get_t()) +
                            " x: " + std::to_string(get_x()) +
                            " y: " + std::to_string(get_y()) +
                            " z: " + std::to_string(get_z()));
    }

} // namespace cherenkov_mc
