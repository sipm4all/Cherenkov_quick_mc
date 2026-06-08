// SPDX-License-Identifier: MIT
//
// cherenkov_mc/physics_vector.h — a four-component (t, x, y, z) vector with a
// configurable Minkowski metric, used for both space-time positions and
// energy-momentum four-vectors.
//
// ROOT-free. Definitions live in src/cherenkov_mc/physics_vector.cxx.
//
#pragma once

#include <array>

namespace cherenkov_mc
{
    class physics_vector
    {
    private:
        std::array<double, 4> _vector;
        std::array<int, 4> _metric;

    public:
        //  Constructors
        physics_vector();
        physics_vector(double t, double x, double y, double z);

        // Overloaded operators
        double &operator[](int index) { return _vector[index]; }
        physics_vector &operator*=(double scale);
        physics_vector &operator/=(double scale) { return *this *= (1. / scale); }
        physics_vector &operator+=(const physics_vector &other);
        physics_vector &operator-=(const physics_vector &other);
        physics_vector operator+(const physics_vector &other);
        physics_vector operator-(const physics_vector &other);
        double operator*(const physics_vector &other);
        physics_vector operator*(double scale);

        // Getter methods
        double get_x() const { return _vector[1]; }
        double get_y() const { return _vector[2]; }
        double get_z() const { return _vector[3]; }
        double get_t() const { return _vector[0]; }
        double get_r() const;
        double get_theta() const;
        double get_phi() const;
        double get_norm() const;
        physics_vector get_direction() const;
        std::array<physics_vector, 3> get_spatial_orth_base();
        std::array<int, 4> get_metric() { return _metric; }

        // Setter methods
        void set_t_x_y_z(double t, double x, double y, double z);
        void set_x_y_z(double x, double y, double z);
        void set_x(double x) { _vector[1] = x; }
        void set_y(double y) { _vector[2] = y; }
        void set_z(double z) { _vector[3] = z; }
        void set_t(double t) { _vector[0] = t; }
        void set_r_theta_phi(double r, double theta, double phi);
        void set_r(double r) { set_r_theta_phi(r, get_theta(), get_phi()); }
        void set_theta(double theta) { set_r_theta_phi(get_r(), theta, get_phi()); }
        void set_phi(double phi) { set_r_theta_phi(get_r(), get_theta(), phi); }
        void set_metric(const std::array<int, 4> &new_metric) { _metric = new_metric; }

        // Vector combination methods
        physics_vector spatial_cross_product(const physics_vector target_2);

        // Vector shift methods
        void rotation_on_an_axis(double angle, physics_vector along_the_vector);

        //  Utility
        void print() const;
    };

} // namespace cherenkov_mc
