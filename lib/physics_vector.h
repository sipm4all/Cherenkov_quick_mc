//  last update: 02/04/2025
//  author: Nicola Rubini - nicola.rubini@bo.infn.it
//  Scope: Cherenkov-related helper functions
//  TODO: -

#pragma once

class physics_vector
{
private:
    std::array<double, 4> _vector;
    std::array<int, 4> _metric;

public:
    //  Contructor
    physics_vector()
    {
        _vector.fill(0.0f);
        _metric = {-1, 1, 1, 1};
    }
    physics_vector(double t, double x, double y, double z)
    {
        _vector = {t, x, y, z};
        _metric = {-1, 1, 1, 1};
    }

    // Getter methods
    double get_x() const { return _vector[1]; }
    double get_y() const { return _vector[2]; }
    double get_z() const { return _vector[3]; }
    double get_t() const { return _vector[0]; }
    double get_r() const { return std::sqrt(get_x() * get_x() + get_y() * get_y() + get_z() * get_z()); }
    double get_theta() const { return std::acos(get_z() / get_r()); }
    double get_phi() const { return std::atan2(get_y(), get_x()); }
    physics_vector get_direction() { return *this * (1. / (get_r())); }
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

    // Overloaded operators
    double &operator[](int index) { return _vector[index]; }
    physics_vector &operator*=(double scale);
    physics_vector &operator*=(const physics_vector &other);
    physics_vector &operator/=(double scale) { return *this *= (1. / scale); }
    physics_vector &operator+=(const physics_vector &other);
    physics_vector &operator-=(const physics_vector &other);
    physics_vector operator+(const physics_vector &other);
    physics_vector operator-(const physics_vector &other);
    double operator*(const physics_vector &other);
    physics_vector operator*(double scale);

    // Vecotr combination methods
    physics_vector spatial_cross_product(const physics_vector target_2) { return physics_vector(0., (*this).get_y() * target_2.get_z() - (*this).get_z() * target_2.get_y(), (*this).get_z() * target_2.get_x() - (*this).get_x() * target_2.get_z(), (*this).get_x() * target_2.get_y() - (*this).get_y() * target_2.get_x()); }
    // Vector shift methods
    void rotation_on_an_axis(double angle, physics_vector along_the_vector);
    physics_vector shift_r_theta_phi(double r, double theta, double phi, physics_vector along_the_base = {0.707107f, 0.707107f, 0.707107f, 0});

    //  Utility
    void print() { cout << "[physics_vector::print] t: " << get_t() << " x: " << get_x() << " y: " << get_y() << " z: " << get_z() << endl; }
};

physics_vector &physics_vector::operator*=(double scale)
{
    for (int i = 0; i < 4; ++i)
        _vector[i] *= scale;
    return *this;
}
physics_vector &physics_vector::operator*=(const physics_vector &other)
{
    for (int i = 0; i < 4; ++i)
        _vector[i] *= _metric[i] * other._vector[i];
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
std::array<physics_vector, 3> physics_vector::get_spatial_orth_base()
{
    //  Define result
    std::array<physics_vector, 3> result;

    //  First base versor is the vector itself
    result[0] = (*this).get_direction();
    result[0].set_t(0);

    //  Define a second perpendicular versor
    if ((fabs(result[0].get_x()) < 1.e-9) && (fabs(result[0].get_y()) < 1.e-9))
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
    // Based on the Rodrigues' rule
    // https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    // Prepare the result vector
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