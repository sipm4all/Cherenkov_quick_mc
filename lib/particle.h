//  last update: 02/04/2025
//  author: Nicola Rubini - nicola.rubini@bo.infn.it
//  Scope: General dynamics and Cherenkov specifcs kinematic helper fuynctions
//  TODO: -

#pragma once

#include "constants.h"
#include "physics_vector.h"

class particle
{
private:
    physics_vector momentum; // (E, px, py, pz)
    physics_vector position; // (t, x, y, z)
    int PDG_id;
    int event_id;
    double mass;
    int sign;      // 1 for particle, -1 for antiparticle
    double charge; // Electrical charge in units of e

public:
    // Constructor
    particle(double px = 0, double py = 0, double pz = 0, double E = 0, double x = 0, double y = 0, double z = 0, double t = 0, double m = 0, int s = 1, double q = 0)
        : momentum{E, px, py, pz}, position{t, x, y, z}, mass(m), sign(s), charge(q) {}

    // Getters
    physics_vector &get_momentum() { return momentum; }
    double get_px() const { return momentum.get_x(); }
    double get_py() const { return momentum.get_y(); }
    double get_pz() const { return momentum.get_z(); }
    double get_E() const { return momentum.get_t(); }
    double get_betax() const { return get_px() / get_E(); }
    double get_betay() const { return get_py() / get_E(); }
    double get_betaz() const { return get_pz() / get_E(); }
    double get_p() const { return std::sqrt(get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz()); }
    double get_beta() const { return get_p() / get_E(); }
    int get_PDG_id() const { return PDG_id; }
    double get_mass() const { return mass; }
    int get_sign() const { return sign; }
    double get_charge() const { return charge; }
    physics_vector &get_position() { return position; }
    double get_x() const { return position.get_x(); }
    double get_y() const { return position.get_y(); }
    double get_z() const { return position.get_z(); }
    double get_t() const { return position.get_t(); }
    double get_event_id() const { return event_id; }

    // Setters
    void set_px(double px) { momentum[1] = px; }
    void set_py(double py) { momentum[2] = py; }
    void set_pz(double pz) { momentum[3] = pz; }
    void set_E(double E) { momentum[0] = E; }
    void set_momentum(double px, double py, double pz);
    void set_momentum(double E, double px, double py, double pz);
    void set_momentum(const std::array<double, 3> &direction, double magnitude);
    void set_momentum(physics_vector momentum) { set_momentum(momentum.get_t(), momentum.get_x(), momentum.get_y(), momentum.get_z()); }
    void set_position(double t, double x, double y, double z) { position = {t, x, y, z}; }
    void set_PDG_id(int id) { PDG_id = id; }
    void set_mass(double m) { mass = m; }
    void set_sign(int s) { sign = s; }
    void set_charge(double q) { charge = q; }
    void set_event_id(int id) { event_id = id; }

    // Methods
    void energy_shift(double energy_loss);
    void step(double time_step);
};

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
    
}