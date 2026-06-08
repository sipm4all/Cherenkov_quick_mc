// SPDX-License-Identifier: MIT
//
// cherenkov_mc/particle.h — a relativistic particle carrying an
// energy-momentum four-vector (E, px, py, pz) and a space-time position
// (t, x, y, z), with helpers for kinematics and straight-line propagation.
//
// ROOT-free. Definitions live in src/cherenkov_mc/particle.cxx.
//
#pragma once

#include <array>

#include "cherenkov_mc/physics_vector.h"

namespace cherenkov_mc
{
    class particle
    {
    private:
        physics_vector momentum; // (E, px, py, pz) GeV/c
        physics_vector position; // (t, x, y, z)    s, cm
        int PDG_id;
        int event_id;
        double mass;
        int sign;      // 1 for particle, -1 for antiparticle
        double charge; // Electrical charge in units of e
        bool absorbed; // Flag for absorption

    public:
        // Constructor
        particle(double px = 0, double py = 0, double pz = 0, double E = 0, double x = 0, double y = 0, double z = 0, double t = 0, double m = 0, int s = 1, double q = 0)
            : momentum{E, px, py, pz}, position{t, x, y, z}, PDG_id(0), event_id(0), mass(m), sign(s), charge(q), absorbed(false) {}

        // Getters
        physics_vector &get_momentum() { return momentum; }
        double get_px() const { return momentum.get_x(); }
        double get_py() const { return momentum.get_y(); }
        double get_pz() const { return momentum.get_z(); }
        double get_E() const { return momentum.get_t(); }
        double get_betax() const { return get_px() / get_E(); }
        double get_betay() const { return get_py() / get_E(); }
        double get_betaz() const { return get_pz() / get_E(); }
        double get_p() const;
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
        int get_event_id() const { return event_id; }
        bool get_absorbed() const { return absorbed; }

        // Setters
        void set_px(double px) { momentum[1] = px; }
        void set_py(double py) { momentum[2] = py; }
        void set_pz(double pz) { momentum[3] = pz; }
        void set_E(double E) { momentum[0] = E; }
        void set_x(double x) { position.set_x(x); }
        void set_y(double y) { position.set_y(y); }
        void set_z(double z) { position.set_z(z); }
        void set_t(double t) { position.set_t(t); }
        void set_position(physics_vector new_position) { position = new_position; }
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
        void set_absorbed(bool abs) { absorbed = abs; }

        // Methods
        void energy_shift(double energy_loss);
        void step(double time_step);
    };

} // namespace cherenkov_mc
