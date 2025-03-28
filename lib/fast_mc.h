//  last update: 21/03/2025
//  author: Nicola Rubini - nicola.rubini@bo.infn.it
//  Scope: Cherenkov fast mc
//  TODO: -

#pragma once

#include "cherenkov.h"

namespace ch_fast_mc
{
    class fast_mc_event;
    class fast_mc_geometry;
    class physics_vector
    {
    private:
        std::array<float, 4> _vector;
        std::array<int, 4> _metric;

    public:
        //  Contructor
        physics_vector()
        {
            _vector.fill(0.0f);
            _metric = {1, 1, 1, -1};
        }
        physics_vector(float x, float y, float z, float t)
        {
            _vector = {x, y, z, t};
            _metric = {1, 1, 1, -1};
        }

        // Getter methods
        float get_x() const { return _vector[1]; }
        float get_y() const { return _vector[2]; }
        float get_z() const { return _vector[3]; }
        float get_t() const { return _vector[0]; }
        float get_r() const { return std::sqrt(_vector[1] * _vector[1] + _vector[2] * _vector[2] + _vector[3] * _vector[3]); }
        float get_theta() const { return std::acos(_vector[3] / get_r()); }
        float get_phi() const { return std::atan2(_vector[2], _vector[1]); }
        physics_vector get_direction() { return (*this) /= (get_r()); }
        std::array<physics_vector, 3> get_spatial_orth_base();
        std::array<int, 4> get_metric() { return _metric; }

        // Setter methods
        void set_x_y_z(float x, float y, float z);
        void set_x(float x) { _vector[1] = x; }
        void set_y(float y) { _vector[2] = y; }
        void set_z(float z) { _vector[3] = z; }
        void set_t(float t) { _vector[0] = t; }
        void set_r_theta_phi(float r, float theta, float phi);
        void set_r(float r) { set_r_theta_phi(r, get_theta(), get_phi()); }
        void set_theta(float theta) { set_r_theta_phi(get_r(), theta, get_phi()); }
        void set_phi(float phi) { set_r_theta_phi(get_r(), get_theta(), phi); }
        void set_metric(const std::array<int, 4> &new_metric) { _metric = new_metric; }

        // Overloaded operators
        float &operator[](int index) { return _vector[index]; }
        physics_vector &operator*=(double scale);
        physics_vector &operator/=(double scale) { return (*this) *= 1. / scale; }
        physics_vector &operator*=(const physics_vector &other);
        // physics_vector operator*(const physics_vector &other) const;
        // physics_vector operator/(const physics_vector &other) const;
        physics_vector operator+(const physics_vector &other) const;
        physics_vector operator-(const physics_vector &other) const;

        // Other methods
        physics_vector spatial_cross_product(physics_vector target_2) { return physics_vector((*this)[2] * target_2[3] - (*this)[3] * target_2[2], (*this)[3] * target_2[1] - (*this)[1] * target_2[3], (*this)[1] * target_2[2] - (*this)[2] * target_2[1], 0.); }
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
    physics_vector physics_vector::operator+(const physics_vector &other) const
    {
        physics_vector result;
        for (int i = 0; i < 4; ++i)
            result._vector[i] = _vector[i] + other._vector[i];
        return result;
    }
    physics_vector physics_vector::operator-(const physics_vector &other) const
    {
        physics_vector result;
        for (int i = 0; i < 4; ++i)
            result._vector[i] = _vector[i] - other._vector[i];
        return result;
    }
    void physics_vector::set_x_y_z(float x, float y, float z)
    {
        set_x(x);
        set_x(y);
        set_x(z);
    }
    void physics_vector::set_r_theta_phi(float r, float theta, float phi)
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
            result[1] = physics_vector(1., 0., 0., 0.);
        else
            result[1] = physics_vector(-result[0].get_y(), result[0].get_x(), 0., 0.);

        //  Use cross product to find last versor
        result[2] = result[1].spatial_cross_product(result[0]);

        //  Impose normalisation
        for (auto &current_vector : result)
            current_vector = current_vector.get_direction();

        // result
        return result;
    }

    class particle
    {
    private:
        physics_vector momentum; // (px, py, pz, E)
        physics_vector position; // (x, y, z, t)
        int PDG_id;
        float mass;
        int sign;     // 1 for particle, -1 for antiparticle
        float charge; // Electrical charge in units of e

    public:
        // Constructor
        particle(float px = 0, float py = 0, float pz = 0, float E = 0, float x = 0, float y = 0, float z = 0, float t = 0, float m = 0, int s = 1, float q = 0)
            : momentum{px, py, pz, E}, position{x, y, z, t}, mass(m), sign(s), charge(q) {}

        // Getters
        physics_vector &get_momentum() { return momentum; }
        float get_px() const { return momentum.get_x(); }
        float get_py() const { return momentum.get_y(); }
        float get_pz() const { return momentum.get_z(); }
        float get_betax() const { return get_px() / get_E(); }
        float get_betay() const { return get_py() / get_E(); }
        float get_betaz() const { return get_pz() / get_E(); }
        float get_p() const { return std::sqrt(get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz()); }
        float get_E() const { return momentum.get_t(); }
        float get_beta() const { return get_p() / get_E(); }
        int get_PDG_id() const { return PDG_id; }
        float get_mass() const { return mass; }
        int get_sign() const { return sign; }
        float get_charge() const { return charge; }
        physics_vector &get_position() { return position; }
        float get_x() const { return position.get_x(); }
        float get_y() const { return position.get_y(); }
        float get_z() const { return position.get_z(); }
        float get_t() const { return position.get_t(); }

        // Setters
        void set_momentum(float px, float py, float pz);
        void set_momentum(float px, float py, float pz, float E);
        void set_momentum(const std::array<float, 3> &direction, float magnitude);
        void set_position(float x, float y, float z, float t) { position = {x, y, z, t}; }
        void set_PDG_id(int id) { PDG_id = id; }
        void set_mass(float m) { mass = m; }
        void set_sign(int s) { sign = s; }
        void set_charge(float q) { charge = q; }

        // Methods
        void energy_shift(float energy_loss);
    };

    void particle::set_momentum(float px, float py, float pz)
    {
        momentum[0] = px;
        momentum[1] = py;
        momentum[2] = pz;
        momentum[3] = std::sqrt(px * px + py * py + pz * pz + mass * mass);
    }
    void particle::set_momentum(float px, float py, float pz, float E)
    {
        momentum = {px, py, pz, E};
        mass = std::sqrt(E * E - (px * px + py * py + pz * pz));
    }
    void particle::set_momentum(const std::array<float, 3> &direction, float magnitude)
    {
        float norm = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
        if (norm > 0)
        {
            momentum[0] = magnitude * direction[0] / norm;
            momentum[1] = magnitude * direction[1] / norm;
            momentum[2] = magnitude * direction[2] / norm;
            momentum[3] = std::sqrt(momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2] + mass * mass);
        }
    }
    void particle::energy_shift(float energy_shift)
    {
        std::array<float, 3> direction = {get_px(), get_py(), get_pz()};
        auto momentum = get_p();
        if ((get_E() + energy_shift) <= mass)
        {
            set_momentum(0., 0., 0.);
            return;
        }
        auto new_momentum = std::sqrt((get_E() + energy_shift) * (get_E() + energy_shift) - mass * mass);
        set_momentum(direction[0] * (new_momentum) / (momentum), direction[1] * (new_momentum) / (momentum), direction[2] * (new_momentum) / (momentum));
    }

    class ch_radiator
    {
    private:
        physics_vector position;
        TGraphErrors refractiveIndexGraph;
        TGraphErrors absorptionLengthGraph;
        TGraphErrors scatteringLengthGraph;
        double refractiveIndex;
        double absorptionLength;
        double scatteringLength;
        double thickness;
        double length;
        double width;

    public:
        //  TODO:
        // Define better thickness, width, length

        // Constructor
        ch_radiator(physics_vector p = physics_vector(0, 0, 0, 0), double t = 1.0, double l = 1.0, double w = 1.0, double ri = 1.0, double a = 1.0, double s = 1.0)
            : position(p), thickness(t), length(l), width(w), refractiveIndex(ri), absorptionLength(a), scatteringLength(s) {}

        // Getters
        const TGraphErrors &get_refractive_index_graph() const { return refractiveIndexGraph; }
        const TGraphErrors &get_absorption_length_graph() const { return absorptionLengthGraph; }
        const TGraphErrors &get_scattering_length_graph() const { return scatteringLengthGraph; }
        double get_refractive_index() const { return refractiveIndex; }
        double get_absorption_length() const { return absorptionLength; }
        double get_scattering_length() const { return scatteringLength; }
        double get_thickness() const { return thickness; }
        double get_length() const { return length; }
        double get_width() const { return width; }
        physics_vector get_position() const { return position; }
        float get_x() const { return position.get_x(); }
        float get_y() const { return position.get_y(); }
        float get_z() const { return position.get_z(); }
        float get_t() const { return position.get_t(); }

        // Setters
        void set_refractive_index_graph(const TGraphErrors &graph) { refractiveIndexGraph = graph; }
        void set_absorption_length_graph(const TGraphErrors &graph) { absorptionLengthGraph = graph; }
        void set_scattering_length_graph(const TGraphErrors &graph) { scatteringLengthGraph = graph; }
        void set_refractive_index(double ri) { refractiveIndex = ri; }
        void set_absorption_length(double a) { absorptionLength = a; }
        void set_scattering_length(double s) { scatteringLength = s; }
        void set_thickness(double t) { thickness = t; }
        void set_length(double l) { length = l; }
        void set_width(double w) { width = w; }
        void set_position(const physics_vector &pos) { position = pos; }

        //  Check the particle is withint the radiator
        bool is_inside_the_radiator(particle target_particle);
    };

    bool ch_radiator::is_inside_the_radiator(particle target_particle)
    {
        if (fabs(target_particle.get_x() - get_x()) > get_length() * 0.5)
            return false;
        if (fabs(target_particle.get_y() - get_y()) > get_width() * 0.5)
            return false;
        if (fabs(target_particle.get_z() - get_z()) > get_thickness() * 0.5)
            return false;
        return true;
    }

    class fast_mc_geometry
    {
    private:
        // List of radiators
        int radiators_counter;
        std::map<int, ch_radiator> radiators; // List of radiators

    public:
        // Constructor
        fast_mc_geometry()
            : radiators_counter(0) {}

        // Getter for all radiators
        int get_n_radiators() { return radiators_counter; }
        const std::map<int, ch_radiator> &get_radiators() const { return radiators; }

        // Setter for all radiators (sets the whole map)
        void set_radiators(const std::map<int, ch_radiator> &new_radiators) { radiators = new_radiators; }

        // Methods
        void add_radiator(ch_radiator target_radiator);
    };

    void fast_mc_geometry::add_radiator(ch_radiator target_radiator)
    {
        radiators[radiators_counter] = target_radiator;
        radiators_counter++;
    }

    class fast_mc_event
    {
    private:
        int particles_counter;
        std::map<int, int> relationships;  // Maps mother particle to daughter particle
        std::map<int, particle> particles; // Maps particle ID to particle object

        // Geometry
        fast_mc_geometry event_geometry;

    public:
        // Constructor
        fast_mc_event() : particles_counter(0) {}

        // Getter for particles_counter
        int get_n_particles() const { return particles_counter; }
        std::map<int, int> &get_relationships() { return relationships; }
        std::map<int, particle> &get_particles() { return particles; }
        void set_particle_by_id(int id, const particle &p) { particles[id] = p; }
        const std::map<int, ch_radiator> &get_radiators() { return event_geometry.get_radiators(); }

        // Setter for particles_counter
        void set_geometry(fast_mc_geometry current_geometry) { event_geometry = current_geometry; }
        void set_relationships(const std::map<int, int> &new_relationships) { relationships = new_relationships; }
        void set_particles(const std::map<int, particle> &new_particles) { particles = new_particles; }
        particle get_id_particle(int id) const;
        void set_relationship(int mother_id, int daughter_id) { relationships[mother_id] = daughter_id; };

        //  Methods
        int add_particle(int mother_id, particle added_particle);
    };

    particle fast_mc_event::get_id_particle(int id) const
    {
        auto it = particles.find(id);
        if (it != particles.end())
            return it->second;
        throw std::out_of_range("Particle ID not found");
    }
    int fast_mc_event::add_particle(int mother_id, particle added_particle)
    {
        particles_counter++;
        particles[particles_counter] = added_particle;
        set_relationship(mother_id, particles_counter);
        return particles_counter;
    }

    class propagator
    {
    private:
        //  Physics limits
        std::array<float, 2> cherenkov_emission_wavelength = {280, 850};

        fast_mc_event &target_event;
        std::map<int, bool> list_of_end_propagation_particles;

    public:
        // Constructors
        propagator(fast_mc_event &event) : target_event(event) {}

        // General methods
        void propagate(std::vector<int> &target_particles, double dt = 1.e-4);
        void run_event();

        // Cherenkov
        std::vector<std::pair<int, particle>> get_photons_prediction(particle target_particle, ch_radiator target_radiator);
    };

    void propagator::propagate(std::vector<int> &target_particles_id, double dt = 1.e-4)
    {
        //  Loop through the target particles
        std::vector<particle *> target_particles;
        for (auto it = target_particles_id.begin(); it != target_particles_id.end();)
        {
            int id = *it;
            auto &particles = target_event.get_particles();

            if (particles.find(id) != particles.end())      // Check if ID exists in particles
                target_particles.push_back(&particles[id]); // Store a pointer to the particle

            it = target_particles_id.erase(it);
        }
        //  Empty particle id vector
        target_particles_id.clear();

        //  Propagate the target particles
        for (auto current_particle : target_particles)
        {
            while (true)
            {
                auto end_of_propagation = true;
                for (auto [id, current_radiator] : target_event.get_radiators())
                {
                    if (!current_radiator.is_inside_the_radiator(*current_particle))
                        continue;
                    end_of_propagation = false;
                    current_particle->set_position(current_particle->get_x() + current_particle->get_betax() * dt, current_particle->get_y() + current_particle->get_betay() * dt, current_particle->get_z() + current_particle->get_betaz() * dt, current_particle->get_t() + dt);
                    if (current_particle->get_mass() == 0)
                        continue;
                    auto ph_prediction = get_photons_prediction(*current_particle, current_radiator);
                    // target_particles_id.push_back();
                }
                if (end_of_propagation)
                    break;
            }
        }
    }
    void propagator::run_event()
    {
        std::vector<int> todo_particles;
        for (auto [id, current_particle] : target_event.get_particles())
            if (id != 0)
                todo_particles.push_back(id);

        while (true)
        {
            propagate(todo_particles);
            if (!todo_particles.size())
                break;
        }
    }
    std::vector<std::pair<int, particle>> propagator::get_photons_prediction(particle target_particle, ch_radiator target_radiator)
    {
        //  Result
        std::vector<std::pair<int, particle>> result;

        //  Average number of photons produced per cm in specified interval of wavelength
        auto avg_ph_emitted = cherenkov::simp_frank_tamm_int(cherenkov_emission_wavelength[0], cherenkov_emission_wavelength[1], target_particle.get_beta(), target_radiator.get_refractive_index());

        //  Decides whether and how many photons to produce
        auto current_photon_probability = gRandom->Uniform(0., 1.);
        auto current_integral_probability = 0.;
        auto n_photons = -1;
        while (true)
        {
            n_photons++;
            current_integral_probability += TMath::PoissonI(n_photons, avg_ph_emitted);
            if (current_photon_probability <= current_integral_probability)
                break;
        }

        //  If no photons are produced return
        if (!n_photons)
            return result;

        //  Generate n photons
        //  Find a base for the momentum vector
        auto direction = target_particle.get_momentum().get_direction();

        //
        return result;
    }
};
