//  last update: 02/04/2025
//  author: Nicola Rubini - nicola.rubini@bo.infn.it
//  Scope: Cherenkov fast mc
//  TODO:   1. In particle you should be able to give the PDG id and set automatically mass, charge, etc.
//          2. Add event time, define a particle internal time (?) To be thought about

#pragma once

#include "probability.h"
#include "particle.h"
#include "physics_vector.h"
#include "cherenkov.h"

namespace ch_fast_mc
{
    class fast_mc_event;
    class fast_mc_geometry;

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
        double get_x() const { return position.get_x(); }
        double get_y() const { return position.get_y(); }
        double get_z() const { return position.get_z(); }
        double get_t() const { return position.get_t(); }

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
        particle get_id_particle(int id) const;
        int get_particle_id(particle target) const;

        // Setter for particles_counter
        void set_geometry(fast_mc_geometry current_geometry) { event_geometry = current_geometry; }
        void set_relationships(const std::map<int, int> &new_relationships) { relationships = new_relationships; }
        void set_particles(const std::map<int, particle> &new_particles) { particles = new_particles; }
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
        double dt;
        std::array<double, 2> cherenkov_emission_wavelength = {280, 850};

        fast_mc_event &target_event;
        std::map<int, bool> list_of_end_propagation_particles;

    public:
        // Constructors
        propagator(fast_mc_event &event) : target_event(event), dt(1.e-12) {}

        // General methods
        void propagate(std::vector<int> &target_particles);
        void run_event();

        //  Photons
        std::vector<std::pair<int, particle>> propagator_photons(particle &target_particle, ch_radiator target_radiator);

        // Cherenkov
        std::vector<int> get_photons_prediction(particle &target_particle, ch_radiator target_radiator);
    };

    void propagator::propagate(std::vector<int> &target_particles_id)
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
                    if (current_particle->get_mass() == 0)
                    {
                        produced_particles = propagator_photons(*current_particle, current_radiator);
                        for (auto id : produced_particles)
                            target_particles_id.push_back(id)
                    }
                    else
                    {
                        current_particle->set_position(current_particle->get_t() + dt, current_particle->get_x() + current_particle->get_betax() * _c0 * 100 * dt, current_particle->get_y() + current_particle->get_betay() * _c0 * 100 * dt, current_particle->get_z() + current_particle->get_betaz() * _c0 * 100 * dt);
                        auto produced_particles = get_photons_prediction(*current_particle, current_radiator);
                        for (auto id : produced_particles)
                            target_particles_id.push_back(id);
                    }
                }
                // break;
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
    std::vector<int> propagator::get_photons_prediction(particle &target_particle, ch_radiator target_radiator)
    {
        //  Result
        std::vector<int> result;

        //  Average number of photons produced per cm in specified interval of wavelength
        auto avg_ph_emitted = dt * target_particle.get_beta() * _c0 * 100 * cherenkov::simp_frank_tamm_int(cherenkov_emission_wavelength[0], cherenkov_emission_wavelength[1], target_particle.get_beta(), target_radiator.get_refractive_index());

        //  Decides whether and how many photons to produce
        auto n_photons = probability::poisson_sampling_Atkinson(avg_ph_emitted);

        //  If no photons are produced return
        if (!n_photons)
            return result;

        //  Generate n photons
        //  Find a base for the momentum vector and the target theta
        auto orth_base = target_particle.get_momentum().get_spatial_orth_base();
        double emission_theta = cherenkov::get_ch_theta(target_particle.get_beta(), target_radiator.get_refractive_index());
        for (auto i_ph = 0; i_ph < n_photons; i_ph++)
        {
            //  Define the new photon momentum
            physics_vector current_photon_momentum = {0., 0., 0., 1.}; // = target_particle.get_momentum().get_direction();
            auto emission_phi = 2 * TMath::Pi() * gRandom->Uniform(0, 1);
            current_photon_momentum.rotation_on_an_axis(emission_theta, orth_base[1]);
            current_photon_momentum.rotation_on_an_axis(emission_phi, orth_base[0]);
            current_photon_momentum.set_t(1.);

            //  Calculate the wavelength of emission
            double emission_wavelength = (cherenkov_emission_wavelength[0]) / (1 - (gRandom->Uniform(0, 1)) * (cherenkov_emission_wavelength[1] - cherenkov_emission_wavelength[0]) / (cherenkov_emission_wavelength[1]));
            double emission_energy = 1.e-9 * 1240. / emission_wavelength;
            current_photon_momentum *= emission_energy;
            particle current_photon;
            current_photon.set_momentum(current_photon_momentum);
            current_photon.set_mass(kPhotonMass);

            result.push_back(target_event.add_particle(target_particle.get_event_id(), current_photon));

            //  Reduce emitting particle energy
            target_particle.energy_shift(-emission_energy);
        }
        //
        return result;
    }
    std::vector<std::pair<int, particle>> propagator::propagator_photons(particle &target_particle, ch_radiator target_radiator)
    {
        auto step_size = photon_interaction::get_step_size();
        
    }
};