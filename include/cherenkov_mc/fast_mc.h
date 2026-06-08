// SPDX-License-Identifier: MIT
//
// cherenkov_mc/fast_mc.h — the fast Monte Carlo itself: radiators, the event
// container (particles + geometry), and the propagator that emits and tracks
// Cherenkov photons.
//
// Declarations live here; definitions in src/cherenkov_mc/fast_mc.cxx.
//
// TODO: derive particle mass / charge from a supplied PDG id; introduce a
//       proper per-event / per-particle time model.
//
#pragma once

#include <array>
#include <map>
#include <vector>

#include <TF1.h>
#include <TGraphErrors.h>

#include "cherenkov_mc/geometry.h"
#include "cherenkov_mc/particle.h"
#include "cherenkov_mc/photon_interactions.h"
#include "cherenkov_mc/physics_vector.h"

namespace cherenkov_mc
{
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
        TF1 *fAbsorptionLength;
        TF1 *fScatteringLength;
        double thickness;
        double length;
        double width;

    public:
        // Constructor
        ch_radiator(physics_vector p = physics_vector(0, 0, 0, 0), double t = 1.0, double l = 1.0, double w = 1.0, double ri = 1.0, double a = 1.0, double s = 1.0)
            : position(p), refractiveIndex(ri), absorptionLength(a), scatteringLength(s),
              fAbsorptionLength((TF1 *)photon_interactions::f_absorption().Clone()),
              fScatteringLength((TF1 *)photon_interactions::f_scattering().Clone()),
              thickness(t), length(l), width(w) {}

        // Getters
        const TGraphErrors &get_refractive_index_graph() const { return refractiveIndexGraph; }
        const TGraphErrors &get_absorption_length_graph() const { return absorptionLengthGraph; }
        const TGraphErrors &get_scattering_length_graph() const { return scatteringLengthGraph; }
        double get_refractive_index() const { return refractiveIndex; }
        double get_absorption_length() const { return absorptionLength; }
        double get_scattering_length() const { return scatteringLength; }
        double get_absorption_length(double lambda) { return fAbsorptionLength->Eval(lambda); }
        double get_scattering_length(double lambda) { return fScatteringLength->Eval(lambda); }
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
        void set_absorption_length(double p0, double p1, double p2, double p3, double p4) { fAbsorptionLength->SetParameters(p0, p1, p2, p3, p4); }
        void set_scattering_length(double p0, double p1) { fScatteringLength->SetParameters(p0, p1); }
        void set_thickness(double t) { thickness = t; }
        void set_length(double l) { length = l; }
        void set_width(double w) { width = w; }
        void set_position(const physics_vector &pos) { position = pos; }

        std::array<double, 2> get_max_step_size(particle target_particle);

        //  Check whether the particle is within the radiator
        bool is_inside_the_radiator(particle target_particle);
    };

    class fast_mc_geometry
    {
    private:
        // List of world containers
        int world_counter;
        std::map<int, geometry *> world_containers; // World box

        // List of radiators
        int radiators_counter;
        std::map<int, ch_radiator> radiators; // List of radiators

    public:
        // Constructor
        fast_mc_geometry()
            : world_counter(0), radiators_counter(0) {}

        // Getter for all radiators
        int get_n_radiators() { return radiators_counter; }
        const std::map<int, ch_radiator> &get_radiators() const { return radiators; }
        const std::map<int, geometry *> &get_world_containers() const { return world_containers; }

        // Setter for all radiators
        void set_radiators(const std::map<int, ch_radiator> &new_radiators) { radiators = new_radiators; }
        void set_world_containers(const std::map<int, geometry *> &new_world_containers) { world_containers = new_world_containers; }

        // Methods
        void add_radiator(ch_radiator target_radiator);
        void add_world_container(geometry *target_geometry);
    };

    class fast_mc_event
    {
    private:
        int particles_counter;
        std::map<int, std::vector<int>> relationships; // Maps a mother particle id to its daughter ids
        std::map<int, particle> particles;             // Maps particle ID to particle object

        // Geometry
        fast_mc_geometry event_geometry;

    public:
        // Constructor
        fast_mc_event() : particles_counter(0) {}

        // Getter for particles_counter
        int get_n_particles() const { return particles_counter; }
        std::map<int, std::vector<int>> &get_relationships() { return relationships; }
        std::map<int, particle> &get_particles() { return particles; }
        void set_particle_by_id(int id, const particle &p) { particles[id] = p; }
        const std::map<int, ch_radiator> &get_radiators() { return event_geometry.get_radiators(); }
        const std::map<int, geometry *> &get_world_containers() { return event_geometry.get_world_containers(); }
        particle get_id_particle(int id) const;

        // Setter for particles_counter
        void set_geometry(fast_mc_geometry current_geometry) { event_geometry = current_geometry; }
        void set_relationships(const std::map<int, std::vector<int>> &new_relationships) { relationships = new_relationships; }
        void set_particles(const std::map<int, particle> &new_particles) { particles = new_particles; }
        void set_relationship(int mother_id, int daughter_id) { relationships[mother_id].push_back(daughter_id); }

        //  Methods
        int add_particle(int mother_id, particle added_particle);
    };

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
        propagator(fast_mc_event &event) : dt(1.e-13), target_event(event) {}

        // General methods
        void propagate(std::vector<int> &target_particles);
        void run_event();

        //  Photons
        std::vector<int> propagator_photons(particle &target_particle, ch_radiator target_radiator, bool &end_of_propagation);

        // Cherenkov
        std::vector<int> get_photons_prediction(particle &target_particle, ch_radiator target_radiator);
    };

} // namespace cherenkov_mc
