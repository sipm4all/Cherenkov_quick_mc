// SPDX-License-Identifier: MIT
#include "cherenkov_mc/fast_mc.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <TMath.h>

#include <mist/hep/globals.h>

#include "cherenkov_mc/cherenkov.h"
#include "cherenkov_mc/constants.h"
#include "cherenkov_mc/probability.h"

namespace cherenkov_mc
{
    // === ch_radiator ===

    bool ch_radiator::is_inside_the_radiator(particle target_particle)
    {
        if (std::fabs(target_particle.get_x() - get_x()) > get_length() * 0.5)
            return false;
        if (std::fabs(target_particle.get_y() - get_y()) > get_width() * 0.5)
            return false;
        if (std::fabs(target_particle.get_z() - get_z()) > get_thickness() * 0.5)
            return false;
        return true;
    }

    std::array<double, 2> ch_radiator::get_max_step_size(particle target_particle)
    {
        //  Get the maximum step size for the particle in the radiator
        auto max_t_x = (get_x() + get_length() * 0.5 - target_particle.get_x()) / (target_particle.get_betax() * _c0 * 100);
        auto max_t_y = (get_y() + get_width() * 0.5 - target_particle.get_y()) / (target_particle.get_betay() * _c0 * 100);
        auto max_t_z = (get_z() + get_thickness() * 0.5 - target_particle.get_z()) / (target_particle.get_betaz() * _c0 * 100);
        auto min_t_x = (-get_x() - get_length() * 0.5 + target_particle.get_x()) / (target_particle.get_betax() * _c0 * 100);
        auto min_t_y = (-get_y() - get_width() * 0.5 + target_particle.get_y()) / (target_particle.get_betay() * _c0 * 100);
        auto min_t_z = (-get_z() - get_thickness() * 0.5 + target_particle.get_z()) / (target_particle.get_betaz() * _c0 * 100);

        // Collect all times in a vector
        std::vector<double> times = {max_t_x, max_t_y, max_t_z, min_t_x, min_t_y, min_t_z};

        // Find the smallest positive time
        float min_positive_time = std::numeric_limits<float>::max();
        for (const auto &time : times)
            if (time > 0 && time < min_positive_time)
                min_positive_time = time;

        // Return the smallest positive time
        return {min_positive_time * target_particle.get_beta() * _c0 * 100, min_positive_time};
    }

    // === fast_mc_geometry ===

    void fast_mc_geometry::add_radiator(ch_radiator target_radiator)
    {
        radiators[radiators_counter] = target_radiator;
        radiators_counter++;
    }

    void fast_mc_geometry::add_world_container(geometry *target_geometry)
    {
        world_containers[world_counter] = target_geometry;
        world_counter++;
    }

    // === fast_mc_event ===

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
        //  Stamp the particle with its own id so later parentage look-ups
        //  (add_particle(mother.get_event_id(), ...)) resolve to a real mother.
        added_particle.set_event_id(particles_counter);
        particles[particles_counter] = added_particle;
        set_relationship(mother_id, particles_counter);
        return particles_counter;
    }

    // === propagator ===

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
                //  Stop once the particle has left every world container
                auto inside_world = false;
                for (auto [id, current_world_container] : target_event.get_world_containers())
                {
                    auto current_position = current_particle->get_position();
                    if (current_world_container->is_inside(current_position))
                    {
                        inside_world = true;
                        break;
                    }
                }
                if (!inside_world)
                    break;

                //  Interact with any radiator the particle currently sits in
                auto end_of_propagation = false;
                auto not_within_radiators = true;
                for (auto [id, current_radiator] : target_event.get_radiators())
                {
                    if (!current_radiator.is_inside_the_radiator(*current_particle))
                        continue;
                    not_within_radiators = false;
                    if (current_particle->get_mass() == 0)
                    {
                        auto produced_particles = propagator_photons(*current_particle, current_radiator, end_of_propagation);
                        for (auto p_id : produced_particles)
                            target_particles_id.push_back(p_id);
                    }
                    else
                    {
                        current_particle->step(dt);
                        auto produced_particles = get_photons_prediction(*current_particle, current_radiator);
                        for (auto p_id : produced_particles)
                            target_particles_id.push_back(p_id);
                    }
                }

                //  A photon that was absorbed (or otherwise terminated) ends here
                if (end_of_propagation)
                    break;

                //  Drift through field-free / vacuum regions of the world
                if (not_within_radiators)
                    current_particle->step(dt);
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

        if (avg_ph_emitted < 0)
            return result;

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
            physics_vector current_photon_momentum = {0., 0., 0., 1.};
            auto emission_phi = 2 * TMath::Pi() * mist::hep::rng().Uniform(0, 1);
            current_photon_momentum.rotation_on_an_axis(emission_theta, orth_base[1]);
            current_photon_momentum.rotation_on_an_axis(emission_phi, orth_base[0]);
            current_photon_momentum.set_t(1.);

            //  Calculate the wavelength of emission
            double emission_wavelength = (cherenkov_emission_wavelength[0]) / (1 - (mist::hep::rng().Uniform(0, 1)) * (cherenkov_emission_wavelength[1] - cherenkov_emission_wavelength[0]) / (cherenkov_emission_wavelength[1]));
            double emission_energy = 1.e-9 * 1240. / emission_wavelength;
            current_photon_momentum *= emission_energy;
            particle current_photon;
            current_photon.set_momentum(current_photon_momentum);
            current_photon.set_position(target_particle.get_position());
            current_photon.set_t(0.);
            current_photon.set_mass(kPhotonMass);

            result.push_back(target_event.add_particle(target_particle.get_event_id(), current_photon));

            //  Reduce emitting particle energy
            target_particle.energy_shift(-emission_energy);
        }
        //
        return result;
    }

    std::vector<int> propagator::propagator_photons(particle &target_particle, ch_radiator target_radiator, bool &end_of_propagation)
    {
        //  Result
        std::vector<int> result;

        //  Check if the photon is absorbed
        if (target_particle.get_absorbed())
        {
            end_of_propagation = true;
            return result;
        }

        //  Get the step size & interaction type
        auto max_step_size = target_radiator.get_max_step_size(target_particle);
        auto step_size = photon_interactions::get_step_size(1. / target_radiator.get_absorption_length(1.e-9 * 1240 / target_particle.get_E()), 1. / target_radiator.get_scattering_length(1.e-9 * 1240 / target_particle.get_E()));

        //  Define if it interacts or not
        if (max_step_size[0] < step_size[0])
        {
            step_size[1] = max_step_size[1];
            target_particle.step(step_size[1] + dt);
            return result;
        }

        //  Define interaction type (sampled on the interaction coefficients, i.e. 1 / length)
        auto interaction_type = photon_interactions::get_interaction(1. / target_radiator.get_absorption_length(1.e-9 * 1240 / target_particle.get_E()), 1. / target_radiator.get_scattering_length(1.e-9 * 1240 / target_particle.get_E()));

        //  Propagate the photon
        target_particle.step(step_size[1]);

        //  Absorption case
        if (interaction_type == photon_interactions::ABSORPTION)
        {
            target_particle.set_absorbed(true);
            end_of_propagation = true;
            return result;
        }

        //  Scattering case: redirect the photon and keep propagating it
        if (interaction_type == photon_interactions::SCATTERING)
        {
            //  Assign through the momentum reference so the (massless) photon
            //  energy/mass bookkeeping is preserved — only the direction changes
            target_particle.get_momentum() = photon_interactions::scattering_interaction(target_particle.get_momentum());
            return result;
        }

        return result;
    }

} // namespace cherenkov_mc
