// SPDX-License-Identifier: MIT
//
// example/run.cxx — end-to-end run of the cherenkov_mc fast Monte Carlo.
//
// Shoots 11.5 GeV electrons through two aerogel radiators, emits and tracks
// Cherenkov photons, applies a SiPM photon-detection efficiency (PDE), and
// writes the resulting spectra to an output ROOT file.
//
// Build:   cmake -B build -DCHERENKOV_MC_BUILD_EXAMPLES=ON && cmake --build build
// Run:     ./build/bin/run [n_events] [data/PDE.root] [out.root]
//
#include <memory>
#include <string>

#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>

#include <mist/hep/globals.h>
#include <mist/logger/logger.h>
#include <mist/logger/progress_bar.h>

#include <cherenkov_mc/cherenkov_mc.h>

int main(int argc, char **argv)
{
    using namespace cherenkov_mc;

    const int n_events = (argc > 1) ? std::stoi(argv[1]) : 100;
    const std::string pde_path = (argc > 2) ? argv[2] : "data/PDE.root";
    const std::string out_path = (argc > 3) ? argv[3] : "run.root";

    //  Load sensor PDE curve
    std::unique_ptr<TFile> pde_file(TFile::Open(pde_path.c_str()));
    if (!pde_file || pde_file->IsZombie())
    {
        mist::logger::error("Cannot open PDE file: " + pde_path);
        return 1;
    }
    auto *pde_13_50 = dynamic_cast<TGraph *>(pde_file->Get("13_50_PDE"));
    if (!pde_13_50)
    {
        mist::logger::error("PDE graph '13_50_PDE' not found in " + pde_path);
        return 1;
    }

    //  Build the geometry: a world box enclosing two aerogel radiators
    fast_mc_geometry full_geometry;

    auto *world = new box(physics_vector(0, 0, 0, 2), std::nullopt,
                          {{"side_x", 10.}, {"side_y", 10.}, {"side_z", 10.}});
    full_geometry.add_world_container(world);

    ch_radiator radiator_AG22_J001(physics_vector(0, 0, 0, 1), 2, 5, 5, 1.0210, 1. / 6., 1. / 100.);
    radiator_AG22_J001.set_absorption_length(0.00753252, 10, 2e6, 1.49469e-05, 8.2902);
    radiator_AG22_J001.set_scattering_length(7.26886e+09, 0);
    full_geometry.add_radiator(radiator_AG22_J001);

    ch_radiator radiator_AG22_J003(physics_vector(0, 0, 0, 3), 2, 5, 5, 1.0207, 1. / 6., 1. / 100.);
    radiator_AG22_J003.set_absorption_length(0.0136385, 10, 2e6, 1.29517e-05, 7.93124);
    radiator_AG22_J003.set_scattering_length(7.62127e+09, 0);
    full_geometry.add_radiator(radiator_AG22_J003);

    //  Output spectra
    TH2F h_hits("h_hits", ";x (cm);y (cm)", 200, -1, 1, 200, -1, 1);
    TH1F h_wvl_emit("h_wvl_emit", ";Photon emitted wavelength (nm);", 700, 200, 900);
    TH1F h_wvl_detd("h_wvl_detd", ";Photon detected wavelength (nm);", 700, 200, 900);
    TH1F h_nph_emit("h_nph_emit", ";# of photons emitted", 300, 0, 300);
    TH1F h_nph_detd("h_nph_detd", ";# of photons detected", 300, 0, 300);

    mist::logger::info("Generating " + std::to_string(n_events) + " events");
    mist::logger::ProgressBar bar("events");

    for (int i_ev = 0; i_ev < n_events; ++i_ev)
    {
        fast_mc_event current_event;
        current_event.set_geometry(full_geometry);
        propagator current_propagator(current_event);

        particle start_electron;
        start_electron.set_mass(kElectronMass);
        start_electron.set_momentum(0, 0, std::sqrt(11.5 * 11.5 - kElectronMass * kElectronMass));
        current_event.add_particle(0, start_electron);

        current_propagator.run_event();

        int n_emit = 0, n_detd = 0;
        for (auto &[id, current_particle] : current_event.get_particles())
        {
            if (current_particle.get_mass() != 0)
                continue;

            const double lambda_nm = 1.e-9 * 1240. / current_particle.get_E();
            ++n_emit;
            h_wvl_emit.Fill(lambda_nm);

            if (current_particle.get_absorbed())
                continue;

            const double detection_prob = 0.01 * pde_13_50->Eval(lambda_nm);
            if (mist::hep::rng().Uniform(0., 1.) > detection_prob)
                continue;

            ++n_detd;
            h_wvl_detd.Fill(lambda_nm);
            h_hits.Fill(current_particle.get_x(), current_particle.get_y());
        }
        h_nph_emit.Fill(n_emit);
        h_nph_detd.Fill(n_detd);
        bar.update(i_ev + 1, n_events);
    }
    bar.finish();

    mist::logger::info("Mean photons emitted/event:  " + std::to_string(h_nph_emit.GetMean()));
    mist::logger::info("Mean photons detected/event: " + std::to_string(h_nph_detd.GetMean()));

    std::unique_ptr<TFile> fout(TFile::Open(out_path.c_str(), "RECREATE"));
    h_hits.Write();
    h_wvl_emit.Write();
    h_wvl_detd.Write();
    h_nph_emit.Write();
    h_nph_detd.Write();
    fout->Close();

    mist::logger::info("Wrote spectra to " + out_path);
    return 0;
}
