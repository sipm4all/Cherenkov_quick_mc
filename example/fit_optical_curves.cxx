// SPDX-License-Identifier: MIT
//
// example/fit_optical_curves.cxx — calibrate the optical models.
//
// Fits the parametrised transmission / scattering-length / absorption-length
// TF1s (cherenkov_mc::photon_interactions) to the measured aerogel curves in a
// "clarity" ROOT file, then reports the fitted parameters — exactly the
// numbers that feed ch_radiator::set_absorption_length / set_scattering_length.
//
// Run:  ./build/bin/fit_optical_curves data/clarity_jap2022_AG22J003.root
//
#include <memory>
#include <string>

#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>

#include <mist/logger/logger.h>

#include <cherenkov_mc/photon_interactions.h>

namespace pi = cherenkov_mc::photon_interactions;

namespace {

void report(const char *what, const TF1 &f)
{
    std::string line = std::string(what) + ":";
    for (int i = 0; i < f.GetNpar(); ++i)
        line += " " + std::string(f.GetParName(i)) + "=" + std::to_string(f.GetParameter(i));
    mist::logger::info(line);
}

} // namespace

int main(int argc, char **argv)
{
    const std::string file = (argc > 1) ? argv[1] : "data/clarity_jap2022_AG22J003.root";

    std::unique_ptr<TFile> fin(TFile::Open(file.c_str()));
    if (!fin || fin->IsZombie())
    {
        mist::logger::error("Cannot open clarity file: " + file);
        return 1;
    }

    auto *gTrans = dynamic_cast<TGraph *>(fin->Get("gTrans"));
    auto *gScatt = dynamic_cast<TGraph *>(fin->Get("gLambdaS"));
    auto *gAbsor = dynamic_cast<TGraph *>(fin->Get("gLambdaA"));
    if (!gTrans || !gScatt || !gAbsor)
    {
        mist::logger::error("Expected graphs gTrans / gLambdaS / gLambdaA not all found in " + file);
        return 1;
    }

    //  Transmission
    pi::f_transmission().SetParameters(1, 1.e10, 1.e10);
    for (int i = 0; i < 4; ++i)
        gTrans->Fit(&pi::f_transmission(), "QN");
    gTrans->Fit(&pi::f_transmission(), "QN");
    report("transmission", pi::f_transmission());

    //  Scattering length
    pi::f_scattering().SetParameters(1, 1, 1);
    for (int i = 0; i < 4; ++i)
        gScatt->Fit(&pi::f_scattering(), "QN");
    gScatt->Fit(&pi::f_scattering(), "QN");
    report("scattering_length", pi::f_scattering());

    //  Absorption length
    pi::f_absorption().SetParameters(0.02, 10, 2e6, 1.e-5, 8);
    for (int i = 0; i < 4; ++i)
        gAbsor->Fit(&pi::f_absorption(), "QN", "", 230, 900);
    gAbsor->Fit(&pi::f_absorption(), "QN", "", 230, 900);
    report("absorption_length", pi::f_absorption());

    return 0;
}
