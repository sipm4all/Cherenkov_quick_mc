struct TreeData
{
    Long64_t time;
    Float_t tset, tint, tlm73, tsipma, tsipmb, tsipmc;
    Float_t ia, ib, ic, va, vb, vc;
    Int_t stable;
    Float_t pa, pb, pc;
};

void quick_look(const char *filename)
{
    // Open the ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie())
    {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get the tree
    TTree *tree = (TTree *)file->Get("tree"); // Change "tree" to the actual name of your TTree
    if (!tree)
    {
        std::cerr << "Error: TTree not found in file!" << std::endl;
        file->Close();
        return;
    }

    // Define structure instance
    TreeData data;

    // Set branch addresses
    tree->SetBranchAddress("time", &data.time);
    tree->SetBranchAddress("tset", &data.tset);
    tree->SetBranchAddress("tint", &data.tint);
    tree->SetBranchAddress("tlm73", &data.tlm73);
    tree->SetBranchAddress("tsipma", &data.tsipma);
    tree->SetBranchAddress("tsipmb", &data.tsipmb);
    tree->SetBranchAddress("tsipmc", &data.tsipmc);
    tree->SetBranchAddress("ia", &data.ia);
    tree->SetBranchAddress("ib", &data.ib);
    tree->SetBranchAddress("ic", &data.ic);
    tree->SetBranchAddress("va", &data.va);
    tree->SetBranchAddress("vb", &data.vb);
    tree->SetBranchAddress("vc", &data.vc);
    tree->SetBranchAddress("stable", &data.stable);
    tree->SetBranchAddress("pa", &data.pa);
    tree->SetBranchAddress("pb", &data.pb);
    tree->SetBranchAddress("pc", &data.pc);

    //  Utility
    std::vector<float> t_sipm_set = {0., 10., 20., 30., 40., 50.};

    //  Output
    std::map<std::string, std::vector<std::array<float, 2>>> data_points;

    // Loop over entries
    Long64_t nEntries = tree->GetEntries();
    auto new_stable = false;
    for (Long64_t i = 0; i < nEntries; i++)
    {
        //  Load entry
        tree->GetEntry(i);

        //  Reject unstable sections
        if (!data.stable)
        {
            new_stable = true;
            continue;
        }

        // Reach a stability plateau
        if (new_stable)
        {
            data_points["tlm73"].push_back({0., 0.});
            data_points["sipmA"].push_back({0., 0.});
            data_points["sipmB"].push_back({0., 0.});
            data_points["sipmC"].push_back({0., 0.});
            data_points["voltA"].push_back({0., 0.});
            data_points["voltB"].push_back({0., 0.});
            data_points["voltC"].push_back({0., 0.});
            data_points["powrA"].push_back({0., 0.});
            data_points["powrB"].push_back({0., 0.});
            data_points["powrC"].push_back({0., 0.});
            new_stable = false;
        }

        // Average on data points
        data_points["tlm73"].back()[0] += data.tlm73;
        data_points["tlm73"].back()[1] += 1;
        data_points["sipmA"].back()[0] += data.tsipma;
        data_points["sipmA"].back()[1] += 1;
        data_points["sipmB"].back()[0] += data.tsipmb;
        data_points["sipmB"].back()[1] += 1;
        data_points["sipmC"].back()[0] += data.tsipmc;
        data_points["sipmC"].back()[1] += 1;
        data_points["voltA"].back()[0] += data.va;
        data_points["voltA"].back()[1] += 1;
        data_points["voltB"].back()[0] += data.vb;
        data_points["voltB"].back()[1] += 1;
        data_points["voltC"].back()[0] += data.vc;
        data_points["voltC"].back()[1] += 1;
        data_points["powrA"].back()[0] += data.pa;
        data_points["powrA"].back()[1] += 1;
        data_points["powrB"].back()[0] += data.pb;
        data_points["powrB"].back()[1] += 1;
        data_points["powrC"].back()[0] += data.pc;
        data_points["powrC"].back()[1] += 1;
    }

    for (auto &[tag, current_points] : data_points)
        for (auto &current_point : current_points)
        {
            auto tot = current_point[0];
            auto npt = current_point[1];
            current_point[0] = tot / npt;
            current_point[1] = std::sqrt(tot / npt) / npt;
        }

    // TGrahpErrors
    TGraphErrors *TsipmA = new TGraphErrors();
    TGraphErrors *TsipmB = new TGraphErrors();
    TGraphErrors *TsipmC = new TGraphErrors();
    TGraphErrors *TsipmA_all = new TGraphErrors();
    TGraphErrors *TsipmB_all = new TGraphErrors();
    TGraphErrors *TsipmC_all = new TGraphErrors();
    TGraphErrors *TsipmA_Btest = new TGraphErrors();
    TGraphErrors *TsipmA_Ctest = new TGraphErrors();
    TGraphErrors *TsipmB_Atest = new TGraphErrors();
    TGraphErrors *TsipmB_Ctest = new TGraphErrors();
    TGraphErrors *TsipmC_Atest = new TGraphErrors();
    TGraphErrors *TsipmC_Btest = new TGraphErrors();
    auto ipnt = -1;
    for (auto target_point : t_sipm_set)
    {
        ipnt++;

        // Sipm A
        auto current_data = data_points["sipmA"][ipnt];
        TsipmA->SetPoint(ipnt, target_point, current_data[0] - data_points["tlm73"][0][0]);
        TsipmA->SetPointError(ipnt, 0, std::sqrt(current_data[1] * current_data[1] + data_points["tlm73"][0][1] * data_points["tlm73"][0][1]));

        // Sipm B
        current_data = data_points["sipmB"][ipnt + t_sipm_set.size()];
        TsipmB->SetPoint(ipnt, target_point, current_data[0] - data_points["tlm73"][0 + t_sipm_set.size()][0]);
        TsipmB->SetPointError(ipnt, 0, std::sqrt(current_data[1] * current_data[1] + data_points["tlm73"][0 + t_sipm_set.size()][1] * data_points["tlm73"][0 + t_sipm_set.size()][1]));

        // Sipm C
        current_data = data_points["sipmC"][ipnt + 2 * t_sipm_set.size()];
        TsipmC->SetPoint(ipnt, target_point, current_data[0] - data_points["tlm73"][0 + 2 * t_sipm_set.size()][0]);
        TsipmC->SetPointError(ipnt, 0, std::sqrt(current_data[1] * current_data[1] + data_points["tlm73"][0 + 2 * t_sipm_set.size()][1] * data_points["tlm73"][0 + 2 * t_sipm_set.size()][1]));

        // Sipm A all on
        current_data = data_points["sipmA"][ipnt + 3 * t_sipm_set.size()];
        TsipmA_all->SetPoint(ipnt, target_point, current_data[0] - data_points["tlm73"][ipnt + 3 * t_sipm_set.size()][0]);
        TsipmA_all->SetPointError(ipnt, 0, std::sqrt(current_data[1] * current_data[1] + data_points["tlm73"][0 + 3 * t_sipm_set.size()][1] * data_points["tlm73"][0 + 3 * t_sipm_set.size()][1]));

        // Sipm B all on
        current_data = data_points["sipmB"][ipnt + 3 * t_sipm_set.size()];
        TsipmB_all->SetPoint(ipnt, target_point, current_data[0] - data_points["tlm73"][ipnt + 3 * t_sipm_set.size()][0]);
        TsipmB_all->SetPointError(ipnt, 0, std::sqrt(current_data[1] * current_data[1] + data_points["tlm73"][0 + 3 * t_sipm_set.size()][1] * data_points["tlm73"][0 + 3 * t_sipm_set.size()][1]));

        // Sipm C all on
        current_data = data_points["sipmC"][ipnt + 3 * t_sipm_set.size()];
        TsipmC_all->SetPoint(ipnt, target_point, current_data[0] - data_points["tlm73"][ipnt + 3 * t_sipm_set.size()][0]);
        TsipmC_all->SetPointError(ipnt, 0, std::sqrt(current_data[1] * current_data[1] + data_points["tlm73"][0 + 3 * t_sipm_set.size()][1] * data_points["tlm73"][0 + 3 * t_sipm_set.size()][1]));

        // Test
        TsipmA_Btest->SetPoint(ipnt, data_points["powrA"][ipnt + 0. * t_sipm_set.size()][0], data_points["sipmB"][ipnt + 0. * t_sipm_set.size()][0] - data_points["tlm73"][ipnt + 0 * t_sipm_set.size()][0]);
        TsipmA_Btest->SetPointError(ipnt, data_points["powrA"][ipnt + 0. * t_sipm_set.size()][1], std::sqrt(data_points["tlm73"][ipnt + 0. * t_sipm_set.size()][1] * data_points["tlm73"][ipnt + 0. * t_sipm_set.size()][1] + data_points["sipmB"][ipnt + 0. * t_sipm_set.size()][1] * data_points["sipmB"][ipnt + 0. * t_sipm_set.size()][1]));

        TsipmA_Ctest->SetPoint(ipnt, data_points["powrA"][ipnt + 0. * t_sipm_set.size()][0], data_points["sipmC"][ipnt + 0. * t_sipm_set.size()][0] - data_points["tlm73"][ipnt + 0 * t_sipm_set.size()][0]);
        TsipmA_Ctest->SetPointError(ipnt, data_points["powrA"][ipnt + 0. * t_sipm_set.size()][1], std::sqrt(data_points["tlm73"][ipnt + 0. * t_sipm_set.size()][1] * data_points["tlm73"][ipnt + 0. * t_sipm_set.size()][1] + data_points["sipmC"][ipnt + 0. * t_sipm_set.size()][1] * data_points["sipmC"][ipnt + 0. * t_sipm_set.size()][1]));

        TsipmB_Atest->SetPoint(ipnt, data_points["powrB"][ipnt + 1. * t_sipm_set.size()][0], data_points["sipmA"][ipnt + 1. * t_sipm_set.size()][0] - data_points["tlm73"][ipnt + 1 * t_sipm_set.size()][0]);
        TsipmB_Atest->SetPointError(ipnt, data_points["powrB"][ipnt + 1. * t_sipm_set.size()][1], std::sqrt(data_points["tlm73"][ipnt + 1. * t_sipm_set.size()][1] * data_points["tlm73"][ipnt + 1. * t_sipm_set.size()][1] + data_points["sipmA"][ipnt + 1. * t_sipm_set.size()][1] * data_points["sipmA"][ipnt + 1. * t_sipm_set.size()][1]));
        TsipmB_Ctest->SetPoint(ipnt, data_points["powrB"][ipnt + 1. * t_sipm_set.size()][0], data_points["sipmC"][ipnt + 1. * t_sipm_set.size()][0] - data_points["tlm73"][ipnt + 1 * t_sipm_set.size()][0]);
        TsipmB_Ctest->SetPointError(ipnt, data_points["powrB"][ipnt + 1. * t_sipm_set.size()][1], std::sqrt(data_points["tlm73"][ipnt + 1. * t_sipm_set.size()][1] * data_points["tlm73"][ipnt + 1. * t_sipm_set.size()][1] + data_points["sipmC"][ipnt + 1. * t_sipm_set.size()][1] * data_points["sipmC"][ipnt + 1. * t_sipm_set.size()][1]));
        TsipmC_Atest->SetPoint(ipnt, data_points["powrC"][ipnt + 2. * t_sipm_set.size()][0], data_points["sipmA"][ipnt + 2. * t_sipm_set.size()][0] - data_points["tlm73"][ipnt + 2 * t_sipm_set.size()][0]);
        TsipmC_Atest->SetPointError(ipnt, data_points["powrC"][ipnt + 2. * t_sipm_set.size()][1], std::sqrt(data_points["tlm73"][ipnt + 2. * t_sipm_set.size()][1] * data_points["tlm73"][ipnt + 2. * t_sipm_set.size()][1] + data_points["sipmA"][ipnt + 2. * t_sipm_set.size()][1] * data_points["sipmA"][ipnt + 2. * t_sipm_set.size()][1]));
        TsipmC_Btest->SetPoint(ipnt, data_points["powrC"][ipnt + 2. * t_sipm_set.size()][0], data_points["sipmB"][ipnt + 2. * t_sipm_set.size()][0] - data_points["tlm73"][ipnt + 2 * t_sipm_set.size()][0]);
        TsipmC_Btest->SetPointError(ipnt, data_points["powrC"][ipnt + 2. * t_sipm_set.size()][1], std::sqrt(data_points["tlm73"][ipnt + 2. * t_sipm_set.size()][1] * data_points["tlm73"][ipnt + 2. * t_sipm_set.size()][1] + data_points["sipmB"][ipnt + 2. * t_sipm_set.size()][1] * data_points["sipmB"][ipnt + 2. * t_sipm_set.size()][1]));
    }

    TF1 *fpol1 = new TF1("fpol1", "pol1(0)");

    TCanvas *cDupmp = new TCanvas();

    TsipmA->Fit(fpol1, "N");
    auto pol1_a = (TF1 *)fpol1->Clone("pol1_a");
    pol1_a->SetLineStyle(kDashed);
    pol1_a->SetLineColor(kRed);
    pol1_a->SetRange(-10, 100);

    TsipmB->Fit(fpol1, "N");
    auto pol1_b = (TF1 *)fpol1->Clone("pol1_b");
    pol1_b->SetLineStyle(kDashed);
    pol1_b->SetLineColor(kBlue);
    pol1_b->SetRange(-10, 100);

    TsipmC->Fit(fpol1, "N");
    auto pol1_c = (TF1 *)fpol1->Clone("pol1_c");
    pol1_c->SetLineStyle(kDashed);
    pol1_c->SetLineColor(kGreen - 1);
    pol1_c->SetRange(-10, 100);

    TCanvas *c_sipm_target_reached = new TCanvas("", "", 750, 750);
    gPad->SetGridy();
    auto hframe = gPad->DrawFrame(0., 0., 60., 60);
    hframe->SetTitle(";Target #deltaT (w.r.t. LM73) (#circC);Reached #deltaT (w.r.t. LM73) (#circC)");

    TsipmA->Draw("SAME PE");
    pol1_a->DrawCopy("SAME");
    TsipmB->Draw("SAME PE");
    pol1_b->DrawCopy("SAME");
    TsipmC->Draw("SAME PE");
    pol1_c->DrawCopy("SAME");

    TsipmA->SetMarkerStyle(20);
    TsipmA->SetMarkerColor(kRed);
    TsipmA->SetLineColor(kRed);
    TsipmB->SetMarkerStyle(20);
    TsipmB->SetMarkerColor(kBlue);
    TsipmB->SetLineColor(kBlue);
    TsipmC->SetMarkerStyle(20);
    TsipmC->SetMarkerColor(kGreen - 1);
    TsipmC->SetLineColor(kGreen - 1);

    TLatex *lLatex = new TLatex();
    lLatex->SetTextSize(0.03);
    lLatex->SetTextColor(kRed);
    lLatex->DrawLatexNDC(0.45, 0.21, Form("Sipm A q: %.2f#pm%.2f m: %.2f#pm%.2f", pol1_a->GetParameter(0), pol1_a->GetParError(0), pol1_a->GetParameter(1), pol1_a->GetParError(1)));
    lLatex->SetTextColor(kBlue);
    lLatex->DrawLatexNDC(0.45, 0.18, Form("Sipm B q: %.2f#pm%.2f m: %.2f#pm%.2f", pol1_b->GetParameter(0), pol1_b->GetParError(0), pol1_b->GetParameter(1), pol1_b->GetParError(1)));
    lLatex->SetTextColor(kGreen - 1);
    lLatex->DrawLatexNDC(0.45, 0.15, Form("Sipm C q: %.2f#pm%.2f m: %.2f#pm%.2f", pol1_c->GetParameter(0), pol1_c->GetParError(0), pol1_c->GetParameter(1), pol1_c->GetParError(1)));

    c_sipm_target_reached->SaveAs("c_sipm_target_reached.png");

    cDupmp->cd();

    TsipmA_all->Fit(fpol1, "N");
    pol1_a = (TF1 *)fpol1->Clone("pol1_a");
    pol1_a->SetLineStyle(kDashed);
    pol1_a->SetLineColor(kRed);
    pol1_a->SetRange(-10, 100);

    TsipmB_all->Fit(fpol1, "N");
    pol1_b = (TF1 *)fpol1->Clone("pol1_b");
    pol1_b->SetLineStyle(kDashed);
    pol1_b->SetLineColor(kBlue);
    pol1_b->SetRange(-10, 100);

    TsipmC_all->Fit(fpol1, "N");
    pol1_c = (TF1 *)fpol1->Clone("pol1_c");
    pol1_c->SetLineStyle(kDashed);
    pol1_c->SetLineColor(kGreen - 1);
    pol1_c->SetRange(-10, 100);

    TCanvas *c_sipm_target_reached_all = new TCanvas("c_sipm_target_reached_all", "", 750, 750);
    gPad->SetGridy();
    hframe->SetMaximum(70);
    hframe->Draw();

    TsipmA_all->Draw("SAME PE");
    pol1_a->DrawCopy("SAME");
    TsipmB_all->Draw("SAME PE");
    pol1_b->DrawCopy("SAME");
    TsipmC_all->Draw("SAME PE");
    pol1_c->DrawCopy("SAME");
    TsipmA_all->SetMarkerStyle(20);
    TsipmA_all->SetMarkerColor(kRed);
    TsipmA_all->SetLineColor(kRed);
    TsipmB_all->SetMarkerStyle(20);
    TsipmB_all->SetMarkerColor(kBlue);
    TsipmB_all->SetLineColor(kBlue);
    TsipmC_all->SetMarkerStyle(20);
    TsipmC_all->SetMarkerColor(kGreen - 1);
    TsipmC_all->SetLineColor(kGreen - 1);
    lLatex->SetTextSize(0.03);
    lLatex->SetTextColor(kRed);
    lLatex->DrawLatexNDC(0.45, 0.21, Form("Sipm A q: %.2f#pm%.2f m: %.2f#pm%.2f", pol1_a->GetParameter(0), pol1_a->GetParError(0), pol1_a->GetParameter(1), pol1_a->GetParError(1)));
    lLatex->SetTextColor(kBlue);
    lLatex->DrawLatexNDC(0.45, 0.18, Form("Sipm B q: %.2f#pm%.2f m: %.2f#pm%.2f", pol1_b->GetParameter(0), pol1_b->GetParError(0), pol1_b->GetParameter(1), pol1_b->GetParError(1)));
    lLatex->SetTextColor(kGreen - 1);
    lLatex->DrawLatexNDC(0.45, 0.15, Form("Sipm C q: %.2f#pm%.2f m: %.2f#pm%.2f", pol1_c->GetParameter(0), pol1_c->GetParError(0), pol1_c->GetParameter(1), pol1_c->GetParError(1)));

    c_sipm_target_reached_all->SaveAs("c_sipm_target_reached_all.png");

    TCanvas *ctest = new TCanvas();
    TGraphErrors *gtest = new TGraphErrors();
    auto i_ter = -1;
    for (auto points : data_points["tlm73"])
    {
        i_ter++;
        gtest->SetPoint(i_ter, i_ter, points[0]);
        gtest->SetPointError(i_ter, 0, points[1]);
    }
    gtest->Draw("ALPE");

    TCanvas *test = new TCanvas("", "", 750, 750);
    gPad->SetGridy();
    auto hframe2 = gPad->DrawFrame(0., -5., 2.5, 10);
    hframe2->SetTitle(";Power (W);Temperature of bystander sensor - LM73 (#circC)");
    hframe2->Draw();
    TsipmA_Btest->Draw("SAME LP");
    TsipmA_Btest->SetMarkerStyle(20);
    TsipmA_Btest->SetMarkerColor(kRed);
    TsipmA_Btest->SetLineColor(kRed);
    TsipmA_Ctest->Draw("SAME LP");
    TsipmA_Ctest->SetMarkerStyle(20);
    TsipmA_Ctest->SetMarkerColor(kBlue);
    TsipmA_Ctest->SetLineColor(kBlue);
    TsipmB_Atest->Draw("SAME LP");
    TsipmB_Atest->SetMarkerStyle(21);
    TsipmB_Atest->SetMarkerColor(kGreen - 1);
    TsipmB_Atest->SetLineColor(kGreen - 1);
    TsipmB_Ctest->Draw("SAME LP");
    TsipmB_Ctest->SetMarkerStyle(21);
    TsipmB_Ctest->SetMarkerColor(kBlue);
    TsipmB_Ctest->SetLineColor(kBlue);
    TsipmB_Ctest->Draw("SAME LP");
    TsipmC_Atest->SetMarkerStyle(22);
    TsipmC_Atest->SetMarkerColor(kGreen - 1);
    TsipmC_Atest->SetLineColor(kGreen - 1);
    TsipmC_Atest->Draw("SAME LP");
    TsipmC_Btest->SetMarkerStyle(22);
    TsipmC_Btest->SetMarkerColor(kRed);
    TsipmC_Btest->SetLineColor(kRed);
    TsipmC_Btest->Draw("SAME LP");

    TLegend *test_legend = new TLegend(0.10, 0.90, 0.45, 0.70);
    test_legend->AddEntry(TsipmA_Btest, "heated A, spect. B", "EP");
    test_legend->AddEntry(TsipmA_Ctest, "heated A, spect. C", "EP");
    test_legend->AddEntry(TsipmB_Atest, "heated B, spect. A", "EP");
    test_legend->AddEntry(TsipmB_Ctest, "heated B, spect. C", "EP");
    test_legend->AddEntry(TsipmC_Atest, "heated C, spect. A", "EP");
    test_legend->AddEntry(TsipmC_Btest, "heated C, spect. B", "EP");
    test_legend->Draw("SAME");

    // Clean up
    file->Close();
}