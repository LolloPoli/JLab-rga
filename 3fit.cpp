#include <iostream>
#include <cmath>
#include <TStyle.h>
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>

void SetStatsBox1(TH1D* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.40);  // Posizione pannello (sinistra)
        stats0->SetX2NDC(0.60);  // Posizione pannello (destra)
        stats0->SetY1NDC(0.48);  // Posizione pannello (basso)
        stats0->SetY2NDC(0.82);  // Posizione pannello (alto)
    }
}

void FourierAnalysis(const char* inputFileName, const char* outputFileName) {
    // Apertura del file ROOT di input
    TFile* inputFile = TFile::Open(inputFileName);
    //if (inputFile) inputFile->ls();
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Errore: impossibile aprire il file " << inputFileName << std::endl;
        return;
    }

    TTree *treeKp = (TTree*)inputFile->Get("Kaon+");
    TTree *treeKm = (TTree*)inputFile->Get("Kaon-");
    TDirectory *Kp_bin = (TDirectory*)inputFile->Get("KaonPlus_Bin");

    TH1D *h_phi_Kp = new TH1D("h_phi_Kp", "Cahn and Boer-Mulders analysis | Kaon+", 60, -TMath::Pi(), TMath::Pi());
    treeKp->Draw("Phi_h >> h_phi_Kp");
    TH1D *h_phi_Kp_Q2 = new TH1D("h_phi_Kp_Q2", "Q^{2} distribution | Kaon+; Q^{2} [GeV^{2}]; counts", 60, 1, 9);
    treeKp->Draw("Q2 >> h_phi_Kp_Q2");
    TH1D *h_phi_Kp_y = new TH1D("h_phi_Kp_y", "y distribution | Kaon+; y; counts", 60, 0.25, 0.75);
    treeKp->Draw("y >> h_phi_Kp_y");
    TH1D *h_phi_Km = new TH1D("h_phi_Km", "Cahn and Boer-Mulders analysis | Kaon-", 60, -TMath::Pi(), TMath::Pi());
    treeKm->Draw("Phi_h >> h_phi_Km");
    // bin
    // z
    TH1D *kp_phi_z1 = (TH1D*)Kp_bin->Get("_Phi_z_1");
    TH1D *kp_phi_z2 = (TH1D*)Kp_bin->Get("_Phi_z_2");
    TH1D *kp_phi_z3 = (TH1D*)Kp_bin->Get("_Phi_z_3");
    TH1D *kp_phi_z4 = (TH1D*)Kp_bin->Get("_Phi_z_4");
    // PhT
    TH1D *kp_phi_PhT1 = (TH1D*)Kp_bin->Get("_Phi_PhT_1");
    TH1D *kp_phi_PhT2 = (TH1D*)Kp_bin->Get("_Phi_PhT_2");
    TH1D *kp_phi_PhT3 = (TH1D*)Kp_bin->Get("_Phi_PhT_3");
    TH1D *kp_phi_PhT4 = (TH1D*)Kp_bin->Get("_Phi_PhT_4");
    // Q^2
    TH1D *kp_phi_Q21 = (TH1D*)Kp_bin->Get("_Phi_Q2_1");
    TH1D *kp_phi_Q22 = (TH1D*)Kp_bin->Get("_Phi_Q2_2");
    TH1D *kp_phi_Q23 = (TH1D*)Kp_bin->Get("_Phi_Q2_3");
    TH1D *kp_phi_Q24 = (TH1D*)Kp_bin->Get("_Phi_Q2_4");
    // xB
    TH1D *kp_phi_xB1 = (TH1D*)Kp_bin->Get("_Phi_xB_1");
    TH1D *kp_phi_xB2 = (TH1D*)Kp_bin->Get("_Phi_xB_2");
    TH1D *kp_phi_xB3 = (TH1D*)Kp_bin->Get("_Phi_xB_3");
    TH1D *kp_phi_xB4 = (TH1D*)Kp_bin->Get("_Phi_xB_4");
    // Q^2 and y for the <k^2_T> extrapolation
    // z
    TH1D *kp_phi_z1_Q2 = (TH1D*)Kp_bin->Get("_Phi_z_1_Q2");
    TH1D *kp_phi_z2_Q2 = (TH1D*)Kp_bin->Get("_Phi_z_2_Q2");
    TH1D *kp_phi_z3_Q2 = (TH1D*)Kp_bin->Get("_Phi_z_3_Q2");
    TH1D *kp_phi_z4_Q2 = (TH1D*)Kp_bin->Get("_Phi_z_4_Q2");
    TH1D *kp_phi_z1_y = (TH1D*)Kp_bin->Get("_Phi_z_1_y");
    TH1D *kp_phi_z2_y = (TH1D*)Kp_bin->Get("_Phi_z_2_y");
    TH1D *kp_phi_z3_y = (TH1D*)Kp_bin->Get("_Phi_z_3_y");
    TH1D *kp_phi_z4_y = (TH1D*)Kp_bin->Get("_Phi_z_4_y");
    // PhT
    TH1D *kp_phi_PhT1_Q2 = (TH1D*)Kp_bin->Get("_Phi_PhT_1_Q2");
    TH1D *kp_phi_PhT2_Q2 = (TH1D*)Kp_bin->Get("_Phi_PhT_2_Q2");
    TH1D *kp_phi_PhT3_Q2 = (TH1D*)Kp_bin->Get("_Phi_PhT_3_Q2");
    TH1D *kp_phi_PhT4_Q2 = (TH1D*)Kp_bin->Get("_Phi_PhT_4_Q2");
    TH1D *kp_phi_PhT1_y = (TH1D*)Kp_bin->Get("_Phi_PhT_1_y");
    TH1D *kp_phi_PhT2_y = (TH1D*)Kp_bin->Get("_Phi_PhT_2_y");
    TH1D *kp_phi_PhT3_y = (TH1D*)Kp_bin->Get("_Phi_PhT_3_y");
    TH1D *kp_phi_PhT4_y = (TH1D*)Kp_bin->Get("_Phi_PhT_4_y");
    // Q^2
    TH1D *kp_phi_Q21_Q2 = (TH1D*)Kp_bin->Get("_Phi_Q2_1_Q2");
    TH1D *kp_phi_Q22_Q2 = (TH1D*)Kp_bin->Get("_Phi_Q2_2_Q2");
    TH1D *kp_phi_Q23_Q2 = (TH1D*)Kp_bin->Get("_Phi_Q2_3_Q2");
    TH1D *kp_phi_Q24_Q2 = (TH1D*)Kp_bin->Get("_Phi_Q2_4_Q2");
    TH1D *kp_phi_Q21_y = (TH1D*)Kp_bin->Get("_Phi_Q2_1_y");
    TH1D *kp_phi_Q22_y = (TH1D*)Kp_bin->Get("_Phi_Q2_2_y");
    TH1D *kp_phi_Q23_y = (TH1D*)Kp_bin->Get("_Phi_Q2_3_y");
    TH1D *kp_phi_Q24_y = (TH1D*)Kp_bin->Get("_Phi_Q2_4_y");
    // xB
    TH1D *kp_phi_xB1_Q2 = (TH1D*)Kp_bin->Get("_Phi_xB_1_Q2");
    TH1D *kp_phi_xB2_Q2 = (TH1D*)Kp_bin->Get("_Phi_xB_2_Q2");
    TH1D *kp_phi_xB3_Q2 = (TH1D*)Kp_bin->Get("_Phi_xB_3_Q2");
    TH1D *kp_phi_xB4_Q2 = (TH1D*)Kp_bin->Get("_Phi_xB_4_Q2");
    TH1D *kp_phi_xB1_y = (TH1D*)Kp_bin->Get("_Phi_xB_1_y");
    TH1D *kp_phi_xB2_y = (TH1D*)Kp_bin->Get("_Phi_xB_2_y");
    TH1D *kp_phi_xB3_y = (TH1D*)Kp_bin->Get("_Phi_xB_3_y");
    TH1D *kp_phi_xB4_y = (TH1D*)Kp_bin->Get("_Phi_xB_4_y");

    if (h_phi_Kp->GetEntries() == 0) {
    std::cerr << "Errore: l'istogramma 'Kaon+ -> Phi_h' è vuoto." << std::endl;
    inputFile->Close();
    return;
    }

    if (h_phi_Km->GetEntries() == 0) {
    std::cerr << "Errore: l'istogramma 'Kaon- -> Phi_h' è vuoto." << std::endl;
    inputFile->Close();
    return;
    }

    // NORMALIZZO L'ISTOGRAMMA
    h_phi_Kp->Scale(1.0 / h_phi_Kp->Integral());
    h_phi_Km->Scale(1.0 / h_phi_Km->Integral());
    kp_phi_z1->Scale(1.0 / kp_phi_z1->Integral());
    kp_phi_z2->Scale(1.0 / kp_phi_z2->Integral());
    kp_phi_z3->Scale(1.0 / kp_phi_z3->Integral());
    kp_phi_z4->Scale(1.0 / kp_phi_z4->Integral());
    kp_phi_PhT1->Scale(1.0 / kp_phi_PhT1->Integral());
    kp_phi_PhT2->Scale(1.0 / kp_phi_PhT2->Integral());
    kp_phi_PhT3->Scale(1.0 / kp_phi_PhT3->Integral());
    kp_phi_PhT4->Scale(1.0 / kp_phi_PhT4->Integral());
    kp_phi_Q21->Scale(1.0 / kp_phi_Q21->Integral());
    kp_phi_Q22->Scale(1.0 / kp_phi_Q22->Integral());
    kp_phi_Q23->Scale(1.0 / kp_phi_Q23->Integral());
    kp_phi_Q24->Scale(1.0 / kp_phi_Q24->Integral());
    kp_phi_xB1->Scale(1.0 / kp_phi_xB1->Integral());
    kp_phi_xB2->Scale(1.0 / kp_phi_xB2->Integral());
    kp_phi_xB3->Scale(1.0 / kp_phi_xB3->Integral());
    kp_phi_xB4->Scale(1.0 / kp_phi_xB4->Integral());

    // Creazione del file ROOT di output
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Errore: impossibile creare il file " << outputFileName << std::endl;
        inputFile->Close();
        return;
    }

    // Creazione del canvas
    TCanvas* c = new TCanvas("Phi_h_Kaon+", "Cahn and Boer-Mulders analysis", 800, 600);
    c->cd();
    c->ToggleEditor();
    // Creazione del fit
    TF1* fitFunc = new TF1("fitFunc", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFunc->SetParameters(0, h_phi_Kp->GetMaximum());
    const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFunc->SetParName(i, paramNames[i]);
    }
    //fitFunc->SetLineColor(kRed);
    h_phi_Kp->GetXaxis()->SetTitle("#Phi_h [Rad]");
    h_phi_Kp->GetYaxis()->SetTitle("Counts");
    h_phi_Kp->Draw("E");
    h_phi_Kp->Fit(fitFunc, "RW");
    fitFunc->Draw("SAME");
    gStyle->SetOptFit(111);
    TPaveStats* stats3 = (TPaveStats*)h_phi_Kp->GetListOfFunctions()->FindObject("stats");
    SetStatsBox1(h_phi_Kp);
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextFont(11);
    double A0_tot = fitFunc->GetParameter(0);
    double A1_tot = fitFunc->GetParameter(1);
    double A2_tot = fitFunc->GetParameter(2);
    double A0_tot_err = fitFunc->GetParError(0);
    double A1_tot_err = fitFunc->GetParError(1);
    double A2_tot_err = fitFunc->GetParError(2);
    double ratio1_tot = A1_tot / A0_tot;
    double ratio2_tot = A2_tot / A0_tot;
    double mean_Q2_tot = h_phi_Kp_Q2->GetMean();
    double mean_y_tot = h_phi_Kp_y->GetMean();
    double mean_Q2_tot_err = h_phi_Kp_Q2->GetMeanError();
    double mean_y_tot_err = h_phi_Kp_y->GetMeanError();
    double num_tot = (1 - mean_y_tot + (mean_y_tot * mean_y_tot*0.5));
    double den_tot = (1 - mean_y_tot);
    double kT_tot = (ratio1_tot * mean_Q2_tot * 0.5)*((1 - mean_y_tot + (mean_y_tot * mean_y_tot*0.5))/(1 - mean_y_tot));
    double ratio1_tot_err = fabs(ratio1_tot) * sqrt(pow(A1_tot_err / A1_tot, 2) + pow(A0_tot_err / A0_tot, 2));
    double ratio2_tot_err = fabs(ratio2_tot) * sqrt(pow(A2_tot_err / A2_tot, 2) + pow(A0_tot_err / A0_tot, 2));
    double dkT_A_tot = 0.5 * mean_Q2_tot * (num_tot / den_tot);
    double dkT_Q_tot = 0.5 * ratio1_tot * (num_tot / den_tot);
    double dkT_y_tot = 0.5 * ratio1_tot * mean_Q2_tot * ((0.5 * mean_y_tot * 2 - 1) * den_tot - num_tot * (-1)) / (den_tot * den_tot);
    double sigma_kT_tot = sqrt(pow(dkT_A_tot * ratio1_tot_err, 2) + pow(dkT_Q_tot * mean_Q2_tot_err, 2) + pow(dkT_y_tot * mean_y_tot_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,z_{1}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_tot, ratio1_tot_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,z_{1}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_tot, ratio2_tot_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{z}_{1} = %.4f #pm %.4f GeV^{2}", kT_tot, sigma_kT_tot));

    // Disegna i residui per controllo
    TCanvas* c2 = new TCanvas("Residui_Kaon+", "Residui del Fit", 800, 400);
    TH1F* h_residui_Kp = (TH1F*)h_phi_Kp->Clone("h_residui_Kp");
    h_residui_Kp->Add(fitFunc, -1);
    h_residui_Kp->SetTitle("Residui del Fit; #Phi_h [rad]; Data - Fit");
    h_residui_Kp->Draw("E");
    TPaveStats* stats4 = (TPaveStats*)h_residui_Kp->GetListOfFunctions()->FindObject("stats");
    if (stats4) {
        stats4->SetX1NDC(0.40);     // Sposta il pannello
        stats4->SetX2NDC(0.60);     
        stats4->SetY1NDC(0.60);     
        stats4->SetY2NDC(0.88);     
    }
    gPad->Update();
    c->Update();
    c->Write(); // Scrive il canvas nel file di output

    // Kaon-

    TCanvas* c3 = new TCanvas("Phi_h_Kaon-", "Cahn and Boer-Mulders analysis", 800, 600);
    c3->ToggleEditor();
    // Creazione del fit
    TF1* fitFunc2 = new TF1("fitFunc2", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Km->GetXaxis()->GetXmin(), h_phi_Km->GetXaxis()->GetXmax());
    fitFunc2->SetParameters(0, h_phi_Kp->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFunc2->SetParName(i, paramNames[i]);
    }
    fitFunc2->SetLineColor(kRed);
    h_phi_Km->GetXaxis()->SetTitle("#Phi_h [Rad]");
    h_phi_Km->GetYaxis()->SetTitle("Counts");
    h_phi_Km->Draw("E");
    h_phi_Km->Fit(fitFunc2, "RW");
    fitFunc2->Draw("SAME");
    gStyle->SetOptFit(111);
    TPaveStats* stats6 = (TPaveStats*)h_phi_Km->GetListOfFunctions()->FindObject("stats");
    if (stats6) {
        stats6->SetX1NDC(0.39);     // Sposta il pannello
        stats6->SetX2NDC(0.61);     
        stats6->SetY1NDC(0.55);     
        stats6->SetY2NDC(0.85);     
    }
    double A0_tot2 = fitFunc2->GetParameter(0);
    double A1_tot2 = fitFunc2->GetParameter(1);
    double A2_tot2 = fitFunc2->GetParameter(2);
    double ratio1_tot2 = A1_tot2 / A0_tot2;
    double ratio2_tot2 = A2_tot2 / A0_tot2;
    latex.DrawLatex(0.15, 0.85, Form("A_{UU}^{cos(#Phi_{h})} = %.4f", ratio1_tot2));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU}^{cos(2#Phi_{h})} = %.4f", ratio2_tot2));
    // Disegna i residui per controllo
    TCanvas* c4 = new TCanvas("Residui_Kaon-", "Residui del Fit", 800, 400);
    TH1F* h_residui_Km = (TH1F*)h_phi_Km->Clone("h_residui_Km");
    h_residui_Km->Add(fitFunc2, -1);
    h_residui_Km->SetTitle("Residui del Fit; #Phi_h [rad]; Data - Fit");
    h_residui_Km->Draw("E");
    TPaveStats* stats5 = (TPaveStats*)h_residui_Km->GetListOfFunctions()->FindObject("stats");
    if (stats5) {
        stats5->SetX1NDC(0.40);     // Sposta il pannello
        stats5->SetX2NDC(0.60);     
        stats5->SetY1NDC(0.60);     
        stats5->SetY2NDC(0.88);     
    }
    gPad->Update();
    c3->Update();
    c3->Write(); // Scrive il canvas nel file di output

    // _____________________________________________________    z    ____________________________________________________________________________
    // z1
    TCanvas* cz1 = new TCanvas("Phi_h_Kaon+_z_1", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cz1->cd();
    cz1->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncz1 = new TF1("fitFuncz1", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncz1->SetParameters(0, kp_phi_z1->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncz1->SetParName(i, paramNames[i]);
    }
    kp_phi_z1->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_z1->GetYaxis()->SetTitle("Counts");
    kp_phi_z1->Draw("E");
    kp_phi_z1->Fit(fitFuncz1, "RW");
    fitFuncz1->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_z1);
    double A0_z1 = fitFuncz1->GetParameter(0);
    double A1_z1 = fitFuncz1->GetParameter(1);
    double A2_z1 = fitFuncz1->GetParameter(2);
    double A0_z1_err = fitFuncz1->GetParError(0);
    double A1_z1_err = fitFuncz1->GetParError(1);
    double A2_z1_err = fitFuncz1->GetParError(2);
    double ratio1_z1 = A1_z1 / A0_z1;
    double ratio2_z1 = A2_z1 / A0_z1;
    double mean_Q2_z1 = kp_phi_z1_Q2->GetMean();
    double mean_y_z1 = kp_phi_z1_y->GetMean();
    double mean_Q2_z1_err = kp_phi_z1_Q2->GetMeanError();
    double mean_y_z1_err = kp_phi_z1_y->GetMeanError();
    double num_z1 = (1 - mean_y_z1 + (mean_y_z1 * mean_y_z1*0.5));
    double den_z1 = (1 - mean_y_z1);
    double kT_z1 = (ratio1_z1 * mean_Q2_z1 * 0.5)*(num_z1/den_z1);
    // calcolo gli errori con le derivate parziali
    double ratio1_z1_err = fabs(ratio1_z1) * sqrt(pow(A1_z1_err / A1_z1, 2) + pow(A0_z1_err / A0_z1, 2));
    double ratio2_z1_err = fabs(ratio2_z1) * sqrt(pow(A2_z1_err / A2_z1, 2) + pow(A0_z1_err / A0_z1, 2));
    double dkT_A_z1 = 0.5 * mean_Q2_z1 * (num_z1 / den_z1);
    double dkT_Q_z1 = 0.5 * ratio1_z1 * (num_z1 / den_z1);
    double dkT_y_z1 = 0.5 * ratio1_z1 * mean_Q2_z1 * ((0.5 * mean_y_z1 * 2 - 1) * den_z1 - num_z1 * (-1)) / (den_z1 * den_z1);
    double sigma_kT_z1 = sqrt(pow(dkT_A_z1 * ratio1_z1_err, 2) + pow(dkT_Q_z1 * mean_Q2_z1_err, 2) + pow(dkT_y_z1 * mean_y_z1_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,z_{1}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_z1, ratio1_z1_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,z_{1}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_z1, ratio2_z1_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{z}_{1} = %.4f #pm %.4f GeV^{2}", kT_z1, sigma_kT_z1));
    gPad->Update();
    cz1->Update();
    cz1->Write();
    // z2
    TCanvas* cz2 = new TCanvas("Phi_h_Kaon+_z_2", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cz2->cd();
    cz2->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncz2 = new TF1("fitFuncz2", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncz2->SetParameters(0, kp_phi_z2->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncz2->SetParName(i, paramNames[i]);
    }
    kp_phi_z2->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_z2->GetYaxis()->SetTitle("Counts");
    kp_phi_z2->Draw("E");
    kp_phi_z2->Fit(fitFuncz2, "RW");
    fitFuncz2->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_z2);
    double A0_z2 = fitFuncz2->GetParameter(0);
    double A1_z2 = fitFuncz2->GetParameter(1);
    double A2_z2 = fitFuncz2->GetParameter(2);
    double A0_z2_err = fitFuncz2->GetParError(0);
    double A1_z2_err = fitFuncz2->GetParError(1);
    double A2_z2_err = fitFuncz2->GetParError(2);
    double ratio1_z2 = A1_z2 / A0_z2;
    double ratio2_z2 = A2_z2 / A0_z2;
    double mean_Q2_z2 = kp_phi_z2_Q2->GetMean();
    double mean_y_z2 = kp_phi_z2_y->GetMean();
    double mean_Q2_z2_err = kp_phi_z2_Q2->GetMeanError();
    double mean_y_z2_err = kp_phi_z2_y->GetMeanError();
    double num_z2 = (1 - mean_y_z2 + (mean_y_z2 * mean_y_z2*0.5));
    double den_z2 = (1 - mean_y_z2);
    double kT_z2 = (ratio1_z2 * mean_Q2_z2 * 0.5)*(num_z2/den_z2);
    // calcolo gli errori con le derivate parziali
    double ratio1_z2_err = fabs(ratio1_z2) * sqrt(pow(A1_z2_err / A1_z2, 2) + pow(A0_z2_err / A0_z2, 2));
    double ratio2_z2_err = fabs(ratio2_z2) * sqrt(pow(A2_z2_err / A2_z2, 2) + pow(A0_z2_err / A0_z2, 2));
    double dkT_A_z2 = 0.5 * mean_Q2_z2 * (num_z2 / den_z2);
    double dkT_Q_z2 = 0.5 * ratio1_z2 * (num_z2 / den_z2);
    double dkT_y_z2 = 0.5 * ratio1_z2 * mean_Q2_z2 * ((0.5 * mean_y_z2 * 2 - 1) * den_z2 - num_z2 * (-1)) / (den_z2 * den_z2);
    double sigma_kT_z2 = sqrt(pow(dkT_A_z2 * ratio1_z2_err, 2) + pow(dkT_Q_z2 * mean_Q2_z2_err, 2) + pow(dkT_y_z2 * mean_y_z2_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,z_{2}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_z2, ratio1_z2_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,z_{2}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_z2, ratio2_z2_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{z}_{2} = %.4f #pm %.4f GeV^{2}", kT_z2, sigma_kT_z2));
    gPad->Update();
    cz2->Update();
    cz2->Write();
    // z3
    TCanvas* cz3 = new TCanvas("Phi_h_Kaon+_z_3", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cz3->cd();
    cz3->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncz3 = new TF1("fitFuncz3", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncz3->SetParameters(0, kp_phi_z3->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncz3->SetParName(i, paramNames[i]);
    }
    kp_phi_z3->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_z3->GetYaxis()->SetTitle("Counts");
    kp_phi_z3->Draw("E");
    kp_phi_z3->Fit(fitFuncz3, "RW");
    fitFuncz3->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_z3);
    double A0_z3 = fitFuncz3->GetParameter(0);
    double A1_z3 = fitFuncz3->GetParameter(1);
    double A2_z3 = fitFuncz3->GetParameter(2);
    double A0_z3_err = fitFuncz3->GetParError(0);
    double A1_z3_err = fitFuncz3->GetParError(1);
    double A2_z3_err = fitFuncz3->GetParError(2);
    double ratio1_z3 = A1_z3 / A0_z3;
    double ratio2_z3 = A2_z3 / A0_z3;
    double mean_Q2_z3 = kp_phi_z3_Q2->GetMean();
    double mean_y_z3 = kp_phi_z3_y->GetMean();
    double mean_Q2_z3_err = kp_phi_z3_Q2->GetMeanError();
    double mean_y_z3_err = kp_phi_z3_y->GetMeanError();
    double num_z3 = (1 - mean_y_z3 + (mean_y_z3 * mean_y_z3*0.5));
    double den_z3 = (1 - mean_y_z3);
    double kT_z3 = (ratio1_z3 * mean_Q2_z3 * 0.5)*(num_z3/den_z3);
    // calcolo gli errori con le derivate parziali
    double ratio1_z3_err = fabs(ratio1_z3) * sqrt(pow(A1_z3_err / A1_z3, 2) + pow(A0_z3_err / A0_z3, 2));
    double ratio2_z3_err = fabs(ratio2_z3) * sqrt(pow(A2_z3_err / A2_z3, 2) + pow(A0_z3_err / A0_z3, 2));
    double dkT_A_z3 = 0.5 * mean_Q2_z3 * (num_z3 / den_z3);
    double dkT_Q_z3 = 0.5 * ratio1_z3 * (num_z3 / den_z3);
    double dkT_y_z3 = 0.5 * ratio1_z3 * mean_Q2_z3 * ((0.5 * mean_y_z3 * 2 - 1) * den_z3 - num_z3 * (-1)) / (den_z3 * den_z3);
    double sigma_kT_z3 = sqrt(pow(dkT_A_z3 * ratio1_z3_err, 2) + pow(dkT_Q_z3 * mean_Q2_z3_err, 2) + pow(dkT_y_z3 * mean_y_z3_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,z_{3}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_z3, ratio1_z3_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,z_{3}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_z3, ratio2_z3_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{z}_{3} = %.4f #pm %.4f GeV^{2}", kT_z3, sigma_kT_z3));
    gPad->Update();
    cz3->Update();
    cz3->Write();
    // z4
    TCanvas* cz4 = new TCanvas("Phi_h_Kaon+_z_4", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cz4->cd();
    cz4->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncz4 = new TF1("fitFuncz4", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncz4->SetParameters(0, kp_phi_z4->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncz4->SetParName(i, paramNames[i]);
    }
    kp_phi_z4->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_z4->GetYaxis()->SetTitle("Counts");
    kp_phi_z4->Draw("E");
    kp_phi_z4->Fit(fitFuncz4, "RW");
    fitFuncz4->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_z4);
    double A0_z4 = fitFuncz4->GetParameter(0);
    double A1_z4 = fitFuncz4->GetParameter(1);
    double A2_z4 = fitFuncz4->GetParameter(2);
    double A0_z4_err = fitFuncz4->GetParError(0);
    double A1_z4_err = fitFuncz4->GetParError(1);
    double A2_z4_err = fitFuncz4->GetParError(2);
    double ratio1_z4 = A1_z4 / A0_z4;
    double ratio2_z4 = A2_z4 / A0_z4;
    double mean_Q2_z4 = kp_phi_z4_Q2->GetMean();
    double mean_y_z4 = kp_phi_z4_y->GetMean();
    double mean_Q2_z4_err = kp_phi_z4_Q2->GetMeanError();
    double mean_y_z4_err = kp_phi_z4_y->GetMeanError();
    double num_z4 = (1 - mean_y_z4 + (mean_y_z4 * mean_y_z4*0.5));
    double den_z4 = (1 - mean_y_z4);
    double kT_z4 = (ratio1_z4 * mean_Q2_z4 * 0.5)*(num_z4/den_z4);
    // calcolo gli errori con le derivate parziali
    double ratio1_z4_err = fabs(ratio1_z4) * sqrt(pow(A1_z4_err / A1_z4, 2) + pow(A0_z4_err / A0_z4, 2));
    double ratio2_z4_err = fabs(ratio2_z4) * sqrt(pow(A2_z4_err / A2_z4, 2) + pow(A0_z4_err / A0_z4, 2));
    double dkT_A_z4 = 0.5 * mean_Q2_z4 * (num_z4 / den_z4);
    double dkT_Q_z4 = 0.5 * ratio1_z4 * (num_z4 / den_z4);
    double dkT_y_z4 = 0.5 * ratio1_z4 * mean_Q2_z4 * ((0.5 * mean_y_z4 * 2 - 1) * den_z4 - num_z4 * (-1)) / (den_z4 * den_z4);
    double sigma_kT_z4 = sqrt(pow(dkT_A_z4 * ratio1_z4_err, 2) + pow(dkT_Q_z4 * mean_Q2_z4_err, 2) + pow(dkT_y_z4 * mean_y_z4_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,z_{4}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_z4, ratio1_z4_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,z_{4}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_z4, ratio2_z4_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{z}_{4} = %.4f #pm %.4f GeV^{2}", kT_z4, sigma_kT_z4));
    gPad->Update();
    cz4->Update();
    cz4->Write();

    // ___________________________________________________________________________________________________________________________________________    
    // _____________________________________________________    PhT    ____________________________________________________________________________
    // PhT1
    TCanvas* cPhT1 = new TCanvas("Phi_h_Kaon+_PhT_1", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cPhT1->cd();
    cPhT1->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncPhT1 = new TF1("fitFuncPhT1", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncPhT1->SetParameters(0, kp_phi_PhT1->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncPhT1->SetParName(i, paramNames[i]);
    }
    kp_phi_PhT1->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_PhT1->GetYaxis()->SetTitle("Counts");
    kp_phi_PhT1->Draw("E");
    kp_phi_PhT1->Fit(fitFuncPhT1, "RW");
    fitFuncPhT1->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_PhT1);
    double A0_PhT1 = fitFuncPhT1->GetParameter(0);
    double A1_PhT1 = fitFuncPhT1->GetParameter(1);
    double A2_PhT1 = fitFuncPhT1->GetParameter(2);
    double A0_PhT1_err = fitFuncPhT1->GetParError(0);
    double A1_PhT1_err = fitFuncPhT1->GetParError(1);
    double A2_PhT1_err = fitFuncPhT1->GetParError(2);
    double ratio1_PhT1 = A1_PhT1 / A0_PhT1;
    double ratio2_PhT1 = A2_PhT1 / A0_PhT1;
    double mean_Q2_PhT1 = kp_phi_PhT1_Q2->GetMean();
    double mean_y_PhT1 = kp_phi_PhT1_y->GetMean();
    double mean_Q2_PhT1_err = kp_phi_PhT1_Q2->GetMeanError();
    double mean_y_PhT1_err = kp_phi_PhT1_y->GetMeanError();
    double num_PhT1 = (1 - mean_y_PhT1 + (mean_y_PhT1 * mean_y_PhT1*0.5));
    double den_PhT1 = (1 - mean_y_PhT1);
    double kT_PhT1 = (ratio1_PhT1 * mean_Q2_PhT1 * 0.5)*(num_PhT1/den_PhT1);
    // calcolo gli errori con le derivate parziali
    double ratio1_PhT1_err = fabs(ratio1_PhT1) * sqrt(pow(A1_PhT1_err / A1_PhT1, 2) + pow(A0_PhT1_err / A0_PhT1, 2));
    double ratio2_PhT1_err = fabs(ratio2_PhT1) * sqrt(pow(A2_PhT1_err / A2_PhT1, 2) + pow(A0_PhT1_err / A0_PhT1, 2));
    double dkT_A_PhT1 = 0.5 * mean_Q2_PhT1 * (num_PhT1 / den_PhT1);
    double dkT_Q_PhT1 = 0.5 * ratio1_PhT1 * (num_PhT1 / den_PhT1);
    double dkT_y_PhT1 = 0.5 * ratio1_PhT1 * mean_Q2_PhT1 * ((0.5 * mean_y_PhT1 * 2 - 1) * den_PhT1 - num_PhT1 * (-1)) / (den_PhT1 * den_PhT1);
    double sigma_kT_PhT1 = sqrt(pow(dkT_A_PhT1 * ratio1_PhT1_err, 2) + pow(dkT_Q_PhT1 * mean_Q2_PhT1_err, 2) + pow(dkT_y_PhT1 * mean_y_PhT1_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,P_{hT,1}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_PhT1, ratio1_PhT1_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,P_{hT,1}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_PhT1, ratio2_PhT1_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{P_{hT}}_{1} = %.4f #pm %.4f GeV^{2}", kT_PhT1, sigma_kT_PhT1));
    gPad->Update();
    cPhT1->Update();
    cPhT1->Write();
    // PhT2
    TCanvas* cPhT2 = new TCanvas("Phi_h_Kaon+_PhT_2", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cPhT2->cd();
    cPhT2->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncPhT2 = new TF1("fitFuncPhT2", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncPhT2->SetParameters(0, kp_phi_PhT2->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncPhT2->SetParName(i, paramNames[i]);
    }
    kp_phi_PhT2->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_PhT2->GetYaxis()->SetTitle("Counts");
    kp_phi_PhT2->Draw("E");
    kp_phi_PhT2->Fit(fitFuncPhT2, "RW");
    fitFuncPhT2->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_PhT2);
    double A0_PhT2 = fitFuncPhT2->GetParameter(0);
    double A1_PhT2 = fitFuncPhT2->GetParameter(1);
    double A2_PhT2 = fitFuncPhT2->GetParameter(2);
    double A0_PhT2_err = fitFuncPhT2->GetParError(0);
    double A1_PhT2_err = fitFuncPhT2->GetParError(1);
    double A2_PhT2_err = fitFuncPhT2->GetParError(2);
    double ratio1_PhT2 = A1_PhT2 / A0_PhT2;
    double ratio2_PhT2 = A2_PhT2 / A0_PhT2;
    double mean_Q2_PhT2 = kp_phi_PhT2_Q2->GetMean();
    double mean_y_PhT2 = kp_phi_PhT2_y->GetMean();
    double mean_Q2_PhT2_err = kp_phi_PhT2_Q2->GetMeanError();
    double mean_y_PhT2_err = kp_phi_PhT2_y->GetMeanError();
    double num_PhT2 = (1 - mean_y_PhT2 + (mean_y_PhT2 * mean_y_PhT2*0.5));
    double den_PhT2 = (1 - mean_y_PhT2);
    double kT_PhT2 = (ratio1_PhT2 * mean_Q2_PhT2 * 0.5)*(num_PhT2/den_PhT2);
    // calcolo gli errori con le derivate parziali
    double ratio1_PhT2_err = fabs(ratio1_PhT2) * sqrt(pow(A1_PhT2_err / A1_PhT2, 2) + pow(A0_PhT2_err / A0_PhT2, 2));
    double ratio2_PhT2_err = fabs(ratio2_PhT2) * sqrt(pow(A2_PhT2_err / A2_PhT2, 2) + pow(A0_PhT2_err / A0_PhT2, 2));
    double dkT_A_PhT2 = 0.5 * mean_Q2_PhT2 * (num_PhT2 / den_PhT2);
    double dkT_Q_PhT2 = 0.5 * ratio1_PhT2 * (num_PhT2 / den_PhT2);
    double dkT_y_PhT2 = 0.5 * ratio1_PhT2 * mean_Q2_PhT2 * ((0.5 * mean_y_PhT2 * 2 - 1) * den_PhT2 - num_PhT2 * (-1)) / (den_PhT2 * den_PhT2);
    double sigma_kT_PhT2 = sqrt(pow(dkT_A_PhT2 * ratio1_PhT2_err, 2) + pow(dkT_Q_PhT2 * mean_Q2_PhT2_err, 2) + pow(dkT_y_PhT2 * mean_y_PhT2_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,P_{hT,2}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_PhT2, ratio1_PhT2_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,P_{hT,2}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_PhT2, ratio2_PhT2_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{P_{hT}}_{2} = %.4f #pm %.4f GeV^{2}", kT_PhT2, sigma_kT_PhT2));
    gPad->Update();
    cPhT2->Update();
    cPhT2->Write();
    // PhT3
    TCanvas* cPhT3 = new TCanvas("Phi_h_Kaon+_PhT_3", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cPhT3->cd();
    cPhT3->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncPhT3 = new TF1("fitFuncPhT3", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncPhT3->SetParameters(0, kp_phi_PhT3->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncPhT3->SetParName(i, paramNames[i]);
    }
    kp_phi_PhT3->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_PhT3->GetYaxis()->SetTitle("Counts");
    kp_phi_PhT3->Draw("E");
    kp_phi_PhT3->Fit(fitFuncPhT3, "RW");
    fitFuncPhT3->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_PhT3);
    double A0_PhT3 = fitFuncPhT3->GetParameter(0);
    double A1_PhT3 = fitFuncPhT3->GetParameter(1);
    double A2_PhT3 = fitFuncPhT3->GetParameter(2);
    double A0_PhT3_err = fitFuncPhT3->GetParError(0);
    double A1_PhT3_err = fitFuncPhT3->GetParError(1);
    double A2_PhT3_err = fitFuncPhT3->GetParError(2);
    double ratio1_PhT3 = A1_PhT3 / A0_PhT3;
    double ratio2_PhT3 = A2_PhT3 / A0_PhT3;
    double mean_Q2_PhT3 = kp_phi_PhT3_Q2->GetMean();
    double mean_y_PhT3 = kp_phi_PhT3_y->GetMean();
    double mean_Q2_PhT3_err = kp_phi_PhT3_Q2->GetMeanError();
    double mean_y_PhT3_err = kp_phi_PhT3_y->GetMeanError();
    double num_PhT3 = (1 - mean_y_PhT3 + (mean_y_PhT3 * mean_y_PhT3*0.5));
    double den_PhT3 = (1 - mean_y_PhT3);
    double kT_PhT3 = (ratio1_PhT3 * mean_Q2_PhT3 * 0.5)*(num_PhT3/den_PhT3);
    // calcolo gli errori con le derivate parziali
    double ratio1_PhT3_err = fabs(ratio1_PhT3) * sqrt(pow(A1_PhT3_err / A1_PhT3, 2) + pow(A0_PhT3_err / A0_PhT3, 2));
    double ratio2_PhT3_err = fabs(ratio2_PhT3) * sqrt(pow(A2_PhT3_err / A2_PhT3, 2) + pow(A0_PhT3_err / A0_PhT3, 2));
    double dkT_A_PhT3 = 0.5 * mean_Q2_PhT3 * (num_PhT3 / den_PhT3);
    double dkT_Q_PhT3 = 0.5 * ratio1_PhT3 * (num_PhT3 / den_PhT3);
    double dkT_y_PhT3 = 0.5 * ratio1_PhT3 * mean_Q2_PhT3 * ((0.5 * mean_y_PhT3 * 2 - 1) * den_PhT3 - num_PhT3 * (-1)) / (den_PhT3 * den_PhT3);
    double sigma_kT_PhT3 = sqrt(pow(dkT_A_PhT3 * ratio1_PhT3_err, 2) + pow(dkT_Q_PhT3 * mean_Q2_PhT3_err, 2) + pow(dkT_y_PhT3 * mean_y_PhT3_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,P_{hT,3}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_PhT3, ratio1_PhT3_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,P_{hT,3}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_PhT3, ratio2_PhT3_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{P_{hT}}_{3} = %.4f #pm %.4f GeV^{2}", kT_PhT3, sigma_kT_PhT3));
    gPad->Update();
    cPhT3->Update();
    cPhT3->Write();
    // PhT4
    TCanvas* cPhT4 = new TCanvas("Phi_h_Kaon+_PhT_4", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cPhT4->cd();
    cPhT4->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncPhT4 = new TF1("fitFuncPhT4", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncPhT4->SetParameters(0, kp_phi_PhT4->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncPhT4->SetParName(i, paramNames[i]);
    }
    kp_phi_PhT4->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_PhT4->GetYaxis()->SetTitle("Counts");
    kp_phi_PhT4->Draw("E");
    kp_phi_PhT4->Fit(fitFuncPhT4, "RW");
    fitFuncPhT4->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_PhT4);
    double A0_PhT4 = fitFuncPhT4->GetParameter(0);
    double A1_PhT4 = fitFuncPhT4->GetParameter(1);
    double A2_PhT4 = fitFuncPhT4->GetParameter(2);
    double A0_PhT4_err = fitFuncPhT4->GetParError(0);
    double A1_PhT4_err = fitFuncPhT4->GetParError(1);
    double A2_PhT4_err = fitFuncPhT4->GetParError(2);
    double ratio1_PhT4 = A1_PhT4 / A0_PhT4;
    double ratio2_PhT4 = A2_PhT4 / A0_PhT4;
    double mean_Q2_PhT4 = kp_phi_PhT4_Q2->GetMean();
    double mean_y_PhT4 = kp_phi_PhT4_y->GetMean();
    double mean_Q2_PhT4_err = kp_phi_PhT4_Q2->GetMeanError();
    double mean_y_PhT4_err = kp_phi_PhT4_y->GetMeanError();
    double num_PhT4 = (1 - mean_y_PhT4 + (mean_y_PhT4 * mean_y_PhT4*0.5));
    double den_PhT4 = (1 - mean_y_PhT4);
    double kT_PhT4 = (ratio1_PhT4 * mean_Q2_PhT4 * 0.5)*(num_PhT4/den_PhT4);
    // calcolo gli errori con le derivate parziali
    double ratio1_PhT4_err = fabs(ratio1_PhT4) * sqrt(pow(A1_PhT4_err / A1_PhT4, 2) + pow(A0_PhT4_err / A0_PhT4, 2));
    double ratio2_PhT4_err = fabs(ratio2_PhT4) * sqrt(pow(A2_PhT4_err / A2_PhT4, 2) + pow(A0_PhT4_err / A0_PhT4, 2));
    double dkT_A_PhT4 = 0.5 * mean_Q2_PhT4 * (num_PhT4 / den_PhT4);
    double dkT_Q_PhT4 = 0.5 * ratio1_PhT4 * (num_PhT4 / den_PhT4);
    double dkT_y_PhT4 = 0.5 * ratio1_PhT4 * mean_Q2_PhT4 * ((0.5 * mean_y_PhT4 * 2 - 1) * den_PhT4 - num_PhT4 * (-1)) / (den_PhT4 * den_PhT4);
    double sigma_kT_PhT4 = sqrt(pow(dkT_A_PhT4 * ratio1_PhT4_err, 2) + pow(dkT_Q_PhT4 * mean_Q2_PhT4_err, 2) + pow(dkT_y_PhT4 * mean_y_PhT4_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,P_{hT,4}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_PhT4, ratio1_PhT4_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,P_{hT,4}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_PhT4, ratio2_PhT4_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{P_{hT}}_{4} = %.4f #pm %.4f GeV^{2}", kT_PhT4, sigma_kT_PhT4));
    gPad->Update();
    cPhT4->Update();
    cPhT4->Write();

    // ___________________________________________________________________________________________________________________________________________    
    // _____________________________________________________    Q2    ____________________________________________________________________________
    // Q21
    TCanvas* cQ21 = new TCanvas("Phi_h_Kaon+_Q2_1", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cQ21->cd();
    cQ21->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncQ21 = new TF1("fitFuncQ21", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncQ21->SetParameters(0, kp_phi_Q21->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncQ21->SetParName(i, paramNames[i]);
    }
    kp_phi_Q21->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_Q21->GetYaxis()->SetTitle("Counts");
    kp_phi_Q21->Draw("E");
    kp_phi_Q21->Fit(fitFuncQ21, "RW");
    fitFuncQ21->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_Q21);
    double A0_Q21 = fitFuncQ21->GetParameter(0);
    double A1_Q21 = fitFuncQ21->GetParameter(1);
    double A2_Q21 = fitFuncQ21->GetParameter(2);
    double A0_Q21_err = fitFuncQ21->GetParError(0);
    double A1_Q21_err = fitFuncQ21->GetParError(1);
    double A2_Q21_err = fitFuncQ21->GetParError(2);
    double ratio1_Q21 = A1_Q21 / A0_Q21;
    double ratio2_Q21 = A2_Q21 / A0_Q21;
    double mean_Q2_Q21 = kp_phi_Q21_Q2->GetMean();
    double mean_y_Q21 = kp_phi_Q21_y->GetMean();
    double mean_Q2_Q21_err = kp_phi_Q21_Q2->GetMeanError();
    double mean_y_Q21_err = kp_phi_Q21_y->GetMeanError();
    double num_Q21 = (1 - mean_y_Q21 + (mean_y_Q21 * mean_y_Q21*0.5));
    double den_Q21 = (1 - mean_y_Q21);
    double kT_Q21 = (ratio1_Q21 * mean_Q2_Q21 * 0.5)*(num_Q21/den_Q21);
    // calcolo gli errori con le derivate parziali
    double ratio1_Q21_err = fabs(ratio1_Q21) * sqrt(pow(A1_Q21_err / A1_Q21, 2) + pow(A0_Q21_err / A0_Q21, 2));
    double ratio2_Q21_err = fabs(ratio2_Q21) * sqrt(pow(A2_Q21_err / A2_Q21, 2) + pow(A0_Q21_err / A0_Q21, 2));
    double dkT_A_Q21 = 0.5 * mean_Q2_Q21 * (num_Q21 / den_Q21);
    double dkT_Q_Q21 = 0.5 * ratio1_Q21 * (num_Q21 / den_Q21);
    double dkT_y_Q21 = 0.5 * ratio1_Q21 * mean_Q2_Q21 * ((0.5 * mean_y_Q21 * 2 - 1) * den_Q21 - num_Q21 * (-1)) / (den_Q21 * den_Q21);
    double sigma_kT_Q21 = sqrt(pow(dkT_A_Q21 * ratio1_Q21_err, 2) + pow(dkT_Q_Q21 * mean_Q2_Q21_err, 2) + pow(dkT_y_Q21 * mean_y_Q21_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,Q^{2}_{1}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_Q21, ratio1_Q21_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,Q^{2}_{1}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_Q21, ratio2_Q21_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{Q^{2}}_{1} = %.4f #pm %.4f GeV^{2}", kT_Q21, sigma_kT_Q21));
    gPad->Update();
    cQ21->Update();
    cQ21->Write();
    // Q22
    TCanvas* cQ22 = new TCanvas("Phi_h_Kaon+_Q2_2", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cQ22->cd();
    cQ22->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncQ22 = new TF1("fitFuncQ22", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncQ22->SetParameters(0, kp_phi_Q22->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncQ22->SetParName(i, paramNames[i]);
    }
    kp_phi_Q22->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_Q22->GetYaxis()->SetTitle("Counts");
    kp_phi_Q22->Draw("E");
    kp_phi_Q22->Fit(fitFuncQ22, "RW");
    fitFuncQ22->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_Q22);
    double A0_Q22 = fitFuncQ22->GetParameter(0);
    double A1_Q22 = fitFuncQ22->GetParameter(1);
    double A2_Q22 = fitFuncQ22->GetParameter(2);
    double A0_Q22_err = fitFuncQ22->GetParError(0);
    double A1_Q22_err = fitFuncQ22->GetParError(1);
    double A2_Q22_err = fitFuncQ22->GetParError(2);
    double ratio1_Q22 = A1_Q22 / A0_Q22;
    double ratio2_Q22 = A2_Q22 / A0_Q22;
    double mean_Q2_Q22 = kp_phi_Q22_Q2->GetMean();
    double mean_y_Q22 = kp_phi_Q22_y->GetMean();
    double mean_Q2_Q22_err = kp_phi_Q22_Q2->GetMeanError();
    double mean_y_Q22_err = kp_phi_Q22_y->GetMeanError();
    double num_Q22 = (1 - mean_y_Q22 + (mean_y_Q22 * mean_y_Q22*0.5));
    double den_Q22 = (1 - mean_y_Q22);
    double kT_Q22 = (ratio1_Q22 * mean_Q2_Q22 * 0.5)*(num_Q22/den_Q22);
    // calcolo gli errori con le derivate parziali
    double ratio1_Q22_err = fabs(ratio1_Q22) * sqrt(pow(A1_Q22_err / A1_Q22, 2) + pow(A0_Q22_err / A0_Q22, 2));
    double ratio2_Q22_err = fabs(ratio2_Q22) * sqrt(pow(A2_Q22_err / A2_Q22, 2) + pow(A0_Q22_err / A0_Q22, 2));
    double dkT_A_Q22 = 0.5 * mean_Q2_Q22 * (num_Q22 / den_Q22);
    double dkT_Q_Q22 = 0.5 * ratio1_Q22 * (num_Q22 / den_Q22);
    double dkT_y_Q22 = 0.5 * ratio1_Q22 * mean_Q2_Q22 * ((0.5 * mean_y_Q22 * 2 - 1) * den_Q22 - num_Q22 * (-1)) / (den_Q22 * den_Q22);
    double sigma_kT_Q22 = sqrt(pow(dkT_A_Q22 * ratio1_Q22_err, 2) + pow(dkT_Q_Q22 * mean_Q2_Q22_err, 2) + pow(dkT_y_Q22 * mean_y_Q22_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,Q^{2}_{2}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_Q22, ratio1_Q22_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,Q^{2}_{2}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_Q22, ratio2_Q22_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{Q^{2}}_{2} = %.4f #pm %.4f GeV^{2}", kT_Q22, sigma_kT_Q22));
    gPad->Update();
    cQ22->Update();
    cQ22->Write();
    // Q23
    TCanvas* cQ23 = new TCanvas("Phi_h_Kaon+_Q2_3", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cQ23->cd();
    cQ23->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncQ23 = new TF1("fitFuncQ23", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncQ23->SetParameters(0, kp_phi_Q23->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncQ23->SetParName(i, paramNames[i]);
    }
    kp_phi_Q23->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_Q23->GetYaxis()->SetTitle("Counts");
    kp_phi_Q23->Draw("E");
    kp_phi_Q23->Fit(fitFuncQ23, "RW");
    fitFuncQ23->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_Q23);
    double A0_Q23 = fitFuncQ23->GetParameter(0);
    double A1_Q23 = fitFuncQ23->GetParameter(1);
    double A2_Q23 = fitFuncQ23->GetParameter(2);
    double A0_Q23_err = fitFuncQ23->GetParError(0);
    double A1_Q23_err = fitFuncQ23->GetParError(1);
    double A2_Q23_err = fitFuncQ23->GetParError(2);
    double ratio1_Q23 = A1_Q23 / A0_Q23;
    double ratio2_Q23 = A2_Q23 / A0_Q23;
    double mean_Q2_Q23 = kp_phi_Q23_Q2->GetMean();
    double mean_y_Q23 = kp_phi_Q23_y->GetMean();
    double mean_Q2_Q23_err = kp_phi_Q23_Q2->GetMeanError();
    double mean_y_Q23_err = kp_phi_Q23_y->GetMeanError();
    double num_Q23 = (1 - mean_y_Q23 + (mean_y_Q23 * mean_y_Q23*0.5));
    double den_Q23 = (1 - mean_y_Q23);
    double kT_Q23 = (ratio1_Q23 * mean_Q2_Q23 * 0.5)*(num_Q23/den_Q23);
    // calcolo gli errori con le derivate parziali
    double ratio1_Q23_err = fabs(ratio1_Q23) * sqrt(pow(A1_Q23_err / A1_Q23, 2) + pow(A0_Q23_err / A0_Q23, 2));
    double ratio2_Q23_err = fabs(ratio2_Q23) * sqrt(pow(A2_Q23_err / A2_Q23, 2) + pow(A0_Q23_err / A0_Q23, 2));
    double dkT_A_Q23 = 0.5 * mean_Q2_Q23 * (num_Q23 / den_Q23);
    double dkT_Q_Q23 = 0.5 * ratio1_Q23 * (num_Q23 / den_Q23);
    double dkT_y_Q23 = 0.5 * ratio1_Q23 * mean_Q2_Q23 * ((0.5 * mean_y_Q23 * 2 - 1) * den_Q23 - num_Q23 * (-1)) / (den_Q23 * den_Q23);
    double sigma_kT_Q23 = sqrt(pow(dkT_A_Q23 * ratio1_Q23_err, 2) + pow(dkT_Q_Q23 * mean_Q2_Q23_err, 2) + pow(dkT_y_Q23 * mean_y_Q23_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,Q^{2}_{3}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_Q23, ratio1_Q23_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,Q^{2}_{3}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_Q23, ratio2_Q23_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{Q^{2}}_{3} = %.4f #pm %.4f GeV^{2}", kT_Q23, sigma_kT_Q23));
    gPad->Update();
    cQ23->Update();
    cQ23->Write();
    // Q24
    TCanvas* cQ24 = new TCanvas("Phi_h_Kaon+_Q2_4", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cQ24->cd();
    cQ24->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncQ24 = new TF1("fitFuncQ24", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncQ24->SetParameters(0, kp_phi_Q24->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncQ24->SetParName(i, paramNames[i]);
    }
    kp_phi_Q24->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_Q24->GetYaxis()->SetTitle("Counts");
    kp_phi_Q24->Draw("E");
    kp_phi_Q24->Fit(fitFuncQ24, "RW");
    fitFuncQ24->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_Q24);
    double A0_Q24 = fitFuncQ24->GetParameter(0);
    double A1_Q24 = fitFuncQ24->GetParameter(1);
    double A2_Q24 = fitFuncQ24->GetParameter(2);
    double A0_Q24_err = fitFuncQ24->GetParError(0);
    double A1_Q24_err = fitFuncQ24->GetParError(1);
    double A2_Q24_err = fitFuncQ24->GetParError(2);
    double ratio1_Q24 = A1_Q24 / A0_Q24;
    double ratio2_Q24 = A2_Q24 / A0_Q24;
    double mean_Q2_Q24 = kp_phi_Q24_Q2->GetMean();
    double mean_y_Q24 = kp_phi_Q24_y->GetMean();
    double mean_Q2_Q24_err = kp_phi_Q24_Q2->GetMeanError();
    double mean_y_Q24_err = kp_phi_Q24_y->GetMeanError();
    double num_Q24 = (1 - mean_y_Q24 + (mean_y_Q24 * mean_y_Q24*0.5));
    double den_Q24 = (1 - mean_y_Q24);
    double kT_Q24 = (ratio1_Q24 * mean_Q2_Q24 * 0.5)*(num_Q24/den_Q24);
    // calcolo gli errori con le derivate parziali
    double ratio1_Q24_err = fabs(ratio1_Q24) * sqrt(pow(A1_Q24_err / A1_Q24, 2) + pow(A0_Q24_err / A0_Q24, 2));
    double ratio2_Q24_err = fabs(ratio2_Q24) * sqrt(pow(A2_Q24_err / A2_Q24, 2) + pow(A0_Q24_err / A0_Q24, 2));
    double dkT_A_Q24 = 0.5 * mean_Q2_Q24 * (num_Q24 / den_Q24);
    double dkT_Q_Q24 = 0.5 * ratio1_Q24 * (num_Q24 / den_Q24);
    double dkT_y_Q24 = 0.5 * ratio1_Q24 * mean_Q2_Q24 * ((0.5 * mean_y_Q24 * 2 - 1) * den_Q24 - num_Q24 * (-1)) / (den_Q24 * den_Q24);
    double sigma_kT_Q24 = sqrt(pow(dkT_A_Q24 * ratio1_Q24_err, 2) + pow(dkT_Q_Q24 * mean_Q2_Q24_err, 2) + pow(dkT_y_Q24 * mean_y_Q24_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,Q^{2}_{4}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_Q24, ratio1_Q24_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,Q^{2}_{4}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_Q24, ratio2_Q24_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{Q^{2}}_{4} = %.4f #pm %.4f GeV^{2}", kT_Q24, sigma_kT_Q24));
    gPad->Update();
    cQ24->Update();
    cQ24->Write();

    // ___________________________________________________________________________________________________________________________________________
    // _____________________________________________________    xB    ____________________________________________________________________________
    // xB1
    TCanvas* cxB1 = new TCanvas("Phi_h_Kaon+_xB_1", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cxB1->cd();
    cxB1->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncxB1 = new TF1("fitFuncxB1", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncxB1->SetParameters(0, kp_phi_xB1->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncxB1->SetParName(i, paramNames[i]);
    }
    kp_phi_xB1->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_xB1->GetYaxis()->SetTitle("Counts");
    kp_phi_xB1->Draw("E");
    kp_phi_xB1->Fit(fitFuncxB1, "RW");
    fitFuncxB1->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_xB1);
    double A0_xB1 = fitFuncxB1->GetParameter(0);
    double A1_xB1 = fitFuncxB1->GetParameter(1);
    double A2_xB1 = fitFuncxB1->GetParameter(2);
    double A0_xB1_err = fitFuncxB1->GetParError(0);
    double A1_xB1_err = fitFuncxB1->GetParError(1);
    double A2_xB1_err = fitFuncxB1->GetParError(2);
    double ratio1_xB1 = A1_xB1 / A0_xB1;
    double ratio2_xB1 = A2_xB1 / A0_xB1;
    double mean_Q2_xB1 = kp_phi_xB1_Q2->GetMean();
    double mean_y_xB1 = kp_phi_xB1_y->GetMean();
    double mean_Q2_xB1_err = kp_phi_xB1_Q2->GetMeanError();
    double mean_y_xB1_err = kp_phi_xB1_y->GetMeanError();
    double num_xB1 = (1 - mean_y_xB1 + (mean_y_xB1 * mean_y_xB1*0.5));
    double den_xB1 = (1 - mean_y_xB1);
    double kT_xB1 = (ratio1_xB1 * mean_Q2_xB1 * 0.5)*(num_xB1/den_xB1);
    // calcolo gli errori con le derivate parziali
    double ratio1_xB1_err = fabs(ratio1_xB1) * sqrt(pow(A1_xB1_err / A1_xB1, 2) + pow(A0_xB1_err / A0_xB1, 2));
    double ratio2_xB1_err = fabs(ratio2_xB1) * sqrt(pow(A2_xB1_err / A2_xB1, 2) + pow(A0_xB1_err / A0_xB1, 2));
    double dkT_A_xB1 = 0.5 * mean_Q2_xB1 * (num_xB1 / den_xB1);
    double dkT_Q_xB1 = 0.5 * ratio1_xB1 * (num_xB1 / den_xB1);
    double dkT_y_xB1 = 0.5 * ratio1_xB1 * mean_Q2_xB1 * ((0.5 * mean_y_xB1 * 2 - 1) * den_xB1 - num_xB1 * (-1)) / (den_xB1 * den_xB1);
    double sigma_kT_xB1 = sqrt(pow(dkT_A_xB1 * ratio1_xB1_err, 2) + pow(dkT_Q_xB1 * mean_Q2_xB1_err, 2) + pow(dkT_y_xB1 * mean_y_xB1_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,x_{B,1}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_xB1, ratio1_xB1_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,x_{B,1}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_xB1, ratio2_xB1_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{x_{B}}_{1} = %.4f #pm %.4f GeV^{2}", kT_xB1, sigma_kT_xB1));
    gPad->Update();
    cxB1->Update();
    cxB1->Write();
    // xB2
    TCanvas* cxB2 = new TCanvas("Phi_h_Kaon+_xB_2", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cxB2->cd();
    cxB2->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncxB2 = new TF1("fitFuncxB2", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncxB2->SetParameters(0, kp_phi_xB2->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncxB2->SetParName(i, paramNames[i]);
    }
    kp_phi_xB2->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_xB2->GetYaxis()->SetTitle("Counts");
    kp_phi_xB2->Draw("E");
    kp_phi_xB2->Fit(fitFuncxB2, "RW");
    fitFuncxB2->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_xB2);
    double A0_xB2 = fitFuncxB2->GetParameter(0);
    double A1_xB2 = fitFuncxB2->GetParameter(1);
    double A2_xB2 = fitFuncxB2->GetParameter(2);
    double A0_xB2_err = fitFuncxB2->GetParError(0);
    double A1_xB2_err = fitFuncxB2->GetParError(1);
    double A2_xB2_err = fitFuncxB2->GetParError(2);
    double ratio1_xB2 = A1_xB2 / A0_xB2;
    double ratio2_xB2 = A2_xB2 / A0_xB2;
    double mean_Q2_xB2 = kp_phi_xB2_Q2->GetMean();
    double mean_y_xB2 = kp_phi_xB2_y->GetMean();
    double mean_Q2_xB2_err = kp_phi_xB2_Q2->GetMeanError();
    double mean_y_xB2_err = kp_phi_xB2_y->GetMeanError();
    double num_xB2 = (1 - mean_y_xB2 + (mean_y_xB2 * mean_y_xB2*0.5));
    double den_xB2 = (1 - mean_y_xB2);
    double kT_xB2 = (ratio1_xB2 * mean_Q2_xB2 * 0.5)*(num_xB2/den_xB2);
    // calcolo gli errori con le derivate parziali
    double ratio1_xB2_err = fabs(ratio1_xB2) * sqrt(pow(A1_xB2_err / A1_xB2, 2) + pow(A0_xB2_err / A0_xB2, 2));
    double ratio2_xB2_err = fabs(ratio2_xB2) * sqrt(pow(A2_xB2_err / A2_xB2, 2) + pow(A0_xB2_err / A0_xB2, 2));
    double dkT_A_xB2 = 0.5 * mean_Q2_xB2 * (num_xB2 / den_xB2);
    double dkT_Q_xB2 = 0.5 * ratio1_xB2 * (num_xB2 / den_xB2);
    double dkT_y_xB2 = 0.5 * ratio1_xB2 * mean_Q2_xB2 * ((0.5 * mean_y_xB2 * 2 - 1) * den_xB2 - num_xB2 * (-1)) / (den_xB2 * den_xB2);
    double sigma_kT_xB2 = sqrt(pow(dkT_A_xB2 * ratio1_xB2_err, 2) + pow(dkT_Q_xB2 * mean_Q2_xB2_err, 2) + pow(dkT_y_xB2 * mean_y_xB2_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,x_{B,2}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_xB2, ratio1_xB2_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,x_{B,2}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_xB2, ratio2_xB2_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{x_{B}}_{2} = %.4f #pm %.4f GeV^{2}", kT_xB2, sigma_kT_xB2));
    gPad->Update();
    cxB2->Update();
    cxB2->Write();
    // xB3
    TCanvas* cxB3 = new TCanvas("Phi_h_Kaon+_xB_3", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cxB3->cd();
    cxB3->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncxB3 = new TF1("fitFuncxB3", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncxB3->SetParameters(0, kp_phi_xB3->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncxB3->SetParName(i, paramNames[i]);
    }
    kp_phi_xB3->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_xB3->GetYaxis()->SetTitle("Counts");
    kp_phi_xB3->Draw("E");
    kp_phi_xB3->Fit(fitFuncxB3, "RW");
    fitFuncxB3->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_xB3);
    double A0_xB3 = fitFuncxB3->GetParameter(0);
    double A1_xB3 = fitFuncxB3->GetParameter(1);
    double A2_xB3 = fitFuncxB3->GetParameter(2);
    double A0_xB3_err = fitFuncxB3->GetParError(0);
    double A1_xB3_err = fitFuncxB3->GetParError(1);
    double A2_xB3_err = fitFuncxB3->GetParError(2);
    double ratio1_xB3 = A1_xB3 / A0_xB3;
    double ratio2_xB3 = A2_xB3 / A0_xB3;
    double mean_Q2_xB3 = kp_phi_xB3_Q2->GetMean();
    double mean_y_xB3 = kp_phi_xB3_y->GetMean();
    double mean_Q2_xB3_err = kp_phi_xB3_Q2->GetMeanError();
    double mean_y_xB3_err = kp_phi_xB3_y->GetMeanError();
    double num_xB3 = (1 - mean_y_xB3 + (mean_y_xB3 * mean_y_xB3*0.5));
    double den_xB3 = (1 - mean_y_xB3);
    double kT_xB3 = (ratio1_xB3 * mean_Q2_xB3 * 0.5)*(num_xB3/den_xB3);
    // calcolo gli errori con le derivate parziali
    double ratio1_xB3_err = fabs(ratio1_xB3) * sqrt(pow(A1_xB3_err / A1_xB3, 2) + pow(A0_xB3_err / A0_xB3, 2));
    double ratio2_xB3_err = fabs(ratio2_xB3) * sqrt(pow(A2_xB3_err / A2_xB3, 2) + pow(A0_xB3_err / A0_xB3, 2));
    double dkT_A_xB3 = 0.5 * mean_Q2_xB3 * (num_xB3 / den_xB3);
    double dkT_Q_xB3 = 0.5 * ratio1_xB3 * (num_xB3 / den_xB3);
    double dkT_y_xB3 = 0.5 * ratio1_xB3 * mean_Q2_xB3 * ((0.5 * mean_y_xB3 * 2 - 1) * den_xB3 - num_xB3 * (-1)) / (den_xB3 * den_xB3);
    double sigma_kT_xB3 = sqrt(pow(dkT_A_xB3 * ratio1_xB3_err, 2) + pow(dkT_Q_xB3 * mean_Q2_xB3_err, 2) + pow(dkT_y_xB3 * mean_y_xB3_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,x_{B,3}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_xB3, ratio1_xB3_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,x_{B,3}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_xB3, ratio2_xB3_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{x_{B}}_{3} = %.4f #pm %.4f GeV^{2}", kT_xB3, sigma_kT_xB3));
    gPad->Update();
    cxB3->Update();
    cxB3->Write();
    // xB4
    TCanvas* cxB4 = new TCanvas("Phi_h_Kaon+_xB_4", "Cahn and Boer-Mulders analysis | Kaon+ | z < 0.4 ", 800, 600);
    cxB4->cd();
    cxB4->ToggleEditor();
    // Creazione del fit
    TF1* fitFuncxB4 = new TF1("fitFuncxB4", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x) + [5]*cos(5*x) + [6]*sin(x) + [7]*sin(2*x) + [8]*sin(3*x) + [9]*sin(4*x) + [10]*sin(5*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFuncxB4->SetParameters(0, kp_phi_xB4->GetMaximum());
    //const char* paramNames[] = {"A0", "A_Cahn", "A_BM", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"};
    for (int i = 0; i < 11; i++) {
        fitFuncxB4->SetParName(i, paramNames[i]);
    }
    kp_phi_xB4->GetXaxis()->SetTitle("#Phi_h [Rad]");
    kp_phi_xB4->GetYaxis()->SetTitle("Counts");
    kp_phi_xB4->Draw("E");
    kp_phi_xB4->Fit(fitFuncxB4, "RW");
    fitFuncxB4->Draw("SAME");
    gStyle->SetOptFit(111);
    SetStatsBox1(kp_phi_xB4);
    double A0_xB4 = fitFuncxB4->GetParameter(0);
    double A1_xB4 = fitFuncxB4->GetParameter(1);
    double A2_xB4 = fitFuncxB4->GetParameter(2);
    double A0_xB4_err = fitFuncxB4->GetParError(0);
    double A1_xB4_err = fitFuncxB4->GetParError(1);
    double A2_xB4_err = fitFuncxB4->GetParError(2);
    double ratio1_xB4 = A1_xB4 / A0_xB4;
    double ratio2_xB4 = A2_xB4 / A0_xB4;
    double mean_Q2_xB4 = kp_phi_xB4_Q2->GetMean();
    double mean_y_xB4 = kp_phi_xB4_y->GetMean();
    double mean_Q2_xB4_err = kp_phi_xB4_Q2->GetMeanError();
    double mean_y_xB4_err = kp_phi_xB4_y->GetMeanError();
    double num_xB4 = (1 - mean_y_xB4 + (mean_y_xB4 * mean_y_xB4*0.5));
    double den_xB4 = (1 - mean_y_xB4);
    double kT_xB4 = (ratio1_xB4 * mean_Q2_xB4 * 0.5)*(num_xB4/den_xB4);
    // calcolo gli errori con le derivate parziali
    double ratio1_xB4_err = fabs(ratio1_xB4) * sqrt(pow(A1_xB4_err / A1_xB4, 2) + pow(A0_xB4_err / A0_xB4, 2));
    double ratio2_xB4_err = fabs(ratio2_xB4) * sqrt(pow(A2_xB4_err / A2_xB4, 2) + pow(A0_xB4_err / A0_xB4, 2));
    double dkT_A_xB4 = 0.5 * mean_Q2_xB4 * (num_xB4 / den_xB4);
    double dkT_Q_xB4 = 0.5 * ratio1_xB4 * (num_xB4 / den_xB4);
    double dkT_y_xB4 = 0.5 * ratio1_xB4 * mean_Q2_xB4 * ((0.5 * mean_y_xB4 * 2 - 1) * den_xB4 - num_xB4 * (-1)) / (den_xB4 * den_xB4);
    double sigma_kT_xB4 = sqrt(pow(dkT_A_xB4 * ratio1_xB4_err, 2) + pow(dkT_Q_xB4 * mean_Q2_xB4_err, 2) + pow(dkT_y_xB4 * mean_y_xB4_err, 2));
    latex.DrawLatex(0.15, 0.85, Form("A_{UU,x_{B,4}}^{cos(#Phi_{h})} = %.4f #pm %.4f", ratio1_xB4, ratio1_xB4_err));
    latex.DrawLatex(0.15, 0.79, Form("A_{UU,x_{B,4}}^{cos(2#Phi_{h})} = %.4f #pm %.4f", ratio2_xB4, ratio2_xB4_err));
    latex.DrawLatex(0.40, 0.85, Form("#LT k_{T}^{2} #GT^{x_{B}}_{4} = %.4f #pm %.4f GeV^{2}", kT_xB4, sigma_kT_xB4));
    gPad->Update();
    cxB4->Update();
    cxB4->Write();

    // _________________________________________________________________________________________________________________________________  

    // PLOT DI <k^2_T>
    double z_range[4] = {0.3, 0.45, 0.575, 0.75};
    double PhT_range[4] = {0.12, 0.375, 0.6, 1};
    double Q2_range[4] = {1.5, 2.5, 3.5, 6};
    double xB_range[4] = {0.11, 0.2, 0.25, 0.5};
    double kT_z_value[4] = {kT_z1, kT_z2, kT_z3, kT_z4};
    double kT_z_value_err[4] = {sigma_kT_z1, sigma_kT_z2, sigma_kT_z3, sigma_kT_z4};
    double kT_PhT_value[4] = {kT_PhT1, kT_PhT2, kT_PhT3, kT_PhT4};
    double kT_PhT_value_err[4] = {sigma_kT_PhT1, sigma_kT_PhT2, sigma_kT_PhT3, sigma_kT_PhT4};
    double kT_Q2_value[4] = {kT_Q21, kT_Q22, kT_Q23, kT_Q24};
    double kT_Q2_value_err[4] = {sigma_kT_Q21, sigma_kT_Q22, sigma_kT_Q23, sigma_kT_Q24};
    double kT_xB_value[4] = {kT_xB1, kT_xB2, kT_xB3, kT_xB4};
    double kT_xB_value_err[4] = {sigma_kT_xB1, sigma_kT_xB2, sigma_kT_xB3, sigma_kT_xB4};
    TGraphErrors* graph_z = new TGraphErrors(4, z_range, kT_z_value, 0, kT_z_value_err);
    TGraphErrors* graph_PhT = new TGraphErrors(4, PhT_range, kT_PhT_value, 0, kT_PhT_value_err);
    TGraphErrors* graph_Q2 = new TGraphErrors(4, Q2_range, kT_Q2_value, 0, kT_Q2_value_err);
    TGraphErrors* graph_xB = new TGraphErrors(4, xB_range, kT_xB_value, 0, kT_xB_value_err);
    // z
    TCanvas* d_z = new TCanvas("K2TvsZ", "#LT k_T^2 #GT vs z", 800, 600);
    graph_z->SetMarkerStyle(20);
    graph_z->SetMarkerSize(1);
    graph_z->SetMarkerColor(kGreen+3);
    graph_z->SetLineColor(kGreen+3);
    graph_z->SetTitle("#LT k_{T}^{2} #GT vs z; z; #LT k_{T}^{2} #GT [GeV^{2}]");
    graph_z->Draw("APL");
    d_z->Update();
    d_z->Write();
    // PhT
    TCanvas* d_PhT = new TCanvas("K2TvsPhT", "#LT k_T^2 #GT vs P_{hT}", 800, 600);
    graph_PhT->SetMarkerStyle(20);
    graph_PhT->SetMarkerSize(1);
    graph_PhT->SetMarkerColor(kRed);
    graph_PhT->SetLineColor(kRed);
    graph_PhT->SetTitle("#LT k_{T}^{2} #GT vs P_{hT}; P_{hT} [GeV]; #LT k_{T}^{2} #GT [GeV^{2}]");
    graph_PhT->Draw("APL");
    d_PhT->Update();
    d_PhT->Write();
    // Q2
    TCanvas* d_Q2 = new TCanvas("K2TvsQ2", "#LT k_T^2 #GT vs Q^{2}", 800, 600);
    graph_Q2->SetMarkerStyle(20);
    graph_Q2->SetMarkerSize(1);
    graph_Q2->SetMarkerColor(kOrange+2);
    graph_Q2->SetLineColor(kOrange+2);
    graph_Q2->SetTitle("#LT k_{T}^{2} #GT vs Q^{2}; Q^{2} [GeV^{2}]; #LT k_{T}^{2} #GT [GeV^{2}]");
    graph_Q2->Draw("APL");
    d_Q2->Update();
    d_Q2->Write();
    // xB
    TCanvas* d_xB = new TCanvas("K2TvsxB", "#LT k_T^2 #GT vs x_{B}", 800, 600);
    graph_xB->SetMarkerStyle(20);
    graph_xB->SetMarkerSize(1);
    graph_xB->SetMarkerColor(kBlue-2);
    graph_xB->SetLineColor(kBlue-2);
    graph_xB->SetTitle("#LT k_{T}^{2} #GT vs x_{B}; x_{B}; #LT k_{T}^{2} #GT [GeV^{2}]");
    graph_xB->Draw("APL");
    d_xB->Update();
    d_xB->Write();

    // _________________________________________________________________________________________________________________________________  

    outputFile->Write();
    outputFile->Close();
    inputFile->Close();

    std::cout << "Analisi completata. Risultati salvati in " << outputFileName << std::endl;
}

int main() {
    const char* inputFileName = "output(+1).root";
    const char* outputFileName = "Fit_output(+1).root";
    FourierAnalysis(inputFileName, outputFileName);
    return 0;
}
