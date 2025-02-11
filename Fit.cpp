#include <iostream>
#include <cmath>
#include <TStyle.h>
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

void SetStatsBox1(TH1D* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.65);  // Posizione pannello (sinistra)
        stats0->SetX2NDC(0.89);  // Posizione pannello (destra)
        stats0->SetY1NDC(0.16);  // Posizione pannello (basso)
        stats0->SetY2NDC(0.42);  // Posizione pannello (alto)
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

    TH1D *h_phi_Kp = new TH1D("h_phi_Kp", "Phi_h distribution | Kaon+", 60, -TMath::Pi(), TMath::Pi());
    treeKp->Draw("Phi_h >> h_phi_Kp");
    TH1D *h_phi_Km = new TH1D("h_phi_Km", "Phi_h distribution | Kaon-", 60, -TMath::Pi(), TMath::Pi());
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
    if (stats3) {
        stats3->SetX1NDC(0.70);     // Sposta il pannello
        stats3->SetX2NDC(0.89);     
        stats3->SetY1NDC(0.16);     
        stats3->SetY2NDC(0.40);     
    }

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
    gPad->Update();
    cxB4->Update();
    cxB4->Write();

    // _________________________________________________________________________________________________________________________________  

    outputFile->Write();
    outputFile->Close();
    inputFile->Close();

    std::cout << "Analisi completata. Risultati salvati in " << outputFileName << std::endl;
}

int main() {
    const char* inputFileName = "output.root";
    const char* outputFileName = "Fit_output.root";
    FourierAnalysis(inputFileName, outputFileName);
    return 0;
}
