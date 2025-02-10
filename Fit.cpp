#include <iostream>
#include <cmath>
#include <TStyle.h>
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

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
    //tree->Print();
    TH1D *h_phi_Kp = new TH1D("h_phi_Kp", "Phi_h distribution | Kaon+", 60, -TMath::Pi(), TMath::Pi());
    treeKp->Draw("Phi_h >> h_phi_Kp");
    TH1D *h_phi_Km = new TH1D("h_phi_Km", "Phi_h distribution | Kaon-", 60, -TMath::Pi(), TMath::Pi());
    treeKm->Draw("Phi_h >> h_phi_Km");

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

    // Creazione del file ROOT di output
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Errore: impossibile creare il file " << outputFileName << std::endl;
        inputFile->Close();
        return;
    }

    TPaveStats* stats = (TPaveStats*)h_phi_Kp->FindObject("stats");
    if (stats) {
        stats->Delete(); // Rimuove il box delle statistiche precedente
    }
    TPaveStats* stats2 = (TPaveStats*)h_phi_Km->FindObject("stats");
    if (stats2) {
        stats2->Delete(); // Rimuove il box delle statistiche precedente
    }

    // Creazione del canvas
    TCanvas* c = new TCanvas("Phi_h_Kaon+", "Cahn and Boer-Mulders analysis", 800, 600);
    c->cd();
    c->ToggleEditor();
    // Creazione del fit
    TF1* fitFunc = new TF1("fitFunc", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x)", h_phi_Kp->GetXaxis()->GetXmin(), h_phi_Kp->GetXaxis()->GetXmax());
    fitFunc->SetParameters(h_phi_Kp->GetMaximum(), -0.01, -0.003, 0.001);
    fitFunc->SetParName(0, "A0");
    fitFunc->SetParName(1, "A_Cahn");
    fitFunc->SetParName(2, "A_BM");
    fitFunc->SetParName(3, "A3");
    fitFunc->SetLineColor(kRed);
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
    TF1* fitFunc2 = new TF1("fitFunc2", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x)", h_phi_Km->GetXaxis()->GetXmin(), h_phi_Km->GetXaxis()->GetXmax());
    fitFunc2->SetParameters(h_phi_Km->GetMaximum(), -0.01, -0.003, 0.001);
    fitFunc2->SetParName(0, "A0");
    fitFunc2->SetParName(1, "A_Cahn");
    fitFunc2->SetParName(2, "A_BM");
    fitFunc2->SetParName(3, "A3");
    fitFunc2->SetLineColor(kRed);
    h_phi_Km->GetXaxis()->SetTitle("#Phi_h [Rad]");
    h_phi_Km->GetYaxis()->SetTitle("Counts");
    h_phi_Km->Draw("E");
    h_phi_Km->Fit(fitFunc2, "RW");
    fitFunc2->Draw("SAME");
    gStyle->SetOptFit(111);
    TPaveStats* stats6 = (TPaveStats*)h_phi_Km->GetListOfFunctions()->FindObject("stats");
    if (stats6) {
        stats6->SetX1NDC(0.70);     // Sposta il pannello
        stats6->SetX2NDC(0.89);     
        stats6->SetY1NDC(0.16);     
        stats6->SetY2NDC(0.40);     
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
