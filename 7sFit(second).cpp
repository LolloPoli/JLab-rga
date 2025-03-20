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

gROOT->SetBatch(kTRUE);
void SetStatsBox1(TH1D* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.80);  // Posizione pannello (sinistra)
        stats0->SetX2NDC(0.98);  // Posizione pannello (destra)
        stats0->SetY1NDC(0.58);  // Posizione pannello (basso)
        stats0->SetY2NDC(0.92);  // Posizione pannello (alto)
    }
}

void SetStatsBox2(TH1D* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.28);  // Posizione pannello (sinistra)
        stats0->SetX2NDC(0.45);  // Posizione pannello (destra)
        stats0->SetY1NDC(0.68);  // Posizione pannello (basso)
        stats0->SetY2NDC(0.88);  // Posizione pannello (alto)
    }
}

void _7sFit(double torus) {
    // Apertura del file ROOT di input
    const char* inputFileName;
    const char* outputFileName;
    if (torus == -1){
        inputFileName = "plot_binned_rga_t-1.root";
        outputFileName = "plot_fit_rga_t-1.root";
    }
    if (torus == +1){
        inputFileName = "plot_binned_rga_t+1.root";
        outputFileName = "plot_fit_rga_t+1.root";
    }
    TFile* inputFile = TFile::Open(inputFileName);
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Errore: impossibile creare il file " << outputFileName << std::endl;
        inputFile->Close();
        return;
    }
    //if (inputFile) inputFile->ls();
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Errore: impossibile aprire il file " << inputFileName << std::endl;
        return;
    }
    
    const int nBin_xQ = 3;
    TDirectory *Kp_hp[3];
    TDirectory *Kp_hm[3];
    TDirectory *KaonBSA[3];
    TCanvas *cAlu[3];
    for(int i = 0; i < nBin_xQ; i++){
        Kp_hp[i] = (TDirectory*)inputFile->Get(Form("KaonPlus_H+_xQ%d", i+1));
        Kp_hm[i] = (TDirectory*)inputFile->Get(Form("KaonPlus_H-_xQ%d", i+1));
        cAlu[i] = (TCanvas*)inputFile->Get(Form("A_LU_z_xQ_bin%d", i+1));
        KaonBSA[i] = outputFile->mkdir(Form("BSA_z_xQ%d", i+1));
    }
    TDirectory* ALU_Dir = outputFile->mkdir("A_LU_z_Binned");

    // Definizione dei bin
    const int nBinz_z = 6;
    const int nBin_Pht = 3;
    double zBins[nBinz_z + 1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9};
    double PtBins[nBin_Pht + 1] = {0.0, 0.25, 0.5, 1.2};
    // VECTOR FOR THE A_LU RESULTS
    vector<vector<vector<double>>> A_LU_values(nBin_xQ,vector<vector<double>>(nBin_Pht, vector<double>(nBinz_z, 0)));
    vector<vector<vector<double>>> A_LU_errors(nBin_xQ,vector<vector<double>>(nBin_Pht, vector<double>(nBinz_z, 0)));
    vector<vector<vector<double>>> A_UU_1_values(nBin_xQ,vector<vector<double>>(nBin_Pht, vector<double>(nBinz_z, 0)));
    vector<vector<vector<double>>> A_UU_1_errors(nBin_xQ,vector<vector<double>>(nBin_Pht, vector<double>(nBinz_z, 0)));
    vector<vector<vector<double>>> A_UU_2_values(nBin_xQ,vector<vector<double>>(nBin_Pht, vector<double>(nBinz_z, 0)));
    vector<vector<vector<double>>> A_UU_2_errors(nBin_xQ,vector<vector<double>>(nBin_Pht, vector<double>(nBinz_z, 0)));

    TH1D* kp_phi_hp[nBin_xQ][nBinz_z][nBin_Pht];
    TH1D* kp_phi_hm[nBin_xQ][nBinz_z][nBin_Pht];

    for (int region = 0; region < nBin_xQ; region++) {
        for (int zbin = 0; zbin < nBinz_z; zbin++) {
            for (int phtbin = 0; phtbin < nBin_Pht; phtbin++) {
                kp_phi_hp[region][zbin][phtbin] = (TH1D*)Kp_hp[region]->Get(Form("_Phi_zBin%d_PhtBin%d_xQ%d_Hp", zbin+1, phtbin+1, region+1));
                kp_phi_hm[region][zbin][phtbin] = (TH1D*)Kp_hm[region]->Get(Form("_Phi_zBin%d_PhtBin%d_xQ%d_Hm", zbin+1, phtbin+1, region+1));
                if (!kp_phi_hp[region][zbin][phtbin]) {
                    std::cerr << "ERROR: Histogram _Phi_zBin" << zbin+1 << "_PhtBin" << phtbin+1 << "_xQ" << region+1 << "_Hp NOT FOUND!" << std::endl;
                }
                if (!kp_phi_hm[region][zbin][phtbin]) {
                    std::cerr << "ERROR: Histogram _Phi_zBin" << zbin+1 << "_PhtBin" << phtbin+1 << "_xQ" << region+1 << "_Hm NOT FOUND!" << std::endl;
                }
            }
        }
    }

    // NORM
    for (int region = 0; region < nBin_xQ; region++) {
        for (int zbin = 0; zbin < nBinz_z; zbin++) {
            for (int phtbin = 0; phtbin < nBin_Pht; phtbin++) {
                if (kp_phi_hp[region][zbin][phtbin] && kp_phi_hp[region][zbin][phtbin]->Integral() > 0)
                kp_phi_hp[region][zbin][phtbin]->Scale(1.0 / kp_phi_hp[region][zbin][phtbin]->Integral());
                if (kp_phi_hm[region][zbin][phtbin] && kp_phi_hm[region][zbin][phtbin]->Integral() > 0)
                kp_phi_hm[region][zbin][phtbin]->Scale(1.0 / kp_phi_hm[region][zbin][phtbin]->Integral());
            }
        }
    }

    // Beam spin asymmetry plot calculation
    TH1D* kp_BSA[nBin_xQ][nBinz_z][nBin_Pht];
    for (int region = 0; region < nBin_xQ; region++) {
        for (int zbin = 0; zbin < nBinz_z; zbin++) {
            for (int phtbin = 0; phtbin < nBin_Pht; phtbin++) {
                if (!kp_phi_hp[region][zbin][phtbin] || !kp_phi_hm[region][zbin][phtbin]) {
                    std::cerr << "ERROR: Skipping BSA calculation for zBin " << zbin+1 << ", PhtBin " << phtbin+1 << ", xQ " << region+1 << " due to missing histograms." << std::endl;
                    continue;
                }
                kp_BSA[region][zbin][phtbin] = (TH1D*)kp_phi_hp[region][zbin][phtbin]->Clone(Form("kp_BSA_z%d_p%d_xQ%d", zbin+1, phtbin+1, region+1));
                kp_BSA[region][zbin][phtbin]->Add(kp_phi_hm[region][zbin][phtbin], -1);
                TH1D* kp_hSum = (TH1D*)kp_phi_hp[region][zbin][phtbin]->Clone();
                kp_hSum->Add(kp_phi_hm[region][zbin][phtbin]);
                kp_BSA[region][zbin][phtbin]->Divide(kp_hSum);
                delete kp_hSum; // Libera memoria
            }
        }
    }
    
    //  CANVAS FOR A_LU
    TCanvas* cBSA[nBin_xQ][nBinz_z][nBin_Pht];
    TF1* fitFuncBSA[nBin_xQ][nBinz_z][nBin_Pht];
    for (int region = 0; region < nBin_xQ; region++) {
        KaonBSA[region]->cd();
        for (int zbin = 0; zbin < nBinz_z; zbin++) {
            for (int phtbin = 0; phtbin < nBin_Pht; phtbin++) {
                // to see if kp_BSA is correctly opened or not
                if (!kp_BSA[region][zbin][phtbin]) {
                    std::cerr << "ERROR: kp_BSA[" << region << "][" << zbin << "][" << phtbin << "] is NULL! Skipping fit." << std::endl;
                    continue; // Salta questo bin se non esiste
                }
                cBSA[region][zbin][phtbin] = new TCanvas(Form("bsa_kp_z%d_p%d_xQ%d", zbin+1, phtbin+1, region+1), "Beam Spin Asymmetry", 800, 600);
                // fitFuncBSA[region][zbin][phtbin] = new TF1(Form("fitFuncBSA_z%d_p%d_xQ%d", zbin+1, phtbin+1, region+1), "[0]*sin(x)", -TMath::Pi(), TMath::Pi());
                fitFuncBSA[region][zbin][phtbin] = new TF1(Form("fitFuncBSA_z%d_p%d_xQ%d", zbin+1, phtbin+1, region+1), "([0]*sin(x)) / ((1 + [1]*cos(x) + [2]*cos(2*x)) + 0.000001)", -TMath::Pi(), TMath::Pi());
                fitFuncBSA[region][zbin][phtbin]->SetParName(0, "A_{LU}");
                fitFuncBSA[region][zbin][phtbin]->SetParName(1, "A_{UU}^{cos#Phi}");
                fitFuncBSA[region][zbin][phtbin]->SetParName(2, "A_{UU}^{cos2#Phi}");
                kp_BSA[region][zbin][phtbin]->SetMarkerStyle(20);
                kp_BSA[region][zbin][phtbin]->SetMarkerSize(1);
                kp_BSA[region][zbin][phtbin]->SetMarkerColor(kBlack);
                kp_BSA[region][zbin][phtbin]->SetLineColor(kBlack);
                kp_BSA[region][zbin][phtbin]->SetTitle(Form("Beam Spin Asymmetry A_{LU} | xQ %d | z Bin %d | Pht Bin %d", region+1, zbin+1, phtbin+1));
                kp_BSA[region][zbin][phtbin]->Draw("E");
                kp_BSA[region][zbin][phtbin]->Fit(fitFuncBSA[region][zbin][phtbin], "RW");
                fitFuncBSA[region][zbin][phtbin]->Draw("SAME");
                gStyle->SetOptFit(1111);
                SetStatsBox2(kp_BSA[region][zbin][phtbin]);
                A_LU_values[region][phtbin][zbin] = fitFuncBSA[region][zbin][phtbin]->GetParameter(0);
                A_LU_errors[region][phtbin][zbin] = fitFuncBSA[region][zbin][phtbin]->GetParError(0);
                A_UU_1_values[region][phtbin][zbin] = fitFuncBSA[region][zbin][phtbin]->GetParameter(1);
                A_UU_1_errors[region][phtbin][zbin] = fitFuncBSA[region][zbin][phtbin]->GetParError(1);
                A_UU_2_values[region][phtbin][zbin] = fitFuncBSA[region][zbin][phtbin]->GetParameter(2);
                A_UU_2_errors[region][phtbin][zbin] = fitFuncBSA[region][zbin][phtbin]->GetParError(2);
                cBSA[region][zbin][phtbin]->Update();
                cBSA[region][zbin][phtbin]->Write();
            }
        }
    }
    
    ALU_Dir->cd();
    TCanvas* cA_LU[nBin_xQ]; 
    TGraphErrors* graphs[nBin_xQ][nBin_Pht];
    TGraphErrors* graphs_Auu1[nBin_xQ][nBin_Pht];
    TGraphErrors* graphs_Auu2[nBin_xQ][nBin_Pht];
    for(int region = 0; region < nBin_xQ; region++){
        cA_LU[region]= new TCanvas(Form("A_LU_z_Binned_xQ%d", region+1), Form("A_LU vs z | xQ bin %d", region+1), 800, 800);
        cA_LU[region]->Divide(1, nBin_Pht);
        for (int i = 0; i < nBin_Pht; i++) {
            graphs[region][i] = new TGraphErrors(nBinz_z);
            graphs_Auu1[region][i] = new TGraphErrors(nBinz_z);
            graphs_Auu2[region][i] = new TGraphErrors(nBinz_z);
            for (int j = 0; j < nBinz_z; j++) {
                graphs[region][i]->SetPoint(j, (zBins[j] + zBins[j+1]) / 2.0, A_LU_values[region][i][j]);
                graphs[region][i]->SetPointError(j, 0, A_LU_errors[region][i][j]);
                graphs_Auu1[region][i]->SetPoint(j, (zBins[j] + zBins[j+1]) / 2.0, A_UU_1_values[region][i][j]);
                graphs_Auu1[region][i]->SetPointError(j, 0, A_UU_1_errors[region][i][j]);
                graphs_Auu2[region][i]->SetPoint(j, (zBins[j] + zBins[j+1]) / 2.0, A_UU_2_values[region][i][j]);
                graphs_Auu2[region][i]->SetPointError(j, 0, A_UU_2_errors[region][i][j]);
            }
            cA_LU[region]->cd(i+1);  // Select the subpad in the canvas
            graphs[region][i]->SetTitle(Form("A_LU vs z (Pht Bin %d) | xQ %d", i+1, region+1));
            graphs_Auu1[region][i]->SetTitle(Form("A_UU vs z (Pht Bin %d) | xQ %d", i+1, region+1));
            graphs[region][i]->GetXaxis()->SetTitle("z");
            graphs_Auu1[region][i]->GetXaxis()->SetTitle("z");
            graphs[region][i]->GetYaxis()->SetTitle("A_LU");
            graphs_Auu1[region][i]->GetYaxis()->SetTitle("A_UU");
            graphs[region][i]->SetMarkerStyle(20);
            graphs[region][i]->SetLineColor(kBlue-2);  // Different colors per xQ region
            graphs[region][i]->SetMarkerColor(kBlue-2);
            graphs_Auu1[region][i]->SetLineColor(kGreen+3); 
            graphs_Auu1[region][i]->SetMarkerColor(kGreen+3);
            graphs_Auu2[region][i]->SetLineColor(kMagenta+3); 
            graphs_Auu2[region][i]->SetMarkerColor(kMagenta+3);
            graphs_Auu1[region][i]->SetMarkerStyle(20);
            graphs_Auu2[region][i]->SetMarkerStyle(20);
            if(torus==-1) graphs[region][i]->GetYaxis()->SetRangeUser(0, 0.06);
            if(torus==+1) graphs[region][i]->GetYaxis()->SetRangeUser(0, 0.1);
            graphs_Auu1[region][i]->GetYaxis()->SetRangeUser(-1, 2);
            graphs[region][i]->GetYaxis()->SetMaxDigits(3);
            graphs[region][i]->Draw("AP");

            gPad->SetGridx();  
            gPad->SetGridy();
        }

        cA_LU[region]->Write();
        cAlu[region]->Write();
    }

    TCanvas* c_3[3];
    for (int region = 0; region < 3; region++) {
        c_3[region] = new TCanvas(Form("A_UU_z_xQ_bin%d_NewFit", region+1), Form("Unbinned A_{UU} vs z for xQ bin %d", region+1), 800, 800);
        c_3[region]->Divide(1, nBin_Pht);
        for (int i = 0; i < nBin_Pht; i++) {
            c_3[region]->cd(i+1);
            graphs_Auu1[region][i]->Draw("AP");
            graphs_Auu2[region][i]->Draw("SAME P");
            //SetStatsBoxGraph(graphs[region][i], nBinz_z);
            TLegend *leg = new TLegend(0.85, 0.65, 0.95, 0.85);
            leg->AddEntry(graphs_Auu1[region][i], "A_{UU} cos#Phi", "ep");
            leg->AddEntry(graphs_Auu2[region][i], "A_{UU} cos2#Phi", "epf");
            leg->Draw();
            gPad->SetGridx();
            gPad->SetGridy();
            gPad->Update();
        }
        c_3[region]->Update();
        c_3[region]->Write();
    }
    
    TCanvas* cCombined[nBin_xQ];

    for (int region = 0; region < nBin_xQ; region++) {
        // Creiamo la nuova canvas combinata
        cCombined[region] = new TCanvas(Form("A_LU_z_Combined_xQ%d", region+1), 
                                        Form("A_LU vs z Combined | xQ bin %d", region+1), 
                                        800, 800);
        cCombined[region]->Divide(1, nBin_Pht);

        for (int i = 0; i < nBin_Pht; i++) {
            cCombined[region]->cd(i+1);  // Selezioniamo il subpad

            // Recuperiamo i grafici dalla prima canvas
            cA_LU[region]->cd(i+1);
            TGraphErrors* graph1 = (TGraphErrors*)gPad->FindObject("Graph");
            
            // Recuperiamo i grafici dalla seconda canvas
            cAlu[region]->cd(i+1);
            TGraphErrors* graph2 = (TGraphErrors*)gPad->FindObject("Graph");

            // Controlliamo se i grafici sono stati trovati
            if (!graph1) {
                std::cout << "⚠️ Warning: Graph1 non trovato in region " << region << ", bin " << i << std::endl;
            }
            if (!graph2) {
                std::cout << "⚠️ Warning: Graph2 non trovato in region " << region << ", bin " << i << std::endl;
            }

            // Disegniamo il primo grafico
            if (graph1) {
                graph1->SetMarkerColor(kRed);
                graph1->SetLineColor(kRed);
                graph1->Draw("AP");
            }

            // Sovrapponiamo il secondo grafico
            if (graph2) {
                graph2->SetMarkerColor(kBlue);
                graph2->SetLineColor(kBlue);
                graph2->Draw("P SAME");  // "SAME" per sovrapporre
            }

            gPad->SetGridx();
            gPad->SetGridy();
            gPad->Modified();  // Forza il ridisegno
            gPad->Update();    // Aggiorna il pad
        }

        cCombined[region]->Update();
        cCombined[region]->Write();  // Salviamo la nuova canvas
    }


    
    //outputFile->Write();
    outputFile->Close();
    inputFile->Close();

    std::cout << "Analyses completed. Results saved on: " << outputFileName << std::endl;

}





