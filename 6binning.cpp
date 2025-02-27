#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

using namespace std;

// Funzione di verosimiglianza negativa da minimizzare
double extractFLU(const double *FLU, const std::vector<double>& phi, const std::vector<double>& bpol, const std::vector<double>& epsilon, const std::vector<int>& helicity) {
    double logLike = 0.0;
    for (size_t i = 0; i < phi.size(); i++) {
        double term = 1 + (helicity[i] * bpol[i]) * FLU[0] * sin(phi[i]);
        // Evitiamo log(0) o log(valori negativi)
        if (term <= 0) term = 1e-10;  
        logLike += -2 * TMath::Log(term);
    }
    return logLike;
}

int main() {
    TFile inFile("out_t-1.root", "READ");
    const char* outputFile = "plot_binned_t-1.root";
    TFile outFile(outputFile, "RECREATE");  // File di output ROOT
    TTree *tree = (TTree*)inFile.Get("Tree");
    TTree *treeEl = (TTree*)inFile.Get("Electron");
    TTree *treeKaonP = (TTree*)inFile.Get("Kaon+");
    TTree *treeKaonPhp = (TTree*)inFile.Get("Kaon+ H+");
    TTree *treeKaonPhm = (TTree*)inFile.Get("Kaon+ H-");
    TTree *treeKaonM = (TTree*)inFile.Get("Kaon-");
    TDirectory* dirKaonp_Hp_xQ[3];
    TDirectory* dirKaonp_Hm_xQ[3];
    for (int region = 0; region < 3; region++) {
        dirKaonp_Hp_xQ[region] = outFile.mkdir(Form("KaonPlus_H+_xQ%d", region+1));
        dirKaonp_Hm_xQ[region] = outFile.mkdir(Form("KaonPlus_H-_xQ%d", region+1));
    }

    double helicity;
    double beta_CMS_Double;
    int _torus;
    TVector3 beta_CMS;
    int electron_Nphe, electron_status, electron_sector;
    double N_up, N_down;
    double electron_PCAL, electron_ECAL, electron_ECIN, electron_CAL_Tot, electron_vz;
    double electron_edge1, electron_edge2, electron_edge3;
    // elettrone
    double electron_px, electron_py, electron_pz, electron_mom, electron_Theta, electron_ThetaDeg, electron_E, electron_W, electron_Q2;
    // gamma
    double gamma_px, gamma_py, gamma_pz, gamma_nu;
    // pion + 
    double pionp_px, pionp_py, pionp_pz;
    // kaone +
    double kaonp_px, kaonp_py, kaonp_pz;
    double kaonp_xF, kaonp_xB, kaonp_Q2, kaonp_z, kaonp_Ph, kaonp_Pt;
    double kaonp_Phi_h, kaonp_Phi, kaonp_Theta, kaonp_y, kaonp_s, kaonp_E, kaonp_W, kaonp_Phi_hDeg;
    double kaonp_Ph_x, kaonp_Ph_y, kaonp_PhT, kaonp_eta, kaonp_M2x, kaonp_chi2pid;
    double kaonp_edge1, kaonp_edge2, kaonp_edge3, kaonp_vz, kaonp_SinPhi;
    double kaonp_rich_Id, kaonp_rich_PID = 0;
    double kaonp_phi_cambio, kaonp_phi_cambio2;
    double kaonp_helicity, kaonp_Phi_Hup, kaonp_Phi_Hdw, kaonp_Pol, kaonp_epsilon;
    std::vector<double> kaonp_Phi_vec, kaonp_Pol_vec, kaonp_eps_vec;
    std::vector<int> kaonp_hel_vec;

    // kaone + Helicity +
    double kaonp_px_hp, kaonp_py_hp, kaonp_pz_hp;
    double kaonp_xF_hp, kaonp_xB_hp, kaonp_Q2_hp, kaonp_z_hp, kaonp_Ph_hp, kaonp_Pt_hp;
    double kaonp_Phi_h_hp, kaonp_Phi_hp, kaonp_Theta_hp, kaonp_y_hp, kaonp_s_hp, kaonp_E_hp, kaonp_W_hp, kaonp_Phi_hDeg_hp;
    double kaonp_Ph_x_hp, kaonp_Ph_y_hp, kaonp_PhT_hp, kaonp_eta_hp, kaonp_M2x_hp, kaonp_chi2pid_hp;
    double kaonp_helicity_hp;
    // kaone + Helicity -
    double kaonp_px_hm, kaonp_py_hm, kaonp_pz_hm;
    double kaonp_xF_hm, kaonp_xB_hm, kaonp_Q2_hm, kaonp_z_hm, kaonp_Ph_hm, kaonp_Pt_hm;
    double kaonp_Phi_h_hm, kaonp_Phi_hm, kaonp_Theta_hm, kaonp_y_hm, kaonp_s_hm, kaonp_E_hm, kaonp_W_hm, kaonp_Phi_hDeg_hm;
    double kaonp_Ph_x_hm, kaonp_Ph_y_hm, kaonp_PhT_hm, kaonp_eta_hm, kaonp_M2x_hm, kaonp_chi2pid_hm;
    double kaonp_helicity_hm;

    // kaone -
    double kaonm_px, kaonm_py, kaonm_pz;
    double kaonm_xF, kaonm_xB, kaonm_Q2, kaonm_z, kaonm_Ph, kaonm_Pt;
    double kaonm_Phi_h, kaonm_Phi, kaonm_Theta, kaonm_y, kaonm_s, kaonm_E, kaonm_W;
    double kaonm_Ph_x, kaonm_Ph_y, kaonm_PhT, kaonm_eta, kaonm_M2x, kaonm_chi2pid;
    double kaonm_edge1, kaonm_edge2, kaonm_edge3, kaonm_vz, kaonm_SinPhi;
    double kaonm_helicity, kaonm_Phi_Hup, kaonm_Phi_Hdw, kaonm_Pol, kaonm_epsilon;
    std::vector<double> kaonm_Phi_vec, kaonm_Pol_vec, kaonm_eps_vec;
    std::vector<int> kaonm_hel_vec;
    // To save all the variables
    // Electron
    treeEl->SetBranchAddress("px", &electron_px);
    treeEl->SetBranchAddress("py", &electron_py);
    treeEl->SetBranchAddress("pz", &electron_pz);
    treeEl->SetBranchAddress("Mom", &electron_mom);
    //treeEl->SetBranchAddress("vz", &electron_vz);
    treeEl->SetBranchAddress("Q2", &electron_Q2);
    //treeEl->SetBranchAddress("W", &electron_W);
    treeEl->SetBranchAddress("Theta", &electron_Theta);
    //treeEl->SetBranchAddress("Nphe", &electron_Nphe);
    //treeEl->SetBranchAddress("status", &electron_status);
    //treeEl->SetBranchAddress("PCAL", &electron_PCAL);
    //treeEl->SetBranchAddress("ECAL", &electron_ECAL);
    //treeEl->SetBranchAddress("ECIN", &electron_ECIN);
    //treeEl->SetBranchAddress("DC_6", &electron_edge1);
    //treeEl->SetBranchAddress("DC_18", &electron_edge2);
    //treeEl->SetBranchAddress("DC_36", &electron_edge3);
    // Kaon +
    treeKaonP->SetBranchAddress("E", &kaonp_E);
    treeKaonP->SetBranchAddress("px", &kaonp_px);
    treeKaonP->SetBranchAddress("py", &kaonp_py);
    treeKaonP->SetBranchAddress("pz", &kaonp_pz);
    treeKaonP->SetBranchAddress("Mom", &kaonp_Ph);
    treeKaonP->SetBranchAddress("W", &kaonp_W);
    treeKaonP->SetBranchAddress("Q2", &kaonp_Q2);
    treeKaonP->SetBranchAddress("xF", &kaonp_xF);
    treeKaonP->SetBranchAddress("xB", &kaonp_xB);
    treeKaonP->SetBranchAddress("y", &kaonp_y);
    treeKaonP->SetBranchAddress("z", &kaonp_z);
    treeKaonP->SetBranchAddress("Pt", &kaonp_Pt);
    treeKaonP->SetBranchAddress("PhT", &kaonp_PhT);
    treeKaonP->SetBranchAddress("Phi_Lab", &kaonp_Phi);
    treeKaonP->SetBranchAddress("Theta_Lab", &kaonp_Theta);
    treeKaonP->SetBranchAddress("Pseudorapidity", &kaonp_eta);
    treeKaonP->SetBranchAddress("Phi_h", &kaonp_Phi_h);
    treeKaonP->SetBranchAddress("Polarization", &kaonp_Pol);
    treeKaonP->SetBranchAddress("epsilon", &kaonp_epsilon);
    //treeKaonP->SetBranchAddress("SinPhi", &kaonp_SinPhi);
    treeKaonP->SetBranchAddress("helicity", &kaonp_helicity);
    treeKaonP->SetBranchAddress("M2_x", &kaonp_M2x);
    treeKaonP->SetBranchAddress("Rich_Id", &kaonp_rich_Id);
    treeKaonP->SetBranchAddress("Rich_PID", &kaonp_rich_PID);
    // Kaon + Helicity positive
    treeKaonPhp->SetBranchAddress("E", &kaonp_E_hp);
    treeKaonPhp->SetBranchAddress("px", &kaonp_px_hp);
    treeKaonPhp->SetBranchAddress("py", &kaonp_py_hp);
    treeKaonPhp->SetBranchAddress("pz", &kaonp_pz_hp);
    treeKaonPhp->SetBranchAddress("Mom", &kaonp_Ph_hp);
    treeKaonPhp->SetBranchAddress("W", &kaonp_W_hp);
    treeKaonPhp->SetBranchAddress("Q2", &kaonp_Q2_hp);
    treeKaonPhp->SetBranchAddress("xF", &kaonp_xF_hp);
    treeKaonPhp->SetBranchAddress("xB", &kaonp_xB_hp);
    treeKaonPhp->SetBranchAddress("y", &kaonp_y_hp);
    treeKaonPhp->SetBranchAddress("z", &kaonp_z_hp);
    treeKaonPhp->SetBranchAddress("Pt", &kaonp_Pt_hp);
    treeKaonPhp->SetBranchAddress("PhT", &kaonp_PhT_hp);
    treeKaonPhp->SetBranchAddress("Phi_Lab", &kaonp_Phi_hp);
    treeKaonPhp->SetBranchAddress("Theta_Lab", &kaonp_Theta_hp);
    treeKaonPhp->SetBranchAddress("Pseudorapidity", &kaonp_eta_hp);
    treeKaonPhp->SetBranchAddress("Phi_h", &kaonp_Phi_h_hp);
    treeKaonPhp->SetBranchAddress("helicity", &kaonp_helicity_hp);
    treeKaonPhp->SetBranchAddress("M2_x", &kaonp_M2x_hp);
    // Kaon + Helicity negative
    treeKaonPhm->SetBranchAddress("E", &kaonp_E_hm);
    treeKaonPhm->SetBranchAddress("px", &kaonp_px_hm);
    treeKaonPhm->SetBranchAddress("py", &kaonp_py_hm);
    treeKaonPhm->SetBranchAddress("pz", &kaonp_pz_hm);
    treeKaonPhm->SetBranchAddress("Mom", &kaonp_Ph_hm);
    treeKaonPhm->SetBranchAddress("W", &kaonp_W_hm);
    treeKaonPhm->SetBranchAddress("Q2", &kaonp_Q2_hm);
    treeKaonPhm->SetBranchAddress("xF", &kaonp_xF_hm);
    treeKaonPhm->SetBranchAddress("xB", &kaonp_xB_hm);
    treeKaonPhm->SetBranchAddress("y", &kaonp_y_hm);
    treeKaonPhm->SetBranchAddress("z", &kaonp_z_hm);
    treeKaonPhm->SetBranchAddress("Pt", &kaonp_Pt_hm);
    treeKaonPhm->SetBranchAddress("PhT", &kaonp_PhT_hm);
    treeKaonPhm->SetBranchAddress("Phi_Lab", &kaonp_Phi_hm);
    treeKaonPhm->SetBranchAddress("Theta_Lab", &kaonp_Theta_hm);
    treeKaonPhm->SetBranchAddress("Pseudorapidity", &kaonp_eta_hm);
    treeKaonPhm->SetBranchAddress("Phi_h", &kaonp_Phi_h_hm);
    treeKaonPhm->SetBranchAddress("helicity", &kaonp_helicity_hm);
    treeKaonPhm->SetBranchAddress("M2_x", &kaonp_M2x_hm);

    // Definizione dei bin
    const int nBinz_z = 6;
    const int nBin_Pht = 3;
    /*
    double x1_min = 0;
    double x1_max = 0.15;
    double x2_min = 0.15;
    double x2_max = 0.35;
    double x3_min = 0.16;
    double x3_max = 1;
    double Q1_min = 1;
    double Q1_max = 2.25;
    double Q2_min = 1.25;
    double Q2_max = 2.5;
    double Q3_min = 2.5;
    double Q3_max = 8;
    */
    double zBins[nBinz_z + 1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9};
    double PtBins[nBin_Pht + 1] = {0.0, 0.25, 0.5, 1.2};
    double xBins[][2] = {{0.05, 0.15}, {0.15, 0.35}, {0.16, 0.67}};
    double Q2Bins[][2] = {{1.0, 2.35}, {1.25, 2.5}, {2.5, 8.0}};
    // VECTOR FOR THE A_LU RESULTS
    vector<vector<vector<double>>> A_LU_values(3, vector<vector<double>>(nBin_Pht, vector<double>(nBinz_z, 0)));
    vector<vector<vector<double>>> A_LU_errors(3, vector<vector<double>>(nBin_Pht, vector<double>(nBinz_z, 0)));

    //
    // UNBINNED MAXIMUM LIKELIHOOD
    Long64_t nEntries_kP = treeKaonP->GetEntries();
    for (int region = 0; region < 3; region++) {
        double x_min = xBins[region][0], x_max = xBins[region][1];
        double Q2_min = Q2Bins[region][0], Q2_max = Q2Bins[region][1];
        // Pre-filtriamo gli eventi della regione xQ
        std::vector<double> phiFiltered, polFiltered, epsilonFiltered, zFiltered, PhtFiltered;
        std::vector<int> helicityFiltered;
        for (Long64_t l = 0; l < nEntries_kP; l++) {
            treeKaonP->GetEntry(l);
            if (kaonp_xB >= x_min && kaonp_xB < x_max && kaonp_Q2 >= Q2_min && kaonp_Q2 < Q2_max) {
                phiFiltered.push_back(kaonp_Phi_h);
                polFiltered.push_back(kaonp_Pol);
                epsilonFiltered.push_back(kaonp_epsilon);
                helicityFiltered.push_back(kaonp_helicity);
                zFiltered.push_back(kaonp_z);
                PhtFiltered.push_back(kaonp_PhT);
            }
        }
        for (int i = 0; i < nBin_Pht; i++) { // Loop su Pht
            for (int j = 0; j < nBinz_z; j++) { // Loop su z
                std::vector<double> phiVector, polVector, epsilonVector;
                std::vector<int> helicityVector;
                // Ora filtriamo in base ai bin di z e Pht
                for (size_t idx = 0; idx < phiFiltered.size(); idx++) {
                    double z_value = zFiltered[idx];    // New vector to store filtered z
                    double Pht_value = PhtFiltered[idx];
                    if (z_value >= zBins[j] && z_value < zBins[j+1] && Pht_value >= PtBins[i] && Pht_value < PtBins[i+1]) {
                        phiVector.push_back(phiFiltered[idx]);
                        polVector.push_back(polFiltered[idx]);
                        epsilonVector.push_back(epsilonFiltered[idx]);
                        helicityVector.push_back(helicityFiltered[idx]);
                    }
                }
                if (phiVector.empty()) continue; // Evitiamo di minimizzare con dati vuoti
                // Creiamo un solo Minimizer
                ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
                minimum->SetMaxFunctionCalls(1000000);
                minimum->SetTolerance(0.001);
                minimum->SetPrintLevel(0);

                ROOT::Math::Functor fbin([&](const double* param) { return extractFLU(param, phiVector, polVector, epsilonVector, helicityVector);}, 1);
                double FLU_MIN = -0.5;
                double FLU_MAX = 0.5;
                minimum->SetLimitedVariable(0, "A_LU", 0.01, 0.0005, FLU_MIN, FLU_MAX);
                minimum->SetFunction(fbin);
                minimum->Minimize();
                A_LU_values[region][i][j] = minimum->X()[0];
                A_LU_errors[region][i][j] = minimum->Errors()[0];
            }
        }
    }

    TGraphErrors* graphs[3][nBin_Pht];
    for (int region = 0; region < 3; region++) {
        for (int i = 0; i < nBin_Pht; i++) {
            graphs[region][i] = new TGraphErrors(nBinz_z);
            for (int j = 0; j < nBinz_z; j++) {
                graphs[region][i]->SetPoint(j, (zBins[j] + zBins[j+1]) / 2.0, A_LU_values[region][i][j]);
                graphs[region][i]->SetPointError(j, 0, A_LU_errors[region][i][j]);
            }
            graphs[region][i]->SetTitle(Form("A_LU vs z | Pht Bin %d | xQ bin %d", i+1, region+1));
            graphs[region][i]->GetXaxis()->SetTitle("z");
            graphs[region][i]->GetYaxis()->SetTitle("A_LU");
            graphs[region][i]->SetMarkerStyle(20);
            graphs[region][i]->SetLineColor(kRed); 
            graphs[region][i]->SetMarkerColor(kRed);
        }
    }
    TCanvas* c[3];
    for (int region = 0; region < 3; region++) {
        c[region] = new TCanvas(Form("A_LU(z)_xQ_bin%d", region+1), Form("A_LU vs z for xQ bin %d", region+1), 800, 800);
        c[region]->Divide(1, nBin_Pht);
        for (int i = 0; i < nBin_Pht; i++) {
            c[region]->cd(i+1);
            graphs[region][i]->Draw("AP");
            gPad->SetGridx();
            gPad->SetGridy();
        }
        c[region]->Update();
        c[region]->Write();
    }
    //

    TH1D *hZBins_Hp[3][nBinz_z], *hZBins_Hm[3][nBinz_z];
    TH1D *hZBins_Phi_Hp[3][nBinz_z][nBin_Pht], *hZBins_Phi_Hm[3][nBinz_z][nBin_Pht];
    TH1D *hPtBins_Hp[3][nBin_Pht], *hPtBins_Hm[3][nBin_Pht];

    for (int region = 0; region < 3; region++) {
        dirKaonp_Hp_xQ[region]->cd();  // Sposta nella directory H+
        for (int i = 0; i < nBinz_z; i++) {
            hZBins_Hp[region][i] = new TH1D(Form("z_Bin%d_xQ%d_Hp", i+1, region+1), 
                                            Form("z Bin %d | xQ%d | H+;z;Counts", i+1, region+1), 
                                            20, zBins[i], zBins[i+1]);
            for(int j = 0; j < nBin_Pht; j++){
                hZBins_Phi_Hp[region][i][j] = new TH1D(Form("_Phi_zBin%d_PhtBin%d_xQ%d_Hp", i+1, j+1, region+1), 
                                                    Form("#Phi_{h} in z Bin %d | P_{hT} Bin %d | xQ Bin %d | H+; #Phi_{h} [rad];Counts", i+1, j+1, region+1), 
                                                    20, -TMath::Pi(), TMath::Pi());
            }
        }

        dirKaonp_Hm_xQ[region]->cd();  // Sposta nella directory H-
        for (int i = 0; i < nBinz_z; i++) {
            hZBins_Hm[region][i] = new TH1D(Form("z_Bin%d_xQ%d_Hm", i+1, region+1), 
                                            Form("z Bin %d | xQ%d | H-;z;Counts", i+1, region+1), 
                                            20, zBins[i], zBins[i+1]);
                for(int j = 0; j < nBin_Pht; j++){
                hZBins_Phi_Hm[region][i][j] = new TH1D(Form("_Phi_zBin%d_PhtBin%d_xQ%d_Hm", i+1, j+1, region+1), 
                                                    Form("#Phi_{h} in z Bin %d | P_{hT} Bin %d| xQ Bin %d | H-; #Phi_{h} [rad];Counts", i+1, j+1, region+1), 
                                                    20, -TMath::Pi(), TMath::Pi());
                }
        }
    }
    
    Long64_t nEntries_kPhp = treeKaonPhp->GetEntries();
    for (Long64_t i = 0; i < nEntries_kPhp; i++) {
        treeKaonPhp->GetEntry(i);
        for (int region = 0; region < 3; region++) {
            if (kaonp_xB_hp >= xBins[region][0] && kaonp_xB_hp < xBins[region][1] &&
                kaonp_Q2_hp >= Q2Bins[region][0] && kaonp_Q2_hp < Q2Bins[region][1]) {
                
                for (int j = 0; j < nBinz_z; j++) {
                    if (kaonp_z_hp >= zBins[j] && kaonp_z_hp < zBins[j+1]) {
                        hZBins_Hp[region][j]->Fill(kaonp_z_hp);
                        if(kaonp_PhT_hp >= PtBins[0] && kaonp_PhT_hp < PtBins[1]) hZBins_Phi_Hp[region][j][0]->Fill(kaonp_Phi_h_hp);
                        if(kaonp_PhT_hp >= PtBins[1] && kaonp_PhT_hp < PtBins[2]) hZBins_Phi_Hp[region][j][1]->Fill(kaonp_Phi_h_hp);
                        if(kaonp_PhT_hp >= PtBins[2] && kaonp_PhT_hp < PtBins[3]) hZBins_Phi_Hp[region][j][2]->Fill(kaonp_Phi_h_hp);
                    }
                }

            }
        }
    }

    Long64_t nEntries_kPhm = treeKaonPhm->GetEntries();
    for (Long64_t i = 0; i < nEntries_kPhm; i++) {
        treeKaonPhm->GetEntry(i);
        for (int region = 0; region < 3; region++) {
            if (kaonp_xB_hm >= xBins[region][0] && kaonp_xB_hm < xBins[region][1] &&
                kaonp_Q2_hm >= Q2Bins[region][0] && kaonp_Q2_hm < Q2Bins[region][1]) {
                
                for (int j = 0; j < nBinz_z; j++) {
                    if (kaonp_z_hm >= zBins[j] && kaonp_z_hm < zBins[j+1]) {
                        hZBins_Hm[region][j]->Fill(kaonp_z_hm);
                        if(kaonp_PhT_hm >= PtBins[0] && kaonp_PhT_hm < PtBins[1]) hZBins_Phi_Hm[region][j][0]->Fill(kaonp_Phi_h_hm);
                        if(kaonp_PhT_hm >= PtBins[1] && kaonp_PhT_hm < PtBins[2]) hZBins_Phi_Hm[region][j][1]->Fill(kaonp_Phi_h_hm);
                        if(kaonp_PhT_hm >= PtBins[2] && kaonp_PhT_hm < PtBins[3]) hZBins_Phi_Hm[region][j][2]->Fill(kaonp_Phi_h_hm);
                    }
                }

            }
        }
    }

    

    outFile.Write();
    outFile.Close();
    inFile.Close();

    cout << "ROOT output file: " << outputFile << endl;
}
