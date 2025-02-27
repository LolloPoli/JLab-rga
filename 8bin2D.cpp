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
    const char* outputFile = "plot_binned_All_t-1.root";
    TFile outFile(outputFile, "RECREATE");  // File di output ROOT
    TTree *tree = (TTree*)inFile.Get("Tree");
    TTree *treeEl = (TTree*)inFile.Get("Electron");
    TTree *treeKaonP = (TTree*)inFile.Get("Kaon+");
    TTree *treeKaonPhp = (TTree*)inFile.Get("Kaon+ H+");
    TTree *treeKaonPhm = (TTree*)inFile.Get("Kaon+ H-");
    TTree *treeKaonM = (TTree*)inFile.Get("Kaon-");
    TDirectory* dirKaonp_Hp_xQ1 = outFile.mkdir("KaonPlus_H+_xQ1");
    TDirectory* dirKaonp_Hm_xQ1 = outFile.mkdir("KaonPlus_H-_xQ1");
    TDirectory* dirKaonp_Hp_xQ2 = outFile.mkdir("KaonPlus_H+_xQ2");
    TDirectory* dirKaonp_Hm_xQ2 = outFile.mkdir("KaonPlus_H-_xQ2");
    TDirectory* dirKaonp_Hp_xQ3 = outFile.mkdir("KaonPlus_H+_xQ3");
    TDirectory* dirKaonp_Hm_xQ3 = outFile.mkdir("KaonPlus_H-_xQ3");

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
    vector<vector<double>> A_LU_values(nBin_Pht, vector<double>(nBinz_z, 0));
    vector<vector<double>> A_LU_errors(nBin_Pht, vector<double>(nBinz_z, 0));

    //
    // UNBINNED MAXIMUM LIKELIHOOD
    for (int i = 0; i < nBin_Pht; i++) { // Loop su Pht
        for (int j = 0; j < nBinz_z; j++) { // Loop su z
            std::vector<double> phiVector, polVector, epsilonVector;
            std::vector<int> helicityVector;
            Long64_t nEntries_kP = treeKaonP->GetEntries();
            for (Long64_t l = 0; l < nEntries_kP; l++){
                treeKaonP->GetEntry(l);
                if (kaonp_z >= zBins[j] && kaonp_z < zBins[j+1] && kaonp_PhT >= PtBins[i] && kaonp_PhT < PtBins[i+1]) {
                    phiVector.push_back(kaonp_Phi_h);
                    polVector.push_back(kaonp_Pol);
                    epsilonVector.push_back(kaonp_epsilon);
                    helicityVector.push_back(kaonp_helicity);
                }
            }

            ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            minimum->SetMaxFunctionCalls(1000000);
            minimum->SetTolerance(0.001);
            minimum->SetPrintLevel(0);

            ROOT::Math::Functor fbin([&](const double* param) { return extractFLU(param, phiVector, polVector, epsilonVector, helicityVector); }, 1);
            double FLU_MIN = -0.5;
            double FLU_MAX = 0.5;
            minimum->SetLimitedVariable(0, "A_LU", 0.01, 0.0005, FLU_MIN, FLU_MAX);
            minimum->SetFunction(fbin);
            minimum->Minimize();
            A_LU_values[i][j] = minimum->X()[0];
            A_LU_errors[i][j] = minimum->Errors()[0];

        }
    }

    TCanvas* c1 = new TCanvas("A_LU(z)_xQ_bin1", "A_LU vs z", 800, 800);
    c1->Divide(1,3);
    TGraphErrors* graphs[nBin_Pht];
    for (int i = 0; i < nBin_Pht; i++) {
        graphs[i] = new TGraphErrors(nBinz_z);
        for (int j = 0; j < nBinz_z; j++) {
            graphs[i]->SetPoint(j, (zBins[j] + zBins[j+1]) / 2.0, A_LU_values[i][j]);
            graphs[i]->SetPointError(j, 0, A_LU_errors[i][j]);
        }
        c1->cd(i+1);
        graphs[i]->SetTitle(Form("A_LU vs z (Pht Bin %d)", i+1));
        graphs[i]->GetXaxis()->SetTitle("z");
        graphs[i]->GetYaxis()->SetTitle("A_LU");
        graphs[i]->SetMarkerStyle(20);
        graphs[i]->SetLineColor(kRed); 
        graphs[i]->SetMarkerColor(kRed);
        graphs[i]->Draw("AP");
        gPad->SetGridx();  
        gPad->SetGridy();
    }
    c1->Update();
    c1->Write();
    //

    dirKaonp_Hp_xQ1->cd();
    // Creazione degli istogrammi per i bin di z
    TH1D *hZBins_hp[nBinz_z];
    TH1D *hZBins_Phi_hp[nBinz_z];
    for (int i = 0; i < nBinz_z; i++) {
        hZBins_hp[i] = new TH1D(Form("z_Bin%d", i+1), Form("z Bin %d | kaon+ | h+;z;Counts", i+1), 20, zBins[i], zBins[i+1]);
        hZBins_Phi_hp[i] = new TH1D(Form("z_Phi_Bin%d", i+1), Form("#Phi_{h} in z Bin %d  | kaon+ | h+; #Phi_{h} [rad];Counts", i+1), 20, -TMath::Pi(), TMath::Pi());
    }
    // Creazione degli istogrammi per i bin di PhT
    TH1D *hPtBins_hp[nBin_Pht];
    TH1D *hPtBins_Phi_hp[nBin_Pht];
    for (int i = 0; i < nBin_Pht; i++) {
        hPtBins_hp[i] = new TH1D(Form("Pht_Bin%d", i+1), Form("P_{hT} Bin %d  | kaon+ | h+;P_{hT} [GeV];Counts", i+1), 20, PtBins[i], PtBins[i+1]);
        hPtBins_Phi_hp[i] = new TH1D(Form("Pht_Phi_Bin%d", i+1), Form("#Phi_{h} in P_{hT} Bin %d  | kaon+ | h+; #Phi_{h} [rad];Counts", i+1), 20, -TMath::Pi(), TMath::Pi());
    }
    //
    //
    //
    TH1D kp_sep1 ("------------ z over PhT -------------------", " running z over 3 bin of P_{hT} |  Kaon+ ;", 60, 1, 8);
    kp_sep1.Draw();
    //
    //
    TH1D *h_Phi_zPt1_hp[nBinz_z];
    TH1D *h_Phi_zPt2_hp[nBinz_z];
    TH1D *h_Phi_zPt3_hp[nBinz_z];
    for (int i = 0; i < nBinz_z; i++) {
        h_Phi_zPt1_hp[i] = new TH1D(Form("_Phi_zBin%d_PhtBin1", i+1), Form("#Phi_{h} in z Bin %d and PhT Bin 1  | kaon+ | h+; #Phi_{h} [rad];Counts", i+1), 20, -TMath::Pi(), TMath::Pi());
        h_Phi_zPt2_hp[i] = new TH1D(Form("_Phi_zBin%d_PhtBin2", i+1), Form("#Phi_{h} in z Bin %d and PhT Bin 1  | kaon+ | h+; #Phi_{h} [rad];Counts", i+1), 20, -TMath::Pi(), TMath::Pi());
        h_Phi_zPt3_hp[i] = new TH1D(Form("_Phi_zBin%d_PhtBin3", i+1), Form("#Phi_{h} in z Bin %d and PhT Bin 1  | kaon+ | h+; #Phi_{h} [rad];Counts", i+1), 20, -TMath::Pi(), TMath::Pi());
    }
    //
    dirKaonp_Hm_xQ1->cd();
    // Creazione degli istogrammi per i bin di z
    TH1D *hZBins_hm[nBinz_z];
    TH1D *hZBins_Phi_hm[nBinz_z];
    for (int i = 0; i < nBinz_z; i++) {
        hZBins_hm[i] = new TH1D(Form("z_Bin%d", i+1), Form("z Bin %d | kaon+ | h-;z;Counts", i+1), 20, zBins[i], zBins[i+1]);
        hZBins_Phi_hm[i] = new TH1D(Form("z_Phi_Bin%d", i+1), Form("#Phi_{h} in z Bin %d  | kaon+ | h-; #Phi_{h} [rad];Counts", i+1), 20, -TMath::Pi(), TMath::Pi());
    }
    // Creazione degli istogrammi per i bin di PhT
    TH1D *hPtBins_hm[nBin_Pht];
    TH1D *hPtBins_Phi_hm[nBin_Pht];
    for (int i = 0; i < nBin_Pht; i++) {
        hPtBins_hm[i] = new TH1D(Form("Pht_Bin%d", i+1), Form("P_{hT} Bin %d  | kaon+ | h-;P_{hT} [GeV];Counts", i+1), 20, PtBins[i], PtBins[i+1]);
        hPtBins_Phi_hm[i] = new TH1D(Form("Pht_Phi_Bin%d", i+1), Form("#Phi_{h} in P_{hT} Bin %d  | kaon+ | h-; #Phi_{h} [rad];Counts", i+1), 20, -TMath::Pi(), TMath::Pi());
    }
    //
    TH1D *h_Phi_zPt1_hm[nBinz_z];
    TH1D *h_Phi_zPt2_hm[nBinz_z];
    TH1D *h_Phi_zPt3_hm[nBinz_z];
    for (int i = 0; i < nBinz_z; i++) {
        h_Phi_zPt1_hm[i] = new TH1D(Form("_Phi_zBin%d_PhtBin1", i+1), Form("#Phi_{h} in z Bin %d and PhT Bin 1  | kaon+ | h-; #Phi_{h} [rad];Counts", i+1), 20, -TMath::Pi(), TMath::Pi());
        h_Phi_zPt2_hm[i] = new TH1D(Form("_Phi_zBin%d_PhtBin2", i+1), Form("#Phi_{h} in z Bin %d and PhT Bin 1  | kaon+ | h-; #Phi_{h} [rad];Counts", i+1), 20, -TMath::Pi(), TMath::Pi());
        h_Phi_zPt3_hm[i] = new TH1D(Form("_Phi_zBin%d_PhtBin3", i+1), Form("#Phi_{h} in z Bin %d and PhT Bin 1  | kaon+ | h-; #Phi_{h} [rad];Counts", i+1), 20, -TMath::Pi(), TMath::Pi());
    }
    //
    /*
    Long64_t nEntries_kP = treeKaonP->GetEntries();
    for (Long64_t i = 0; i < nEntries_kP; i++){
        treeKaonP->GetEntry(i);
        kaonp_Phi_vec.push_back(kaonp_Phi_h);
        kaonp_Pol_vec.push_back(kaonp_Pol);
        kaonp_eps_vec.push_back(kaonp_epsilon);
        kaonp_hel_vec.push_back(kaonp_helicity);
    }
    */
    // combination H+
    Long64_t nEntries_kPhp = treeKaonPhp->GetEntries();
    for (Long64_t i = 0; i < nEntries_kPhp; i++) {
        treeKaonPhp->GetEntry(i);
        for (int j = 0; j < nBinz_z; j++) {
            if (kaonp_z_hp >= zBins[j] && kaonp_z_hp < zBins[j+1]) {
                hZBins_hp[j]->Fill(kaonp_z_hp);
                hZBins_Phi_hp[j]->Fill(kaonp_Phi_h_hp);
                if(kaonp_PhT_hp >= PtBins[0] && kaonp_PhT_hp < PtBins[1]) h_Phi_zPt1_hp[j]->Fill(kaonp_Phi_h_hp);
                if(kaonp_PhT_hp >= PtBins[1] && kaonp_PhT_hp < PtBins[2]) h_Phi_zPt2_hp[j]->Fill(kaonp_Phi_h_hp);
                if(kaonp_PhT_hp >= PtBins[2] && kaonp_PhT_hp < PtBins[3]) h_Phi_zPt3_hp[j]->Fill(kaonp_Phi_h_hp);
            }
        }
        for (int j = 0; j < nBin_Pht; j++) {
            if (kaonp_PhT_hp >= PtBins[j] && kaonp_PhT_hp < PtBins[j+1]) {
                hPtBins_hp[j]->Fill(kaonp_PhT_hp);
                hPtBins_Phi_hp[j]->Fill(kaonp_Phi_h_hp);
            }
        }
    }
    // combination H-
    Long64_t nEntries_kPhm = treeKaonPhm->GetEntries();
    for (Long64_t i = 0; i < nEntries_kPhm; i++) {
        treeKaonPhm->GetEntry(i);
        for (int j = 0; j < nBinz_z; j++) {
            if (kaonp_z_hm >= zBins[j] && kaonp_z_hm < zBins[j+1]) {
                hZBins_hm[j]->Fill(kaonp_z_hm);
                hZBins_Phi_hm[j]->Fill(kaonp_Phi_h_hm);
                if(kaonp_PhT_hm >= PtBins[0] && kaonp_PhT_hm < PtBins[1]) h_Phi_zPt1_hm[j]->Fill(kaonp_Phi_h_hm);
                if(kaonp_PhT_hm >= PtBins[1] && kaonp_PhT_hm < PtBins[2]) h_Phi_zPt2_hm[j]->Fill(kaonp_Phi_h_hm);
                if(kaonp_PhT_hm >= PtBins[2] && kaonp_PhT_hm < PtBins[3]) h_Phi_zPt3_hm[j]->Fill(kaonp_Phi_h_hm);
            }
        }
        for (int j = 0; j < nBin_Pht; j++) {
            if (kaonp_PhT_hm >= PtBins[j] && kaonp_PhT_hm < PtBins[j+1]) {
                hPtBins_hm[j]->Fill(kaonp_PhT_hm);
                hPtBins_Phi_hm[j]->Fill(kaonp_Phi_h_hm);
            }
        }
    }


    /*
    // Minimizer
    ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimum->SetMaxFunctionCalls(1000000);
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(1);

    ROOT::Math::Functor fbin([&](const double* param) { return extractFLU(param, kaonp_Phi_vec, kaonp_Pol_vec, kaonp_eps_vec, kaonp_hel_vec); }, 1);
    double FLU_MIN = -0.5;
    double FLU_MAX = 0.5;
    minimum->SetLimitedVariable(0, "A_LU", 0.01, 0.0005, FLU_MIN, FLU_MAX);
    minimum->SetFunction(fbin);
    minimum->Minimize();
    const double result = minimum->X()[0];
    const double error = minimum->Errors()[0];
    */

    outFile.Write();
    outFile.Close();
    inFile.Close();

    cout << "ROOT output file: " << outputFile << endl;
}
