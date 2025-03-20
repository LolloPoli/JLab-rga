#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <iostream>

using namespace std;
namespace fs = std::filesystem;
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

void SetStatsBoxGraph(TGraphErrors* graph, int nPoints) {
    TPaveText* statsBox = new TPaveText(0.80, 0.58, 0.98, 0.92, "NDC");
    statsBox->SetFillColor(0);  // Sfondo trasparente
    statsBox->SetTextAlign(12); // Allineamento del testo a sinistra
    // Calcola statistiche base (media e RMS su asse Y)
    double mean = 0, rms = 0, sum = 0, sum2 = 0;
    double x, y;
    for (int i = 0; i < nPoints; i++) {
        graph->GetPoint(i, x, y);
        sum += y;
        sum2 += y * y;
    }
    if (nPoints > 0) {
        mean = sum / nPoints;
        rms = sqrt(sum2 / nPoints - mean * mean);
    }
    // Aggiungi testi al box
    statsBox->AddText(Form("Entries: %d", nPoints));
    //statsBox->AddText(Form("Mean: %.4f", mean));
    //statsBox->AddText(Form("RMS: %.4f", rms));
    // Disegna il box
    statsBox->Draw();
}

// Funzione di verosimiglianza negativa da minimizzare
double extractFLU(const double *FLU, const std::vector<double>& phi, const std::vector<double>& bpol, const std::vector<double>& epsilon, const std::vector<int>& helicity) {
    double logLike = 0.0;
    for (size_t i = 0; i < phi.size(); i++) {
        // deve tendere a 1 siccome è una probabilità normalizzata 
        double term = 1 + (helicity[i] * bpol[i]) * FLU[0] * sin(phi[i]);
        // Evitiamo log(0) o log(valori negativi)
        if (term <= 0) term = 1e-10;  
        // lavoriamo con il negativo della log-verosimiglianza siccome minuit minimizza (cercassimo il massimo no)
        logLike += -2 * TMath::Log(term);
    }
    return logLike;
}

double extractFLU_2(const double *FLU, const std::vector<double>& phi, const std::vector<double>& bpol, const std::vector<double>& epsilon, const std::vector<int>& helicity) {
    double logLike = 0.0;
    for (size_t i = 0; i < phi.size(); i++) {
        double numerator = (helicity[i] * bpol[i]) * FLU[0] * sin(phi[i]);
        double denominator = 1 + FLU[1] * cos(phi[i]) + FLU[2] * cos(2 * phi[i]);
        double term = 1 + numerator / denominator;
        //double denominator = FLU[1] * cos(phi[i]) + FLU[2] * cos(2 * phi[i]);
        //double denominator = FLU[1]*sin(2*phi[i]) + FLU[2] * cos(phi[i]);
        //double term = 1 + numerator + denominator;

        // Evitiamo log(0) o valori negativi
        if (denominator <= 1e-8) denominator = 1e-8;
        if (term <= 1e-10) term = 1e-10;

        logLike += -2 * TMath::Log(term);
    }
    return logLike;
}

void _6binning(double torus, bool rich_yes) {
    string inputDir;
    if(rich_yes) inputDir = "output_fall2018_torus-1_rich"; 
    else{
        if(torus == -1) {
            inputDir = "output_fall2018_torus-1"; 
        } else if(torus == +1){
            inputDir = "output_fall2018_torus+1";
        }
    }
    TChain chainTree("Tree");
    TChain chainElectron("Electron");
    TChain chainKaonP("Kaon+");
    TChain chainKaonPhp("Kaon+ H+");
    TChain chainKaonPhm("Kaon+ H-");
    TChain chainKaonM("Kaon-");
    int fileCount = 0;
    for (const auto &entry : fs::directory_iterator(inputDir)) {
        if (entry.path().extension() == ".root") {
            string filePath = entry.path().string();
            chainTree.Add(Form("%s/Tree", filePath.c_str()));
            chainElectron.Add(Form("%s/Electron", filePath.c_str()));
            chainKaonP.Add(Form("%s/Kaon+", filePath.c_str()));
            chainKaonPhp.Add(Form("%s/Kaon+ H+", filePath.c_str()));
            chainKaonPhm.Add(Form("%s/Kaon+ H-", filePath.c_str()));
            chainKaonM.Add(Form("%s/Kaon-", filePath.c_str()));
            fileCount++;
        }
    }
    if (fileCount == 0) {
        cerr << "Nessun file .root trovato in " << inputDir << endl;
        return;
    }

    const char* outputFile; 
    if(rich_yes)outputFile = "plot_binned_rga_t-1_rich.root";
    else{
        if(torus == -1){
            outputFile = "plot_binned_rga_t-1.root";
        } else if(torus == +1){
            outputFile = "plot_binned_rga_t+1.root";
        }
    }
    TFile outFile(outputFile, "RECREATE");  // File di output ROOT

    TDirectory* dirKaonp_Hp_xQ[3];
    TDirectory* dirKaonp_Hm_xQ[3];
    for (int region = 0; region < 3; region++) {
        dirKaonp_Hp_xQ[region] = outFile.mkdir(Form("KaonPlus_H+_xQ%d", region+1));
        dirKaonp_Hm_xQ[region] = outFile.mkdir(Form("KaonPlus_H-_xQ%d", region+1));
    }
    TDirectory* dirBin_z = outFile.mkdir("z Binning Population");
    TDirectory* dirBin_PhT = outFile.mkdir("PhT Binning Population");
    TDirectory* dirAsym = outFile.mkdir("A_LU & A_UU");

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
    chainElectron.SetBranchAddress("px", &electron_px);
    chainElectron.SetBranchAddress("py", &electron_py);
    chainElectron.SetBranchAddress("pz", &electron_pz);
    chainElectron.SetBranchAddress("Mom", &electron_mom);
    //chainElectron.SetBranchAddress("vz", &electron_vz);
    chainElectron.SetBranchAddress("Q2", &electron_Q2);
    //chainElectron.SetBranchAddress("W", &electron_W);
    chainElectron.SetBranchAddress("Theta", &electron_Theta);
    // Kaon +
    chainKaonP.SetBranchAddress("E", &kaonp_E);
    chainKaonP.SetBranchAddress("px", &kaonp_px);
    chainKaonP.SetBranchAddress("py", &kaonp_py);
    chainKaonP.SetBranchAddress("pz", &kaonp_pz);
    chainKaonP.SetBranchAddress("Mom", &kaonp_Ph);
    chainKaonP.SetBranchAddress("W", &kaonp_W);
    chainKaonP.SetBranchAddress("Q2", &kaonp_Q2);
    chainKaonP.SetBranchAddress("xF", &kaonp_xF);
    chainKaonP.SetBranchAddress("xB", &kaonp_xB);
    chainKaonP.SetBranchAddress("y", &kaonp_y);
    chainKaonP.SetBranchAddress("z", &kaonp_z);
    chainKaonP.SetBranchAddress("Pt", &kaonp_Pt);
    chainKaonP.SetBranchAddress("PhT", &kaonp_PhT);
    chainKaonP.SetBranchAddress("Phi_Lab", &kaonp_Phi);
    chainKaonP.SetBranchAddress("Theta_Lab", &kaonp_Theta);
    chainKaonP.SetBranchAddress("Pseudorapidity", &kaonp_eta);
    chainKaonP.SetBranchAddress("Phi_h", &kaonp_Phi_h);
    chainKaonP.SetBranchAddress("Polarization", &kaonp_Pol);
    chainKaonP.SetBranchAddress("epsilon", &kaonp_epsilon);
    //chainKaonP.SetBranchAddress("SinPhi", &kaonp_SinPhi);
    chainKaonP.SetBranchAddress("helicity", &kaonp_helicity);
    chainKaonP.SetBranchAddress("M2_x", &kaonp_M2x);
    chainKaonP.SetBranchAddress("Rich_Id", &kaonp_rich_Id);
    chainKaonP.SetBranchAddress("Rich_PID", &kaonp_rich_PID);
    // Kaon + Helicity positive
    chainKaonPhp.SetBranchAddress("E", &kaonp_E_hp);
    chainKaonPhp.SetBranchAddress("px", &kaonp_px_hp);
    chainKaonPhp.SetBranchAddress("py", &kaonp_py_hp);
    chainKaonPhp.SetBranchAddress("pz", &kaonp_pz_hp);
    chainKaonPhp.SetBranchAddress("Mom", &kaonp_Ph_hp);
    chainKaonPhp.SetBranchAddress("W", &kaonp_W_hp);
    chainKaonPhp.SetBranchAddress("Q2", &kaonp_Q2_hp);
    chainKaonPhp.SetBranchAddress("xF", &kaonp_xF_hp);
    chainKaonPhp.SetBranchAddress("xB", &kaonp_xB_hp);
    chainKaonPhp.SetBranchAddress("y", &kaonp_y_hp);
    chainKaonPhp.SetBranchAddress("z", &kaonp_z_hp);
    chainKaonPhp.SetBranchAddress("Pt", &kaonp_Pt_hp);
    chainKaonPhp.SetBranchAddress("PhT", &kaonp_PhT_hp);
    chainKaonPhp.SetBranchAddress("Phi_Lab", &kaonp_Phi_hp);
    chainKaonPhp.SetBranchAddress("Theta_Lab", &kaonp_Theta_hp);
    chainKaonPhp.SetBranchAddress("Pseudorapidity", &kaonp_eta_hp);
    chainKaonPhp.SetBranchAddress("Phi_h", &kaonp_Phi_h_hp);
    chainKaonPhp.SetBranchAddress("helicity", &kaonp_helicity_hp);
    chainKaonPhp.SetBranchAddress("M2_x", &kaonp_M2x_hp);
    // Kaon + Helicity negative
    chainKaonPhm.SetBranchAddress("E", &kaonp_E_hm);
    chainKaonPhm.SetBranchAddress("px", &kaonp_px_hm);
    chainKaonPhm.SetBranchAddress("py", &kaonp_py_hm);
    chainKaonPhm.SetBranchAddress("pz", &kaonp_pz_hm);
    chainKaonPhm.SetBranchAddress("Mom", &kaonp_Ph_hm);
    chainKaonPhm.SetBranchAddress("W", &kaonp_W_hm);
    chainKaonPhm.SetBranchAddress("Q2", &kaonp_Q2_hm);
    chainKaonPhm.SetBranchAddress("xF", &kaonp_xF_hm);
    chainKaonPhm.SetBranchAddress("xB", &kaonp_xB_hm);
    chainKaonPhm.SetBranchAddress("y", &kaonp_y_hm);
    chainKaonPhm.SetBranchAddress("z", &kaonp_z_hm);
    chainKaonPhm.SetBranchAddress("Pt", &kaonp_Pt_hm);
    chainKaonPhm.SetBranchAddress("PhT", &kaonp_PhT_hm);
    chainKaonPhm.SetBranchAddress("Phi_Lab", &kaonp_Phi_hm);
    chainKaonPhm.SetBranchAddress("Theta_Lab", &kaonp_Theta_hm);
    chainKaonPhm.SetBranchAddress("Pseudorapidity", &kaonp_eta_hm);
    chainKaonPhm.SetBranchAddress("Phi_h", &kaonp_Phi_h_hm);
    chainKaonPhm.SetBranchAddress("helicity", &kaonp_helicity_hm);
    chainKaonPhm.SetBranchAddress("M2_x", &kaonp_M2x_hm);

    TH1D *hist_count_xQ2 = new TH1D("hist_count_xQ2_bin", "Number of events in the 3 bin of x_{B} and Q^{2}; Bin; counts", 3, 0.5, 3.5);
    // Definizione dei bin
    const int nBin_z = 6;
    const int nBin_Pht = 3;
    const int nBin_z_2 = 3;
    const int nBin_Pht_2 = 9;
    double zBins[nBin_z + 1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9};
    double PtBins_2[nBin_Pht_2 + 1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1};
    //double PtBins[nBin_Pht + 1] = {0.0, 0.25, 0.5, 1.2};
    double PtBins[nBin_Pht + 1];
    double zBins_2[nBin_z_2 + 1];
    if(rich_yes){
        PtBins[0] = 0.0, zBins_2[0] = 0.2;
        PtBins[1] = 0.25, zBins_2[1] = 0.35;
        PtBins[2] = 0.4, zBins_2[2] = 0.45;
        PtBins[3] = 1.2, zBins_2[3] = 0.9;
    }
    else{
        if(torus == -1){
            PtBins[0] = 0.0, zBins_2[0] = 0.2;
            PtBins[1] = 0.25, zBins_2[1] = 0.4;
            PtBins[2] = 0.5, zBins_2[2] = 0.55;
            PtBins[3] = 1.2, zBins_2[3] = 0.9;
        } else if(torus == +1){
            PtBins[0] = 0.0;
            PtBins[1] = 0.35;
            PtBins[2] = 0.55;
            PtBins[3] = 1.2;
        }
    }
    std::vector<std::vector<double>> xBins;
    std::vector<std::vector<double>> Q2Bins;
    if(torus == -1){
        xBins = {{0.05, 0.18}, {0.18, 0.44}, {0.18, 0.7}};
        Q2Bins = {{1.0, 2.715}, {1.25, 3}, {3, 9.0}};
    } else if(torus == +1){
        xBins = {{0.05, 0.13}, {0.13, 0.4}, {0.15, 0.7}};
        Q2Bins = {{1.0, 2.2}, {1, 2.4}, {2.4, 9.0}};
    }
    TH1D *hist_count_z[3][3]; // bin di xQ e PhT
    TH1D *hist_count_PhT[3][3];
    for (int region = 0; region < 3; region++) {
        for (int i = 0; i < 3; i++) {
            hist_count_z[region][i] = new TH1D(Form("hist_count_xQ%d_Pht%d", region+1, i+1),
                Form("Bin population | xQ^{2} bin %d | P_{hT} bin %d;z;counts", region+1, i+1), nBin_z, zBins);  // 6 bin di z tra 0.2 e 0.9
            hist_count_PhT[region][i] = new TH1D(Form("hist_count_xQ%d_z%d", region+1, i+1),
                Form("Bin population | xQ^{2} bin %d | z bin %d;P_{hT};counts", region+1, i+1), nBin_Pht_2, PtBins_2);  
        }
    }
    
    // VECTOR FOR THE A_LU RESULTS
    // z
    vector<vector<vector<double>>> A_LU_values(3, vector<vector<double>>(nBin_Pht, vector<double>(nBin_z, 0)));
    vector<vector<vector<double>>> A_LU_errors(3, vector<vector<double>>(nBin_Pht, vector<double>(nBin_z, 0)));
    vector<vector<vector<double>>> A_LU_values_2(3, vector<vector<double>>(nBin_Pht, vector<double>(nBin_z, 0)));
    vector<vector<vector<double>>> A_LU_errors_2(3, vector<vector<double>>(nBin_Pht, vector<double>(nBin_z, 0)));
    vector<vector<vector<double>>> B_values_2(3, vector<vector<double>>(nBin_Pht, vector<double>(nBin_z, 0)));
    vector<vector<vector<double>>> B_errors_2(3, vector<vector<double>>(nBin_Pht, vector<double>(nBin_z, 0)));
    vector<vector<vector<double>>> C_values_2(3, vector<vector<double>>(nBin_Pht, vector<double>(nBin_z, 0)));
    vector<vector<vector<double>>> C_errors_2(3, vector<vector<double>>(nBin_Pht, vector<double>(nBin_z, 0)));
    // PhT
    vector<vector<vector<double>>> A_LU_PhT_values(3, vector<vector<double>>(nBin_z_2, vector<double>(nBin_Pht_2, 0)));
    vector<vector<vector<double>>> A_LU_PhT_errors(3, vector<vector<double>>(nBin_z_2, vector<double>(nBin_Pht_2, 0)));
    vector<vector<vector<double>>> A_LU_PhT_values_2(3, vector<vector<double>>(nBin_z_2, vector<double>(nBin_Pht_2, 0)));
    vector<vector<vector<double>>> A_LU_PhT_errors_2(3, vector<vector<double>>(nBin_z_2, vector<double>(nBin_Pht_2, 0)));
    vector<vector<vector<double>>> B_PhT_values_2(3, vector<vector<double>>(nBin_z_2, vector<double>(nBin_Pht_2, 0)));
    vector<vector<vector<double>>> B_PhT_errors_2(3, vector<vector<double>>(nBin_z_2, vector<double>(nBin_Pht_2, 0)));
    vector<vector<vector<double>>> C_PhT_values_2(3, vector<vector<double>>(nBin_z_2, vector<double>(nBin_Pht_2, 0)));
    vector<vector<vector<double>>> C_PhT_errors_2(3, vector<vector<double>>(nBin_z_2, vector<double>(nBin_Pht_2, 0)));

    //
    // UNBINNED MAXIMUM LIKELIHOOD
    Long64_t nEntries_kP = chainKaonP.GetEntries();
    for (int region = 0; region < 3; region++) {
        double x_min, Q2_min, x_max, Q2_max;
        x_min = xBins[region][0], x_max = xBins[region][1];
        Q2_min = Q2Bins[region][0], Q2_max = Q2Bins[region][1];
        // Pre-filtriamo gli eventi della regione xQ
        std::vector<double> phiFiltered, polFiltered, epsilonFiltered, zFiltered, PhtFiltered;
        std::vector<int> helicityFiltered;
        int eventBin_xQ2_count = 0;
        for (Long64_t l = 0; l < nEntries_kP; l++) {
            chainKaonP.GetEntry(l);
            if (kaonp_xB >= x_min && kaonp_xB < x_max && kaonp_Q2 >= Q2_min && kaonp_Q2 < Q2_max) {
                phiFiltered.push_back(kaonp_Phi_h);
                polFiltered.push_back(kaonp_Pol);
                epsilonFiltered.push_back(kaonp_epsilon);
                helicityFiltered.push_back(kaonp_helicity);
                zFiltered.push_back(kaonp_z);
                PhtFiltered.push_back(kaonp_PhT);
                eventBin_xQ2_count++;
                for (int i = 0; i < 3; i++) {  // Loop su Pht
                    if (kaonp_PhT >= PtBins[i] && kaonp_PhT < PtBins[i+1]) {
                        hist_count_z[region][i]->Fill(kaonp_z);
                    }
                }
            }
        }
        hist_count_xQ2->SetBinContent(region + 1, eventBin_xQ2_count);
        for (int i = 0; i < nBin_Pht; i++) { // Loop su Pht
            for (int j = 0; j < nBin_z; j++) { // Loop su z
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
                // Minimizer, Minuit2 metodo usato in particle physics, Migrad è l'algoritmo che usa il gradiente per cercare il minimo
                ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
                minimum->SetMaxFunctionCalls(1000000); // iterazioni
                minimum->SetTolerance(0.0005);         // tolleranza per la convergenza
                minimum->SetPrintLevel(0);             // non voglio output
                
                // [&] cattura per riferimento -> prende i vettori dei dati senza doverli copiare
                ROOT::Math::Functor fbin([&](const double* param) { return extractFLU(param, phiVector, polVector, epsilonVector, helicityVector);}, 1);
                double FLU_MIN = -0.5, FLU_MAX = 0.5;
                minimum->SetLimitedVariable(0, "A_LU", 0.01, 0.0005, FLU_MIN, FLU_MAX);
                minimum->SetFunction(fbin);
                minimum->Minimize();
                A_LU_values[region][i][j] = minimum->X()[0];
                A_LU_errors[region][i][j] = minimum->Errors()[0];

                ROOT::Math::Functor fbin2([&](const double* param) { return extractFLU_2(param, phiVector, polVector, epsilonVector, helicityVector);}, 3);
                double B_MIN = -0.8, B_MAX = 0.8;
                double C_MIN = -0.8, C_MAX = 0.8;
                minimum->SetLimitedVariable(0, "A_LU", 0.01, 0.005, FLU_MIN, FLU_MAX);
                minimum->SetLimitedVariable(1, "B", 0.01, 0.005, B_MIN, B_MAX);
                minimum->SetLimitedVariable(2, "C", 0.01, 0.005, C_MIN, C_MAX);
                minimum->SetFunction(fbin2);
                minimum->Minimize();
                A_LU_values_2[region][i][j] = minimum->X()[0];
                A_LU_errors_2[region][i][j] = minimum->Errors()[0];
                B_values_2[region][i][j] = minimum->X()[1];
                B_errors_2[region][i][j] = minimum->Errors()[1];
                C_values_2[region][i][j] = minimum->X()[2];
                C_errors_2[region][i][j] = minimum->Errors()[2];
            }
        }
    }
    
    hist_count_xQ2->GetYaxis()->SetRangeUser(0, 3e6);
    hist_count_xQ2->Write();
    //
    dirAsym->cd();
    TGraphErrors* graphs[3][nBin_Pht];
    for (int region = 0; region < 3; region++) {
        for (int i = 0; i < nBin_Pht; i++) {
            graphs[region][i] = new TGraphErrors(nBin_z);
            for (int j = 0; j < nBin_z; j++) {
                graphs[region][i]->SetPoint(j, (zBins[j] + zBins[j+1]) / 2.0, A_LU_values[region][i][j]);
                graphs[region][i]->SetPointError(j, 0, A_LU_errors[region][i][j]);
            }
            graphs[region][i]->SetTitle(Form("Unbinned A_{LU} vs z | Pht Bin %d | xQ bin %d", i+1, region+1));
            graphs[region][i]->GetXaxis()->SetTitle("z");
            graphs[region][i]->GetYaxis()->SetTitle("A_LU");
            graphs[region][i]->SetMarkerStyle(20);
            graphs[region][i]->SetLineColor(kRed+1); 
            graphs[region][i]->SetMarkerColor(kRed+1);
            graphs[region][i]->GetYaxis()->SetRangeUser(0, 0.06);
            graphs[region][i]->GetYaxis()->SetMaxDigits(3);
            SetStatsBoxGraph(graphs[region][i], nBin_z);
        }
    }
    TCanvas* c[3];
    TCanvas* cline1 = new TCanvas("------------- z -------------", "", 800,800);
    cline1->Update();
    cline1->Write();
    for (int region = 0; region < 3; region++) {
        c[region] = new TCanvas(Form("A_LU_z_xQ_bin%d", region+1), Form("Unbinned A_{LU} vs z for xQ bin %d", region+1), 800, 800);
        c[region]->Divide(1, nBin_Pht);
        for (int i = 0; i < nBin_Pht; i++) {
            c[region]->cd(i+1);
            graphs[region][i]->Draw("AP");
            //SetStatsBoxGraph(graphs[region][i], nBin_z);
            gPad->SetGridx();
            gPad->SetGridy();
            gPad->Update();
        }
        c[region]->Update();
        c[region]->Write();
    }
    // more defined fit
    TGraphErrors* graphs_2[3][nBin_Pht];
    TGraphErrors* graphs_Auu1[3][nBin_Pht];
    TGraphErrors* graphs_Auu2[3][nBin_Pht];
    for (int region = 0; region < 3; region++) {
        for (int i = 0; i < nBin_Pht; i++) {
            graphs_2[region][i] = new TGraphErrors(nBin_z);
            graphs_Auu1[region][i] = new TGraphErrors(nBin_z);
            graphs_Auu2[region][i] = new TGraphErrors(nBin_z);
            for (int j = 0; j < nBin_z; j++) {
                graphs_2[region][i]->SetPoint(j, (zBins[j] + zBins[j+1]) / 2.0, A_LU_values_2[region][i][j]);
                graphs_2[region][i]->SetPointError(j, 0, A_LU_errors_2[region][i][j]);
                graphs_Auu1[region][i]->SetPoint(j, (zBins[j] + zBins[j+1]) / 2.0, B_values_2[region][i][j]);
                graphs_Auu1[region][i]->SetPointError(j, 0, B_errors_2[region][i][j]);
                graphs_Auu2[region][i]->SetPoint(j, (zBins[j] + zBins[j+1]) / 2.0, C_values_2[region][i][j]);
                graphs_Auu2[region][i]->SetPointError(j, 0, C_errors_2[region][i][j]);
            }
            graphs_2[region][i]->SetTitle(Form("Unbinned A_{LU} vs z | P_{hT} Bin %d | xQ bin %d", i+1, region+1));
            graphs_2[region][i]->GetXaxis()->SetTitle("z");
            graphs_2[region][i]->GetYaxis()->SetTitle("A_LU");
            graphs_2[region][i]->SetMarkerStyle(20);
            graphs_2[region][i]->SetLineColor(kBlue-2); 
            graphs_2[region][i]->SetMarkerColor(kBlue-2);
            graphs_2[region][i]->GetYaxis()->SetRangeUser(0, 0.06);
            graphs_2[region][i]->GetYaxis()->SetMaxDigits(3);
            // A_UU
            graphs_Auu1[region][i]->SetTitle(Form("Unbinned A_{UU} vs z | P_{hT} Bin %d | xQ bin %d", i+1, region+1));
            graphs_Auu1[region][i]->GetXaxis()->SetTitle("z");
            graphs_Auu1[region][i]->GetYaxis()->SetTitle("A_UU");
            graphs_Auu1[region][i]->SetLineColor(kGreen+3); 
            graphs_Auu1[region][i]->SetMarkerColor(kGreen+3);
            graphs_Auu1[region][i]->SetMarkerStyle(20);
            graphs_Auu1[region][i]->GetYaxis()->SetRangeUser(-0.5, 0.5);
            graphs_Auu2[region][i]->GetYaxis()->SetRangeUser(-1, 1);
            graphs_Auu2[region][i]->SetMarkerStyle(20);
            graphs_Auu2[region][i]->SetLineColor(kMagenta+3); 
            graphs_Auu2[region][i]->SetMarkerColor(kMagenta+3);
            
        }
    }
    TCanvas* c_2[3];
    for (int region = 0; region < 3; region++) {
        c_2[region] = new TCanvas(Form("A_LU_z_xQ_bin%d_NewFit", region+1), Form("Unbinned A_{LU} vs z for xQ bin %d", region+1), 800, 800);
        c_2[region]->Divide(1, nBin_Pht);
        for (int i = 0; i < nBin_Pht; i++) {
            c_2[region]->cd(i+1);
            graphs_2[region][i]->Draw("AP");
            graphs[region][i]->Draw("SAME P");
            //SetStatsBoxGraph(graphs[region][i], nBin_z);
            TLegend *leg = new TLegend(0.85, 0.65, 0.95, 0.85);
            leg->AddEntry(graphs[region][i], "A_{LU} fit_1", "ep");
            leg->AddEntry(graphs_2[region][i], "A_{LU} fit_2", "epf");
            leg->Draw();
            gPad->SetGridx();
            gPad->SetGridy();
            gPad->Update();
        }
        c_2[region]->Update();
        c_2[region]->Write();
    }
    // canvas per Auu
    TCanvas* c_3[3];
    for (int region = 0; region < 3; region++) {
        c_3[region] = new TCanvas(Form("A_UU_z_xQ_bin%d_NewFit", region+1), Form("Unbinned A_{UU} vs z for xQ bin %d", region+1), 800, 800);
        c_3[region]->Divide(1, nBin_Pht);
        for (int i = 0; i < nBin_Pht; i++) {
            c_3[region]->cd(i+1);
            graphs_Auu1[region][i]->Draw("AP");
            graphs_Auu2[region][i]->Draw("SAME P");
            //SetStatsBoxGraph(graphs[region][i], nBin_z);
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
    //
    // UNBINNED MAXIMUM LIKELIHOOD
    //Long64_t nEntries_kP = chainKaonP.GetEntries();
    // _______________________________________________________________ PhT _______________________________________________________________________________

    for (int region = 0; region < 3; region++) {
        double x_min, Q2_min, x_max, Q2_max;
        x_min = xBins[region][0], x_max = xBins[region][1];
        Q2_min = Q2Bins[region][0], Q2_max = Q2Bins[region][1];
        // Pre-filtriamo gli eventi della regione xQ
        std::vector<double> phiFiltered, polFiltered, epsilonFiltered, zFiltered, PhtFiltered;
        std::vector<int> helicityFiltered;
        //int eventBin_xQ2_count = 0;
        for (Long64_t l = 0; l < nEntries_kP; l++) {
            chainKaonP.GetEntry(l);
            if (kaonp_xB >= x_min && kaonp_xB < x_max && kaonp_Q2 >= Q2_min && kaonp_Q2 < Q2_max) {
                phiFiltered.push_back(kaonp_Phi_h);
                polFiltered.push_back(kaonp_Pol);
                epsilonFiltered.push_back(kaonp_epsilon);
                helicityFiltered.push_back(kaonp_helicity);
                zFiltered.push_back(kaonp_z);
                PhtFiltered.push_back(kaonp_PhT);
                //eventBin_xQ2_count++;
                for (int i = 0; i < 3; i++) {  // Loop su Pht
                    if (kaonp_z >= zBins_2[i] && kaonp_z < zBins_2[i+1]) {
                        hist_count_PhT[region][i]->Fill(kaonp_PhT);
                    }
                }
            }
        }
        //hist_count_xQ2->SetBinContent(region + 1, eventBin_xQ2_count);
        for (int i = 0; i < nBin_z_2; i++) { // Loop su z
            for (int j = 0; j < nBin_Pht_2; j++) { // Loop su PhT
                std::vector<double> phiVector, polVector, epsilonVector;
                std::vector<int> helicityVector;
                // Ora filtriamo in base ai bin di z e Pht
                for (size_t idx = 0; idx < phiFiltered.size(); idx++) {
                    double z_value = zFiltered[idx];    // New vector to store filtered z
                    double Pht_value = PhtFiltered[idx];
                    if (Pht_value >= PtBins_2[j] && Pht_value < PtBins_2[j+1] && z_value >= zBins_2[i] && z_value < zBins_2[i+1]) {
                        phiVector.push_back(phiFiltered[idx]);
                        polVector.push_back(polFiltered[idx]);
                        epsilonVector.push_back(epsilonFiltered[idx]);
                        helicityVector.push_back(helicityFiltered[idx]);
                    }
                }
                if (phiVector.empty()) continue; // Evitiamo di minimizzare con dati vuoti
                // Minimizer, Minuit2 metodo usato in particle physics, Migrad è l'algoritmo che usa il gradiente per cercare il minimo
                ROOT::Math::Minimizer* minimum2 = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
                minimum2->SetMaxFunctionCalls(1000000); // iterazioni
                minimum2->SetTolerance(0.0005);         // tolleranza per la convergenza
                minimum2->SetPrintLevel(0);             // non voglio output
                
                // [&] cattura per riferimento -> prende i vettori dei dati senza doverli copiare
                ROOT::Math::Functor fbin_p([&](const double* param) { return extractFLU(param, phiVector, polVector, epsilonVector, helicityVector);}, 1);
                double FLU_MIN = -0.5, FLU_MAX = 0.5;
                minimum2->SetLimitedVariable(0, "A_LU", 0.01, 0.0005, FLU_MIN, FLU_MAX);
                minimum2->SetFunction(fbin_p);
                minimum2->Minimize();
                A_LU_PhT_values[region][i][j] = minimum2->X()[0];
                A_LU_PhT_errors[region][i][j] = minimum2->Errors()[0];

                ROOT::Math::Functor fbin2_p([&](const double* param) { return extractFLU_2(param, phiVector, polVector, epsilonVector, helicityVector);}, 3);
                double B_MIN = -0.8, B_MAX = 0.8;
                double C_MIN = -0.8, C_MAX = 0.8;
                minimum2->SetLimitedVariable(0, "A_LU", 0.01, 0.005, FLU_MIN, FLU_MAX);
                minimum2->SetLimitedVariable(1, "B", 0.01, 0.005, B_MIN, B_MAX);
                minimum2->SetLimitedVariable(2, "C", 0.01, 0.005, C_MIN, C_MAX);
                minimum2->SetFunction(fbin2_p);
                minimum2->Minimize();
                A_LU_PhT_values_2[region][i][j] = minimum2->X()[0];
                A_LU_PhT_errors_2[region][i][j] = minimum2->Errors()[0];
                B_PhT_values_2[region][i][j] = minimum2->X()[1];
                B_PhT_errors_2[region][i][j] = minimum2->Errors()[1];
                C_PhT_values_2[region][i][j] = minimum2->X()[2];
                C_PhT_errors_2[region][i][j] = minimum2->Errors()[2];
            }
        }
    }
    dirAsym->cd();
    //hist_count_xQ2->GetYaxis()->SetRangeUser(0, 3e6);
    //hist_count_xQ2->Draw("HIST");
    //
    TGraphErrors* graphs_PhT[3][nBin_z_2];
    for (int region = 0; region < 3; region++) {
        for (int i = 0; i < nBin_z_2; i++) {
            graphs_PhT[region][i] = new TGraphErrors(nBin_Pht_2);
            for (int j = 0; j < nBin_Pht_2; j++) {
                graphs_PhT[region][i]->SetPoint(j, (PtBins_2[j] + PtBins_2[j+1]) / 2.0, A_LU_PhT_values[region][i][j]);
                graphs_PhT[region][i]->SetPointError(j, 0, A_LU_PhT_errors[region][i][j]);
            }
            graphs_PhT[region][i]->SetTitle(Form("Unbinned A_{LU} vs P_{hT} | z Bin %d | xQ bin %d", i+1, region+1));
            graphs_PhT[region][i]->GetXaxis()->SetTitle("P_{hT}");
            graphs_PhT[region][i]->GetYaxis()->SetTitle("A_LU");
            graphs_PhT[region][i]->SetMarkerStyle(20);
            graphs_PhT[region][i]->SetLineColor(kRed+1); 
            graphs_PhT[region][i]->SetMarkerColor(kRed+1);
            graphs_PhT[region][i]->GetYaxis()->SetRangeUser(0, 0.06);
            graphs_PhT[region][i]->GetYaxis()->SetMaxDigits(3);
            SetStatsBoxGraph(graphs_PhT[region][i], nBin_Pht_2);
        }
    }
    TCanvas* c_PhT[3];
    TCanvas* cline2 = new TCanvas("------------ PhT ------------", "", 800,800);
    cline2->Update();
    cline2->Write();
    for (int region = 0; region < 3; region++) {
        c_PhT[region] = new TCanvas(Form("A_LU_PhT_xQ_bin%d", region+1), Form("Unbinned A_{LU} vs P_{hT} for xQ bin %d", region+1), 800, 800);
        c_PhT[region]->Divide(1, nBin_z_2);
        for (int i = 0; i < nBin_z_2; i++) {
            c_PhT[region]->cd(i+1);
            graphs_PhT[region][i]->Draw("AP");
            //SetStatsBoxGraph(graphs_PhT[region][i], nBin_z);
            gPad->SetGridx();
            gPad->SetGridy();
            gPad->Update();
        }
        c_PhT[region]->Update();
        c_PhT[region]->Write();
    }
    
    // more defined fit
    TGraphErrors* graphs_PhT_2[3][nBin_z_2];
    TGraphErrors* graphs_PhT_Auu1[3][nBin_z_2];
    TGraphErrors* graphs_PhT_Auu2[3][nBin_z_2];
    for (int region = 0; region < 3; region++) {
        for (int i = 0; i < nBin_z_2; i++) {
            graphs_PhT_2[region][i] = new TGraphErrors(nBin_Pht_2);
            graphs_PhT_Auu1[region][i] = new TGraphErrors(nBin_Pht_2);
            graphs_PhT_Auu2[region][i] = new TGraphErrors(nBin_Pht_2);
            for (int j = 0; j < nBin_Pht_2; j++) {
                graphs_PhT_2[region][i]->SetPoint(j, (PtBins_2[j] + PtBins_2[j+1]) / 2.0, A_LU_PhT_values_2[region][i][j]);
                graphs_PhT_2[region][i]->SetPointError(j, 0, A_LU_PhT_errors_2[region][i][j]);
                graphs_PhT_Auu1[region][i]->SetPoint(j, (PtBins_2[j] + PtBins_2[j+1]) / 2.0, B_PhT_values_2[region][i][j]);
                graphs_PhT_Auu1[region][i]->SetPointError(j, 0, B_PhT_errors_2[region][i][j]);
                graphs_PhT_Auu2[region][i]->SetPoint(j, (PtBins_2[j] + PtBins_2[j+1]) / 2.0, C_PhT_values_2[region][i][j]);
                graphs_PhT_Auu2[region][i]->SetPointError(j, 0, C_PhT_errors_2[region][i][j]);
            }
            graphs_PhT_2[region][i]->SetTitle(Form("Unbinned A_{LU} vs P_{hT} | z Bin %d | xQ bin %d", i+1, region+1));
            graphs_PhT_2[region][i]->GetXaxis()->SetTitle("P_{hT}");
            graphs_PhT_2[region][i]->GetYaxis()->SetTitle("A_LU");
            graphs_PhT_2[region][i]->SetMarkerStyle(20);
            graphs_PhT_2[region][i]->SetLineColor(kBlue-2); 
            graphs_PhT_2[region][i]->SetMarkerColor(kBlue-2);
            graphs_PhT_2[region][i]->GetYaxis()->SetRangeUser(0, 0.06);
            graphs_PhT_2[region][i]->GetYaxis()->SetMaxDigits(3);
            // A_UU
            graphs_PhT_Auu1[region][i]->SetTitle(Form("Unbinned A_{UU} vs P_{hT} | z Bin %d | xQ bin %d", i+1, region+1));
            graphs_PhT_Auu1[region][i]->GetXaxis()->SetTitle("P_{hT}");
            graphs_PhT_Auu1[region][i]->GetYaxis()->SetTitle("A_UU");
            graphs_PhT_Auu1[region][i]->SetLineColor(kGreen+3); 
            graphs_PhT_Auu1[region][i]->SetMarkerColor(kGreen+3);
            graphs_PhT_Auu1[region][i]->SetMarkerStyle(20);
            graphs_PhT_Auu1[region][i]->GetYaxis()->SetRangeUser(-0.5, 0.5);
            graphs_PhT_Auu2[region][i]->GetYaxis()->SetRangeUser(-1, 1);
            graphs_PhT_Auu2[region][i]->SetMarkerStyle(20);
            graphs_PhT_Auu2[region][i]->SetLineColor(kMagenta+3); 
            graphs_PhT_Auu2[region][i]->SetMarkerColor(kMagenta+3);
        }
    }
    TCanvas* c_PhT_2[3];
    for (int region = 0; region < 3; region++) {
        c_PhT_2[region] = new TCanvas(Form("A_LU_PhT_xQ_bin%d_NewFit", region+1), Form("Unbinned A_{LU} vs P_{hT} for xQ bin %d", region+1), 800, 800);
        c_PhT_2[region]->Divide(1, nBin_z_2);
        for (int i = 0; i < nBin_z_2; i++) {
            c_PhT_2[region]->cd(i+1);
            graphs_PhT_2[region][i]->Draw("AP");
            graphs_PhT[region][i]->Draw("SAME P");
            //SetStatsBoxGraph(graphs_PhT[region][i], nBin_z);
            TLegend *leg = new TLegend(0.85, 0.65, 0.95, 0.85);
            leg->AddEntry(graphs_PhT[region][i], "A_{LU} fit_1", "ep");
            leg->AddEntry(graphs_PhT_2[region][i], "A_{LU} fit_2", "epf");
            leg->Draw();
            gPad->SetGridx();
            gPad->SetGridy();
            gPad->Update();
        }
        c_PhT_2[region]->Update();
        c_PhT_2[region]->Write();
    }
    // canvas per Auu
    TCanvas* c_PhT_3[3];
    for (int region = 0; region < 3; region++) {
        c_PhT_3[region] = new TCanvas(Form("A_UU_PhT_xQ_bin%d_NewFit", region+1), Form("Unbinned A_{UU} vs P_{hT} for xQ bin %d", region+1), 800, 800);
        c_PhT_3[region]->Divide(1, nBin_z_2);
        for (int i = 0; i < nBin_z_2; i++) {
            c_PhT_3[region]->cd(i+1);
            graphs_PhT_Auu1[region][i]->Draw("AP");
            graphs_PhT_Auu2[region][i]->Draw("SAME P");
            //SetStatsBoxGraph(graphs_PhT[region][i], nBin_z);
            TLegend *leg = new TLegend(0.85, 0.65, 0.95, 0.85);
            leg->AddEntry(graphs_PhT_Auu1[region][i], "A_{UU} cos#Phi", "ep");
            leg->AddEntry(graphs_PhT_Auu2[region][i], "A_{UU} cos2#Phi", "epf");
            leg->Draw();
            gPad->SetGridx();
            gPad->SetGridy();
            gPad->Update();
        }
        c_PhT_3[region]->Update();
        c_PhT_3[region]->Write();
    }
    
    //___________________________________________________________________________________________________________________________________________________
    
    TH1D *hZBins_Hp[3][nBin_z], *hZBins_Hm[3][nBin_z];
    TH1D *hZBins_Phi_Hp[3][nBin_z][nBin_Pht], *hZBins_Phi_Hm[3][nBin_z][nBin_Pht];
    TH1D *hPtBins_Hp[3][nBin_Pht], *hPtBins_Hm[3][nBin_Pht];

    for (int region = 0; region < 3; region++) {
        for (int i = 0; i < nBin_z; i++) {
            hZBins_Hp[region][i] = new TH1D(Form("z_Bin%d_xQ%d_Hp", i+1, region+1), 
                                            Form("z Bin %d | xQ%d | H+;z;Counts", i+1, region+1), 
                                            20, zBins[i], zBins[i+1]);
            for(int j = 0; j < nBin_Pht; j++){
                hZBins_Phi_Hp[region][i][j] = new TH1D(Form("_Phi_zBin%d_PhtBin%d_xQ%d_Hp", i+1, j+1, region+1), 
                                                    Form("#Phi_{h} in z Bin %d | P_{hT} Bin %d | xQ Bin %d | H+; #Phi_{h} [rad];Counts", i+1, j+1, region+1), 
                                                    20, -TMath::Pi(), TMath::Pi());
            }
        }
        
        for (int i = 0; i < nBin_z; i++) {
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
    
    Long64_t nEntries_kPhp = chainKaonPhp.GetEntries();
    for (Long64_t i = 0; i < nEntries_kPhp; i++) {
        chainKaonPhp.GetEntry(i);
        for (int region = 0; region < 3; region++) {
            if (kaonp_xB_hp >= xBins[region][0] && kaonp_xB_hp < xBins[region][1] &&
                kaonp_Q2_hp >= Q2Bins[region][0] && kaonp_Q2_hp < Q2Bins[region][1]) {
                
                for (int j = 0; j < nBin_z; j++) {
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

    Long64_t nEntries_kPhm = chainKaonPhm.GetEntries();
    for (Long64_t i = 0; i < nEntries_kPhm; i++) {
        chainKaonPhm.GetEntry(i);
        for (int region = 0; region < 3; region++) {
            if (kaonp_xB_hm >= xBins[region][0] && kaonp_xB_hm < xBins[region][1] &&
                kaonp_Q2_hm >= Q2Bins[region][0] && kaonp_Q2_hm < Q2Bins[region][1]) {
                for (int j = 0; j < nBin_z; j++) {
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
    
    for (int region = 0; region < 3; region++) {
        for (int i = 0; i < nBin_z; i++) {
            dirKaonp_Hp_xQ[region]->cd();
            hZBins_Hp[region][i]->Write();
            dirKaonp_Hm_xQ[region]->cd(); 
            hZBins_Hm[region][i]->Write();
    
            for (int j = 0; j < nBin_Pht; j++) {
                dirKaonp_Hp_xQ[region]->cd();
                hZBins_Phi_Hp[region][i][j]->Write();
                dirKaonp_Hm_xQ[region]->cd(); 
                hZBins_Phi_Hm[region][i][j]->Write();
            }
        }
    }
    
    for(int region = 0; region < 3; region++){
        dirBin_PhT->cd();
        for(int i = 0; i < nBin_z_2; i++){
            //hist_count_PhT[region][i]->Draw("hist");
            hist_count_PhT[region][i]->Write();
        }
        dirBin_z->cd();
        for(int i = 0; i < nBin_Pht; i++){
            //hist_count_PhT[region][i]->Draw("hist");
            hist_count_z[region][i]->Write();
        }
    }

    //outFile.Write();
    outFile.Close();
    //inFile.Close();

    cout << "ROOT output file: " << outputFile << endl;
}
