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
// to download the data
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rga/torus+1/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rga/output_fall2018_torus+1
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rga/torus-1/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rga/output_fall2018_torus-1
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rga/torus-1/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rga/output_fall2018_torus-1_rich
// 
void _5plot(double torus, bool rich_yes) {
    // Istanza di TDatabasePDG per accedere alle proprietà delle particelle, non va in conflitto con clas12database, sono indipendenti
    auto db2 = TDatabasePDG::Instance();
    // Variables
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
    double kaonp_helicity, kaonp_Phi_Hup, kaonp_Phi_Hdw;

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
    double kaonm_helicity, kaonm_Phi_Hup, kaonm_Phi_Hdw;
    // open file
    //TFile inFile("out_t-1.root", "READ");

    //string inputDir = "output_fall2018_torus+1"; 
    string inputDir;
    if(rich_yes) inputDir = "output_fall2018_torus-1_rich"; 
    else{
        if (torus == -1) {
            inputDir = "output_fall2018_torus-1"; 
        } else if (torus == +1){
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
    // creo un output root 
    //const char* outputFile = "plot_rga_t+1.root";
    const char* outputFile; 
    if(rich_yes)outputFile = "plot_rga_t-1_rich.root";
    else{
        if (torus == -1){
            outputFile = "plot_rga_t-1.root";
        } else if (torus == +1){
            outputFile = "plot_rga_t+1.root";
        }
    }
    TFile outFile(outputFile, "RECREATE");  // File di output ROOT
    // Creiamo un TChain per ogni TTree

    TDirectory* dirKaonp = outFile.mkdir("KaonPlus");
    TDirectory* dirKaonp_Hp = outFile.mkdir("KaonPlus_H+");
    TDirectory* dirKaonp_Hm = outFile.mkdir("KaonPlus_H-");
    //TDirectory* dirKaonp_Bin = outFile.mkdir("KaonPlus_Bin");
    TDirectory* dirKaonm = outFile.mkdir("KaonMinus");
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
    //chainElectron.SetBranchAddress("Nphe", &electron_Nphe);
    //chainElectron.SetBranchAddress("status", &electron_status);
    //chainElectron.SetBranchAddress("PCAL", &electron_PCAL);
    //chainElectron.SetBranchAddress("ECAL", &electron_ECAL);
    //chainElectron.SetBranchAddress("ECIN", &electron_ECIN);
    //chainElectron.SetBranchAddress("DC_6", &electron_edge1);
    //chainElectron.SetBranchAddress("DC_18", &electron_edge2);
    //chainElectron.SetBranchAddress("DC_36", &electron_edge3);
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
    // Kaon -
    chainKaonM.SetBranchAddress("E", &kaonm_E);
    chainKaonM.SetBranchAddress("px", &kaonm_px);
    chainKaonM.SetBranchAddress("py", &kaonm_py);
    chainKaonM.SetBranchAddress("pz", &kaonm_pz);
    chainKaonM.SetBranchAddress("Mom", &kaonm_Ph);
    chainKaonM.SetBranchAddress("W", &kaonm_W);
    chainKaonM.SetBranchAddress("Q2", &kaonm_Q2);
    chainKaonM.SetBranchAddress("xF", &kaonm_xF);
    chainKaonM.SetBranchAddress("xB", &kaonm_xB);
    chainKaonM.SetBranchAddress("y", &kaonm_y);
    chainKaonM.SetBranchAddress("z", &kaonm_z);
    chainKaonM.SetBranchAddress("Pt", &kaonm_Pt);
    chainKaonM.SetBranchAddress("PhT", &kaonm_PhT);
    chainKaonM.SetBranchAddress("Phi_Lab", &kaonm_Phi);
    chainKaonM.SetBranchAddress("Theta_Lab", &kaonm_Theta);
    chainKaonM.SetBranchAddress("Pseudorapidity", &kaonm_eta);
    chainKaonM.SetBranchAddress("Phi_h", &kaonm_Phi_h);
    //chainKaonM.SetBranchAddress("SinPhi", &kaonm_SinPhi);
    chainKaonM.SetBranchAddress("helicity", &kaonm_helicity);

    // KAON+ PLOT
    dirKaonp->cd();
    // Mom
    TH2D kp_MomVsPhT ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  Kaon+  |; P_{hT} [GeV]; Mom [GeV]", 200, 0, 1, 200, 0.9, 7);
    TH2D kp_MomVsXb ("_MomVsXb", "Correlation Mom vs x_{B}  |  Kaon+  |; x_{B}; Mom [GeV]", 200, 0, 0.8, 200, 0.9, 7);
    TH2D kp_MomVsXf ("_MomVsXf", "Correlation Mom vs x_F  |  Kaon+  |; x_F; Mom [GeV]", 200, 0, 0.6, 200, 0.9, 7);
    TH2D kp_MomVsZ ("_MomVsZ", "Correlation Mom vs Z  |  Kaon+  |; z; Mom [GeV]", 200, 0.2, 0.9, 200, 0.9, 7);
    TH2D kp_MomVsY ("_MomVsY", "Correlation Mom vs Y  |  Kaon+  |; y; Mom [GeV]", 200, 0.2, 0.75, 200, 0.9, 7);
    TH2D kp_MomVsEta ("_MomVsEta", "Correlation Mom vs Eta  |  Kaon+  |; Eta; Mom [GeV]", 200, 1.5, 3.0, 200, 0.9, 7);
    TH2D kp_MomVsTheta ("_MomVsTheta", "Correlation Mom vs Theta  |  Kaon+  |; Mom [GeV]; Theta [Rad]", 200, 0.9, 7, 200, 0.09, 0.4);
    TH2D kp_MomVsPhi_h ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  Kaon+  |; #Phi_{h} [Rad]; Mom [GeV]", 200, -TMath::Pi(), TMath::Pi(), 200, 0.9, 7);
    TH2D kp_MomVsM2x ("_MomVsM2x", "correlation Mom vs M2x | Kaon+ |; Mom [GeV]; M2x[GeV^{2}]", 200, 1, 7, 200, 2.45, 11);
    // Q2
    TH2D kp_Q2VsXb ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  Kaon+  |; x_{B}; Q^{2} [GeV^{2}]", 200, 0, 0.8, 200, 1, 8);
    TH2D kp_Q2VsXf ("_Q2VsXf", "Correlation Q^{2} vs x_F  |  Kaon+  |; x_F; Q^{2} [GeV^{2}]", 200, 0, 0.6, 200, 1, 7);
    TH2D kp_Q2VsMom ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  Kaon+  |; Mom [GeV]; Q^{2} [GeV^{2}]", 200, 0.9, 6, 200, 1, 7);
    TH2D kp_Q2VsPhT ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  Kaon+  |; P_{hT} [GeV]; Q^{2} [GeV^{2}]", 200, 0, 1.2, 200, 1, 6);
    TH2D kp_Q2VsZ ("_Q2VsZ", "Correlation Q^{2} vs Z  |  Kaon+  |; z; Q^{2} [GeV^{2}]", 200, 0.2, 0.9, 200, 1, 6);
    TH2D kp_Q2VsY ("_Q2VsY", "Correlation Q^{2} vs Y  |  Kaon+  |; y; Q^{2} [GeV^{2}]", 200, 0.2, 0.75, 200, 0.9, 7);
    TH2D kp_Q2VsEta ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  Kaon+  |; Eta; Q^{2} [GeV^{2}]", 200, 1.5, 3.0, 200, 1, 6);
    TH2D kp_Q2VsPhi_h ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  Kaon+  |; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", 200, -TMath::Pi(), TMath::Pi(), 200, 1, 6);
    // PhT
    TH2D kp_PhTvsZ ("_PhTvsZ", "Correlation P_{hT} vs Z  |  Kaon+  |; z; P_{hT} [GeV]", 200, 0.2, 1, 200, 0, 1.2);
    TH2D kp_PhTvsXb ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  Kaon+  |; x_{B}; P_{hT} [GeV]", 200, 0, 0.8, 200, 0, 1.2);
    TH2D kp_PhTvsEta ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  Kaon+  |; Eta; P_{hT} [GeV]", 200, 1.5, 3.0, 200, 0, 1.2);
    TH2D kp_PhTvsPhi_h ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  Kaon+  |; #Phi_{h} [Rad]; P_{hT} [GeV]", 200, -TMath::Pi(), TMath::Pi(), 200, 0, 1.2);
    // Z
    TH2D kp_zVsXb ("_zVsXb", "Correlation Z vs x_{B}  |  Kaon+  |; x_{B}; z", 200, 0, 0.8, 200, 0.2, 0.9);
    TH2D kp_zVsXf ("_zVsXf", "Correlation Z vs x_F  |  Kaon+  |; x_F; z", 200, 0, 0.6, 200, 0.2, 0.9);
    TH2D kp_zVsEta ("_zVsEta", "Correlation Z vs Eta  |  Kaon+  |; Eta; z", 200, 1.5, 3.0, 200, 0.2, 0.9);
    TH2D kp_zVsPhi_h ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  Kaon+  |; #Phi_{h} [Rad]; z", 200, -TMath::Pi(), TMath::Pi(), 200, 0.2, 0.9);
    //
    TH2D kp_xBvsY ("_xBvsY", "Correlation y vs x_{B}  |  Kaon+  |; x_{B}; y", 200, 0, 0.8, 200, 0.25, 0.75);
    // Angles
    TH2D kp_ThetaVsPhi_h ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  Kaon+  |; #Phi_{h} [Rad]; Theta [Rad]", 200, -TMath::Pi(), TMath::Pi(), 200, 0.1, 0.4);
    //TH2D el_ThetaVsMom ("el_ThetaVsMom", "Correlation Theta vs Mom  |  electron  |; Mom [GeV]; Theta [Deg]", 200, 1.5, 9, 200, 5, 40);
    //TH1D kp_Phi_h ("_Phi_h", "#Phi_{h} distribution | Kaon+ ; #Phi_{h} [Deg]; counts", 200, 0, 360);
    //TH1D kp_Phi_Hp ("_Phi_H+", "#Phi_{h} distribution with positive helicity | Kaon+ ; #Phi_{h} [Deg]; counts", 200, -TMath::Pi(), TMath::Pi());
    //TH1D kp_Phi_Hm ("_Phi_H-", "#Phi_{h} distribution with negative helicity | Kaon+ ; #Phi_{h} [Deg]; counts", 200, -TMath::Pi(), TMath::Pi());

    // KAON+ PLOT Helicity+
    dirKaonp_Hp->cd();
    // Mom
    TH2D kp_MomVsPhT_Hp ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  Kaon+  | Helicity > 0; P_{hT} [GeV]; Mom [GeV]", 200, 0, 1, 200, 0.9, 7);
    TH2D kp_MomVsXb_Hp ("_MomVsXb", "Correlation Mom vs x_{B}  |  Kaon+  | Helicity > 0; x_{B}; Mom [GeV]", 200, 0, 0.8, 200, 0.9, 7);
    TH2D kp_MomVsXf_Hp ("_MomVsXf", "Correlation Mom vs x_F  |  Kaon+  | Helicity > 0; x_F; Mom [GeV]", 200, 0, 0.6, 200, 0.9, 7);
    TH2D kp_MomVsZ_Hp ("_MomVsZ", "Correlation Mom vs Z  |  Kaon+  | Helicity > 0; z; Mom [GeV]", 200, 0.2, 0.9, 200, 0.9, 7);
    TH2D kp_MomVsY_Hp ("_MomVsY", "Correlation Mom vs Y  |  Kaon+  | Helicity > 0; y; Mom [GeV]", 200, 0.2, 0.75, 200, 0.9, 7);
    TH2D kp_MomVsEta_Hp ("_MomVsEta", "Correlation Mom vs Eta  |  Kaon+  | Helicity > 0; Eta; Mom [GeV]", 200, 1.5, 3.0, 200, 0.9, 7);
    TH2D kp_MomVsTheta_Hp ("_MomVsTheta", "Correlation Mom vs Theta  |  Kaon+  | Helicity > 0; Mom [GeV]; Theta [Rad]", 200, 0.9, 7, 200, 0.09, 0.4);
    TH2D kp_MomVsPhi_h_Hp ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  Kaon+  | Helicity > 0; #Phi_{h} [Rad]; Mom [GeV]", 200, -TMath::Pi(), TMath::Pi(), 200, 0.9, 7);
    TH2D kp_MomVsM2x_Hp ("_MomVsM2x", "correlation Mom vs M2x | Kaon+ |; Mom [GeV]; M2x[GeV^{2}]", 200, 1, 7, 200, 2.45, 11);
    // Q2
    TH2D kp_Q2VsXb_Hp ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  Kaon+  | Helicity > 0; x_{B}; Q^{2} [GeV^{2}]", 200, 0, 0.8, 200, 1, 8);
    TH2D kp_Q2VsXf_Hp ("_Q2VsXf", "Correlation Q^{2} vs x_F  |  Kaon+  | Helicity > 0; x_F; Q^{2} [GeV^{2}]", 200, 0, 0.6, 200, 1, 7);
    TH2D kp_Q2VsMom_Hp ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  Kaon+  | Helicity > 0; Mom [GeV]; Q^{2} [GeV^{2}]", 200, 0.9, 6, 200, 1, 7);
    TH2D kp_Q2VsPhT_Hp ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  Kaon+  | Helicity > 0; P_{hT} [GeV]; Q^{2} [GeV^{2}]", 200, 0, 1.2, 200, 1, 6);
    TH2D kp_Q2VsZ_Hp ("_Q2VsZ", "Correlation Q^{2} vs Z  |  Kaon+  | Helicity > 0; z; Q^{2} [GeV^{2}]", 200, 0.2, 0.9, 200, 1, 6);
    TH2D kp_Q2VsY_Hp ("_Q2VsY", "Correlation Q^{2} vs Y  |  Kaon+  | Helicity > 0; y; Q^{2} [GeV^{2}]", 200, 0.2, 0.75, 200, 0.9, 7);
    TH2D kp_Q2VsEta_Hp ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  Kaon+  | Helicity > 0; Eta; Q^{2} [GeV^{2}]", 200, 1.5, 3.0, 200, 1, 6);
    TH2D kp_Q2VsPhi_h_Hp ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  Kaon+  | Helicity > 0; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", 200, -TMath::Pi(), TMath::Pi(), 200, 1, 6);
    // PhT
    TH2D kp_PhTvsZ_Hp ("_PhTvsZ", "Correlation P_{hT} vs Z  |  Kaon+  | Helicity > 0; z; P_{hT} [GeV]", 200, 0.2, 1, 200, 0, 1.2);
    TH2D kp_PhTvsXb_Hp ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  Kaon+  | Helicity > 0; x_{B}; P_{hT} [GeV]", 200, 0, 0.8, 200, 0, 1.2);
    TH2D kp_PhTvsEta_Hp ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  Kaon+  | Helicity > 0; Eta; P_{hT} [GeV]", 200, 1.5, 3.0, 200, 0, 1.2);
    TH2D kp_PhTvsPhi_h_Hp ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  Kaon+  | Helicity > 0; #Phi_{h} [Rad]; P_{hT} [GeV]", 200, -TMath::Pi(), TMath::Pi(), 200, 0, 1.2);
    // Z
    TH2D kp_zVsXb_Hp ("_zVsXb", "Correlation Z vs x_{B}  |  Kaon+  | Helicity > 0; x_{B}; z", 200, 0, 0.8, 200, 0.2, 0.9);
    TH2D kp_zVsXf_Hp ("_zVsXf", "Correlation Z vs x_F  |  Kaon+  | Helicity > 0; x_F; z", 200, 0, 0.6, 200, 0.2, 0.9);
    TH2D kp_zVsEta_Hp ("_zVsEta", "Correlation Z vs Eta  |  Kaon+  | Helicity > 0; Eta; z", 200, 1.5, 3.0, 200, 0.2, 0.9);
    TH2D kp_zVsPhi_h_Hp ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  Kaon+  | Helicity > 0; #Phi_{h} [Rad]; z", 200, -TMath::Pi(), TMath::Pi(), 200, 0.2, 0.9);
    //
    TH2D kp_xBvsY_Hp ("_xBvsY", "Correlation y vs x_{B}  |  Kaon+  | Helicity > 0; x_{B}; y", 200, 0, 0.8, 200, 0.25, 0.75);
    //
    TH2D kp_ThetaVsPhi_h_Hp ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  Kaon+  | Helicity > 0; #Phi_{h} [Rad]; Theta [Rad]", 200, -TMath::Pi(), TMath::Pi(), 200, 0.1, 0.4);
    

    // KAON+ PLOT Helicity-
    dirKaonp_Hm->cd();
    // Mom
    TH2D kp_MomVsPhT_Hm ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  Kaon+  | Helicity < 0; P_{hT} [GeV]; Mom [GeV]", 200, 0, 1, 200, 0.9, 7);
    TH2D kp_MomVsXb_Hm ("_MomVsXb", "Correlation Mom vs x_{B}  |  Kaon+  | Helicity < 0; x_{B}; Mom [GeV]", 200, 0, 0.8, 200, 0.9, 7);
    TH2D kp_MomVsXf_Hm ("_MomVsXf", "Correlation Mom vs x_F  |  Kaon+  | Helicity < 0; x_F; Mom [GeV]", 200, 0, 0.6, 200, 0.9, 7);
    TH2D kp_MomVsZ_Hm ("_MomVsZ", "Correlation Mom vs Z  |  Kaon+  | Helicity < 0; z; Mom [GeV]", 200, 0.2, 0.9, 200, 0.9, 7);
    TH2D kp_MomVsY_Hm ("_MomVsY", "Correlation Mom vs Y  |  Kaon+  | Helicity < 0; y; Mom [GeV]", 200, 0.2, 0.75, 200, 0.9, 7);
    TH2D kp_MomVsEta_Hm ("_MomVsEta", "Correlation Mom vs Eta  |  Kaon+  | Helicity < 0; Eta; Mom [GeV]", 200, 1.5, 3.0, 200, 0.9, 7);
    TH2D kp_MomVsTheta_Hm ("_MomVsTheta", "Correlation Mom vs Theta  |  Kaon+  | Helicity < 0; Mom [GeV]; Theta [Rad]", 200, 0.9, 7, 200, 0.09, 0.4);
    TH2D kp_MomVsPhi_h_Hm ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  Kaon+  | Helicity < 0; #Phi_{h} [Rad]; Mom [GeV]", 200, -TMath::Pi(), TMath::Pi(), 200, 0.9, 7);
    TH2D kp_MomVsM2x_Hm ("_MomVsM2x", "correlation Mom vs M2x | Kaon+ |; Mom [GeV]; M2x[GeV^{2}]", 200, 1, 7, 200, 2.45, 11);
    // Q2
    TH2D kp_Q2VsXb_Hm ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  Kaon+  | Helicity < 0; x_{B}; Q^{2} [GeV^{2}]", 200, 0, 0.8, 200, 1, 8);
    TH2D kp_Q2VsXf_Hm ("_Q2VsXf", "Correlation Q^{2} vs x_F  |  Kaon+  | Helicity < 0; x_F; Q^{2} [GeV^{2}]", 200, 0, 0.6, 200, 1, 7);
    TH2D kp_Q2VsMom_Hm ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  Kaon+  | Helicity < 0; Mom [GeV]; Q^{2} [GeV^{2}]", 200, 0.9, 6, 200, 1, 7);
    TH2D kp_Q2VsPhT_Hm ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  Kaon+  | Helicity < 0; P_{hT} [GeV]; Q^{2} [GeV^{2}]", 200, 0, 1.2, 200, 1, 6);
    TH2D kp_Q2VsZ_Hm ("_Q2VsZ", "Correlation Q^{2} vs Z  |  Kaon+  | Helicity < 0; z; Q^{2} [GeV^{2}]", 200, 0.2, 0.9, 200, 1, 6);
    TH2D kp_Q2VsY_Hm ("_Q2VsY", "Correlation Q^{2} vs Y  |  Kaon+  | Helicity < 0; y; Q^{2} [GeV^{2}]", 200, 0.2, 0.75, 200, 0.9, 7);
    TH2D kp_Q2VsEta_Hm ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  Kaon+  | Helicity < 0; Eta; Q^{2} [GeV^{2}]", 200, 1.5, 3.0, 200, 1, 6);
    TH2D kp_Q2VsPhi_h_Hm ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  Kaon+  | Helicity < 0; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", 200, -TMath::Pi(), TMath::Pi(), 200, 1, 6);
    // PhT
    TH2D kp_PhTvsZ_Hm ("_PhTvsZ", "Correlation P_{hT} vs Z  |  Kaon+  | Helicity < 0; z; P_{hT} [GeV]", 200, 0.2, 1, 200, 0, 1.2);
    TH2D kp_PhTvsXb_Hm ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  Kaon+  | Helicity < 0; x_{B}; P_{hT} [GeV]", 200, 0, 0.8, 200, 0, 1.2);
    TH2D kp_PhTvsEta_Hm ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  Kaon+  | Helicity < 0; Eta; P_{hT} [GeV]", 200, 1.5, 3.0, 200, 0, 1.2);
    TH2D kp_PhTvsPhi_h_Hm ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  Kaon+  | Helicity < 0; #Phi_{h} [Rad]; P_{hT} [GeV]", 200, -TMath::Pi(), TMath::Pi(), 200, 0, 1.2);
    // Z
    TH2D kp_zVsXb_Hm ("_zVsXb", "Correlation Z vs x_{B}  |  Kaon+  | Helicity < 0; x_{B}; z", 200, 0, 0.8, 200, 0.2, 0.9);
    TH2D kp_zVsXf_Hm ("_zVsXf", "Correlation Z vs x_F  |  Kaon+  | Helicity < 0; x_F; z", 200, 0, 0.6, 200, 0.2, 0.9);
    TH2D kp_zVsEta_Hm ("_zVsEta", "Correlation Z vs Eta  |  Kaon+  | Helicity < 0; Eta; z", 200, 1.5, 3.0, 200, 0.2, 0.9);
    TH2D kp_zVsPhi_h_Hm ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  Kaon+  | Helicity < 0; #Phi_{h} [Rad]; z", 200, -TMath::Pi(), TMath::Pi(), 200, 0.2, 0.9);
    //
    TH2D kp_xBvsY_Hm ("_xBvsY", "Correlation y vs x_{B}  |  Kaon+  | Helicity < 0; x_{B}; y", 200, 0, 0.8, 200, 0.25, 0.75);
    //
    TH2D kp_ThetaVsPhi_h_Hm ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  Kaon+  | Helicity < 0; #Phi_{h} [Rad]; Theta [Rad]", 200, -TMath::Pi(), TMath::Pi(), 200, 0.1, 0.4);
    
    // KAON- PLOT
    dirKaonm->cd();
    // Mom
    TH2D km_MomVsPhT ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  Kaon-  |; P_{hT} [GeV]; Mom[GeV]", 200, 0, 1, 200, 0.9, 6);
    TH2D km_MomVsXb ("_MomVsXb", "Correlation Mom vs x_{B}  |  Kaon-  |; x_{B}; Mom[GeV]", 200, 0, 0.8, 200, 0.9, 6);
    TH2D km_MomVsXf ("_MomVsXf", "Correlation Mom vs x_F  |  Kaon-  |; x_F; Mom[GeV]", 200, 0, 0.6, 200, 0.9, 6);
    TH2D km_MomVsZ ("_MomVsZ", "Correlation Mom vs Z  |  Kaon-  |; z; Mom[GeV]", 200, 0.2, 0.9, 200, 0.9, 6);
    TH2D km_MomVsY ("_MomVsY", "Correlation Mom vs Y  |  Kaon-  |; y; Mom[GeV]", 200, 0.2, 0.75, 200, 0.9, 6);
    TH2D km_MomVsEta ("_MomVsEta", "Correlation Mom vs Eta  |  Kaon-  |; Eta; Mom[GeV]", 200, 1.5, 3, 200, 0.9, 6);
    TH2D km_MomVsTheta ("_MomVsTheta", "Correlation Mom vs Theta  |  Kaon-  |; Mom [GeV]; Theta [Rad]", 200, 0.9, 6, 200, 0.09, 0.4);
    TH2D km_MomVsPhi_h ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  Kaon-  |; #Phi_{h} [Rad]; Mom[GeV]", 200, -TMath::Pi(), TMath::Pi(), 200, 0.9, 6);
    // Q2
    TH2D km_Q2VsXb ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  Kaon-  |; x_{B}; Q^{2} [GeV^{2}]", 200, 0, 0.8, 200, 1, 6);
    TH2D km_Q2VsXf ("_Q2VsXf", "Correlation Q^{2} vs x_F  |  Kaon-  |; x_F; Q^{2} [GeV^{2}]", 200, 0, 0.6, 200, 1, 6);
    TH2D km_Q2VsMom ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  Kaon-  |; Mom [GeV]; Q^{2} [GeV^{2}]", 200, 0.9, 6, 200, 1, 6);
    TH2D km_Q2VsPhT ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  Kaon-  |; P_{hT} [GeV]; Q^{2} [GeV^{2}]", 200, 0, 1.1, 200, 1, 6);
    TH2D km_Q2VsZ ("_Q2VsZ", "Correlation Q^{2} vs Z  |  Kaon-  |; z; Q^{2} [GeV^{2}]", 200, 0.2, 0.9, 200, 1, 6);
    TH2D km_Q2VsY ("_Q2VsY", "Correlation Q^{2} vs Y  |  Kaon-  |; y; Q^{2} [GeV^{2}]", 200, 0.2, 0.75, 200, 0.9, 6);
    TH2D km_Q2VsEta ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  Kaon-  |; Eta; Q^{2} [GeV^{2}]", 200, 1.5, 3, 200, 1, 6);
    TH2D km_Q2VsPhi_h ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  Kaon-  |; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", 200, -TMath::Pi(), TMath::Pi(), 200, 1, 6);
    // PhT
    TH2D km_PhTvsZ ("_PhTvsZ", "Correlation P_{hT} vs Z  |  Kaon-  |; z; P_{hT} [GeV]", 200, 0.2, 1, 200, 0, 1.1);
    TH2D km_PhTvsXb ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  Kaon-  |; x_{B}; P_{hT} [GeV]", 200, 0, 0.8, 200, 0, 1.1);
    TH2D km_PhTvsEta ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  Kaon-  |; Eta; P_{hT} [GeV]", 200, 1.5, 3, 200, 0, 1.1);
    TH2D km_PhTvsPhi_h ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  Kaon-  |; #Phi_{h} [Rad]; P_{hT} [GeV]", 200, -TMath::Pi(), TMath::Pi(), 200, 0, 1.1);
    // Z
    TH2D km_zVsXb ("_zVsXb", "Correlation Z vs x_{B}  |  Kaon-  |; x_{B}; z", 200, 0, 0.8, 200, 0.2, 0.9);
    TH2D km_zVsXf ("_zVsXf", "Correlation Z vs x_F  |  Kaon-  |; x_F; z", 200, 0, 0.6, 200, 0.2, 0.9);
    TH2D km_zVsEta ("_zVsEta", "Correlation Z vs Eta  |  Kaon-  |; Eta; z", 200, 1.5, 3, 200, 0.2, 0.9);
    TH2D km_zVsPhi_h ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  Kaon-  |; #Phi_{h} [Rad]; z", 200, -TMath::Pi(), TMath::Pi(), 200, 0.2, 0.9);
    //
    TH2D km_xBvsY ("_xBvsY", "Correlation y vs x_{B}  |  Kaon-  |; x_{B}; y", 200, 0, 0.8, 200, 0.25, 0.75);
    // Angles
    TH2D km_ThetaVsPhi_h ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  Kaon-  |; #Phi_{h} [Rad]; Theta [Rad]", 200, -TMath::Pi(), TMath::Pi(), 200, 0.1, 0.4);
    //TH1D km_Phi_Hp ("_Phi_H+", "#Phi_{h} distribution with positive helicity | Kaon- ; #Phi_{h} [Deg]; counts", 200, -TMath::Pi(), TMath::Pi());
    //TH1D km_Phi_Hm ("_Phi_H-", "#Phi_{h} distribution with negative helicity | Kaon- ; #Phi_{h} [Deg]; counts", 200, -TMath::Pi(), TMath::Pi());

    // ORA RIEMPI I GRAFICI
    // Kaon+
    Long64_t nEntries_kp = chainKaonP.GetEntries();
    for (Long64_t i = 0; i < nEntries_kp; i++) {
        // Mom
        chainKaonP.GetEntry(i);
        kp_MomVsPhT.Fill(kaonp_PhT, kaonp_Ph);
        kp_MomVsEta.Fill(kaonp_eta, kaonp_Ph);
        kp_MomVsM2x.Fill(kaonp_Ph, kaonp_M2x);
        kp_MomVsPhi_h.Fill(kaonp_Phi_h, kaonp_Ph);
        kp_MomVsTheta.Fill(kaonp_Ph, kaonp_Theta);
        kp_MomVsXb.Fill(kaonp_xB, kaonp_Ph);
        kp_MomVsXf.Fill(kaonp_xF, kaonp_Ph);
        kp_MomVsY.Fill(kaonp_y, kaonp_Ph);
        kp_MomVsZ.Fill(kaonp_z, kaonp_Ph);
        // Q2
        kp_Q2VsEta.Fill(kaonp_eta, kaonp_Q2);
        kp_Q2VsMom.Fill(kaonp_Ph, kaonp_Q2);
        kp_Q2VsPhi_h.Fill(kaonp_Phi_h, kaonp_Q2);
        kp_Q2VsPhT.Fill(kaonp_PhT, kaonp_Q2);
        kp_Q2VsXb.Fill(kaonp_xB, kaonp_Q2);
        kp_Q2VsXf.Fill(kaonp_xF, kaonp_Q2);
        kp_Q2VsY.Fill(kaonp_y, kaonp_Q2);
        kp_Q2VsZ.Fill(kaonp_z, kaonp_Q2);
        // PhT
        kp_PhTvsEta.Fill(kaonp_eta, kaonp_PhT);
        kp_PhTvsXb.Fill(kaonp_xB, kaonp_PhT);
        kp_PhTvsZ.Fill(kaonp_z, kaonp_PhT);
        kp_PhTvsPhi_h.Fill(kaonp_Phi_h, kaonp_PhT);
        // z
        kp_zVsEta.Fill(kaonp_eta, kaonp_z);
        kp_zVsPhi_h.Fill(kaonp_Phi_h, kaonp_z);
        kp_zVsXb.Fill(kaonp_xB, kaonp_z);
        kp_zVsXf.Fill(kaonp_xF, kaonp_z);
        //
        kp_xBvsY.Fill(kaonp_xB, kaonp_y);
        kp_ThetaVsPhi_h.Fill(kaonp_Phi_h, kaonp_Theta);
    }

    // Kaon+ H+
    Long64_t nEntries_kp_Hp = chainKaonPhp.GetEntries();
    for (Long64_t i = 0; i < nEntries_kp_Hp; i++) {
        // Mom
        chainKaonPhp.GetEntry(i);
        kp_MomVsPhT_Hp.Fill(kaonp_PhT_hp, kaonp_Ph_hp);
        kp_MomVsEta_Hp.Fill(kaonp_eta_hp, kaonp_Ph_hp);
        kp_MomVsM2x_Hp.Fill(kaonp_Ph_hp, kaonp_M2x_hp);
        kp_MomVsPhi_h_Hp.Fill(kaonp_Phi_h_hp, kaonp_Ph_hp);
        kp_MomVsTheta_Hp.Fill(kaonp_Ph_hp, kaonp_Theta_hp);
        kp_MomVsXb_Hp.Fill(kaonp_xB_hp, kaonp_Ph_hp);
        kp_MomVsXf_Hp.Fill(kaonp_xF_hp, kaonp_Ph_hp);
        kp_MomVsY_Hp.Fill(kaonp_y_hp, kaonp_Ph_hp);
        kp_MomVsZ_Hp.Fill(kaonp_z_hp, kaonp_Ph_hp);
        // Q2
        kp_Q2VsEta_Hp.Fill(kaonp_eta_hp, kaonp_Q2_hp);
        kp_Q2VsMom_Hp.Fill(kaonp_Ph_hp, kaonp_Q2_hp);
        kp_Q2VsPhi_h_Hp.Fill(kaonp_Phi_h_hp, kaonp_Q2_hp);
        kp_Q2VsPhT_Hp.Fill(kaonp_PhT_hp, kaonp_Q2_hp);
        kp_Q2VsXb_Hp.Fill(kaonp_xB_hp, kaonp_Q2_hp);
        kp_Q2VsXf_Hp.Fill(kaonp_xF_hp, kaonp_Q2_hp);
        kp_Q2VsY_Hp.Fill(kaonp_y_hp, kaonp_Q2_hp);
        kp_Q2VsZ_Hp.Fill(kaonp_z_hp, kaonp_Q2_hp);
        // PhT
        kp_PhTvsEta_Hp.Fill(kaonp_eta_hp, kaonp_PhT_hp);
        kp_PhTvsXb_Hp.Fill(kaonp_xB_hp, kaonp_PhT_hp);
        kp_PhTvsZ_Hp.Fill(kaonp_z_hp, kaonp_PhT_hp);
        kp_PhTvsPhi_h_Hp.Fill(kaonp_Phi_h_hp, kaonp_PhT_hp);
        // z
        kp_zVsEta_Hp.Fill(kaonp_eta_hp, kaonp_z_hp);
        kp_zVsPhi_h_Hp.Fill(kaonp_Phi_h_hp, kaonp_z_hp);
        kp_zVsXb_Hp.Fill(kaonp_xB_hp, kaonp_z_hp);
        kp_zVsXf_Hp.Fill(kaonp_xF_hp, kaonp_z_hp);
        //
        kp_xBvsY_Hp.Fill(kaonp_xB_hp, kaonp_y_hp);
        kp_ThetaVsPhi_h_Hp.Fill(kaonp_Phi_h_hp, kaonp_Theta_hp);
    }

    // Kaon+ H-
    Long64_t nEntries_kp_Hm = chainKaonPhm.GetEntries();
    for (Long64_t i = 0; i < nEntries_kp_Hm; i++) {
        // Mom
        chainKaonPhm.GetEntry(i);
        kp_MomVsPhT_Hm.Fill(kaonp_PhT_hm, kaonp_Ph_hm);
        kp_MomVsEta_Hm.Fill(kaonp_eta_hm, kaonp_Ph_hm);
        kp_MomVsM2x_Hm.Fill(kaonp_Ph_hm, kaonp_M2x_hm);
        kp_MomVsPhi_h_Hm.Fill(kaonp_Phi_h_hm, kaonp_Ph_hm);
        kp_MomVsTheta_Hm.Fill(kaonp_Ph_hm, kaonp_Theta_hm);
        kp_MomVsXb_Hm.Fill(kaonp_xB_hm, kaonp_Ph_hm);
        kp_MomVsXf_Hm.Fill(kaonp_xF_hm, kaonp_Ph_hm);
        kp_MomVsY_Hm.Fill(kaonp_y_hm, kaonp_Ph_hm);
        kp_MomVsZ_Hm.Fill(kaonp_z_hm, kaonp_Ph_hm);
        // Q2
        kp_Q2VsEta_Hm.Fill(kaonp_eta_hm, kaonp_Q2_hm);
        kp_Q2VsMom_Hm.Fill(kaonp_Ph_hm, kaonp_Q2_hm);
        kp_Q2VsPhi_h_Hm.Fill(kaonp_Phi_h_hm, kaonp_Q2_hm);
        kp_Q2VsPhT_Hm.Fill(kaonp_PhT_hm, kaonp_Q2_hm);
        kp_Q2VsXb_Hm.Fill(kaonp_xB_hm, kaonp_Q2_hm);
        kp_Q2VsXf_Hm.Fill(kaonp_xF_hm, kaonp_Q2_hm);
        kp_Q2VsY_Hm.Fill(kaonp_y_hm, kaonp_Q2_hm);
        kp_Q2VsZ_Hm.Fill(kaonp_z_hm, kaonp_Q2_hm);
        // PhT
        kp_PhTvsEta_Hm.Fill(kaonp_eta_hm, kaonp_PhT_hm);
        kp_PhTvsXb_Hm.Fill(kaonp_xB_hm, kaonp_PhT_hm);
        kp_PhTvsZ_Hm.Fill(kaonp_z_hm, kaonp_PhT_hm);
        kp_PhTvsPhi_h_Hm.Fill(kaonp_Phi_h_hm, kaonp_PhT_hm);
        // z
        kp_zVsEta_Hm.Fill(kaonp_eta_hm, kaonp_z_hm);
        kp_zVsPhi_h_Hm.Fill(kaonp_Phi_h_hm, kaonp_z_hm);
        kp_zVsXb_Hm.Fill(kaonp_xB_hm, kaonp_z_hm);
        kp_zVsXf_Hm.Fill(kaonp_xF_hm, kaonp_z_hm);
        //
        kp_xBvsY_Hm.Fill(kaonp_xB_hm, kaonp_y_hm);
        kp_ThetaVsPhi_h_Hm.Fill(kaonp_Phi_h_hm, kaonp_Theta_hm);
    }

    // Kaon-
    Long64_t nEntries_km = chainKaonM.GetEntries();
    for (Long64_t i = 0; i < nEntries_km; i++) {
        chainKaonM.GetEntry(i);
        km_MomVsPhT.Fill(kaonm_PhT, kaonm_Ph);
        km_MomVsEta.Fill(kaonm_eta, kaonm_Ph);
        //km_MomVsM2x.Fill(kaonm_M2x, kaonm_Ph);
        km_MomVsPhi_h.Fill(kaonm_Phi_h, kaonm_Ph);
        km_MomVsTheta.Fill(kaonm_Ph, kaonm_Theta);
        km_MomVsXb.Fill(kaonm_xB, kaonm_Ph);
        km_MomVsXf.Fill(kaonm_xF, kaonm_Ph);
        km_MomVsY.Fill(kaonm_y, kaonm_Ph);
        km_MomVsZ.Fill(kaonm_z, kaonm_Ph);
        // Q2
        km_Q2VsEta.Fill(kaonm_eta, kaonm_Q2);
        km_Q2VsMom.Fill(kaonm_Ph, kaonm_Q2);
        km_Q2VsPhi_h.Fill(kaonm_Phi_h, kaonm_Q2);
        km_Q2VsPhT.Fill(kaonm_PhT, kaonm_Q2);
        km_Q2VsXb.Fill(kaonm_xB, kaonm_Q2);
        km_Q2VsXf.Fill(kaonm_xF, kaonm_Q2);
        km_Q2VsY.Fill(kaonm_y, kaonm_Q2);
        km_Q2VsZ.Fill(kaonm_z, kaonm_Q2);
        // PhT
        km_PhTvsEta.Fill(kaonm_eta, kaonm_PhT);
        km_PhTvsXb.Fill(kaonm_xB, kaonm_PhT);
        km_PhTvsZ.Fill(kaonm_z, kaonm_PhT);
        km_PhTvsPhi_h.Fill(kaonm_Phi_h, kaonm_PhT);
        // z
        km_zVsEta.Fill(kaonm_eta, kaonm_z);
        km_zVsPhi_h.Fill(kaonm_Phi_h, kaonm_z);
        km_zVsXb.Fill(kaonm_xB, kaonm_z);
        km_zVsXf.Fill(kaonm_xF, kaonm_z);
        //
        km_xBvsY.Fill(kaonm_xB, kaonm_y);
        km_ThetaVsPhi_h.Fill(kaonm_Phi_h, kaonm_Theta);
    }

    // Binning plot
    dirKaonp->cd();
    // Q2 vs xB
    TCanvas *c1 = new TCanvas("Q2_vs_xB_Bin", "Q^{2} vs x_{B} bin", 800, 600);
    kp_Q2VsXb.Draw("COLZ");
    std::vector<double> xB1, Q21, xB2, Q22, xB3, Q23;
    if(torus == -1){
        xB1= {0.05, 0.18, 0.18, 0.05, 0.05};  
        Q21= {1.0, 1.0, 2.715, 2.715, 1.0};
        xB2= {0.18, 0.45, 0.45, 0.18, 0.18};  
        Q22= {1.25, 1.25, 3, 3, 1.25};;
        xB3= {0.18, 0.7, 0.7, 0.18, 0.18};  
        Q23= {3, 3, 8, 8, 3};
    } else if(torus == +1){
        xB1= {0.05, 0.13, 0.13, 0.05, 0.05};  
        Q21= {1.0, 1.0, 2.2, 2.2, 1.0};
        xB2= {0.13, 0.4, 0.4, 0.13, 0.13};  
        Q22= {1, 1, 2.4, 2.4, 1};;
        xB3= {0.15, 0.7, 0.7, 0.15, 0.15};  
        Q23= {2.4, 2.4, 8, 8, 2.4};
    }
    TPolyLine *rect1 = new TPolyLine(5, xB1.data(), Q21.data());
    //rect1->SetLineWidth(2);
    rect1->Draw("same");
    TPolyLine *rect2 = new TPolyLine(5, xB2.data(), Q22.data());
    //rect2->SetLineWidth(2);
    rect2->Draw("same");
    TPolyLine *rect3 = new TPolyLine(5, xB3.data(), Q23.data());
    //rect3->SetLineWidth(2);
    rect3->Draw("same");
    c1->Update();
    c1->Write();
    //
    // z vs PhT
    TCanvas *c2 = new TCanvas("z_vs_PhT_Bin_z", "z vs P_{hT} bin", 800, 600);
    kp_PhTvsZ.Draw("COLZ");
    std::vector<TPolyLine*> gridLines;
    // Definizione dei bin
    double zBins_1[7] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9};
    double pBins_1[4];  // Dichiarato fuori per essere visibile ovunque

    if (rich_yes) {
        double tempBins[] = {0.0, 0.25, 0.4, 1.2};
        std::copy(std::begin(tempBins), std::end(tempBins), std::begin(pBins_1));
    } else {
        if (torus == -1) {
            double tempBins[] = {0.0, 0.25, 0.5, 1.2};
            std::copy(std::begin(tempBins), std::end(tempBins), std::begin(pBins_1));
        } else if (torus == +1) {
            double tempBins[] = {0.0, 0.35, 0.55, 1.2};
            std::copy(std::begin(tempBins), std::end(tempBins), std::begin(pBins_1));
        }
    }

    // Disegna il bordo del rettangolo (area totale del grafico)
    double z_border[] = {zBins_1[0], zBins_1[6], zBins_1[6], zBins_1[0], zBins_1[0]};
    double p_border[] = {pBins_1[0], pBins_1[0], pBins_1[3], pBins_1[3], pBins_1[0]};
    gridLines.push_back(new TPolyLine(5, z_border, p_border));

    // Linee verticali (divisori di z)
    for (size_t i = 1; i < 6; i++) {  // Escludiamo i bordi estremi
        double z_temp[] = {zBins_1[i], zBins_1[i]};
        double p_temp[] = {pBins_1[0], pBins_1[3]};
        gridLines.push_back(new TPolyLine(2, z_temp, p_temp));
    }

    // Linee orizzontali (divisori di p)
    for (size_t i = 1; i < 3; i++) {  // Escludiamo i bordi estremi
        double z_temp[] = {zBins_1[0], zBins_1[6]};
        double p_temp[] = {pBins_1[i], pBins_1[i]};
        gridLines.push_back(new TPolyLine(2, z_temp, p_temp));
    }

    // Disegna tutte le linee sulla canvas
    for (auto line : gridLines) {
        line->Draw("same");
    }
    // Disegna tutte le linee sulla canvas
    for (auto line : gridLines) {
        line->Draw("same");
    }
    c2->Update();
    c2->Write();
    //
    // PhT vs z
    double PtBins_2[10] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1};
    double zBins_2[4];
    if(rich_yes){ 
        double tempBins[] = {0.2, 0.35, 0.45, 0.9};
        std::copy(std::begin(tempBins), std::end(tempBins), std::begin(zBins_2));
    } else {
        double tempBins[] = {0.2, 0.4, 0.55, 0.9};
        std::copy(std::begin(tempBins), std::end(tempBins), std::begin(zBins_2));
    }
    TCanvas *c3= new TCanvas("z_vs_PhT_Bin_PhT", "z vs P_{hT} bin", 800, 600);
    kp_PhTvsZ.Draw("COLZ");  // Disegna l'istogramma di base
    std::vector<TPolyLine*> gridLines2;
    // Disegna il bordo del rettangolo (area totale del grafico)
    double z_border_2[] = {zBins_2[0], zBins_2[3], zBins_2[3], zBins_2[0], zBins_2[0]};
    double pt_border_2[] = {PtBins_2[0], PtBins_2[0], PtBins_2[9], PtBins_2[9], PtBins_2[0]};
    gridLines2.push_back(new TPolyLine(5, z_border_2, pt_border_2));
    // Linee verticali (divisori di z)
    for (int i = 1; i < 3; i++) {  // Escludiamo i bordi estremi
        double z_temp[] = {zBins_2[i], zBins_2[i]};
        double pt_temp[] = {PtBins_2[0], PtBins_2[9]};
        gridLines2.push_back(new TPolyLine(2, z_temp, pt_temp));
    }
    // Linee orizzontali (divisori di Pt)
    for (int i = 1; i < 9; i++) {  // Escludiamo i bordi estremi
        double z_temp[] = {zBins_2[0], zBins_2[3]};
        double pt_temp[] = {PtBins_2[i], PtBins_2[i]};
        gridLines2.push_back(new TPolyLine(2, z_temp, pt_temp));
    }
    // Disegna tutte le linee sulla canvas
    for (auto line : gridLines2) {
        line->Draw("same");
    }
    c3->Update();
    c3->Write();

    outFile.Write();
    outFile.Close();
    //chain.Close();

    cout << "ROOT output file: " << outputFile << endl;
}
