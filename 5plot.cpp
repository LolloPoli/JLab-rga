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
#include "clas12reader.h"
#include <iostream>
#include "HipoChain.h"
#include "ccdb_reader.h"
#include "rcdb_reader.h"
#include "rich.h"
#include "particle.h"
#include "particle_detector.h"

using namespace clas12;
using namespace std;


int main() {
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
    TFile inFile("out_t-1.root", "READ");
    // creo un output root 
    const char* outputFile = "plot_t-1.root";
    TFile outFile(outputFile, "RECREATE");  // File di output ROOT
    TTree *tree = (TTree*)inFile.Get("Tree");
    TTree *treeEl = (TTree*)inFile.Get("Electron");
    TTree *treeKaonP = (TTree*)inFile.Get("Kaon+");
    TTree *treeKaonPhp = (TTree*)inFile.Get("Kaon+ H+");
    TTree *treeKaonPhm = (TTree*)inFile.Get("Kaon+ H-");
    TTree *treeKaonM = (TTree*)inFile.Get("Kaon-");
    TDirectory* dirKaonp = outFile.mkdir("KaonPlus");
    TDirectory* dirKaonp_Hp = outFile.mkdir("KaonPlus_H+");
    TDirectory* dirKaonp_Hm = outFile.mkdir("KaonPlus_H-");
    //TDirectory* dirKaonp_Bin = outFile.mkdir("KaonPlus_Bin");
    TDirectory* dirKaonm = outFile.mkdir("KaonMinus");
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
    // Kaon -
    treeKaonM->SetBranchAddress("E", &kaonm_E);
    treeKaonM->SetBranchAddress("px", &kaonm_px);
    treeKaonM->SetBranchAddress("py", &kaonm_py);
    treeKaonM->SetBranchAddress("pz", &kaonm_pz);
    treeKaonM->SetBranchAddress("Mom", &kaonm_Ph);
    treeKaonM->SetBranchAddress("W", &kaonm_W);
    treeKaonM->SetBranchAddress("Q2", &kaonm_Q2);
    treeKaonM->SetBranchAddress("xF", &kaonm_xF);
    treeKaonM->SetBranchAddress("xB", &kaonm_xB);
    treeKaonM->SetBranchAddress("y", &kaonm_y);
    treeKaonM->SetBranchAddress("z", &kaonm_z);
    treeKaonM->SetBranchAddress("Pt", &kaonm_Pt);
    treeKaonM->SetBranchAddress("PhT", &kaonm_PhT);
    treeKaonM->SetBranchAddress("Phi_Lab", &kaonm_Phi);
    treeKaonM->SetBranchAddress("Theta_Lab", &kaonm_Theta);
    treeKaonM->SetBranchAddress("Pseudorapidity", &kaonm_eta);
    treeKaonM->SetBranchAddress("Phi_h", &kaonm_Phi_h);
    //treeKaonM->SetBranchAddress("SinPhi", &kaonm_SinPhi);
    treeKaonM->SetBranchAddress("helicity", &kaonm_helicity);

    // KAON+ PLOT
    dirKaonp->cd();
    // Mom
    TH2D kp_MomVsPhT ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  Kaon+  |; P_{hT} [GeV]; Mom [GeV]", 80, 0, 1, 80, 0.9, 7);
    TH2D kp_MomVsXb ("_MomVsXb", "Correlation Mom vs x_{B}  |  Kaon+  |; x_{B}; Mom [GeV]", 80, 0, 0.8, 80, 0.9, 7);
    TH2D kp_MomVsXf ("_MomVsXf", "Correlation Mom vs x_F  |  Kaon+  |; x_F; Mom [GeV]", 80, 0, 0.6, 80, 0.9, 7);
    TH2D kp_MomVsZ ("_MomVsZ", "Correlation Mom vs Z  |  Kaon+  |; z; Mom [GeV]", 80, 0.2, 0.9, 80, 0.9, 7);
    TH2D kp_MomVsY ("_MomVsY", "Correlation Mom vs Y  |  Kaon+  |; y; Mom [GeV]", 80, 0.2, 0.75, 80, 0.9, 7);
    TH2D kp_MomVsEta ("_MomVsEta", "Correlation Mom vs Eta  |  Kaon+  |; Eta; Mom [GeV]", 80, 1.5, 3.0, 80, 0.9, 7);
    TH2D kp_MomVsTheta ("_MomVsTheta", "Correlation Mom vs Theta  |  Kaon+  |; Mom [GeV]; Theta [Rad]", 80, 0.9, 7, 80, 0.09, 0.4);
    TH2D kp_MomVsPhi_h ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  Kaon+  |; #Phi_{h} [Rad]; Mom [GeV]", 80, -TMath::Pi(), TMath::Pi(), 80, 0.9, 7);
    TH2D kp_MomVsM2x ("_MomVsM2x", "correlation Mom vs M2x | Kaon+ |; Mom [GeV]; M2x[GeV^{2}]", 60, 1, 7, 60, 2.45, 11);
    // Q2
    TH2D kp_Q2VsXb ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  Kaon+  |; x_{B}; Q^{2} [GeV^{2}]", 80, 0, 0.8, 80, 1, 8);
    TH2D kp_Q2VsXf ("_Q2VsXf", "Correlation Q^{2} vs x_F  |  Kaon+  |; x_F; Q^{2} [GeV^{2}]", 80, 0, 0.6, 80, 1, 7);
    TH2D kp_Q2VsMom ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  Kaon+  |; Mom [GeV]; Q^{2} [GeV^{2}]", 80, 0.9, 6, 80, 1, 7);
    TH2D kp_Q2VsPhT ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  Kaon+  |; P_{hT} [GeV]; Q^{2} [GeV^{2}]", 80, 0, 1.2, 80, 1, 6);
    TH2D kp_Q2VsZ ("_Q2VsZ", "Correlation Q^{2} vs Z  |  Kaon+  |; z; Q^{2} [GeV^{2}]", 80, 0.2, 0.9, 80, 1, 6);
    TH2D kp_Q2VsY ("_Q2VsY", "Correlation Q^{2} vs Y  |  Kaon+  |; y; Q^{2} [GeV^{2}]", 80, 0.2, 0.75, 80, 0.9, 7);
    TH2D kp_Q2VsEta ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  Kaon+  |; Eta; Q^{2} [GeV^{2}]", 80, 1.5, 3.0, 80, 1, 6);
    TH2D kp_Q2VsPhi_h ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  Kaon+  |; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", 80, -TMath::Pi(), TMath::Pi(), 80, 1, 6);
    // PhT
    TH2D kp_PhTvsZ ("_PhTvsZ", "Correlation P_{hT} vs Z  |  Kaon+  |; z; P_{hT} [GeV]", 80, 0.2, 1, 80, 0, 1.2);
    TH2D kp_PhTvsXb ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  Kaon+  |; x_{B}; P_{hT} [GeV]", 80, 0, 0.8, 80, 0, 1.2);
    TH2D kp_PhTvsEta ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  Kaon+  |; Eta; P_{hT} [GeV]", 80, 1.5, 3.0, 80, 0, 1.2);
    TH2D kp_PhTvsPhi_h ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  Kaon+  |; #Phi_{h} [Rad]; P_{hT} [GeV]", 80, -TMath::Pi(), TMath::Pi(), 80, 0, 1.2);
    // Z
    TH2D kp_zVsXb ("_zVsXb", "Correlation Z vs x_{B}  |  Kaon+  |; x_{B}; z", 80, 0, 0.8, 80, 0.2, 0.9);
    TH2D kp_zVsXf ("_zVsXf", "Correlation Z vs x_F  |  Kaon+  |; x_F; z", 80, 0, 0.6, 80, 0.2, 0.9);
    TH2D kp_zVsEta ("_zVsEta", "Correlation Z vs Eta  |  Kaon+  |; Eta; z", 80, 1.5, 3.0, 80, 0.2, 0.9);
    TH2D kp_zVsPhi_h ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  Kaon+  |; #Phi_{h} [Rad]; z", 80, -TMath::Pi(), TMath::Pi(), 80, 0.2, 0.9);
    // Angles
    TH2D kp_ThetaVsPhi_h ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  Kaon+  |; #Phi_{h} [Rad]; Theta [Rad]", 80, -TMath::Pi(), TMath::Pi(), 80, 0, 0.2*TMath::Pi());
    //TH2D el_ThetaVsMom ("el_ThetaVsMom", "Correlation Theta vs Mom  |  electron  |; Mom [GeV]; Theta [Deg]", 80, 1.5, 9, 80, 5, 40);
    //TH1D kp_Phi_h ("_Phi_h", "#Phi_{h} distribution | Kaon+ ; #Phi_{h} [Deg]; counts", 60, 0, 360);
    //TH1D kp_Phi_Hp ("_Phi_H+", "#Phi_{h} distribution with positive helicity | Kaon+ ; #Phi_{h} [Deg]; counts", 60, -TMath::Pi(), TMath::Pi());
    //TH1D kp_Phi_Hm ("_Phi_H-", "#Phi_{h} distribution with negative helicity | Kaon+ ; #Phi_{h} [Deg]; counts", 60, -TMath::Pi(), TMath::Pi());

    // KAON+ PLOT Helicity+
    dirKaonp_Hp->cd();
    // Mom
    TH2D kp_MomVsPhT_Hp ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  Kaon+  | Helicity > 0; P_{hT} [GeV]; Mom [GeV]", 80, 0, 1, 80, 0.9, 7);
    TH2D kp_MomVsXb_Hp ("_MomVsXb", "Correlation Mom vs x_{B}  |  Kaon+  | Helicity > 0; x_{B}; Mom [GeV]", 80, 0, 0.8, 80, 0.9, 7);
    TH2D kp_MomVsXf_Hp ("_MomVsXf", "Correlation Mom vs x_F  |  Kaon+  | Helicity > 0; x_F; Mom [GeV]", 80, 0, 0.6, 80, 0.9, 7);
    TH2D kp_MomVsZ_Hp ("_MomVsZ", "Correlation Mom vs Z  |  Kaon+  | Helicity > 0; z; Mom [GeV]", 80, 0.2, 0.9, 80, 0.9, 7);
    TH2D kp_MomVsY_Hp ("_MomVsY", "Correlation Mom vs Y  |  Kaon+  | Helicity > 0; y; Mom [GeV]", 80, 0.2, 0.75, 80, 0.9, 7);
    TH2D kp_MomVsEta_Hp ("_MomVsEta", "Correlation Mom vs Eta  |  Kaon+  | Helicity > 0; Eta; Mom [GeV]", 80, 1.5, 3.0, 80, 0.9, 7);
    TH2D kp_MomVsTheta_Hp ("_MomVsTheta", "Correlation Mom vs Theta  |  Kaon+  | Helicity > 0; Mom [GeV]; Theta [Rad]", 80, 0.9, 7, 80, 0.09, 0.4);
    TH2D kp_MomVsPhi_h_Hp ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  Kaon+  | Helicity > 0; #Phi_{h} [Rad]; Mom [GeV]", 80, -TMath::Pi(), TMath::Pi(), 80, 0.9, 7);
    TH2D kp_MomVsM2x_Hp ("_MomVsM2x", "correlation Mom vs M2x | Kaon+ |; Mom [GeV]; M2x[GeV^{2}]", 60, 1, 7, 60, 2.45, 11);
    // Q2
    TH2D kp_Q2VsXb_Hp ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  Kaon+  | Helicity > 0; x_{B}; Q^{2} [GeV^{2}]", 80, 0, 0.8, 80, 1, 8);
    TH2D kp_Q2VsXf_Hp ("_Q2VsXf", "Correlation Q^{2} vs x_F  |  Kaon+  | Helicity > 0; x_F; Q^{2} [GeV^{2}]", 80, 0, 0.6, 80, 1, 7);
    TH2D kp_Q2VsMom_Hp ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  Kaon+  | Helicity > 0; Mom [GeV]; Q^{2} [GeV^{2}]", 80, 0.9, 6, 80, 1, 7);
    TH2D kp_Q2VsPhT_Hp ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  Kaon+  | Helicity > 0; P_{hT} [GeV]; Q^{2} [GeV^{2}]", 80, 0, 1.2, 80, 1, 6);
    TH2D kp_Q2VsZ_Hp ("_Q2VsZ", "Correlation Q^{2} vs Z  |  Kaon+  | Helicity > 0; z; Q^{2} [GeV^{2}]", 80, 0.2, 0.9, 80, 1, 6);
    TH2D kp_Q2VsY_Hp ("_Q2VsY", "Correlation Q^{2} vs Y  |  Kaon+  | Helicity > 0; y; Q^{2} [GeV^{2}]", 80, 0.2, 0.75, 80, 0.9, 7);
    TH2D kp_Q2VsEta_Hp ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  Kaon+  | Helicity > 0; Eta; Q^{2} [GeV^{2}]", 80, 1.5, 3.0, 80, 1, 6);
    TH2D kp_Q2VsPhi_h_Hp ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  Kaon+  | Helicity > 0; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", 80, -TMath::Pi(), TMath::Pi(), 80, 1, 6);
    // PhT
    TH2D kp_PhTvsZ_Hp ("_PhTvsZ", "Correlation P_{hT} vs Z  |  Kaon+  | Helicity > 0; z; P_{hT} [GeV]", 80, 0.2, 1, 80, 0, 1.2);
    TH2D kp_PhTvsXb_Hp ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  Kaon+  | Helicity > 0; x_{B}; P_{hT} [GeV]", 80, 0, 0.8, 80, 0, 1.2);
    TH2D kp_PhTvsEta_Hp ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  Kaon+  | Helicity > 0; Eta; P_{hT} [GeV]", 80, 1.5, 3.0, 80, 0, 1.2);
    TH2D kp_PhTvsPhi_h_Hp ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  Kaon+  | Helicity > 0; #Phi_{h} [Rad]; P_{hT} [GeV]", 80, -TMath::Pi(), TMath::Pi(), 80, 0, 1.2);
    // Z
    TH2D kp_zVsXb_Hp ("_zVsXb", "Correlation Z vs x_{B}  |  Kaon+  | Helicity > 0; x_{B}; z", 80, 0, 0.8, 80, 0.2, 0.9);
    TH2D kp_zVsXf_Hp ("_zVsXf", "Correlation Z vs x_F  |  Kaon+  | Helicity > 0; x_F; z", 80, 0, 0.6, 80, 0.2, 0.9);
    TH2D kp_zVsEta_Hp ("_zVsEta", "Correlation Z vs Eta  |  Kaon+  | Helicity > 0; Eta; z", 80, 1.5, 3.0, 80, 0.2, 0.9);
    TH2D kp_zVsPhi_h_Hp ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  Kaon+  | Helicity > 0; #Phi_{h} [Rad]; z", 80, -TMath::Pi(), TMath::Pi(), 80, 0.2, 0.9);
    //
    TH2D kp_ThetaVsPhi_h_Hp ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  Kaon+  | Helicity > 0; #Phi_{h} [Rad]; Theta [Rad]", 80, -TMath::Pi(), TMath::Pi(), 80, 0, 0.2*TMath::Pi());
    

    // KAON+ PLOT Helicity-
    dirKaonp_Hm->cd();
    // Mom
    TH2D kp_MomVsPhT_Hm ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  Kaon+  | Helicity < 0; P_{hT} [GeV]; Mom [GeV]", 80, 0, 1, 80, 0.9, 7);
    TH2D kp_MomVsXb_Hm ("_MomVsXb", "Correlation Mom vs x_{B}  |  Kaon+  | Helicity < 0; x_{B}; Mom [GeV]", 80, 0, 0.8, 80, 0.9, 7);
    TH2D kp_MomVsXf_Hm ("_MomVsXf", "Correlation Mom vs x_F  |  Kaon+  | Helicity < 0; x_F; Mom [GeV]", 80, 0, 0.6, 80, 0.9, 7);
    TH2D kp_MomVsZ_Hm ("_MomVsZ", "Correlation Mom vs Z  |  Kaon+  | Helicity < 0; z; Mom [GeV]", 80, 0.2, 0.9, 80, 0.9, 7);
    TH2D kp_MomVsY_Hm ("_MomVsY", "Correlation Mom vs Y  |  Kaon+  | Helicity < 0; y; Mom [GeV]", 80, 0.2, 0.75, 80, 0.9, 7);
    TH2D kp_MomVsEta_Hm ("_MomVsEta", "Correlation Mom vs Eta  |  Kaon+  | Helicity < 0; Eta; Mom [GeV]", 80, 1.5, 3.0, 80, 0.9, 7);
    TH2D kp_MomVsTheta_Hm ("_MomVsTheta", "Correlation Mom vs Theta  |  Kaon+  | Helicity < 0; Mom [GeV]; Theta [Rad]", 80, 0.9, 7, 80, 0.09, 0.4);
    TH2D kp_MomVsPhi_h_Hm ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  Kaon+  | Helicity < 0; #Phi_{h} [Rad]; Mom [GeV]", 80, -TMath::Pi(), TMath::Pi(), 80, 0.9, 7);
    TH2D kp_MomVsM2x_Hm ("_MomVsM2x", "correlation Mom vs M2x | Kaon+ |; Mom [GeV]; M2x[GeV^{2}]", 60, 1, 7, 60, 2.45, 11);
    // Q2
    TH2D kp_Q2VsXb_Hm ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  Kaon+  | Helicity < 0; x_{B}; Q^{2} [GeV^{2}]", 80, 0, 0.8, 80, 1, 8);
    TH2D kp_Q2VsXf_Hm ("_Q2VsXf", "Correlation Q^{2} vs x_F  |  Kaon+  | Helicity < 0; x_F; Q^{2} [GeV^{2}]", 80, 0, 0.6, 80, 1, 7);
    TH2D kp_Q2VsMom_Hm ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  Kaon+  | Helicity < 0; Mom [GeV]; Q^{2} [GeV^{2}]", 80, 0.9, 6, 80, 1, 7);
    TH2D kp_Q2VsPhT_Hm ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  Kaon+  | Helicity < 0; P_{hT} [GeV]; Q^{2} [GeV^{2}]", 80, 0, 1.2, 80, 1, 6);
    TH2D kp_Q2VsZ_Hm ("_Q2VsZ", "Correlation Q^{2} vs Z  |  Kaon+  | Helicity < 0; z; Q^{2} [GeV^{2}]", 80, 0.2, 0.9, 80, 1, 6);
    TH2D kp_Q2VsY_Hm ("_Q2VsY", "Correlation Q^{2} vs Y  |  Kaon+  | Helicity < 0; y; Q^{2} [GeV^{2}]", 80, 0.2, 0.75, 80, 0.9, 7);
    TH2D kp_Q2VsEta_Hm ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  Kaon+  | Helicity < 0; Eta; Q^{2} [GeV^{2}]", 80, 1.5, 3.0, 80, 1, 6);
    TH2D kp_Q2VsPhi_h_Hm ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  Kaon+  | Helicity < 0; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", 80, -TMath::Pi(), TMath::Pi(), 80, 1, 6);
    // PhT
    TH2D kp_PhTvsZ_Hm ("_PhTvsZ", "Correlation P_{hT} vs Z  |  Kaon+  | Helicity < 0; z; P_{hT} [GeV]", 80, 0.2, 1, 80, 0, 1.2);
    TH2D kp_PhTvsXb_Hm ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  Kaon+  | Helicity < 0; x_{B}; P_{hT} [GeV]", 80, 0, 0.8, 80, 0, 1.2);
    TH2D kp_PhTvsEta_Hm ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  Kaon+  | Helicity < 0; Eta; P_{hT} [GeV]", 80, 1.5, 3.0, 80, 0, 1.2);
    TH2D kp_PhTvsPhi_h_Hm ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  Kaon+  | Helicity < 0; #Phi_{h} [Rad]; P_{hT} [GeV]", 80, -TMath::Pi(), TMath::Pi(), 80, 0, 1.2);
    // Z
    TH2D kp_zVsXb_Hm ("_zVsXb", "Correlation Z vs x_{B}  |  Kaon+  | Helicity < 0; x_{B}; z", 80, 0, 0.8, 80, 0.2, 0.9);
    TH2D kp_zVsXf_Hm ("_zVsXf", "Correlation Z vs x_F  |  Kaon+  | Helicity < 0; x_F; z", 80, 0, 0.6, 80, 0.2, 0.9);
    TH2D kp_zVsEta_Hm ("_zVsEta", "Correlation Z vs Eta  |  Kaon+  | Helicity < 0; Eta; z", 80, 1.5, 3.0, 80, 0.2, 0.9);
    TH2D kp_zVsPhi_h_Hm ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  Kaon+  | Helicity < 0; #Phi_{h} [Rad]; z", 80, -TMath::Pi(), TMath::Pi(), 80, 0.2, 0.9);
    //
    TH2D kp_ThetaVsPhi_h_Hm ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  Kaon+  | Helicity < 0; #Phi_{h} [Rad]; Theta [Rad]", 80, -TMath::Pi(), TMath::Pi(), 80, 0, 0.2*TMath::Pi());
    
    // KAON- PLOT
    dirKaonm->cd();
    // Mom
    TH2D km_MomVsPhT ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  Kaon-  |; P_{hT} [GeV]; Mom[GeV]", 60, 0, 1, 60, 0.9, 6);
    TH2D km_MomVsXb ("_MomVsXb", "Correlation Mom vs x_{B}  |  Kaon-  |; x_{B}; Mom[GeV]", 60, 0, 0.8, 60, 0.9, 6);
    TH2D km_MomVsXf ("_MomVsXf", "Correlation Mom vs x_F  |  Kaon-  |; x_F; Mom[GeV]", 60, 0, 0.6, 60, 0.9, 6);
    TH2D km_MomVsZ ("_MomVsZ", "Correlation Mom vs Z  |  Kaon-  |; z; Mom[GeV]", 60, 0.2, 0.9, 60, 0.9, 6);
    TH2D km_MomVsY ("_MomVsY", "Correlation Mom vs Y  |  Kaon-  |; y; Mom[GeV]", 60, 0.2, 0.75, 60, 0.9, 6);
    TH2D km_MomVsEta ("_MomVsEta", "Correlation Mom vs Eta  |  Kaon-  |; Eta; Mom[GeV]", 60, 1.5, 3, 60, 0.9, 6);
    TH2D km_MomVsTheta ("_MomVsTheta", "Correlation Mom vs Theta  |  Kaon-  |; Mom [GeV]; Theta [Rad]", 80, 0.9, 6, 80, 0.09, 0.4);
    TH2D km_MomVsPhi_h ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  Kaon-  |; #Phi_{h} [Rad]; Mom[GeV]", 60, -TMath::Pi(), TMath::Pi(), 60, 0.9, 6);
    // Q2
    TH2D km_Q2VsXb ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  Kaon-  |; x_{B}; Q^{2} [GeV^{2}]", 60, 0, 0.8, 60, 1, 6);
    TH2D km_Q2VsXf ("_Q2VsXf", "Correlation Q^{2} vs x_F  |  Kaon-  |; x_F; Q^{2} [GeV^{2}]", 60, 0, 0.6, 60, 1, 6);
    TH2D km_Q2VsMom ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  Kaon-  |; Mom [GeV]; Q^{2} [GeV^{2}]", 60, 0.9, 6, 60, 1, 6);
    TH2D km_Q2VsPhT ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  Kaon-  |; P_{hT} [GeV]; Q^{2} [GeV^{2}]", 60, 0, 1.1, 60, 1, 6);
    TH2D km_Q2VsZ ("_Q2VsZ", "Correlation Q^{2} vs Z  |  Kaon-  |; z; Q^{2} [GeV^{2}]", 60, 0.2, 0.9, 60, 1, 6);
    TH2D km_Q2VsY ("_Q2VsY", "Correlation Q^{2} vs Y  |  Kaon-  |; y; Q^{2} [GeV^{2}]", 80, 0.2, 0.75, 80, 0.9, 6);
    TH2D km_Q2VsEta ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  Kaon-  |; Eta; Q^{2} [GeV^{2}]", 60, 1.5, 3, 60, 1, 6);
    TH2D km_Q2VsPhi_h ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  Kaon-  |; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", 60, -TMath::Pi(), TMath::Pi(), 60, 1, 6);
    // PhT
    TH2D km_PhTvsZ ("_PhTvsZ", "Correlation P_{hT} vs Z  |  Kaon-  |; z; P_{hT} [GeV]", 60, 0.2, 1, 60, 0, 1.1);
    TH2D km_PhTvsXb ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  Kaon-  |; x_{B}; P_{hT} [GeV]", 60, 0, 0.8, 60, 0, 1.1);
    TH2D km_PhTvsEta ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  Kaon-  |; Eta; P_{hT} [GeV]", 60, 1.5, 3, 60, 0, 1.1);
    TH2D km_PhTvsPhi_h ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  Kaon-  |; #Phi_{h} [Rad]; P_{hT} [GeV]", 60, -TMath::Pi(), TMath::Pi(), 60, 0, 1.1);
    // Z
    TH2D km_zVsXb ("_zVsXb", "Correlation Z vs x_{B}  |  Kaon-  |; x_{B}; z", 60, 0, 0.8, 60, 0.2, 0.9);
    TH2D km_zVsXf ("_zVsXf", "Correlation Z vs x_F  |  Kaon-  |; x_F; z", 60, 0, 0.6, 60, 0.2, 0.9);
    TH2D km_zVsEta ("_zVsEta", "Correlation Z vs Eta  |  Kaon-  |; Eta; z", 60, 1.5, 3, 60, 0.2, 0.9);
    TH2D km_zVsPhi_h ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  Kaon-  |; #Phi_{h} [Rad]; z", 60, -TMath::Pi(), TMath::Pi(), 60, 0.2, 0.9);
    // Angles
    TH2D km_ThetaVsPhi_h ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  Kaon-  |; #Phi_{h} [Rad]; Theta [Rad]", 60, -TMath::Pi(), TMath::Pi(), 60, 0, 0.2*TMath::Pi());
    //TH1D km_Phi_Hp ("_Phi_H+", "#Phi_{h} distribution with positive helicity | Kaon- ; #Phi_{h} [Deg]; counts", 60, -TMath::Pi(), TMath::Pi());
    //TH1D km_Phi_Hm ("_Phi_H-", "#Phi_{h} distribution with negative helicity | Kaon- ; #Phi_{h} [Deg]; counts", 60, -TMath::Pi(), TMath::Pi());

    // ORA RIEMPI I GRAFICI
    // Kaon+
    Long64_t nEntries_kp = treeKaonP->GetEntries();
    for (Long64_t i = 0; i < nEntries_kp; i++) {
        // Mom
        treeKaonP->GetEntry(i);
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
        kp_ThetaVsPhi_h.Fill(kaonp_Phi_h, kaonp_Theta);
    }

    // Kaon+ H+
    Long64_t nEntries_kp_Hp = treeKaonPhp->GetEntries();
    for (Long64_t i = 0; i < nEntries_kp_Hp; i++) {
        // Mom
        treeKaonPhp->GetEntry(i);
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
        kp_ThetaVsPhi_h_Hp.Fill(kaonp_Phi_h_hp, kaonp_Theta_hp);
    }

    // Kaon+ H-
    Long64_t nEntries_kp_Hm = treeKaonPhm->GetEntries();
    for (Long64_t i = 0; i < nEntries_kp_Hm; i++) {
        // Mom
        treeKaonPhm->GetEntry(i);
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
        kp_ThetaVsPhi_h_Hm.Fill(kaonp_Phi_h_hm, kaonp_Theta_hm);
    }

    // Kaon-
    Long64_t nEntries_km = treeKaonM->GetEntries();
    for (Long64_t i = 0; i < nEntries_km; i++) {
        treeKaonM->GetEntry(i);
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
        km_ThetaVsPhi_h.Fill(kaonm_Phi_h, kaonm_Theta);
    }

    // Binning plot
    dirKaonp->cd();
    // Q2 vs xB
    TCanvas *c1 = new TCanvas("Q2_vs_xB_Bin", "Q^{2} vs x_{B} bin", 800, 600);
    kp_Q2VsXb.Draw("COLZ");
    double xB1[] = {0.05, 0.15, 0.15, 0.05, 0.05};  
    double Q21[] = {1.0, 1.0, 2.35, 2.35, 1.0};
    double xB2[] = {0.15, 0.35, 0.35, 0.15, 0.15};  
    double Q22[] = {1.25, 1.25, 2.5, 2.5, 1.25};
    double xB3[] = {0.16, 0.67, 0.67, 0.16, 0.16};  
    double Q23[] = {2.5, 2.5, 8, 8, 2.5};
    TPolyLine *rect1 = new TPolyLine(5, xB1, Q21);
    //rect1->SetLineWidth(2);
    rect1->Draw("same");
    TPolyLine *rect2 = new TPolyLine(5, xB2, Q22);
    //rect2->SetLineWidth(2);
    rect2->Draw("same");
    TPolyLine *rect3 = new TPolyLine(5, xB3, Q23);
    //rect3->SetLineWidth(2);
    rect3->Draw("same");
    c1->Update();
    c1->Write();
    // z vs PhT
    TCanvas *c2 = new TCanvas("z_vs_PhT", "z vs P_{hT} bin", 800, 600);
    kp_PhTvsZ.Draw("COLZ");
    double z1[] = {0.2, 0.9, 0.9, 0.2, 0.2};
    double p1[] = {0.0, 0.0, 1.1, 1.1, 0.0};
    double z2[] = {0.3, 0.3};
    double p2[] = {0.0, 1.1};
    double z3[] = {0.4, 0.4};
    double p3[] = {0.0, 1.1};
    double z4[] = {0.5, 0.5};
    double p4[] = {0.0, 1.1};
    double z5[] = {0.6, 0.6};
    double p5[] = {0.0, 1.1};
    double z6[] = {0.7, 0.7};
    double p6[] = {0.0, 1.1};
    double z7[] = {0.2, 0.9};
    double p7[] = {0.25, 0.25};
    double z8[] = {0.2, 0.9};
    double p8[] = {0.5, 0.5};
    TPolyLine *a1 = new TPolyLine(5, z1, p1);
    TPolyLine *a2 = new TPolyLine(2, z2, p2);
    TPolyLine *a3 = new TPolyLine(2, z3, p3);
    TPolyLine *a4 = new TPolyLine(2, z4, p4);
    TPolyLine *a5 = new TPolyLine(2, z5, p5);
    TPolyLine *a6 = new TPolyLine(2, z6, p6);
    TPolyLine *a7 = new TPolyLine(2, z7, p7);
    TPolyLine *a8 = new TPolyLine(2, z8, p8);
    //a1->SetLineWidth(2);
    a1->Draw("same");
    a2->Draw("same");
    a3->Draw("same");
    a4->Draw("same");
    a5->Draw("same");
    a6->Draw("same");
    a7->Draw("same");
    a8->Draw("same");
    c2->Update();
    c2->Write();
    // PhT vs z
    

    outFile.Write();
    outFile.Close();
    inFile.Close();

    cout << "ROOT output file: " << outputFile << endl;
}
