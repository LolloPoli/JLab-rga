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

// creo la funzione per i TLorenztVector
// region_part è una classe di clas12 che rappresenta una particella includendo le info relative ai dati del rivelatore
//_ptr dovrebbe essere un alias intelligente (smart pointer), il quale gestisce la memoria automaticamente
void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), p4.M());
}

// copiata e incollata da github non so se vunzioni
void FTElectronCorrection(TLorentzVector& el){
	auto p_electron_pre_correction = el;
	double p_mag = p_electron_pre_correction.P();
	double p_mag_corrected = (p_mag-0.03689+0.1412*p_mag-0.04316*pow(p_mag,2)+0.007046*pow(p_mag,3)-0.0004055*pow(p_mag,4))/p_mag;
	el.SetXYZM(el.X()*p_mag_corrected,el.Y()*p_mag_corrected,el.Z()*p_mag_corrected,0.00051099891);
}

bool InvertHelicity(int runNumber){
  if(runNumber == 11) return false; //SIMULATION
  if(runNumber == 0) return false; //SIMULATION
  if(runNumber >= 5032 && runNumber <= 5666) return true; //RGA FALL2018
  if(runNumber >= 6616 && runNumber <= 6783) return true; //RGA SPRING2019
  if(runNumber >= 6120 && runNumber <= 6604) return true; //RGB SPRING2019
  if(runNumber >= 11093 && runNumber <= 11283) return false; //RGB FALL2019
  if(runNumber >= 11284 && runNumber <= 11300) return true; //RGB FALL2019
  if(runNumber >= 11323 && runNumber <= 11571) return false; //RGB SPRING2020
  cout <<Form("Run %d not included in invertHelicity method\n",runNumber);
  return false;
}
// Electron PID
// B1
bool ElectronPID_ForwardDetector(int status){
    return (status > -4000 && status <= 2000);
}
// B2
bool ElectronPID_nPhe(int nphe){
    return nphe > 2;
}
// B3
bool ElectronPID_PCAL(double pcal){
    return pcal > 0.07;
}
// B4
bool ElectronPID_CalSFcut(int sector, int runnum, double p, double cal_energy){
    // calorimeter sampling fraction cut
    // getSector(), Nrun, electron mom, electron dep E in every cal, ECAL, PCAL, ECIN? or idk the itter o pre shower
  double scale = 3.5; // how many std away from mean to cut on
  // Common calculation for mean and std
  double mean = 0.0;
  double std = 0.0;
  std::vector<std::vector<double>> e_cal_sampl_mu;
  std::vector<std::vector<double>> e_cal_sampl_sigma;

  if ((runnum == 11) || ((runnum >= 5032 && runnum <= 5666) || (runnum >= 6616 && runnum <= 6783))) { // RGA
    e_cal_sampl_mu = {
      {0.2531, 0.2550, 0.2514, 0.2494, 0.2528, 0.2521},
      {-0.6502, -0.7472, -0.7674, -0.4913, -0.3988, -0.703},
      {4.939, 5.350, 5.102, 6.440, 6.149, 4.957}
    };

    e_cal_sampl_sigma = {
      {0.002726, 0.004157, 0.00522, 0.005398, 0.008453, 0.006553},
      {1.062, 0.859, 0.5564, 0.6576, 0.3242, 0.4423},
      {-4.089, -3.318, -2.078, -2.565, -0.8223, -1.274}
    };
  } else if (runnum >= 6120 && runnum <= 6604) { // RGB
    e_cal_sampl_mu = {
      {0.2520, 0.2520, 0.2479, 0.2444, 0.2463, 0.2478},
      {-0.8615, -0.8524, -0.6848, -0.5521, -0.5775, -0.7327},
      {5.596, 6.522, 5.752, 5.278, 6.430, 5.795}
    };

    e_cal_sampl_sigma = {
      {-0.02963, -0.1058, -0.05087, -0.04524, -0.02951, -0.01769},
      {20.4, 129.3, 0.6191, 0.6817, 20.84, 8.44},
      {-41.44, -101.6, -2.673, -2.606, -42.67, -21.73}
    };
  } else if ((runnum >= 11323 && runnum <= 11571) || (runnum >= 11093 && runnum <= 11300)) {
    // RGB winter 2020 // (also using this for RGB fall 2019, but it should be updated! TODO)
    e_cal_sampl_mu = {
      {0.2433, 0.2421, 0.2415, 0.2486, 0.2419, 0.2447},
      {-0.8052, -1.0495, -1.1747, -0.5170, -0.6840, -0.9022},
      {5.2750, 4.4886, 4.4935, 5.9044, 5.6716, 4.9288}
    };
    e_cal_sampl_sigma = {
      {0.0120, 0.0164, 0.0120, 0.0108, 0.0147, 0.0077},
      {0.1794, 0.1519, 0.1379, 0.1838, 0.0494, 0.3509},
      {-0.0695, 0.1553, 0.3300, 0.4330, 1.1032, -0.7996}
    };
  }

  // Calculation of mean and std
  mean = e_cal_sampl_mu[0][sector] + (e_cal_sampl_mu[1][sector] / 1000) * (p - e_cal_sampl_mu[2][sector]) * (p - e_cal_sampl_mu[2][sector]);
  std = e_cal_sampl_sigma[0][sector] + e_cal_sampl_sigma[1][sector] / (10 * (p - e_cal_sampl_sigma[2][sector]));

  // Return result
  return ((cal_energy / p) > (mean - scale * std)) && ((cal_energy / p) < (mean + scale * std));
}
// B5
bool ElectronPID_DiagonalCut(clas12::region_part_ptr pa, TLorentzVector *p){
  bool response = true;
  double einner = pa->cal(ECIN)->getEnergy();
  //double einner = pa->cal(ECIN)->getEnergy()+pa->cal(ECOUT)->getEnergy();
  if(p->P() > 4.5){
    double xxx = einner/p->P();
    double yyy = pa->cal(PCAL)->getEnergy()/p->P();
    if(xxx + yyy < (0.2) )response = false;
    //cout <<Form("Diagonal cut check: %lf + %lf = %lf > 0.2) ? %d",xxx,yyy,xxx+yyy,response);
    //cin.get();
  }
  return response;
}
//B6
bool ElectronPID_VertexCut(double vz){
    return (vz > -8 && vz < 3);
}
//
// HADRON PID 
// C2
// they say 1.2 GeV but is very high, so I must try with 2M_k as lower limit, which will remove the largest fraction of resonances
bool HadronPID_Reson(double mom, double m2x){
  if(mom > 1.1){
    if(m2x > 2.32) return true;
    else return false;
  }  
  else return false;
  /*
  bool signal = false;
  if(mom >= 1){
    if(mom < 1.2){
      if(m2x <= 1.8 || m2x >= 5.8) signal = true;
    }
    if(mom >= 1.2 && mom < 1.3){
      if(m2x <= 1.8 || m2x >= 4.2) signal = true;
    } 
    if(mom >= 1.3 && mom <= 1.4){
      // stay avray from Lambda (1115, 1405, 1520, 1600, 1670, 1690, 1710, 1800), 1520 and 1670 are the two strongest resonances
      // Sigma (1385, 1660, 1670)
      // N* (1535, 1650)
      if(m2x <= 1.8 || m2x >= 3.8) signal = true; 
    }
    if(mom > 1.4 && mom < 1.5){
      if(m2x < 2 || m2x > 2.6) signal = true;
    }
    if(mom >= 1.5) signal = true;
  }
  return signal;
  */
}

bool HadronPID_Q2(double q2){
  return q2 > 1;
}

bool KinematicPID_xF(double xF){
  return xF > 0;
}

bool KinematicPID_W(double w){
  return w > 2;
}

bool HadronPID_DC(double edge1, double edge2, double edge3, int torus){
  if(torus == -1){  // INBENDING HADRON 
    if(edge1 > 3 && edge2 > 3 && edge3 > 7) return true;
    else return false;
  }
  else if(torus == 1){  // OUTBENDING HADRON  
    if(edge1 > 3 && edge2 > 3 && edge3 > 9) return true;
    else return false;
  }
  else return false;
}

bool ElectronPID_DC(double edge1, double edge2, double edge3, int torus){
  if(torus == -1){  // INBENDING ELECTRON 
    if(edge1 > 5 && edge2 > 5 && edge3 > 10) return true;
    else return false;
  }
  else if(torus == 1){  // OUTBENDING ELECTRON
    if(edge1 > 3 && edge2 > 3 && edge3 > 10) return true;
    else return false;
  }
  else return false;
}

bool KinematicPID_Vertex(double vz, int torus, int pid){
  bool signal = false;
  if(torus == -1){  // INBENDING
    if(pid > 20){
      if(-10 < vz && vz < 2.5) signal = true;
    } 
    if(pid < 20){
      if(-8 < vz && vz < 3) signal = true;
    }
  }
  if(torus == 1){   // OUTBENDING
    if(pid > 20){
      if(-8 < vz && vz < 3) signal = true;
    } 
    if(pid < 20){
      if(-10 < vz && vz < 2.5) signal = true;
    }
  }
  return signal;
}

bool HadronPID_Chi2Pid(double chi2){
  return (std::abs(chi2) < 3);
}

bool HadronPID_z(double z){
  return z > 0.2;
}

float BeamPolarization(Int_t run, Bool_t v) { //nella booleana metti true se vuoi la pol e false se vuoi l'errore
  /* RGA */
  if     (run>= 5032 && run<= 5332) return v? 0.8592 : 0.01290; //rga_fa18, before Wien angle change
  else if(run>= 5333 && run<= 5666) return v? 0.8922 : 0.02509; //rga_fa18, after Wien angle change
  else if(run>= 6616 && run<= 6783) return v? 0.8453 : 0.01474; //rga_sp19 https://logbooks.jlab.org/entry/3677077
  /* RGB */
  else if(run>= 6142 && run<= 6149) return v? 0.81132 : 0.01505; //rgb_sp19
  else if(run>= 6150 && run<= 6188) return v? 0.82137 : 0.01491; //https://clasweb.jlab.org/wiki/images/c/ca/Moller_Runs_15Jan.pdf
  else if(run>= 6189 && run<= 6260) return v? 0.83598 : 0.01475;
  else if(run>= 6261 && run<= 6339) return v? 0.80770 : 0.01449;
  else if(run>= 6340 && run<= 6342) return v? 0.85536 : 0.01484;
  else if(run>= 6344 && run<= 6399) return v? 0.87038 : 0.01474;
  else if(run>= 6420 && run<= 6476) return v? 0.88214 : 0.01502;
  else if(run>= 6479 && run<= 6532) return v? 0.86580 : 0.01460;
  else if(run>= 6533 && run<= 6603) return v? 0.87887 : 0.01454;
  else if(run>=11013 && run<=11309) return v? 0.84983 : 0.02929; //rgb_fa19
  else if(run>=11323 && run<=11334) return v? 0.87135 : 0.01464; //rgb_wi20
  else if(run>=11335 && run<=11387) return v? 0.85048 : 0.01530;
  else if(run>=11389 && run<=11571) return v? 0.84262 : 0.01494; //NOTE: table last updated 1/15/2020, but run ended on 1/30
  /* MC */
  else if(run==11) return v? 0.86 : 0.0; //MC
  else {
    fprintf(stderr,"ERROR: RundepPolarization unknown for run %d\n",run);
    return 0.0;
  }
}

void start() {
    //const string filePath = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_005032.hipo";
    clas12root::HipoChain chain;
    // > 5423 -> torus +1
    // < 5419 -> torus -1
    /*
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_005160.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_005162.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_005163.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_005164.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_005165.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_005166.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_005167.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_005168.hipo");
    */
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/nSidis_005423.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/nSidis_005424.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/nSidis_005425.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/nSidis_005426.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/nSidis_005428.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/nSidis_005429.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/nSidis_005430.hipo");
    chain.Add("/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/nSidis_005431.hipo");
    // Istanza di TDatabasePDG per accedere alle proprietà delle particelle, non va in conflitto con clas12database, sono indipendenti
    auto db2 = TDatabasePDG::Instance();
    // Definizione di alcuni TLorentzVector utili
    TLorentzVector beam(0, 0, 10.6, 10.6); // Fascio con energia 10.6 GeV
    TLorentzVector target(0, 0, 0, db2->GetParticle(2212)->Mass());  // Bersaglio 
    TLorentzVector el(0, 0, 0, db2->GetParticle(11)->Mass());        // Elettrone
    TLorentzVector pr(0, 0, 0, db2->GetParticle(2212)->Mass());      // Protone
    TLorentzVector gm(0, 0, 0, db2->GetParticle(22)->Mass());        // Photon
    TLorentzVector pip(0,0,0,db2->GetParticle(211)->Mass());         // Pion+
    TLorentzVector pim(0,0,0,db2->GetParticle(-211)->Mass());        // Pion-
    TLorentzVector kp(0,0,0,db2->GetParticle(321)->Mass());          // Kaon+
    TLorentzVector km(0,0,0,db2->GetParticle(-321)->Mass());         // Kaon-
    TLorentzVector Lab = beam + target;

    // variabili
    TVector3 yAxis_Lab (0,1,0);
    // info taglio
    unsigned int helicity;
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
    double kaonp_Phi_h, kaonp_Phi, kaonp_Theta, kaonp_y, kaonp_s, kaonp_E, kaonp_W;
    double kaonp_Ph_x, kaonp_Ph_y, kaonp_PhT, kaonp_eta, kaonp_M2x, kaonp_chi2pid;
    double kaonp_edge1, kaonp_edge2, kaonp_edge3, kaonp_vz, kaonp_SinPhi;
    double kaonp_rich_Id, kaonp_rich_PID = 0;
    // kaone -
    double kaonm_px, kaonm_py, kaonm_pz;
    double kaonm_xF, kaonm_xB, kaonm_Q2, kaonm_z, kaonm_Ph, kaonm_Pt;
    double kaonm_Phi_h, kaonm_Phi, kaonm_Theta, kaonm_y, kaonm_s, kaonm_E, kaonm_W;
    double kaonm_Ph_x, kaonm_Ph_y, kaonm_PhT, kaonm_eta, kaonm_M2x, kaonm_chi2pid;
    double kaonm_edge1, kaonm_edge2, kaonm_edge3, kaonm_vz, kaonm_SinPhi;
    // BINNING
    const int nZBins = 3;
    const int nPtBins = 3;
    double zBins[nZBins + 1] = {0.2, 0.4, 0.6, 1.0};  // Binning in z
    double PtBins[nPtBins + 1] = {0.0, 0.25, 0.5, 1.2}; // Binning in Pt
    TH1D* h_phi[nZBins][nPtBins];
    //clas12reader c12(filePath);
    //chain.SetRcdbData("rcdb.root");
    auto config_c12=chain.GetC12Reader();
    auto& c12=chain.C12ref();
    clas12databases::SetRCDBRemoteConnection();
    chain.WriteRcdbData("rcdb.root");
    //clas12databases::SetRCDBRootConnection("rcdb.root");
    clas12databases db;
    c12->connectDataBases(&db);
    //c12.showBanks();
    // creo un output root 
    TFile outFile("output.root", "RECREATE");  // File di output ROOT
    TDirectory* dirKaonp = outFile.mkdir("KaonPlus");
    TDirectory* dirKaonp_Bin = outFile.mkdir("KaonPlus_Bin");
    TDirectory* dirKaonm = outFile.mkdir("KaonMinus");
    TTree treeEl("Electron", "a");
    TTree treePh("Photon", "");
    TTree treePionP("Pion+", "");
    TTree treeKaonP("Kaon+", "");
    TTree treeKaonM("Kaon-", "");
    // salvataggio dei tree
    // Electron
    treeEl.Branch("px", &electron_px, "px/D");
    treeEl.Branch("py", &electron_py, "py/D");
    treeEl.Branch("pz", &electron_pz, "pz/D");
    treeEl.Branch("Mom", &electron_mom, "Mom/D");
    treeEl.Branch("vz", &electron_vz, "vz/D");
    treeEl.Branch("Q2", &electron_Q2, "Q2/D");
    treeEl.Branch("W", &electron_W, "W/D");
    treeEl.Branch("Theta", &electron_Theta, "Theta/D");
    treeEl.Branch("Nphe", &electron_Nphe, "Nphe/I");
    treeEl.Branch("status", &electron_status, "status/I");
    treeEl.Branch("PCAL", &electron_PCAL, "PCAL/D");
    treeEl.Branch("ECAL", &electron_ECAL, "ECAL/D");
    treeEl.Branch("ECIN", &electron_ECIN, "ECIN/D");
    treeEl.Branch("DC_6", &electron_edge1, "DC_6/D");
    treeEl.Branch("DC_18", &electron_edge2, "DC_18/D");
    treeEl.Branch("DC_36", &electron_edge3, "DC_36/D");
    // Gamma
    treePh.Branch("px", &gamma_px, "px/D");
    treePh.Branch("py", &gamma_py, "py/D");
    treePh.Branch("pz", &gamma_pz, "pz/D");
    // Pion+
    treePionP.Branch("px", &pionp_px, "px/D");
    treePionP.Branch("py", &pionp_py, "py/D");
    treePionP.Branch("pz", &pionp_pz, "pz/D");
    // Kaon +
    treeKaonP.Branch("E", &kaonp_E, "E/D");
    treeKaonP.Branch("px", &kaonp_px, "px/D");
    treeKaonP.Branch("py", &kaonp_py, "py/D");
    treeKaonP.Branch("pz", &kaonp_pz, "pz/D");
    treeKaonP.Branch("Mom", &kaonp_Ph, "Mom/D");
    treeKaonP.Branch("W", &kaonp_W, "W/D");
    treeKaonP.Branch("Q2", &kaonp_Q2, "Q2/D");
    treeKaonP.Branch("xF", &kaonp_xF, "xF/D");
    treeKaonP.Branch("xB", &kaonp_xB, "xB/D");
    treeKaonP.Branch("y", &kaonp_y, "y/D");
    treeKaonP.Branch("z", &kaonp_z, "z/D");
    treeKaonP.Branch("Pt", &kaonp_Pt, "Pt/D");
    treeKaonP.Branch("PhT", &kaonp_PhT, "PhT/D");
    treeKaonP.Branch("Phi_Lab", &kaonp_Phi, "Phi_Lab/D");
    treeKaonP.Branch("Theta_Lab", &kaonp_Theta, "Theta_Lab/D");
    treeKaonP.Branch("Pseudorapidity", &kaonp_eta, "Pseudorapidity/D");
    treeKaonP.Branch("Phi_h", &kaonp_Phi_h, "Phi_h/D");
    treeKaonP.Branch("SinPhi", &kaonp_SinPhi, "SinPhi/D");
    treeKaonP.Branch("M2_x", &kaonp_M2x, "M2_x/D");
    treeKaonP.Branch("Rich_Id", &kaonp_rich_Id, "Rich_Id/D");
    treeKaonP.Branch("Rich_PID", &kaonp_rich_PID, "Rich_PID/D");

    // Kaon -
    treeKaonM.Branch("E", &kaonm_E, "E/D");
    treeKaonM.Branch("px", &kaonm_px, "px/D");
    treeKaonM.Branch("py", &kaonm_py, "py/D");
    treeKaonM.Branch("pz", &kaonm_pz, "pz/D");
    treeKaonM.Branch("Mom", &kaonm_Ph, "Mom/D");
    treeKaonM.Branch("W", &kaonm_W, "W/D");
    treeKaonM.Branch("Q2", &kaonm_Q2, "Q2/D");
    treeKaonM.Branch("xF", &kaonm_xF, "xF/D");
    treeKaonM.Branch("xB", &kaonm_xB, "xB/D");
    treeKaonM.Branch("y", &kaonm_y, "y/D");
    treeKaonM.Branch("z", &kaonm_z, "z/D");
    treeKaonM.Branch("Pt", &kaonm_Pt, "Pt/D");
    treeKaonM.Branch("PhT", &kaonm_PhT, "PhT/D");
    treeKaonM.Branch("Phi_Lab", &kaonm_Phi, "Phi_Lab/D");
    treeKaonM.Branch("Theta_Lab", &kaonm_Theta, "Theta_Lab/D");
    treeKaonM.Branch("Pseudorapidity", &kaonm_eta, "Pseudorapidity/D");
    treeKaonM.Branch("Phi_h", &kaonm_Phi_h, "Phi_h/D");
    treeKaonM.Branch("SinPhi", &kaonm_SinPhi, "SinPhi/D");

    // KAON+ PLOT
    dirKaonp->cd();
    // Mom
    TH2D kp_MomVsPhT ("_MomVsPhT", "Correlation Mom vs PhT  |  Kaon+  |; PhT [GeV]; Mom [GeV]", 80, 0, 1, 80, 0.9, 7);
    TH2D kp_MomVsXb ("_MomVsXb", "Correlation Mom vs x_B  |  Kaon+  |; x_B; Mom [GeV]", 80, 0, 1, 80, 0.9, 7);
    TH2D kp_MomVsXf ("_MomVsXf", "Correlation Mom vs x_F  |  Kaon+  |; x_F; Mom [GeV]", 80, 0, 0.6, 80, 0.9, 7);
    TH2D kp_MomVsZ ("_MomVsZ", "Correlation Mom vs Z  |  Kaon+  |; z; Mom [GeV]", 80, 0.2, 0.9, 80, 0.9, 7);
    TH2D kp_MomVsY ("_MomVsY", "Correlation Mom vs Y  |  Kaon+  |; y; Mom [GeV]", 80, 0.2, 0.75, 80, 0.9, 7);
    TH2D kp_MomVsEta ("_MomVsEta", "Correlation Mom vs Eta  |  Kaon+  |; Eta; Mom [GeV]", 80, 1.5, 3.0, 80, 0.9, 7);
    TH2D kp_MomVsTheta ("_MomVsTheta", "Correlation Mom vs Theta  |  Kaon+  |; Theta [Rad]; Mom [GeV]", 80, 0, 0.2*TMath::Pi(), 80, 0.9, 7);
    TH2D kp_MomVsPhi_h ("_MomVsPhi_h", "Correlation Mom vs Phi_h  |  Kaon+  |; Phi_h [Rad]; Mom [GeV]", 80, -TMath::Pi(), TMath::Pi(), 80, 0.9, 7);
    // Q2
    TH2D kp_Q2VsXb ("_Q2VsXb", "Correlation Q^2 vs x_B  |  Kaon+  |; x_B; Q^2 [GeV^2]", 80, 0, 1, 80, 1, 8);
    TH2D kp_Q2VsXf ("_Q2VsXf", "Correlation Q^2 vs x_F  |  Kaon+  |; x_F; Q^2 [GeV^2]", 80, 0, 0.6, 80, 1, 7);
    TH2D kp_Q2VsMom ("_Q2VsMom", "Correlation Q^2 vs Mom  |  Kaon+  |; Mom [GeV]; Q^2 [GeV^2]", 80, 0.9, 6, 80, 1, 7);
    TH2D kp_Q2VsPhT ("_Q2VsPhT", "Correlation Q^2 vs PhT  |  Kaon+  |; PhT [GeV]; Q^2 [GeV^2]", 80, 0, 1.2, 80, 1, 6);
    TH2D kp_Q2VsZ ("_Q2VsZ", "Correlation Q^2 vs Z  |  Kaon+  |; z; Q^2 [GeV^2]", 80, 0.2, 0.9, 80, 1, 6);
    TH2D kp_Q2VsY ("_Q2VsY", "Correlation Q^2 vs Y  |  Kaon+  |; y; Q^2 [GeV^2]", 80, 0.2, 0.75, 80, 0.9, 7);
    TH2D kp_Q2VsEta ("_Q2VsEta", "Correlation Q^2 vs Eta  |  Kaon+  |; Eta; Q^2 [GeV^2]", 80, 1.5, 3.0, 80, 1, 6);
    TH2D kp_Q2VsPhi_h ("_Q2VsPhi_h", "Correlation Q^2 vs Phi_h  |  Kaon+  |; Phi_h [Rad]; Q^2 [GeV^2]", 80, -TMath::Pi(), TMath::Pi(), 80, 1, 6);
    // PhT
    TH2D kp_PhTvsZ ("_PhTvsZ", "Correlation PhT vs Z  |  Kaon+  |; z; PhT [GeV]", 80, 0.2, 1, 80, 0, 1.2);
    TH2D kp_PhTvsXb ("_PhTvsXb", "Correlation PhT vs Xb  |  Kaon+  |; x_B; PhT [GeV]", 80, 0, 1, 80, 0, 1.2);
    TH2D kp_PhTvsEta ("_PhTvsEta", "Correlation PhT vs Eta  |  Kaon+  |; Eta; PhT [GeV]", 80, 1.5, 3.0, 80, 0, 1.2);
    TH2D kp_PhTvsPhi_h ("_PhTvsPhi_h", "Correlation PhT vs Phi_h  |  Kaon+  |; Phi_h [Rad]; PhT [GeV]", 80, -TMath::Pi(), TMath::Pi(), 80, 0, 1.2);
    // Z
    TH2D kp_zVsXb ("_zVsXb", "Correlation Z vs x_B  |  Kaon+  |; x_B; z", 80, 0, 1, 80, 0.2, 0.9);
    TH2D kp_zVsXf ("_zVsXf", "Correlation Z vs x_F  |  Kaon+  |; x_F; z", 80, 0, 0.6, 80, 0.2, 0.9);
    TH2D kp_zVsEta ("_zVsEta", "Correlation Z vs Eta  |  Kaon+  |; Eta; z", 80, 1.5, 3.0, 80, 0.2, 0.9);
    TH2D kp_zVsPhi_h ("_zVsPhi_h", "Correlation Z vs Phi_h  |  Kaon+  |; Phi_h [Rad]; z", 80, -TMath::Pi(), TMath::Pi(), 80, 0.2, 0.9);
    // Angles
    TH2D kp_ThetaVsPhi_h ("_ThetaVsPhi_h", "Correlation Theta vs Phi_h  |  Kaon+  |; Phi_h [Rad]; Theta [Rad]", 80, -TMath::Pi(), TMath::Pi(), 80, 0, 0.2*TMath::Pi());
    TH2D el_ThetaVsMom ("el_ThetaVsMom", "Correlation Theta vs Mom  |  electron  |; Mom [GeV]; Theta [Deg]", 80, 1.5, 9, 80, 5, 40);
    TH2D kp_MomVsM2x ("_MomVsM2x", "correlation Mom vs M2x | Kaon+ |; Mom [GeV]; M2x[GeV^2]", 60, 1, 7, 60, 2.3, 11);
    // KAON- PLOT
    dirKaonm->cd();
    // Mom
    TH2D km_MomVsPhT ("_MomVsPhT", "Correlation Mom vs PhT  |  Kaon-  |; PhT [GeV]; Mom[GeV]", 60, 0, 1, 60, 0.9, 6);
    TH2D km_MomVsXb ("_MomVsXb", "Correlation Mom vs x_B  |  Kaon-  |; x_B; Mom[GeV]", 60, 0, 1, 60, 0.9, 6);
    TH2D km_MomVsXf ("_MomVsXf", "Correlation Mom vs x_F  |  Kaon-  |; x_F; Mom[GeV]", 60, 0, 0.6, 60, 0.9, 6);
    TH2D km_MomVsZ ("_MomVsZ", "Correlation Mom vs Z  |  Kaon-  |; z; Mom[GeV]", 60, 0.2, 0.9, 60, 0.9, 6);
    TH2D km_MomVsY ("_MomVsY", "Correlation Mom vs Y  |  Kaon-  |; y; Mom[GeV]", 60, 0.2, 0.75, 60, 0.9, 6);
    TH2D km_MomVsEta ("_MomVsEta", "Correlation Mom vs Eta  |  Kaon-  |; Eta; Mom[GeV]", 60, 1.5, 3, 60, 0.9, 6);
    TH2D km_MomVsTheta ("_MomVsTheta", "Correlation Mom vs Theta  |  Kaon-  |; Theta [Rad]; Mom[GeV]", 60, 0, 0.2*TMath::Pi(), 60, 0.9, 6);
    TH2D km_MomVsPhi_h ("_MomVsPhi_h", "Correlation Mom vs Phi_h  |  Kaon-  |; Phi_h [Rad]; Mom[GeV]", 60, -TMath::Pi(), TMath::Pi(), 60, 0.9, 6);
    // Q2
    TH2D km_Q2VsXb ("_Q2VsXb", "Correlation Q^2 vs x_B  |  Kaon-  |; x_B; Q^2 [GeV^2]", 60, 0, 1, 60, 1, 6);
    TH2D km_Q2VsXf ("_Q2VsXf", "Correlation Q^2 vs x_F  |  Kaon-  |; x_F; Q^2 [GeV^2]", 60, 0, 0.6, 60, 1, 6);
    TH2D km_Q2VsMom ("_Q2VsMom", "Correlation Q^2 vs Mom  |  Kaon-  |; Mom [GeV]; Q^2 [GeV^2]", 60, 0.9, 6, 60, 1, 6);
    TH2D km_Q2VsPhT ("_Q2VsPhT", "Correlation Q^2 vs PhT  |  Kaon-  |; PhT [GeV]; Q^2 [GeV^2]", 60, 0, 1.1, 60, 1, 6);
    TH2D km_Q2VsZ ("_Q2VsZ", "Correlation Q^2 vs Z  |  Kaon-  |; z; Q^2 [GeV^2]", 60, 0.2, 0.9, 60, 1, 6);
    TH2D km_Q2VsY ("_Q2VsY", "Correlation Q^2 vs Y  |  Kaon-  |; y; Q^2 [GeV^2]", 80, 0.2, 0.75, 80, 0.9, 6);
    TH2D km_Q2VsEta ("_Q2VsEta", "Correlation Q^2 vs Eta  |  Kaon-  |; Eta; Q^2 [GeV^2]", 60, 1.5, 3, 60, 1, 6);
    TH2D km_Q2VsPhi_h ("_Q2VsPhi_h", "Correlation Q^2 vs Phi_h  |  Kaon-  |; Phi_h [Rad]; Q^2 [GeV^2]", 60, -TMath::Pi(), TMath::Pi(), 60, 1, 6);
    // PhT
    TH2D km_PhTvsZ ("_PhTvsZ", "Correlation PhT vs Z  |  Kaon-  |; z; PhT [GeV]", 60, 0.2, 1, 60, 0, 1.1);
    TH2D km_PhTvsXb ("_PhTvsXb", "Correlation PhT vs Xb  |  Kaon-  |; x_B; PhT [GeV]", 60, 0, 1, 60, 0, 1.1);
    TH2D km_PhTvsEta ("_PhTvsEta", "Correlation PhT vs Eta  |  Kaon-  |; Eta; PhT [GeV]", 60, 1.5, 3, 60, 0, 1.1);
    TH2D km_PhTvsPhi_h ("_PhTvsPhi_h", "Correlation PhT vs Phi  |  Kaon-  |; Phi_h [Rad]; PhT [GeV]", 60, -TMath::Pi(), TMath::Pi(), 60, 0, 1.1);
    // Z
    TH2D km_zVsXb ("_zVsXb", "Correlation Z vs x_B  |  Kaon-  |; x_B; z", 60, 0, 1, 60, 0.2, 0.9);
    TH2D km_zVsXf ("_zVsXf", "Correlation Z vs x_F  |  Kaon-  |; x_F; z", 60, 0, 0.6, 60, 0.2, 0.9);
    TH2D km_zVsEta ("_zVsEta", "Correlation Z vs Eta  |  Kaon-  |; Eta; z", 60, 1.5, 3, 60, 0.2, 0.9);
    TH2D km_zVsPhi_h ("_zVsPhi_h", "Correlation Z vs Phi_h  |  Kaon-  |; Phi_h [Rad]; z", 60, -TMath::Pi(), TMath::Pi(), 60, 0.2, 0.9);
    // Angles
    TH2D km_ThetaVsPhi_h ("_ThetaVsPhi_h", "Correlation Theta vs Phi_h  |  Kaon-  |; Phi_h [Rad]; Theta [Rad]", 60, -TMath::Pi(), TMath::Pi(), 60, 0, 0.2*TMath::Pi());

    // KAON + BIN
    dirKaonp_Bin->cd();
    // Bin 1
    TH1D kp_Phi_z1 ("_Phi_z_1", "Phi_h (z < 0.4)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_PhT1 ("_Phi_PhT_1", "Phi_h (PhT < 0.25 GeV)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_Q21 ("_Phi_Q2_1", "Phi_h (Q^2 < 2 GeV^2)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_xB1 ("_Phi_xB_1", "Phi_h (x_B < 0.15)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    // Bin 2
    TH1D kp_Phi_z2 ("_Phi_z_2", "Phi_h (0.4 < z < 0.5)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_PhT2 ("_Phi_PhT_2", "Phi_h (0.25 < PhT < 0.5 GeV)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_Q22 ("_Phi_Q2_2", "Phi_h (2 < Q^2 < 3 GeV^2)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_xB2 ("_Phi_xB_2", "Phi_h (0.15 < x_B < 0.25)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    // Bin 3
    TH1D kp_Phi_z3 ("_Phi_z_3", "Phi_h (0.5 < z < 0.65)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_PhT3 ("_Phi_PhT_3", "Phi_h (0.5 < PhT < 0.7 GeV)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_Q23 ("_Phi_Q2_3", "Phi_h (3 < Q^2 < 4 GeV^2)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_xB3 ("_Phi_xB_3", "Phi_h (0.25 < x_B < 0.35)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    // Bin 4
    TH1D kp_Phi_z4 ("_Phi_z_4", "Phi_h (z > 0.65)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_PhT4 ("_Phi_PhT_4", "Phi_h (PhT > 0.7 GeV)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_Q24 ("_Phi_Q2_4", "Phi_h (Q^2 > 4 GeV^2)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());
    TH1D kp_Phi_xB4 ("_Phi_xB_4", "Phi_h (x_B > 0.35)  |  Kaon+ ; #Phi_h [rad]; Counts", 60, -TMath::Pi(), TMath::Pi());

    // controlla se è disponibile l'istanza Run Conditions DataBase e accede ai metadati di RCDB
    if (config_c12->rcdb()) {
        cout << "rcdb trovato correttamente" << endl;
        auto& rcdbData = config_c12->rcdb()->current();
        //cout << "Event count: " << rcdbData.event_count << endl;
        //cout << "Beam energy: " << rcdbData.beam_energy << " GeV" << endl;
        //cout << "Beam current: " << rcdbData.beam_current << " μA" << endl;
    }
    // Calibration Constants DataBase, e electron_sf contiene informazioni delle frazioni energetiche degli elettroni nei calorimentri
    // copia e incollata così, è una tabella di calibrazione

    if (config_c12->ccdb()) {
        cout << "ccdb trovato correttamente" << endl;
        auto& ccdbElSF = config_c12->ccdb()->requestTableDoubles("/calibration/eb/electron_sf");
        //cout << "CCDB Electron Sampling Fraction (first two rows):" << endl;
        //cout << "  Row 0, Col 5: " << ccdbElSF[0][5] << endl;
        //cout << "  Row 1, Col 5: " << ccdbElSF[1][5] << endl;
    }
    
    auto elRegion=FD;
    gBenchmark->Start("db");
    size_t eventCount = 0;
    const size_t maxEvents = 1e8;
    while (chain.Next()){
        eventCount++;
        if (eventCount >= maxEvents) {
            break;  
        }
        auto electrons = c12->getByID(11);  // Recupera tutte le particelle con PID = 11 
        int e_PID = 11;
        auto protons = c12->getByID(2212); 
        auto pionp = c12->getByID(211);
        auto kaonp = c12->getByID(321);
        int kp_PID = 321;
        auto kaonm = c12->getByID(-321);
        int km_PID = -321;
        auto virtual_gamma = c12->getByID(22);
        auto N_run = c12->runconfig()->getRun();
        auto torus = c12->runconfig()->getTorus();
        if(InvertHelicity(N_run)) helicity = -c12->event()->getHelicity();
        else helicity = c12->event()->getHelicity();
        for (auto& e : electrons) {
            SetLorentzVector(el, e);  // filled with the TLorentzVector that I want to fill and with the particle required
            TLorentzVector Elec = el;
            electron_px = e->getPx();
            electron_py = e->getPy();
            electron_pz = e->getPz();
            electron_E = el.E(); 
            // Math.pow(Math.pow(particle_mass(2212),2)+2*particle_mass(2212)*nu - Q2, 0.5);
            TLorentzVector q = beam - el;           // Virtual photon
            electron_Q2 = -q.M2();
            gamma_nu = beam.E() - el.E();
            electron_W = pow(pow(pr.M(),2)+2*pr.M()*gamma_nu - electron_Q2, 0.5);
            electron_vz = e->par()->getVz();
            electron_mom = el.P();
            electron_Nphe = e->che(HTCC)->getNphe();  //number of photoelectron
            electron_status = e->getStatus();      
            electron_sector = e->getSector();
            electron_PCAL = e->cal(PCAL)->getEnergy();    
            electron_ECAL = e->cal(ECAL)->getEnergy();
            electron_ECIN = e->cal(ECIN)->getEnergy();
            electron_CAL_Tot = electron_ECAL + electron_PCAL + electron_ECIN;
            electron_edge1 = e->traj(6,6)->getEdge();
            electron_edge2 = e->traj(6,18)->getEdge();
            electron_edge3 = e->traj(6,36)->getEdge();
            // primo numero, detector, secondo il layer, 6 = Drift Chamber
            //cal(PCAL)
            electron_Theta = el.Theta();
            electron_ThetaDeg = (electron_Theta * 180)/(TMath::Pi());
            // SCATTERED ELECTRON PID
            if(ElectronPID_ForwardDetector(electron_status) && ElectronPID_ForwardDetector(electron_status) && ElectronPID_PCAL(electron_PCAL) &&
            ElectronPID_CalSFcut(electron_sector, N_run, electron_mom, electron_CAL_Tot) && ElectronPID_DiagonalCut(e, &Elec) && 
            KinematicPID_W(electron_W) && ElectronPID_VertexCut(electron_vz) && ElectronPID_DC(electron_edge1, electron_edge2, electron_edge3, torus) &&
            KinematicPID_Vertex(electron_vz, torus, e_PID)){
                el_ThetaVsMom.Fill(electron_mom, electron_ThetaDeg);
                //treeEl.Fill();
                /*
                for (auto& p : protons) {
                    SetLorentzVector(pr, p);   
                }
                // vitual photon
                for(auto& gam : virtual_gamma) {
                    SetLorentzVector(gm, gam);
                    gamma_px = gam->getPx();
                    gamma_py = gam->getPy();
                    gamma_pz = gam->getPz();
                    treePh.Fill();
                }
                // Pion+
                for(auto& Pip : pionp) {
                    SetLorentzVector(pip, Pip);
                    pionp_px = Pip->getPx();
                    pionp_py = Pip->getPy();
                    pionp_pz = Pip->getPz();
                    treePionP.Fill();
                }
                */
                // Kaon+
                for(auto& Kp : kaonp) {
                    SetLorentzVector(kp, Kp);
                    TLorentzVector P_h = kp;
                    TVector3 qVec = q.Vect();
                    kaonp_s = (beam + target).M2();         // center of mass energy
                    //kaonp_s = (2*0.938*10.6 + 0.938*0.938);
                    double nu = beam.E() - el.E();          // Virtual photon energy
                    beta_CMS = Lab.BoostVector();
                    TLorentzVector Ph_CMS = P_h;
                    Ph_CMS.Boost(-beta_CMS);
                    kaonp_y = target.Dot(q) / target.Dot(beam);
                    TLorentzVector M2x = (Lab - el - kp);
                    kaonp_M2x = M2x.M2();
                    kaonp_chi2pid = Kp->getChi2Pid();
                    if(kaonp_y <= 0.75){
                        kaonp_Ph = P_h.P();
                        kaonp_E = P_h.E();
                        kaonp_W = sqrt((kaonp_E + electron_E)*(kaonp_E + electron_E) - (kaonp_Ph + electron_mom)*(kaonp_Ph + electron_mom));
                        kaonp_px = Kp->getPx();
                        kaonp_py = Kp->getPy();
                        kaonp_pz = Kp->getPz();
                        kaonp_Q2 = -q.M2();                              // Q^2
                        kaonp_xB = kaonp_Q2 / (2 * target.Dot(q));       // x_B     
                        kaonp_z = target.Dot(P_h) / target.Dot(q);     
                        kaonp_Pt = P_h.Perp(q.Vect());                 
                        kaonp_Phi = P_h.Phi();
                        kaonp_eta = P_h.PseudoRapidity();
                        kaonp_Theta = P_h.Theta();
                        kaonp_xF = (2 * Ph_CMS.Pz()) / sqrt(kaonp_s);         
                        // Trento Convention - Scattering Plane Axis calculation
                        TVector3 zAxis = q.Vect().Unit();           // virtual photon direction
                        TVector3 l_vect = beam.Vect();              // spatial component of the beam lepton
                        TVector3 sl_vect = el.Vect();               // spatial component of the scattered lepton           
                        TVector3 yAxis = l_vect.Cross(sl_vect).Unit();
                        TVector3 xAxis = yAxis.Cross(zAxis);        // x-Axis of the scattering plane
                        TVector3 Ph_T = P_h.Vect() - (P_h.Vect() * zAxis) * zAxis; 
                        kaonp_PhT = Ph_T.Mag();                     // just a check to see if the calculated transverse momenta are correct
                        kaonp_Ph_x = Ph_T.Dot(xAxis);
                        kaonp_Ph_y = Ph_T.Dot(yAxis);
                        kaonp_edge1 = Kp->traj(6,6)->getEdge();
                        kaonp_edge2 = Kp->traj(6,18)->getEdge();
                        kaonp_edge3 = Kp->traj(6,36)->getEdge();
                        kaonp_vz = Kp->par()->getVz();
                        kaonp_SinPhi = TMath::Sin(kaonp_Phi_h);
                        kaonp_Phi_h = TMath::ATan2(kaonp_Ph_y, kaonp_Ph_x); 
                        kaonp_rich_Id = Kp->rich()->getId();
                        if(Kp->rich()->getBest_PID() == 321) kaonp_rich_PID = Kp->rich()->getBest_PID();
                        if(HadronPID_Reson(kaonp_Ph, kaonp_M2x) && HadronPID_Q2(kaonp_Q2) && KinematicPID_xF(kaonp_xF) &&
                        HadronPID_DC(kaonp_edge1, kaonp_edge2, kaonp_edge3, torus) && KinematicPID_Vertex(kaonp_vz, torus, kp_PID) &&
                        HadronPID_Chi2Pid(kaonp_chi2pid) && HadronPID_z(kaonp_z)){
                            // normalization of Phi_h inside -pi,pi
                            /*
                            if(kaonp_Phi_h > TMath::Pi()){
                            kaonp_Phi_h -= 2*TMath::Pi();
                            }
                            else if(kaonp_Phi_h < -TMath::Pi()){
                                kaonp_Phi_h += 2*TMath::Pi();
                            }
                            */
                            // Histograms
                            // z
                            if(kaonp_z < 0.4) kp_Phi_z1.Fill(kaonp_Phi_h);
                            else if(kaonp_z >= 0.4 && kaonp_z <0.5) kp_Phi_z2.Fill(kaonp_Phi_h);
                            else if(kaonp_z >= 0.5 && kaonp_z <0.65) kp_Phi_z3.Fill(kaonp_Phi_h);
                            else if(kaonp_z >= 0.65) kp_Phi_z4.Fill(kaonp_Phi_h);
                            // PhT
                            if(kaonp_PhT < 0.25) kp_Phi_PhT1.Fill(kaonp_Phi_h);
                            else if(kaonp_PhT >= 0.25 && kaonp_PhT <0.5) kp_Phi_PhT2.Fill(kaonp_Phi_h);
                            else if(kaonp_PhT >= 0.5 && kaonp_PhT <0.7) kp_Phi_PhT3.Fill(kaonp_Phi_h);
                            else if(kaonp_PhT >= 0.7) kp_Phi_PhT4.Fill(kaonp_Phi_h);
                            // Q2
                            if(kaonp_Q2 < 2) kp_Phi_Q21.Fill(kaonp_Phi_h);
                            else if(kaonp_Q2 >= 2 && kaonp_Q2 < 3) kp_Phi_Q22.Fill(kaonp_Phi_h);
                            else if(kaonp_Q2 >= 3 && kaonp_Q2 < 4) kp_Phi_Q23.Fill(kaonp_Phi_h);
                            else if(kaonp_Q2 >= 4) kp_Phi_Q24.Fill(kaonp_Phi_h);
                            // xB
                            if(kaonp_xB < 0.15) kp_Phi_xB1.Fill(kaonp_Phi_h);
                            else if(kaonp_xB >= 0.15 && kaonp_xB <0.25) kp_Phi_xB2.Fill(kaonp_Phi_h);
                            else if(kaonp_xB >= 0.25 && kaonp_xB <0.35) kp_Phi_xB3.Fill(kaonp_Phi_h);
                            else if(kaonp_xB >= 0.35) kp_Phi_xB4.Fill(kaonp_Phi_h);
                            // Mom
                            kp_MomVsPhT.Fill(kaonp_PhT, kaonp_Ph);
                            kp_MomVsXb.Fill(kaonp_xB, kaonp_Ph);
                            kp_MomVsXf.Fill(kaonp_xF, kaonp_Ph);
                            kp_MomVsPhi_h.Fill(kaonp_Phi_h, kaonp_Ph);
                            kp_MomVsEta.Fill(kaonp_eta, kaonp_Ph);
                            kp_MomVsZ.Fill(kaonp_z, kaonp_Ph);
                            kp_MomVsY.Fill(kaonp_y, kaonp_Ph);
                            kp_MomVsTheta.Fill(kaonp_Theta, kaonp_Ph);
                            // Q^2
                            kp_Q2VsEta.Fill(kaonp_eta, kaonp_Q2);
                            kp_Q2VsXf.Fill(kaonp_xF, kaonp_Q2);
                            kp_Q2VsMom.Fill(kaonp_Ph, kaonp_Q2);
                            kp_Q2VsPhi_h.Fill(kaonp_Phi_h, kaonp_Q2);
                            kp_Q2VsPhT.Fill(kaonp_PhT, kaonp_Q2);
                            kp_Q2VsXb.Fill(kaonp_xB, kaonp_Q2);
                            kp_Q2VsZ.Fill(kaonp_z, kaonp_Q2);
                            kp_Q2VsY.Fill(kaonp_y, kaonp_Q2);
                            // PhT
                            kp_PhTvsEta.Fill(kaonp_eta, kaonp_PhT);
                            kp_PhTvsXb.Fill(kaonp_xB, kaonp_PhT);
                            kp_PhTvsZ.Fill(kaonp_z, kaonp_PhT);
                            kp_PhTvsPhi_h.Fill(kaonp_Phi_h, kaonp_PhT);
                            // z
                            kp_zVsEta.Fill(kaonp_eta, kaonp_z);
                            kp_zVsXf.Fill(kaonp_xF, kaonp_z);
                            kp_zVsPhi_h.Fill(kaonp_Phi_h, kaonp_z);
                            kp_zVsXb.Fill(kaonp_xB, kaonp_z);
                            kp_MomVsM2x.Fill(kaonp_Ph, kaonp_M2x);
                            // Angles
                            kp_ThetaVsPhi_h.Fill(kaonp_Phi_h, kaonp_Theta);
                            treeKaonP.Fill();  // Riempie solo i valori di z nel range desiderato
                            if(helicity < 0) N_down += 1;
                            if(helicity > 0) N_up += 1;
                        }
                    }
                }
                // Kaon-
                for(auto& Km : kaonm) {
                    SetLorentzVector(km, Km);
                    TLorentzVector P_h = km;
                    TLorentzVector q = beam - el;           // Virtual photon
                    TVector3 qVec = q.Vect();
                    kaonm_s = (beam + target).M2();         // center of mass energy
                    //beta_CMS = ((beam.P() + target.P()))/(beam.E() + target.E());
                    //kaonm_s = (2*0.938*10.6 + 0.938*0.938);
                    double nu = beam.E() - el.E();          // Virtual photon energy
                    kaonm_y = target.Dot(q) / target.Dot(beam);
                    beta_CMS = Lab.BoostVector();
                    TLorentzVector Ph_CMS = P_h;
                    Ph_CMS.Boost(-beta_CMS);
                    TLorentzVector M2x = (Lab - el - km);
                    kaonm_chi2pid = Km->getChi2Pid();
                    kaonm_M2x = M2x.M2();
                    if(kaonm_y <= 0.75){
                        kaonm_Ph = P_h.P();
                        kaonm_E = P_h.E();
                        kaonm_W = P_h.M();
                        kaonm_px = Km->getPx();
                        kaonm_py = Km->getPy();
                        kaonm_pz = Km->getPz();
                        kaonm_Q2 = -q.M2();                              // Q^2
                        kaonm_xB = kaonm_Q2 / (2 * target.Dot(q));       // x_B     
                        kaonm_z = target.Dot(P_h) / target.Dot(q);       
                        kaonm_Pt = P_h.Perp(q.Vect());                 
                        kaonm_Phi = P_h.Phi();
                        kaonm_eta = P_h.PseudoRapidity();
                        kaonm_Theta = P_h.Theta();
                        kaonm_xF = (2 * Ph_CMS.Pz()) / sqrt(kaonm_s);         
                        // Trento Convention - Scattering Plane Axis calculation
                        TVector3 zAxis = q.Vect().Unit();           // virtual photon direction
                        TVector3 l_vect = beam.Vect();              // spatial component of the beam lepton
                        TVector3 sl_vect = el.Vect();               // spatial component of the scattered lepton           
                        TVector3 yAxis = l_vect.Cross(sl_vect).Unit();
                        TVector3 xAxis = yAxis.Cross(zAxis);        // x-Axis of the scattering plane
                        TVector3 Ph_T = P_h.Vect() - (P_h.Vect() * zAxis) * zAxis; 
                        kaonm_PhT = Ph_T.Mag();                     // just a check to see if the calculated transverse momenta are correct
                        kaonm_Ph_x = Ph_T.Dot(xAxis);
                        kaonm_Ph_y = Ph_T.Dot(yAxis);
                        kaonm_Phi_h = TMath::ATan2(kaonm_Ph_y, kaonm_Ph_x); 
                        kaonm_edge1 = Km->traj(6,6)->getEdge();
                        kaonm_edge2 = Km->traj(6,18)->getEdge();
                        kaonm_edge3 = Km->traj(6,36)->getEdge();
                        kaonm_vz = Km->par()->getVz();
                        kaonm_SinPhi = TMath::Sin(kaonm_Phi_h);
                        if(HadronPID_Reson(kaonm_Ph, kaonm_M2x) && HadronPID_Q2(kaonm_Q2) && KinematicPID_xF(kaonm_xF) && 
                        HadronPID_DC(kaonm_edge1, kaonm_edge2, kaonm_edge3, torus) && KinematicPID_Vertex(kaonm_vz, torus, km_PID) &&
                        HadronPID_Chi2Pid(kaonm_chi2pid) && HadronPID_z(kaonm_z)){
                            // normalization of Phi_h inside -pi,pi
                            /*
                            if(kaonm_Phi_h > TMath::Pi()){
                              kaonm_Phi_h -= 2*TMath::Pi();
                            }
                            else if(kaonm_Phi_h < -TMath::Pi()){
                              kaonm_Phi_h += 2*TMath::Pi();
                            }
                            */
                            // Histograms
                            // Mom
                            km_MomVsPhT.Fill(kaonm_PhT, kaonm_Ph);
                            km_MomVsXb.Fill(kaonm_xB, kaonm_Ph);
                            km_MomVsXf.Fill(kaonm_xF, kaonm_Ph);
                            km_MomVsPhi_h.Fill(kaonm_Phi_h, kaonm_Ph);
                            km_MomVsEta.Fill(kaonm_eta, kaonm_Ph);
                            km_MomVsZ.Fill(kaonm_z, kaonm_Ph);
                            km_MomVsY.Fill(kaonm_y, kaonm_Ph);
                            km_MomVsTheta.Fill(kaonm_Theta, kaonm_Ph);
                            // Q^2
                            km_Q2VsEta.Fill(kaonm_eta, kaonm_Q2);
                            km_Q2VsXf.Fill(kaonm_xF, kaonm_Q2);
                            km_Q2VsMom.Fill(kaonm_Ph, kaonm_Q2);
                            km_Q2VsPhi_h.Fill(kaonm_Phi_h, kaonm_Q2);
                            km_Q2VsPhT.Fill(kaonm_PhT, kaonm_Q2);
                            km_Q2VsXb.Fill(kaonm_xB, kaonm_Q2);
                            km_Q2VsZ.Fill(kaonm_z, kaonm_Q2);
                            km_Q2VsY.Fill(kaonm_y, kaonm_Q2);
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
                            // Angles
                            km_ThetaVsPhi_h.Fill(kaonm_Phi_h, kaonm_Theta);
                            treeKaonM.Fill();  // Riempie solo i valori di z nel range desiderato
                            if(helicity < 0) N_down += 1;
                            if(helicity > 0) N_up += 1;
                        }
                    }
                }
            }
        }
    }
    //treeKaonP.Draw("Mom/PhT", "", "colz");
    //kp_MomVsPhT.Draw("colz");
    outFile.Write();
    outFile.Close();
    gBenchmark->Stop("db");
    gBenchmark->Print("db");
    double A_LU = (N_up - N_down) / (N_up + N_down);
    //cout << "N_up = " << N_up << " --- N_down = " << N_down << endl;
    //cout << "A_LU = " << A_LU << endl;
    //cout << "Processed " << eventCount << " events from file: " << filePath << endl;
    cout << "ROOT output file: output.root" << endl;
}
