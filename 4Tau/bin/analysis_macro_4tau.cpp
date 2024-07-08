#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <string> 
#include <map>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "HtoAA/Utilities/interface/Config.h"
#include "HtoAA/Utilities/interface/RoccoR.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "HtoAA/Utilities/interface/json.h"
#include "HtoAA/Utilities/interface/Jets.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "HtoAA/Utilities/interface/PileUp.h"
#include "HtoAA/Utilities/interface/functions.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "CondFormats/BTauObjects/interface/BTagEntry.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "TBranch.h"

using namespace std;

float getEffectiveArea(float eta) {
  float effArea =  0.1440;
  float absEta = fabs(eta);
  if (absEta<1.0) effArea = 0.1440;
  else if (absEta < 1.4790) effArea = 0.1562;
  else if (absEta < 2.0) effArea = 0.1032;
  else if (absEta < 2.2) effArea = 0.0859;
  else if (absEta < 2.3) effArea = 0.1116;
  else if (absEta < 2.4) effArea = 0.1321;
  else if (absEta < 5.0) effArea = 0.1654;
  return effArea;
} 

bool passedFilters(std::map<string,int> * flags, std::vector<TString> met_filters) {
  bool passed = true; 
  for (std::map<string,int>::iterator it = flags->begin(); it != flags->end(); ++it) {
    TString Flag(it->first);
    for (auto met_filter : met_filters) {
      if (met_filter==Flag && it->second == 0) {
	passed = false;
	break;
      }
    }
  }
  return passed;

}

float impactParameter(float * v0, float * v, float * p, float &dxy, float &dz) {

  float mod = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);

  float n[3];
  n[0] = p[0]/mod;
  n[1] = p[1]/mod;
  n[2] = p[2]/mod;

  float t = 
    (v0[0]-v[0])*n[0] + 
    (v0[1]-v[1])*n[1] +
    (v0[2]-v[2])*n[2];

  float vert[3];
  vert[0] = v[0] + n[0]*t;
  vert[1] = v[1] + n[1]*t;
  vert[2] = v[2] + n[2]*t;



  float ip = TMath::Sqrt(
			 (vert[0]-v0[0])*(vert[0]-v0[0])+
			 (vert[1]-v0[1])*(vert[1]-v0[1])+
			 (vert[2]-v0[2])*(vert[2]-v0[2])
			 );


  t = (v0[0]-v[0])*n[0] + (v0[1]-v[1])*n[1];
  vert[0] = v[0] + n[0]*t;
  vert[1] = v[1] + n[1]*t;
  vert[2] = v[2] + n[2]*t;

  dxy = TMath::Sqrt(
		    (vert[0]-v0[0])*(vert[0]-v0[0])+
		    (vert[1]-v0[1])*(vert[1]-v0[1])
		    );

  dz = abs(vert[2]-v0[2]);

  return ip;

}

int main(int argc, char * argv[]) {
  
  if (argc<2) {
    std::cout << "Usage of the program : Hto4TausAnalysis [file_list]" << std::endl;
    std::cout << "file_list : file list of RooT files to be processed. To run on Data the string has to include the string \"Data\"." << std::endl;
    exit(1);
  }


  // **** configuration
  Config cfg(argv[1]);

  const bool isData = cfg.get<bool>("IsData");
  const int era = cfg.get<int>("Era");
  const bool applyHiggsPtWeight = cfg.get<bool>("ApplyHiggsPtWeight");

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
  const float etaMuonLowCut  = cfg.get<float>("etaMuonLowCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");

  // topological cuts
  const float dRMuonsCut   = cfg.get<float>("dRMuonsCut");
  const bool sameSign      = cfg.get<bool>("SameSignMuons");

  // track selection
  const float dRIsoMuon       = cfg.get<float>("dRIsoMuon");
  const float ptTrkLooseCut   = cfg.get<float>("ptTrkLooseCut");
  const float ptTrkCut        = cfg.get<float>("ptTrkCut");
  const float etaTrkCut       = cfg.get<float>("etaTrkCut");
  const float dxyTrkLooseCut  = cfg.get<float>("dxyTrkLooseCut");
  const float dxyTrkCut       = cfg.get<float>("dxyTrkCut");
  const float dzTrkLooseCut   = cfg.get<float>("dzTrkLooseCut");
  const float dzTrkCut        = cfg.get<float>("dzTrkCut");

  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");
  const string jsonFile = cfg.get<string>("jsonFile");

  // trigger
  const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");
  const string dimuonTriggerName = cfg.get<string>("DiMuonTriggerName");
  const string muonHighPtFilterName = cfg.get<string>("MuonHighPtFilterName");
  const string muonLowPtFilterName1 = cfg.get<string>("MuonLowPtFilterName1");
  const string muonLowPtFilterName2 = cfg.get<string>("MuonLowPtFilterName2");
  const string dimuonDzFilterName = cfg.get<string>("DimuonDzFilterName");
  const string dimuonSameSignFilterName = cfg.get<string>("DimuonSameSignFilterName");
  // trigger matching
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 
  const float effDzSS        = cfg.get<float>("effDzSS");
  const unsigned int numberOfMuons = cfg.get<unsigned int>("NumberOfMuons");

  // track id/iso 
  const float trackIsoSF = cfg.get<float>("TrackIsolationSF"); 
  const float InterceptErr = cfg.get<float>("InterceptErr");
  const float SlopeErr = cfg.get<float>("SlopeErr");
  const float Correlation = cfg.get<float>("Correlation");

  // btagging
  const bool ApplyBTagVeto = cfg.get<bool>("ApplyBTagVeto");
  const string bTagAlgorithm = cfg.get<string>("BTagAlgorithm");
  const string bTagDiscriminator1 = cfg.get<string>("BTagDiscriminator1");
  const string bTagDiscriminator2 = cfg.get<string>("BTagDiscriminator2");
  const string bTagDiscriminator3 = cfg.get<string>("BTagDiscriminator3");
  const float btagCut = cfg.get<float>("BTagCut");
  const float bjetEta = cfg.get<float>("BJetEta");
  const float bjetPt = cfg.get<float>("BJetPt");
  //  const bool applyHEM = cfg.get<bool>("ApplyHEM");

  // ***********************************
  // ************* HEM issue ***********
  // ***********************************

  //bool applyHEM = era==2018; 
  bool applyHEM = false;
  float etaMinJetHEM = -3.2;
  float etaMaxJetHEM = -1.2;
  float phiMinJetHEM = -1.77;
  float phiMaxJetHEM = -0.67;
  float ptJetHEM = 30.0;
  float dphiHEM = 0.3;

  float etaMinEleHEM = -3.0;
  float etaMaxEleHEM = -1.4;
  float phiMinEleHEM = -1.57;
  float phiMaxEleHEM = -0.87;
  float ptEleHEM = 30.0;

  float dxyElectronCut = 0.05;
  float dzElectronCut = 0.1;
  float isoElectronCut = 0.15;

  unsigned int runHEM = 319077;
  float weightHEM = 0.35;

  // ***** Prefiring **********  
  bool applyPrefire = era==2016 || era==2017;
  // bool applyPrefire = false;

  TString BTagAlgorithm(bTagAlgorithm);
  TString BTagDiscriminator1(bTagDiscriminator1);
  TString BTagDiscriminator2(bTagDiscriminator2);
  TString BTagDiscriminator3(bTagDiscriminator3);

  TString DiMuonTriggerName(dimuonTriggerName);
  TString MuonHighPtFilterName(muonHighPtFilterName);
  TString MuonLowPtFilterName1(muonLowPtFilterName1);
  TString MuonLowPtFilterName2(muonLowPtFilterName2);
  TString DiMuonDzFilterName(dimuonDzFilterName);
  TString DiMuonSameSignFilterName(dimuonSameSignFilterName);

  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");

  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);

  //  const string Muon17TriggerFile = cfg.get<string>("Muon17TriggerEff");
  //  const string Muon8TriggerFile = cfg.get<string>("Muon8TriggerEff");
  const string Muon17TriggerFile = cfg.get<string>("MuonHighPtTriggerEff");
  const string Muon8TriggerFile = cfg.get<string>("MuonLowPtTriggerEff");

  // Higgs pt reweighting
  const string higgsPtFileName = cfg.get<string>("HiggsPtFileName");
  TString HiggsPtFileName(higgsPtFileName);
  const bool isVH = cfg.get<bool>("IsVH");

  const string correctionsWorkspaceFileName = cfg.get<string>("CorrectionWorkspaceFileName");
  TString CorrectionsWorkspaceFileName(correctionsWorkspaceFileName);

  const string roccorFileName = cfg.get<string>("RochesterCorrections");

  // ********** end of configuration *******************
  string cmsswBase = (getenv ("CMSSW_BASE"));

  std::ifstream fileList(argv[2]);
  std::ifstream fileListX(argv[2]);

  // event info
  ULong64_t event_nr;
  unsigned int event_run;
  unsigned int event_luminosityblock;

  // prefiring weight
  float prefiringweight;
  float prefiringweightup;
  float prefiringweightdown;
  float rho;

  // electrons
  unsigned int electron_count;
  float electron_pt[100];
  float electron_eta[100];
  float electron_phi[100];
  float electron_dxy[100];
  float electron_dz[100];
  float electron_mva_wp80_noIso_Fall17_v2[100];
  bool electron_pass_conversion[100];
  UChar_t electron_nmissinginnerhits[100]; 
  float electron_superclusterEta[100];
  float electron_r03_sumChargedHadronPt[100];
  float electron_r03_sumNeutralHadronEt[100];
  float electron_r03_sumPhotonEt[100];

  // tracks 
  UInt_t track_count;
  int track_ID[1000];
  float track_px[1000];
  float track_py[1000];
  float track_pz[1000];
  float track_pt[1000];
  float track_eta[1000];
  float track_phi[1000];
  float track_charge[1000];
  float track_mass[1000];
  float track_dxy[1000];
  float track_dxyerr[1000];
  float track_dz[1000];
  float track_dzerr[1000];
  bool track_highPurity[1000];
  // muons
  UInt_t muon_count;
  UInt_t muon_nMuonStations[1000];
  UInt_t muon_nMuonHits[1000];
  UInt_t muon_nPixelHits[1000];
  UInt_t muon_nTrackerHits[1000];
  float muon_px[1000];
  float muon_py[1000];
  float muon_pz[1000];
  float muon_pt[1000];
  float muon_px_uncorr[1000];
  float muon_py_uncorr[1000];
  float muon_pz_uncorr[1000];
  float muon_pt_uncorr[1000];
  float muon_eta[1000];
  float muon_phi[1000];
  float muon_pterror[1000];
  float muon_chi2[1000];
  float muon_ndof[1000];
  float muon_charge[1000];
  float muon_dxy[1000];
  float muon_dxyerr[1000];
  float muon_dz[1000];
  float muon_dzerr[1000];
  float muon_chargedHadIso[1000];
  float muon_neutralHadIso[1000];
  float muon_photonIso[1000];
  float muon_puIso[1000];
  float muon_r03_sumChargedHadronPt[1000];
  float muon_r03_sumChargedParticlePt[1000];
  float muon_r04_sumChargedHadronPt[1000];
  float muon_r04_sumChargedParticlePt[1000];

  bool muon_isPF[1000];
  bool muon_isGlobal[1000];
  bool muon_isTracker[1000];
  bool muon_isTight[1000];
  bool muon_isLoose[1000];
  bool muon_isMedium[1000];
  bool muon_isICHEP[1000];

  UInt_t genparticles_count;
  Float_t genparticles_e[1000];
  Float_t genparticles_px[1000];
  Float_t genparticles_py[1000];
  Float_t genparticles_pz[1000];
  Float_t genparticles_vx[1000];
  Float_t genparticles_vy[1000];
  Float_t genparticles_vz[1000];
  Int_t genparticles_pdgid[1000];
  Int_t genparticles_status[1000];
  UInt_t genparticles_info[1000];

  UInt_t   pfjet_count;
  Float_t  pfjet_e[200];     //[pfjet_count]
  Float_t  pfjet_pt[200];    //[pfjet_count]
  Float_t  pfjet_px[200];    //[pfjet_count]
  Float_t  pfjet_py[200];    //[pfjet_count]
  Float_t  pfjet_eta[200];   //[pfjet_count]
  Float_t  pfjet_phi[200];   //[pfjet_count]
  Float_t  pfjet_neutralhadronicenergy[200];   //[pfjet_count]
  Float_t  pfjet_chargedhadronicenergy[200];   //[pfjet_count]
  Float_t  pfjet_neutralemenergy[200];   //[pfjet_count]
  Float_t  pfjet_chargedemenergy[200];   //[pfjet_count]
  Float_t  pfjet_muonenergy[200];   //[pfjet_count]
  Float_t  pfjet_chargedmuonenergy[200];   //[pfjet_count]
  UInt_t   pfjet_chargedmulti[200];   //[pfjet_count]
  UInt_t   pfjet_neutralmulti[200];   //[pfjet_count]
  UInt_t   pfjet_chargedhadronmulti[200];   //[pfjet_count]
  Int_t    pfjet_flavour[200];   //[pfjet_count]
  Float_t  pfjet_btag[200][10];   //[pfjet_count]
  Float_t  pfjet_jecUncertainty[200]; //[pfjet_count]

  float genweight;

  float metx;
  float mety;
  float met;
  float metphi;
  
  float primvertex_x;
  float primvertex_y;
  float primvertex_z;

   // Trigger
  unsigned int trigobject_count;
  float trigobject_px[1000];
  float trigobject_py[1000];
  float trigobject_pz[1000];
  float trigobject_pt[1000];
  float  trigobject_eta[1000];
  float trigobject_phi[1000];
  bool trigobject_filters[1000][200];

  float numtruepileupinteractions;

  //unsigned int iLeadingPosTrig = 0;
  //vector<bool> trigobject_filter; trigobject_filter.clear();

  std::map<std::string, int> * hltriggerresults = new std::map<std::string, int>() ;
  std::map<std::string, int> * hltriggerprescales = new std::map<std::string, int>() ;
  std::vector<std::string>   * hltfilters = new std::vector<std::string>();
  std::vector<std::string>   * btagdiscriminators = new std::vector<std::string>();
  std::map<std::string, int> * flags = new std::map<std::string, int>();

  std::string rootFileName(argv[2]);
  
  std::string chainName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");
  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;
  if (TStrName.Contains("Signal")) {
    std::cout << "=============================" << std::endl;
    std::cout << "=== Running on Signal MC ====" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;
  }

  TString FullName = TStrName;      
  
  TFile * file = new TFile(FullName+TString(".root"),"recreate");

  file->cd("");
  
   // Muons
  TH1D * muonCountH = new TH1D("muonCountH","",11,-0.5,10.5);
  TH1D * nGoodMuonsH = new TH1D("nGoodMuonsH","",11,-0.5,10.5);
  TH1D * puWeightH = new TH1D("puWeightH","",250,0,5);
  TH1D * triggerWeightH = new TH1D("triggerWeightH","",100,0,2);

  // histograms after dimuon selection
  TH1D * ptLeadingMuH = new TH1D("ptLeadingMuH","",400,0,400);
  TH1D * ptTrailingMuH = new TH1D("ptTrailingMuH","",400,0,400);

  TH1D * ptLeadingIsoMuH = new TH1D("ptLeadingIsoMuH","",400,0,400);
  TH1D * ptTrailingIsoMuH = new TH1D("ptTrailingIsoMuH","",400,0,400);

  TH1D * ptLeadingNonIsoMuH = new TH1D("ptLeadingNonIsoMuH","",400,0,400);
  TH1D * ptTrailingNonIsoMuH = new TH1D("ptTrailingNonIsoMuH","",400,0,400);

  TH1D * etaLeadingMuH = new TH1D("etaLeadingMuH","",48,-2.4,2.4);
  TH1D * etaTrailingMuH = new TH1D("etaTrailingMuH","",48,-2.4,2.4);
  
  TH1D * dimuonMassH = new TH1D("dimuonMassH","",500,0,500);
  TH1D * dimuonMassIsoH = new TH1D("dimuonMassIsoH","",500,0,500);
  TH1D * dimuonMassNonIsoH = new TH1D("dimuonMassNonIsoH","",500,0,500);

  TH1D * nTracksLeadingMuH = new TH1D("nTracksLeadingMuH","",21,-0.5,20.5);
  TH1D * nTracksTrailingMuH = new TH1D("nTracksTrailingMuH","",21,-0.5,20.5);
  TH1D * nSigTracksLeadingMuH = new TH1D("nSigTracksLeadingMuH","",21,-0.5,20.5);
  TH1D * nSigTracksTrailingMuH = new TH1D("nSigTracksTrailingMuH","",21,-0.5,20.5);
  TH1D * nSoftTracksLeadingMuH = new TH1D("nSoftTracksLeadingMuH","",21,-0.5,20.5);
  TH1D * nSoftTracksTrailingMuH = new TH1D("nSoftTracksTrailingMuH","",21,-0.5,20.5);

  TH1D * ptRatioTrkLeadingMuH = new TH1D("ptRatioTrkLeadingMuH","",500,0.,5.);
  TH1D * ptRatioTrkTrailingMuH = new TH1D("ptRatioTrkTrailingMuH","",500,0.,5.);
  TH1D * ptRatioTrkMuH = new TH1D("ptRatioTrkMuH","",500,0.,5.);
  
  TH1D * ptRatioHadLeadingMuH = new TH1D("ptRatioHadLeadingMuH","",500,0.,5.);
  TH1D * ptRatioHadTrailingMuH = new TH1D("ptRatioHadTrailingMuH","",500,0.,5.);
  TH1D * ptRatioHadMuH = new TH1D("ptRatioHadMuH","",500,0.,5.);

  TH1D * IsoTrkLeadingMuH = new TH1D("IsoTrkLeadingMuH","",500,0.,5.);
  TH1D * IsoTrkTrailingMuH = new TH1D("IsoTrkTrailingMuH","",500,0.,5.);
  TH1D * IsoTrkMuH = new TH1D("IsoTrkMuH","",500,0.,5.);
  
  TH1D * IsoHadLeadingMuH = new TH1D("IsoHadLeadingMuH","",500,0.,5.);
  TH1D * IsoHadTrailingMuH = new TH1D("IsoHadTrailingMuH","",500,0.,5.);
  TH1D * IsoHadMuH = new TH1D("IsoHadMuH","",500,0.,5.);
  
  // isolated muons
  TH1D * ptTrackLeadingMuH = new TH1D("ptTrackLeadingMuH","",100,0,100);
  TH1D * etaTrackLeadingMuH = new TH1D("etaTrackLeadingMuH","",48,-2.4,2.4);
  TH1D * dxyTrackLeadingMuH = new TH1D("dxyTrackLeadingMuH","",200,-0.5,0.5);
  TH1D * dzTrackLeadingMuH = new TH1D("dzTrackLeadingMuH","",200,-1,1);

  TH1D * ptTrackTrailingMuH = new TH1D("ptTrackTrailingMuH","",100,0,100);
  TH1D * etaTrackTrailingMuH = new TH1D("etaTrackTrailingMuH","",48,-2.4,2.4);
  TH1D * dxyTrackTrailingMuH = new TH1D("dxyTrackTrailingMuH","",200,-0.5,0.5);
  TH1D * dzTrackTrailingMuH = new TH1D("dzTrackTrailingMuH","",200,-1,1);

  TH1D * ptTrackN1H = new TH1D("ptTrackN1H","",100,0,100);
  TH1D * etaTrackN1H = new TH1D("etaTrackN1H","",48,-2.4,2.4);
  TH1D * dxyTrackN1H = new TH1D("dxyTrackN1H","",200,-0.5,0.5);
  TH1D * dzTrackN1H = new TH1D("dzTrackN1H","",200,-1,1);

  TH1D * ptTrackH = new TH1D("ptTrackH","",100,0,100);
  TH1D * etaTrackH = new TH1D("etaTrackH","",48,-2.4,2.4);
  TH1D * dxyTrackH = new TH1D("dxyTrackH","",200,-0.5,0.5);
  TH1D * dzTrackH = new TH1D("dzTrackH","",200,-1,1);

  // Signal region
  TH1D * InvMassLeadingH = new TH1D("InvMassLeadingH","",100,0.,20.);
  TH1D * InvMassTrailingH = new TH1D("InvMassTrailingH","",100,0.,20.);
  TH1D * InvMassH = new TH1D("InvMassH","",100,0.,20.);
  TH2D * InvMass2DH = new TH2D("InvMass2DH","",100,0.,20.,100,0.,20.);

  // btag variations ->

  TH1D * InvMassH_btagUp = new TH1D("InvMassH_btagUp","",100,0.,20.);
  TH2D * InvMass2DH_btagUp = new TH2D("InvMass2DH_btagUp","",100,0.,20.,100,0.,20.);

  TH1D * InvMassH_btagDown = new TH1D("InvMassH_btagDown","",100,0.,20.);
  TH2D * InvMass2DH_btagDown = new TH2D("InvMass2DH_btagDown","",100,0.,20.,100,0.,20.);

  TH1D * InvMassH_mistagUp = new TH1D("InvMassH_mistagUp","",100,0.,20.);
  TH2D * InvMass2DH_mistagUp = new TH2D("InvMass2DH_mistagUp","",100,0.,20.,100,0.,20.);

  TH1D * InvMassH_mistagDown = new TH1D("InvMassH_mistagDown","",100,0.,20.);
  TH2D * InvMass2DH_mistagDown = new TH2D("InvMass2DH_mistagDown","",100,0.,20.,100,0.,20.);

  TH1D * InvMassH_jesUp = new TH1D("InvMassH_jesUp","",100,0.,20.);
  TH2D * InvMass2DH_jesUp = new TH2D("InvMass2DH_jesUp","",100,0.,20.,100,0.,20.);

  TH1D * InvMassH_jesDown = new TH1D("InvMassH_jesDown","",100,0.,20.);
  TH2D * InvMass2DH_jesDown = new TH2D("InvMass2DH_jesDown","",100,0.,20.,100,0.,20.);

  // trk iso SF variations
  TH1D * InvMassH_trkIsoUp = new TH1D("InvMassH_trkIsoUp","",100,0.,20.);
  TH2D * InvMass2DH_trkIsoUp = new TH2D("InvMass2DH_trkIsoUp","",100,0.,20.,100,0.,20.);
  
  TH1D * InvMassH_trkIsoDown = new TH1D("InvMassH_trkIsoDown","",100,0.,20.);
  TH2D * InvMass2DH_trkIsoDown = new TH2D("InvMass2DH_trkIsoDown","",100,0.,20.,100,0.,20.);

  // prefiring variations
  TH1D * InvMassH_prefireUp = new TH1D("InvMassH_prefireUp","",100,0.,20.);
  TH2D * InvMass2DH_prefireUp = new TH2D("InvMass2DH_prefireUp","",100,0.,20.,100,0.,20.);
  
  TH1D * InvMassH_prefireDown = new TH1D("InvMassH_prefireDown","",100,0.,20.);
  TH2D * InvMass2DH_prefireDown = new TH2D("InvMass2DH_prefireDown","",100,0.,20.,100,0.,20.);
  

  //
  
  TH1D * MetSelH = new TH1D("MetH","",400,0.,400.);
  TH1D * mTtotSelH = new TH1D("mTtotSelH","",400,0.,400.);

  TH1D * ptLeadingMuSelH = new TH1D("ptLeadingMuSelH","",100,0,100);
  TH1D * ptLeadingMuTrkSelH = new TH1D("ptLeadingMuTrkSelH","",100,0,100);
  TH1D * ptTrailingMuSelH = new TH1D("ptTrailingMuSelH","",100,0,100);
  TH1D * ptTrailingMuTrkSelH = new TH1D("ptTrailingMuTrkSelH","",100,0,100);

  TH1D * ptMuSelH = new TH1D("ptMuSelH","",100,0,100);
  TH1D * ptTrkSelH = new TH1D("ptTrkSelH","",100,0,100);

  TH1D * etaLeadingMuSelH = new TH1D("etaLeadingMuSelH","",48,-2.4,2.4);
  TH1D * etaLeadingMuTrkSelH = new TH1D("etaLeadingMuTrkSelH","",48,-2.4,2.4);
  TH1D * etaTrailingMuSelH = new TH1D("etaTrailingMuSelH","",48,-2.4,2.4);
  TH1D * etaTrailingMuTrkSelH = new TH1D("etaTrailingMuTrkSelH","",48,-2.4,2.4);

  TH1D * etaMuSelH = new TH1D("etaMuSelH","",48,-2.4,2.4);
  TH1D * etaTrkSelH = new TH1D("etaTrkSelH","",48,-2.4,2.4);

  TH1D * dRMuTrkSelH = new TH1D("dRMuTrkSelH","",50,0.,0.5);

  TH1D * dimuonMassSelH = new TH1D("dimuonMassSelH","",500,0,500);
  TH1D * invMass2Mu2TrkSelH = new TH1D("invMass2Mu2TrkSelH","",500,0,500);

  // Counters
  TH1D * counter_InputEventsH=new TH1D("counter_InputEventsH","",1,0.,2.);
  TH1D * counter_MuonSizeGTE2H=new TH1D("counter_MuonSizeGTE2H","",1,0.,2.);
  TH1D * counter_MuonKinematicsH=new TH1D("counter_MuonKinematicsH","",1,0.,2.);         
  TH1D * counter_nMuTrackSigH=new TH1D("counter_nMuTrackSigH","",1,0.,2.);
  TH1D * counter_FinalEventsH=new TH1D("counter_FinalEventsH","",1,0.,2.);         
  TH1D * counter_ControlEventsH=new TH1D("counter_ControlEventsH","",1,0.,2.);         
  TH1D * counter_ControlXEventsH=new TH1D("counter_ControlXEventsH","",1,0.,2.);         
  TH1D * counter_ControlYEventsH=new TH1D("counter_ControlYEventsH","",1,0.,2.);        

  // BTag counters
  TH1D * counter_btagCorrUp = new TH1D("counter_btagCorrUp","",1,0.,2.);
  TH1D * counter_btagCorrDown = new TH1D("counter_btagCorrDown","",1,0.,2.);

  TH1D * counter_btagUncorrUp = new TH1D("counter_btagUncorrUp","",1,0.,2.);
  TH1D * counter_btagUncorrDown = new TH1D("counter_btagUncorrDown","",1,0.,2.);

  TH1D * counter_mistagCorrUp = new TH1D("counter_mistagCorrUp","",1,0.,2.);
  TH1D * counter_mistagCorrDown = new TH1D("counter_mistagCorrDown","",1,0.,2.);

  TH1D * counter_mistagUncorrUp = new TH1D("counter_mistagUncorrUp","",1,0.,2.);
  TH1D * counter_mistagUncorrDown = new TH1D("counter_mistagUncorrDown","",1,0.,2.);

  TH1D * counter_prefireUp = new TH1D("counter_prefireUp","",1,0.,2.);
  TH1D * counter_prefireDown = new TH1D("counter_prefireDown","",1,0.,2.);



  TH1D * counter_btagH = new TH1D("counter_btagH","",1,0.,2.);
  TH1D * counter_btag_jesUpH = new TH1D("counter_btag_jesUpH","",1,0.,2.);
  TH1D * counter_btag_jesDownH = new TH1D("counter_btag_jesDownH","",1,0.,2.);

  TH1D * histWeightsH = new TH1D("histWeightsH","",1,0.,2.);
  TH1D * histWeightsSingleMuH = new TH1D("histWeightsSingleMuH","",1,0.,2.);
  TH1D * histWeightsTripleMuH = new TH1D("histWeightsTripleMuH","",1,0.,2.);
  TH1D * histWeightsDoubleMuSSH = new TH1D("histWeightsDoubleMuSSH","",1,0.,2.);
  TH1D * histWeightsAllTriggersH = new TH1D("histWeightsAllTriggersH","",1,0.,2.);

  // Background studies
  // N23 
  TH1D * InvMassN23leadingH = new TH1D("InvMassN23leadingH","",100,0.,20.);
  TH1D * InvMassN23trailingH = new TH1D("InvMassN23trailingH","",100,0.,20.);
  TH1D * InvMassN23H = new TH1D("InvMassN23H","",100,0.,20.);
  TH1D * InvMassN45H = new TH1D("InvMassN45H","",100,0.,20.);

  TH1D * ptMuN23H = new TH1D("ptMuN23H","",100,0,100);
  TH1D * ptTrkN23H = new TH1D("ptTrkN23H","",100,0,100);

  TH1D * etaMuN23H = new TH1D("etaMuN23H","",48,-2.4,2.4);
  TH1D * etaTrkN23H = new TH1D("etaTrkN23H","",48,-2.4,2.4);

  TH1D * dRMuTrkN23H = new TH1D("dRMuTrkN23H","",50,0.,0.5);

  // N1trk, N23trk 

  TH1D * InvMassHardestNtrk23leadingH = new TH1D("InvMassHardestNtrk23leadingH","",100,0.,20.);
  TH1D * InvMassHardestNtrk23trailingH = new TH1D("InvMassHardestNtrk23trailingH","",100,0.,20.);
  TH1D * InvMassHardestNtrk23H = new TH1D("InvMassHardestNtrk23H","",100,0.,20.);

  TH1D * InvMassSoftestNtrk23leadingH = new TH1D("InvMassSoftestNtrk23leadingH","",100,0.,20.);
  TH1D * InvMassSoftestNtrk23trailingH = new TH1D("InvMassSoftestNtrk23trailingH","",100,0.,20.);
  TH1D * InvMassSoftestNtrk23H = new TH1D("InvMassSoftestNtrk23H","",100,0.,20.);

  TH1D * InvMassHardestNtrk1leadingH = new TH1D("InvMassHardestNtrk1leadingH","",100,0.,20.);
  TH1D * InvMassHardestNtrk1trailingH = new TH1D("InvMassHardestNtrk1trailingH","",100,0.,20.);
  TH1D * InvMassHardestNtrk1H = new TH1D("InvMassHardestNtrk1H","",100,0.,20.);

  TH1D * InvMassSoftestNtrk1leadingH = new TH1D("InvMassSoftestNtrk1leadingH","",100,0.,20.);
  TH1D * InvMassSoftestNtrk1trailingH = new TH1D("InvMassSoftestNtrk1trailingH","",100,0.,20.);
  TH1D * InvMassSoftestNtrk1H = new TH1D("InvMassSoftestNtrk1H","",100,0.,20.);

  // Correlation Plots
  TH1D * InvMassTrackPlusMuon1D_ControlH = new TH1D("InvMassTrackPlusMuon1D_ControlH","",100,0.,20.); 
  TH2D * InvMassTrackPlusMuon2D_ControlH = new TH2D("InvMassTrackPlusMuon2D_ControlH","",100,0.,20.,100,0.,20.);

  TH1D * InvMassTrackPlusMuon1D_ControlXH = new TH1D("InvMassTrackPlusMuon1D_ControlXH","",100,0.,20.); 
  TH2D * InvMassTrackPlusMuon2D_ControlXH = new TH2D("InvMassTrackPlusMuon2D_ControlXH","",100,0.,20.,100,0.,20.);
   
  TH1D * InvMassTrackPlusMuon1D_ControlYH = new TH1D("InvMassTrackPlusMuon1D_ControlYH","",100,0.,20.); 
  TH2D * InvMassTrackPlusMuon2D_ControlYH = new TH2D("InvMassTrackPlusMuon2D_ControlYH","",100,0.,20.,100,0.,20.);

  // Monte Carlo information
  TH1D * deltaRMuonPionH = new TH1D("deltaRMuonPionH","",400,0,4);
  TH1D * pionPtH = new TH1D("pionPtH","",200,0,200);
  TH1D * muonPtH = new TH1D("muonPtH","",200,0,200);

  // generator variables
  float higgsPt;

  float genmu_Pt;
  float genmu_Eta;
  float genmu_P;
  float genmu_Q;

  float gentrk_Pt;
  float gentrk_Eta;
  float gentrk_P;
  float gentrk_Q;

  float genmutrk_Pt;
  float genmutrk_Eta;
  float genmutrk_P;
  float genmutrk_DR;
  float genmutrk_M;

  TTree * higgsTree = new TTree("higgsTree","");
  higgsTree->Branch("HiggsPt",&higgsPt,"HiggsPt/F");

  float genmu1_Pt;
  float genmu1_Eta;
  float genmu1_P;
  float genmu1_Q;
  unsigned int genmu1_ntrk;
  unsigned int genmu1_npart;
  unsigned int genmu1_ntrkIp;
  unsigned int genmu1_npartIp;

  float gentrk1_Pt;
  float gentrk1_Eta;
  float gentrk1_P;
  float gentrk1_Q;

  float genmutrk1_Pt;
  float genmutrk1_Eta;
  float genmutrk1_P;
  float genmutrk1_DR;
  float genmutrk1_M;
  
  float genmu2_Pt;
  float genmu2_Eta;
  float genmu2_P;
  float genmu2_Q;
  unsigned int genmu2_ntrk;
  unsigned int genmu2_npart;
  unsigned int genmu2_ntrkIp;
  unsigned int genmu2_npartIp;
  
  float gentrk2_Pt;
  float gentrk2_Eta;
  float gentrk2_P;
  float gentrk2_Q;

  float genmutrk2_Pt;
  float genmutrk2_Eta;
  float genmutrk2_P;
  float genmutrk2_DR;
  float genmutrk2_M;
  
  float dimuons_DR;
  float dimuons_Pt;
  float dimuons_P;

  unsigned int npart;
  unsigned int npartIp;
  unsigned int ntrk;
  unsigned int ntrkIp;

  float npu;
  float xvertex;
  float yvertex;
  float zvertex;

  float vgen_x;
  float vgen_y;
  float vgen_z;

  float puweight;

  TTree * dimuonsTree = new TTree("dimuonsTree","");

  dimuonsTree->Branch("xvertex",&xvertex,"xvertex/F");
  dimuonsTree->Branch("yvertex",&yvertex,"yvertex/F");
  dimuonsTree->Branch("zvertex",&zvertex,"zvertex/F");

  dimuonsTree->Branch("xvertex_gen",&vgen_x,"xvertex_gen/F");
  dimuonsTree->Branch("yvertex_gen",&vgen_y,"yvertex_gen/F");
  dimuonsTree->Branch("zvertex_gen",&vgen_z,"zvertex_gen/F");

  dimuonsTree->Branch("puweight",&puweight,"puweight/F");

  dimuonsTree->Branch("genmu1_P",&genmu1_P,"genmu1_P/F");
  dimuonsTree->Branch("genmu1_Pt",&genmu1_Pt,"genmu1_Pt/F");
  dimuonsTree->Branch("genmu1_Eta",&genmu1_Eta,"genmu1_Eta/F");
  dimuonsTree->Branch("genmu1_Q",&genmu1_Q,"genmu1_Q/F");
  dimuonsTree->Branch("genmu1_ntrk",&genmu1_ntrk,"genmu1_ntrkIp/I");
  dimuonsTree->Branch("genmu1_npart",&genmu1_npart,"genmu1_npartIp/I");
  dimuonsTree->Branch("genmu1_ntrkIp",&genmu1_ntrkIp,"genmu1_ntrkIp/I");
  dimuonsTree->Branch("genmu1_npartIp",&genmu1_npartIp,"genmu1_npartIp/I");

  dimuonsTree->Branch("gentrk1_P",&gentrk1_P,"gentrk1_P/F");
  dimuonsTree->Branch("gentrk1_Pt",&gentrk1_Pt,"gentrk1_Pt/F");
  dimuonsTree->Branch("gentrk1_Eta",&gentrk1_Eta,"gentrk1_Eta/F");
  dimuonsTree->Branch("gentrk1_Q",&gentrk1_Q,"gentrk1_Q/F");

  dimuonsTree->Branch("genmutrk1_P",&genmutrk1_P,"genmutrk1_P/F");
  dimuonsTree->Branch("genmutrk1_Pt",&genmutrk1_Pt,"genmutrk1_Pt/F");
  dimuonsTree->Branch("genmutrk1_Eta",&genmutrk1_Eta,"genmutrk1_Eta/F");
  dimuonsTree->Branch("genmutrk1_DR",&genmutrk1_DR,"genmutrk1_DR/F");
  dimuonsTree->Branch("genmutrk1_M",&genmutrk1_M,"genmutrk1_M/F");

  dimuonsTree->Branch("genmu2_P",&genmu2_P,"genmu2_P/F");
  dimuonsTree->Branch("genmu2_Pt",&genmu2_Pt,"genmu2_Pt/F");
  dimuonsTree->Branch("genmu2_Eta",&genmu2_Eta,"genmu2_Eta/F");
  dimuonsTree->Branch("genmu2_Q",&genmu2_Q,"genmu2_Q/F");
  dimuonsTree->Branch("genmu2_ntrk",&genmu2_ntrk,"genmu2_ntrk/i");
  dimuonsTree->Branch("genmu2_npart",&genmu2_npart,"genmu2_npart/i");
  dimuonsTree->Branch("genmu2_ntrkIp",&genmu2_ntrkIp,"genmu2_ntrkIp/i");
  dimuonsTree->Branch("genmu2_npartIp",&genmu2_npartIp,"genmu2_npartIp/i");

  dimuonsTree->Branch("gentrk2_P",&gentrk2_P,"gentrk2_P/F");
  dimuonsTree->Branch("gentrk2_Pt",&gentrk2_Pt,"gentrk2_Pt/F");
  dimuonsTree->Branch("gentrk2_Eta",&gentrk2_Eta,"gentrk2_Eta/F");
  dimuonsTree->Branch("gentrk2_Q",&gentrk2_Q,"gentrk2_Q/F");

  dimuonsTree->Branch("genmutrk2_P",&genmutrk2_P,"genmutrk2_P/F");
  dimuonsTree->Branch("genmutrk2_Pt",&genmutrk2_Pt,"genmutrk2_Pt/F");
  dimuonsTree->Branch("genmutrk2_Eta",&genmutrk2_Eta,"genmutrk2_Eta/F");
  dimuonsTree->Branch("genmutrk2_DR",&genmutrk2_DR,"genmutrk2_DR/F");
  dimuonsTree->Branch("genmutrk2_M",&genmutrk2_M,"genmutrk2_M/F");

  dimuonsTree->Branch("dimuons_DR",&dimuons_DR,"dimuons_DR/F");
  dimuonsTree->Branch("dimuons_Pt",&dimuons_Pt,"dimuons_Pt/F");
  dimuonsTree->Branch("dimuons_P",&dimuons_P,"dimuons_P/F");

  dimuonsTree->Branch("npart",&npart,"npart/i");
  dimuonsTree->Branch("npartIp",&npartIp,"npartIp/i");

  dimuonsTree->Branch("ntrk",&ntrk,"ntrk/i");
  dimuonsTree->Branch("ntrkIp",&ntrkIp,"ntrkIp/i");
  dimuonsTree->Branch("npu",&npu,"npu/F");

  TTree * tuple = new TTree("tuple","selected events");
  ULong64_t t_event;
  unsigned int t_run;

  unsigned int t_nmu;
  float t_mupt[100];
  float t_mueta[100];
  float t_muphi[100];
  float t_muiso[100];
  float t_mudxy[100];
  float t_mudz[100];
  float t_muq[100];
  unsigned int t_mu1_index;
  unsigned int t_mu2_index;

  float t_mu1_mutrkmass;
  float t_mu1_trkpt;
  float t_mu1_trketa;
  float t_mu1_trkphi;

  float t_mu2_mutrkmass;
  float t_mu2_trkpt;
  float t_mu2_trketa;
  float t_mu2_trkphi;

  unsigned int t_mu1_nsoft;
  unsigned int t_mu2_nsoft;
  unsigned int t_njets;
  unsigned int t_njetspt20;
  float t_met;
  float t_metphi;
  bool t_metfilters;
  bool t_badmufilter;

  tuple->Branch("event",&t_event,"event/l");
  tuple->Branch("run",&t_run,"run/i");
  tuple->Branch("nmuons",&t_nmu,"nmuons/i");
  tuple->Branch("ptmuon",t_mupt,"ptmuon[nmuons]/F");
  tuple->Branch("etamuon",t_mueta,"etamuon[nmuons]/F");
  tuple->Branch("phimuon",t_muphi,"phimuon[nmuons]/F");
  tuple->Branch("isomuon",t_muiso,"isomuon[nmuons]/F");
  tuple->Branch("dxymuon",t_mudxy,"dxymuon[nmuons]/F");
  tuple->Branch("dzmuon",t_mudz,"dzmuon[nmuons]/F");
  tuple->Branch("qmuon",t_muq,"qmuon[nmuons]/F");

  tuple->Branch("mu1index",&t_mu1_index,"mu1index/i");
  tuple->Branch("mu1trkmass",&t_mu1_mutrkmass,"mu1trkmass/F");
  tuple->Branch("mu1trkpt",&t_mu1_trkpt,"mu1trkpt/F");
  tuple->Branch("mu1trketa",&t_mu1_trketa,"mu1trketa/F");
  tuple->Branch("mu1trkphi",&t_mu1_trkphi,"mu1trkphi/F");
  tuple->Branch("mu1nsoft",&t_mu1_nsoft,"mu1nsoft/i");

  tuple->Branch("mu2index",&t_mu2_index,"mu2index/i");
  tuple->Branch("mu2trkmass",&t_mu2_mutrkmass,"mu2trkmass/F");
  tuple->Branch("mu2trkpt",&t_mu2_trkpt,"mu2trkpt/F");
  tuple->Branch("mu2trketa",&t_mu2_trketa,"mu2trketa/F");
  tuple->Branch("mu2trkphi",&t_mu2_trkphi,"mu2trkphi/F");
  tuple->Branch("mu2nsoft",&t_mu2_nsoft,"mu2nsoft/i");
  
  tuple->Branch("njets",&t_njets,"njets/i");
  tuple->Branch("njetspt20",&t_njetspt20,"njetspt20/i");
  tuple->Branch("met",&t_met,"met/F");
  tuple->Branch("metphi",&t_metphi,"metphi/F");
  tuple->Branch("metfilters",&t_metfilters,"metfilters/O");
  tuple->Branch("badmufilter",&t_badmufilter,"badmufilter/O");


  // Run-lumi selector
  string fullPathToJsonFile = cmsswBase + "/src/HtoAA/data/json/" + jsonFile;
  std::vector<Period> periods;  
  if (isData) { // read the good runs 
    std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
    if (inputFileStream.fail() ) {
      std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
      std::cout << "please check" << std::endl;
      std::cout << "quitting program" << std::endl;
      exit(-1);
    }
    
    for(std::string s; std::getline(inputFileStream, s); ) {
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }
  }

  // PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/HtoAA/data/PileUpDistrib/"+PileUpDataFile,"read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/HtoAA/data/PileUpDistrib/"+PileUpMCFile, "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  // Higgs reweighting 
  TFile * higgsPtFile = NULL;new TFile(TString(cmsswBase)+"/src/HtoAA/data/"+HiggsPtFileName);
  TH1D * higgsPtH = NULL;
  TH1D * higgsPt_WPlusH = NULL;
  TH1D * higgsPt_WMinusH = NULL;
  TH1D * higgsPt_ZH = NULL;
  if (applyHiggsPtWeight) { 
    std::cout << "ApplyHiggsPtWeight = " << applyHiggsPtWeight << std::endl;
    TString fullpath_HiggsPtFile = TString(cmsswBase)+"/src/HtoAA/data/"+HiggsPtFileName;
    higgsPtFile = new TFile(fullpath_HiggsPtFile);
    if (higgsPtFile->IsZombie()) {
      std::cout << fullpath_HiggsPtFile << "  not found" << std::endl;
      exit(-1);
    }
    if (isVH) {
      std::cout << "IsVH = " << isVH << std::endl;
      higgsPt_WPlusH = (TH1D*)higgsPtFile->Get("kfactor_WplusH");
      higgsPt_WMinusH = (TH1D*)higgsPtFile->Get("kfactor_WminusH");
      higgsPt_ZH = (TH1D*)higgsPtFile->Get("kfactor_ZH");
      if (higgsPt_WPlusH==NULL) {
	std::cout << "histogram kfactor_WplusH is not found in file " << fullpath_HiggsPtFile << std::endl;
	exit(-1);
      }
      if (higgsPt_WMinusH==NULL) {
	std::cout << "histogram kfactor_WminusH is not found in file " << fullpath_HiggsPtFile << std::endl;
	exit(-1);
      }
      if (higgsPt_ZH==NULL) {
	std::cout << "histogram kfactor_ZH is not found in file " << fullpath_HiggsPtFile << std::endl;
	exit(-1);
      }
    }
    else {
      higgsPtH = (TH1D*)higgsPtFile->Get("kfactor");
      if (higgsPtH==NULL) {
	std::cout << "histogram kfactor is not found in file " << fullpath_HiggsPtFile << std::endl;
	exit(-1);	
      }
    }
  }

  const bool applyBTagSF = cfg.get<bool>("ApplyBTagSF");

  // BTag SF file
  const string BtagSfFile = cfg.get<string>("BtagSfFile");
  BTagCalibration calib;
  BTagCalibrationReader reader_B;
  BTagCalibrationReader reader_C;
  BTagCalibrationReader reader_Light;
  if (applyBTagSF && !isData) {
    calib = BTagCalibration(bTagAlgorithm, BtagSfFile, true);
    reader_B = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central",
				     {"up","down","up_correlated","down_correlated","up_uncorrelated","down_uncorrelated"});
    reader_C = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central",
				     {"up","down","up_correlated","down_correlated","up_uncorrelated","down_uncorrelated"});
    reader_Light = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central",
					 {"up","down","up_correlated","down_correlated","up_uncorrelated","down_uncorrelated"});
    reader_B.load(calib, BTagEntry::FLAV_B, "comb");
    reader_C.load(calib, BTagEntry::FLAV_C, "comb");
    reader_Light.load(calib, BTagEntry::FLAV_UDSG, "incl");
  }


  // BTAG efficiency for various flavours ->
  TString fileBtagEff = (TString)cfg.get<string>("BtagMCeffFile");
  TFile *fileTagging  = new TFile(fileBtagEff);
  TH2F  *tagEff_B     = 0;
  TH2F  *tagEff_C     = 0;
  TH2F  *tagEff_Light = 0;
  TRandom3 *rand = new TRandom3();
  if (applyBTagSF && !isData) {
    tagEff_B     = (TH2F*)fileTagging->Get("btag_eff_b");
    tagEff_C     = (TH2F*)fileTagging->Get("btag_eff_c");
    tagEff_Light = (TH2F*)fileTagging->Get("btag_eff_oth");
  }
  float MaxBJetPt = 1000.;
  float MinBJetPt = 20.;
  float MaxBJetEta = 2.4;
  float MinBJetEta = 0.0;

  // MET filters for data
  std::vector<TString> metfilters;
  if (isData) {
    if (era==2016) {
      metfilters = {
	"Flag_HBHENoiseFilter",
	"Flag_HBHENoiseIsoFilter",
	"Flag_globalSuperTightHalo2016Filter",
	"Flag_EcalDeadCellTriggerPrimitiveFilter",
	"Flag_goodVertices",
	"Flag_BadPFMuonFilter",
	"Flag_eeBadScFilter"      
      };
    }
    if (era==2017) {    
      metfilters = {
	"Flag_HBHENoiseFilter",
	"Flag_HBHENoiseIsoFilter",
	"Flag_globalSuperTightHalo2016Filter",
	"Flag_EcalDeadCellTriggerPrimitiveFilter",
	"Flag_goodVertices",
	"Flag_BadPFMuonFilter",
	"ecalBadCalibReducedMINIAODFilter",
	"Flag_eeBadScFilter"
      };
    }
    if (era==2018) {
      metfilters = {
	"Flag_HBHENoiseFilter",
	"Flag_HBHENoiseIsoFilter",
	"Flag_globalSuperTightHalo2016Filter",
	"Flag_EcalDeadCellTriggerPrimitiveFilter",
	"Flag_goodVertices",
	"Flag_BadPFMuonFilter",
	"ecalBadCalibReducedMINIAODFilter",
	"Flag_eeBadScFilter"
      };
    }
  }
  else {
    if (era==2016) {
      metfilters = {
	"Flag_HBHENoiseFilter",
	"Flag_HBHENoiseIsoFilter",
	"Flag_globalSuperTightHalo2016Filter",
	"Flag_EcalDeadCellTriggerPrimitiveFilter",
	"Flag_goodVertices",
	"Flag_BadPFMuonFilter"
      };
    }
    if (era==2017) {    
      metfilters = {
	"Flag_HBHENoiseFilter",
	"Flag_HBHENoiseIsoFilter",
	"Flag_globalSuperTightHalo2016Filter",
	"Flag_EcalDeadCellTriggerPrimitiveFilter",
	"Flag_goodVertices",
	"Flag_BadPFMuonFilter",
	"ecalBadCalibReducedMINIAODFilter"
      };
    }
    if (era==2018) {
      metfilters = {
	"Flag_HBHENoiseFilter",
	"Flag_HBHENoiseIsoFilter",
	"Flag_globalSuperTightHalo2016Filter",
	"Flag_EcalDeadCellTriggerPrimitiveFilter",
	"Flag_goodVertices",
	"Flag_BadPFMuonFilter",
	"ecalBadCalibReducedMINIAODFilter"
      };
    }
  }

  std::vector<TString> badmufilter = {
    "Flag_BadPFMuonFilter"
  };

  // Trigger efficiencies
  ScaleFactor * SF_muon17 = new ScaleFactor();
  SF_muon17->init_ScaleFactor(TString(cmsswBase)+"/src/HtoAA/data/"+TString(Muon17TriggerFile));
  ScaleFactor * SF_muon8 = new ScaleFactor();
  SF_muon8->init_ScaleFactor(TString(cmsswBase)+"/src/HtoAA/data/"+TString(Muon8TriggerFile));

  // Rochester corrections
  std::string RochesterFileName = cmsswBase+"/src/HtoAA/data/Roccor/"+roccorFileName;
  RoccoR rc; 
  TString RoccorFileName = TString(roccorFileName);
  std::cout << std::endl;
  bool applyRoccoR = false;
  if (RoccorFileName=="None"||RoccorFileName=="none") {
    std::cout << "No Rochester corrections are applied" << std::endl;
  }
  else {
    applyRoccoR = true;
    std::cout << "Rochester corrections : " << RoccorFileName << std::endl;
    rc = RoccoR(RochesterFileName);
  }

  // isolation efficiency for era 2017
  TH1D * Mu17Iso_data = NULL;
  TH1D * Mu8Iso_data = NULL;
  TH1D * Mu17Iso_mc = NULL;
  TH1D * Mu8Iso_mc = NULL;
  float muonIsoMax = 10.;
  if (era==2017) {
    TString pathMuon17 = TString(cmsswBase)+"/src/HtoAA/data/"+TString(Muon17TriggerFile);
    TString pathMuon8  = TString(cmsswBase)+"/src/HtoAA/data/"+TString(Muon8TriggerFile);
    TFile * fileMu17 = new TFile(pathMuon17);
    TFile * fileMu8 = new TFile(pathMuon8);

    Mu17Iso_data = (TH1D*)fileMu17->Get("ZMass_Iso_Data");
    if (Mu17Iso_data==NULL) {
      std::cout << "histogram ZMass_Iso_Data is not found in file " << std::endl;
      std::cout << pathMuon17 << std::endl;
      exit(-1);
    }
    Mu17Iso_mc = (TH1D*)fileMu17->Get("ZMass_Iso_MC");
    if (Mu17Iso_mc==NULL) {
      std::cout << "histogram ZMass_Iso_MC is not found in file " << std::endl;
      std::cout << pathMuon17 << std::endl;
      exit(-1);
    }

    Mu8Iso_data = (TH1D*)fileMu8->Get("ZMass_Iso_Data");
    if (Mu8Iso_data==NULL) {
      std::cout << "histogram ZMass_Iso_Data is not found in file " << std::endl;
      std::cout << pathMuon8 << std::endl;
      exit(-1);
    }
    Mu8Iso_mc = (TH1D*)fileMu8->Get("ZMass_Iso_MC");
    if (Mu8Iso_mc==NULL) {
      std::cout << "histogram ZMass_Iso_MC is not found in file " << std::endl;
      std::cout << pathMuon8 << std::endl;
      exit(-1);
    }
    int nbins = Mu17Iso_data->GetNbinsX();
    muonIsoMax = Mu17Iso_data->GetBinLowEdge(nbins+1);
  }

  // Correction workspace
  //  TString corrWorkspaceFileName = TString(cmsswBase)+"/src/HtoAA/data/"+CorrectionsWorkspaceFileName;
  TFile * correctionWorkSpaceFile = new TFile(CorrectionsWorkspaceFileName);
  RooWorkspace *correctionWS = (RooWorkspace*)correctionWorkSpaceFile->Get("w");

  TString filen;
  int iFiles = 0;
  int events = 0;
  int counterFiles = 0;
  while (fileListX >> filen)
    counterFiles++;

  while (fileList >> filen) {
   iFiles++;
   cout << "file " << iFiles << " out of " << counterFiles << " : " << filen << endl;
   
   TString FileName = TString(filen);
   TFile * file_ = TFile::Open(FileName);

   if (file_==NULL) continue;
   
   if (file_->GetListOfKeys()->GetSize() == 0)
     continue; 

   if (file_->GetEND() > file_->GetSize())
     continue; 

   if (file_->GetSeekKeys()<=file_->GetEND()-file_->GetSize())
     continue;

   if (file_->IsZombie()) {
     cout << "cannot open file " << FileName << std::endl;
     continue;
   }


   // sum of weights (needed for normalization)
   TTree * _inittree = (TTree*)file_->Get(TString(initNtupleName));
   if (_inittree!=NULL) {
     Float_t Genweight;
     if (!isData)
       _inittree->SetBranchAddress("genweight",&Genweight);
     Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
     std::cout << "Number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
     for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
       _inittree->GetEntry(iEntry);
       if (isData)
	 histWeightsH->Fill(1.,1.);
       else
	 histWeightsH->Fill(1.,Genweight);
     }
   }


   TTree * tree_ = (TTree*)file_->Get(TString(chainName));
   
   if (tree_==NULL) continue;

   tree_->SetMaxVirtualSize(3000000);
   // event info
   tree_->SetBranchAddress("event_nr", &event_nr);
   tree_->SetBranchAddress("event_run", &event_run);
   tree_->SetBranchAddress("event_luminosityblock", &event_luminosityblock);

   // Prefiring weight and rho (density of diffuse noise)
   tree_->SetBranchAddress("rho",&rho);
   tree_->SetBranchAddress("prefiringweight",&prefiringweight);
   tree_->SetBranchAddress("prefiringweightup",&prefiringweightup);
   tree_->SetBranchAddress("prefiringweightdown",&prefiringweightdown);

   // Primary vertex
   tree_->SetBranchAddress("primvertex_x",&primvertex_x);
   tree_->SetBranchAddress("primvertex_y",&primvertex_y);
   tree_->SetBranchAddress("primvertex_z",&primvertex_z);

   // Muons
   tree_->SetBranchAddress("muon_count", &muon_count);
   tree_->SetBranchAddress("muon_nMuonStations", muon_nMuonStations);
   tree_->SetBranchAddress("muon_nMuonHits", muon_nMuonHits);
   tree_->SetBranchAddress("muon_nPixelHits", muon_nPixelHits);
   tree_->SetBranchAddress("muon_nTrackerHits", muon_nTrackerHits);
   tree_->SetBranchAddress("muon_px", muon_px);
   tree_->SetBranchAddress("muon_py", muon_py);
   tree_->SetBranchAddress("muon_pz", muon_pz);
   tree_->SetBranchAddress("muon_pt", muon_pt);
   tree_->SetBranchAddress("muon_eta", muon_eta);
   tree_->SetBranchAddress("muon_phi", muon_phi);
   tree_->SetBranchAddress("muon_pterror", muon_pterror);
   tree_->SetBranchAddress("muon_chi2", muon_chi2);
   tree_->SetBranchAddress("muon_ndof", muon_ndof);
   tree_->SetBranchAddress("muon_charge", muon_charge);
   tree_->SetBranchAddress("muon_dxy", muon_dxy);
   tree_->SetBranchAddress("muon_dxyerr", muon_dxyerr);
   tree_->SetBranchAddress("muon_dz", muon_dz);
   tree_->SetBranchAddress("muon_dzerr", muon_dzerr);
   tree_->SetBranchAddress("muon_r03_sumChargedHadronPt",muon_r03_sumChargedHadronPt);
   tree_->SetBranchAddress("muon_r03_sumChargedParticlePt",muon_r03_sumChargedParticlePt);
   tree_->SetBranchAddress("muon_r04_sumChargedHadronPt",muon_r04_sumChargedHadronPt);
   tree_->SetBranchAddress("muon_r04_sumChargedParticlePt",muon_r04_sumChargedParticlePt);
   tree_->SetBranchAddress("muon_chargedHadIso", muon_chargedHadIso);
   tree_->SetBranchAddress("muon_neutralHadIso", muon_neutralHadIso);
   tree_->SetBranchAddress("muon_photonIso", muon_photonIso);
   tree_->SetBranchAddress("muon_puIso", muon_puIso);
   tree_->SetBranchAddress("muon_isMedium", muon_isMedium);
   //   tree_->SetBranchAddress("muon_isICHEP", muon_isICHEP);

   // Electrons 
   tree_->SetBranchAddress("electron_count",&electron_count);
   tree_->SetBranchAddress("electron_pt",electron_pt);
   tree_->SetBranchAddress("electron_eta",electron_eta);
   tree_->SetBranchAddress("electron_phi",electron_phi);
   tree_->SetBranchAddress("electron_dxy",electron_dxy);
   tree_->SetBranchAddress("electron_dz",electron_dz);
   tree_->SetBranchAddress("electron_mva_wp80_noIso_Fall17_v2",electron_mva_wp80_noIso_Fall17_v2);
   tree_->SetBranchAddress("electron_pass_conversion",electron_pass_conversion);
   tree_->SetBranchAddress("electron_nmissinginnerhits",electron_nmissinginnerhits);
   tree_->SetBranchAddress("electron_superclusterEta",electron_superclusterEta);
   tree_->SetBranchAddress("electron_r03_sumChargedHadronPt",electron_r03_sumChargedHadronPt);
   tree_->SetBranchAddress("electron_r03_sumNeutralHadronEt",electron_r03_sumNeutralHadronEt);
   tree_->SetBranchAddress("electron_r03_sumPhotonEt",electron_r03_sumPhotonEt);

   // MET
   tree_->SetBranchAddress("pfmetcorr_ex", &metx);
   tree_->SetBranchAddress("pfmetcorr_ey", &mety);
   tree_->SetBranchAddress("pfmetcorr_pt", &met);
   tree_->SetBranchAddress("pfmetcorr_phi",&metphi);
 

   // Tracks
   tree_->SetBranchAddress("track_count", &track_count);
   tree_->SetBranchAddress("track_ID", track_ID);
   tree_->SetBranchAddress("track_px", track_px);
   tree_->SetBranchAddress("track_py", track_py);
   tree_->SetBranchAddress("track_pz", track_pz);
   tree_->SetBranchAddress("track_pt", track_pt);
   tree_->SetBranchAddress("track_eta", track_eta);
   tree_->SetBranchAddress("track_phi", track_phi);
   tree_->SetBranchAddress("track_mass", track_mass);
   tree_->SetBranchAddress("track_charge", track_charge);
   tree_->SetBranchAddress("track_dxy", track_dxy);
   tree_->SetBranchAddress("track_dxyerr", track_dxyerr);
   tree_->SetBranchAddress("track_dz", track_dz);
   tree_->SetBranchAddress("track_dzerr",track_dzerr);
   tree_->SetBranchAddress("track_highPurity", track_highPurity);

   // trigger objects
   tree_->SetBranchAddress("trigobject_count", &trigobject_count);
   tree_->SetBranchAddress("trigobject_px", trigobject_px);
   tree_->SetBranchAddress("trigobject_py", trigobject_py);
   tree_->SetBranchAddress("trigobject_pz", trigobject_pz);
   tree_->SetBranchAddress("trigobject_pt", trigobject_pt);
   tree_->SetBranchAddress("trigobject_eta", trigobject_eta);
   tree_->SetBranchAddress("trigobject_phi", trigobject_phi);
   tree_->SetBranchAddress("trigobject_filters",trigobject_filters);

   // Additional trigger objects
   tree_->SetBranchAddress("run_hltfilters",&hltfilters);
   tree_->SetBranchAddress("run_btagdiscriminators", &btagdiscriminators);
   tree_->SetBranchAddress("hltriggerresults",&hltriggerresults);
   tree_->SetBranchAddress("hltriggerprescales",&hltriggerprescales);
   tree_->SetBranchAddress("flags",&flags);

   tree_->SetBranchAddress("numtruepileupinteractions",&numtruepileupinteractions);

   tree_->SetBranchAddress("pfjet_count",&pfjet_count);
   tree_->SetBranchAddress("pfjet_e",pfjet_e);
   tree_->SetBranchAddress("pfjet_pt",pfjet_pt);
   tree_->SetBranchAddress("pfjet_px",pfjet_px);
   tree_->SetBranchAddress("pfjet_py",pfjet_py);
   tree_->SetBranchAddress("pfjet_eta",pfjet_eta);
   tree_->SetBranchAddress("pfjet_phi",pfjet_phi);
   tree_->SetBranchAddress("pfjet_neutralhadronicenergy",pfjet_neutralhadronicenergy);
   tree_->SetBranchAddress("pfjet_chargedhadronicenergy",pfjet_chargedhadronicenergy);
   tree_->SetBranchAddress("pfjet_neutralemenergy",pfjet_neutralemenergy);
   tree_->SetBranchAddress("pfjet_chargedemenergy",pfjet_chargedemenergy);
   tree_->SetBranchAddress("pfjet_muonenergy",pfjet_muonenergy);
   tree_->SetBranchAddress("pfjet_chargedmuonenergy",pfjet_chargedmuonenergy);
   tree_->SetBranchAddress("pfjet_chargedmulti",pfjet_chargedmulti);
   tree_->SetBranchAddress("pfjet_neutralmulti",pfjet_neutralmulti);
   tree_->SetBranchAddress("pfjet_flavour",pfjet_flavour);
   tree_->SetBranchAddress("pfjet_btag",pfjet_btag);
   tree_->SetBranchAddress("pfjet_jecUncertainty",pfjet_jecUncertainty);

   if (!isData) {
     tree_->SetBranchAddress("genweight",&genweight);
     tree_->SetBranchAddress("genparticles_count", &genparticles_count);
     tree_->SetBranchAddress("genparticles_e", genparticles_e);
     tree_->SetBranchAddress("genparticles_px", genparticles_px);
     tree_->SetBranchAddress("genparticles_py", genparticles_py);
     tree_->SetBranchAddress("genparticles_pz", genparticles_pz);
     tree_->SetBranchAddress("genparticles_vx", genparticles_vx);
     tree_->SetBranchAddress("genparticles_vy", genparticles_vy);
     tree_->SetBranchAddress("genparticles_vz", genparticles_vz);
     tree_->SetBranchAddress("genparticles_pdgid", genparticles_pdgid);
     tree_->SetBranchAddress("genparticles_status", genparticles_status);
     tree_->SetBranchAddress("genparticles_info", genparticles_info);
   }   

   int numberOfCandidates = tree_->GetEntries();

   std::cout << "number of events = " << numberOfCandidates << std::endl;
   
   TRandom3 rand;

   for (int iCand=0; iCand<numberOfCandidates; iCand++) {

     //     std::cout << "Ok 0" << std::endl;
     
     tree_->GetEntry(iCand);

     events++;
     if (events%10000==0) cout << "   processed events : " << events << endl;

     float weight = 1;
     if (!isData) {
       weight *= genweight;
     }

     puweight = 1.0;
     npu = 0.0;
     if (!isData) {
       puweight = float(PUofficial->get_PUweight(double(numtruepileupinteractions)));
       //       std::cout << "n(true interactions) = " << numtruepileupinteractions << "   :  PU weight = " << puweight << std::endl; 
       npu = numtruepileupinteractions;
     }

     // **************************************************
     // ********* generator studies **********************
     // **************************************************

     //     std::vector<unsigned int> posPion; posPion.clear();
     //     std::vector<unsigned int> negPion; negPion.clear();
     //     std::vector<unsigned int> posMuon; posMuon.clear();
     //     std::vector<unsigned int> negMuon; negMuon.clear();
     //     std::vector<unsigned int> allPions; allPions.clear();
     std::vector<unsigned int> genMuons; genMuons.clear();
     std::vector<unsigned int> genPions; genPions.clear();
     std::vector<unsigned int> gen_muons; gen_muons.clear(); // for Rochester corrections 

     vgen_x = 0;
     vgen_y = 0;
     vgen_z = 0;

     xvertex = primvertex_x;
     yvertex = primvertex_y;
     zvertex = primvertex_z;

     bool HiggsFound = false;
     unsigned int higgsIndex = 0;
     bool isWplus = false;
     bool isWminus = false;
     bool isZ = false;
     if (!isData) {
       // 
       // std::cout << "Generated particles = " << genparticles_count << std::endl;
       //       
       npart = 0;
       ntrk = 0;
       npartIp = 0;
       ntrkIp = 0;
       for (unsigned int iT=0; iT<track_count; ++iT) {
	 float chargeTrk = TMath::Abs(track_charge[iT]);
	 if (chargeTrk>0.5) {
	   if (track_pt[iT]>ptTrkLooseCut&&TMath::Abs(track_eta[iT])<etaTrkCut) {
	     ntrk++;
	     /*
	       std::cout << "track : " << ntrk << "  pdgId = " << track_ID[iT] 
	       << "   pT = " << track_pt[iT] 
	       << "   eta = " << track_eta[iT]
	       << "   dxy = " << track_dxy[iT]
	       << "   dz  = " << track_dz[iT] 
	       << "   purity = " << track_highPurity[iT] << std::endl;
	     */
	     if (abs(track_dxy[iT])<dxyTrkLooseCut&&TMath::Abs(track_dz[iT])<dzTrkLooseCut)
	       ntrkIp++;
	   }
	 }
       }
       // std::cout << std::endl;

       for (unsigned int iP=0; iP<genparticles_count; ++iP) {
	 if (genparticles_status[iP]==1&&
	     (genparticles_pdgid[iP]==13||genparticles_pdgid[iP]==-13)) {
	   gen_muons.push_back(iP);
	 }
	 if (genparticles_status[iP]==1&&(genparticles_info[iP]==12||genparticles_info[iP]==1)) {
	   if (genparticles_pdgid[iP]==13||genparticles_pdgid[iP]==-13) genMuons.push_back(iP);
	   if (genparticles_pdgid[iP]==211||genparticles_pdgid[iP]==-211) { 
	     TLorentzVector partLV; partLV.SetXYZT(genparticles_px[iP],
                                                   genparticles_py[iP],
                                                   genparticles_pz[iP],
                                                   genparticles_e[iP]);
	     float partEta = TMath::Abs(partLV.Eta());
	     float partPt = partLV.Pt();
	     if (partPt>ptTrkCut&&partEta<etaTrkCut)
	       genPions.push_back(iP);
	   }
	 }
	 if (genparticles_pdgid[iP]==24) isWplus = true;
	 if (genparticles_pdgid[iP]==-24) isWminus = true;
	 if (genparticles_pdgid[iP]==23) isZ = true;
	 if (genparticles_pdgid[iP]==35) {
	   higgsIndex = iP;
	   HiggsFound = true;
	   TLorentzVector higgsLV; higgsLV.SetXYZT(genparticles_px[iP],
						   genparticles_py[iP],
						   genparticles_pz[iP],
						   genparticles_e[iP]);
	   vgen_x = genparticles_vx[iP];
	   vgen_y = genparticles_vy[iP];
	   vgen_z = genparticles_vz[iP];
	   
	 }
       }

       for (unsigned int iP=0; iP<genparticles_count; ++iP) {
	 if (genparticles_status[iP]==1) {
	   int pdgid = TMath::Abs(genparticles_pdgid[iP]);
	   if (pdgid==11||pdgid==13||pdgid==211) {
	     TLorentzVector partLV;
	     partLV.SetXYZT(genparticles_px[iP],
			    genparticles_py[iP],
			    genparticles_pz[iP],
			    genparticles_e[iP]);
	     
	     if (partLV.Pt()>ptTrkLooseCut&&TMath::Abs(partLV.Eta())<etaTrkCut) {	       
	       //	   if (partLV.Pt()>0.0&&TMath::Abs(partLV.Eta())<5.0) {	       
	       npart++;
	       float v0[3];
	       float v[3];
	       float ppart[3];
	       ppart[0] = partLV.Px();
	       ppart[1] = partLV.Py();
	       ppart[2] = partLV.Pz();
	       v0[0] = vgen_x;
	       v0[1] = vgen_y;
	       v0[2] = vgen_z;
	       v[0] = genparticles_vx[iP];
	       v[1] = genparticles_vy[iP];
	       v[2] = genparticles_vz[iP];
	       float dxypart = 0;
	       float dzpart  = 0;
	       float ip = impactParameter(v0, v, ppart, dxypart, dzpart);
	       //	     std::cout << npart << "  pdgId = " << genparticles_pdgid[iP] 
	       //		       << "   pT = " << partLV.Pt() 
	       //		       << "   eta = " << partLV.Eta()
	       //		       << "   dxy = " << dxypart
	       //		       << "   dz  = " << dzpart << std::endl;
	       //	     if (dxypart<dxyTrkLooseCut&&dzpart<dzTrkLooseCut) {
	       if (ip<dzTrkLooseCut) {
		 npartIp++;
	       }
	     }
	   }
	 }
       }
       //       std::cout << std::endl;
       //       if (posMuon.size()==2||negMuon.size()==2) {
       //	 std::cout << "H->aa->4tau : " << std::endl;
       //	 std::cout << "Number of mu-   : " << negMuon.size() << std::endl;
       //	 std::cout << "Number of mu+   : " << posMuon.size() << std::endl;
       //	 std::cout << "Number of pion- : " << negPion.size() << std::endl;
       //	 std::cout << "Number of pion+ : " << posPion.size() << std::endl;
       //	 std::cout << std::endl;
       //       }

       //       std::vector<unsigned int> genPions; genPions.clear();

       /*
       int qmuon = 0;
       if (posMuon.size()==2&&negMuon.size()==0&&negPion.size()==2&&posPion.size()==0) {
	 genMuons = posMuon;
	 genPions = negPion;
	 qmuon = 1;
       }

       if (posMuon.size()==0&&negMuon.size()==2&&negPion.size()==0&&posPion.size()==2) {
	 genMuons = negMuon;
	 genPions = posPion;
	 qmuon = -1;
       }

       for (unsigned int iM=0; iM<negMuon.size(); ++iM)
	 genMuons.push_back(negMuon.at(iM));

       for (unsigned int iM=0; iM<posMuon.size(); ++iM)
	 genMuons.push_back(posMuon.at(iM));
       */

       //       if (genMuons.size()==2&&genPions.size()==2) {
       if (genMuons.size()==2&&genPions.size()==2) {
	 unsigned int muonIndex1 = genMuons.at(0);
	 TLorentzVector muonLV1; muonLV1.SetXYZT(genparticles_px[muonIndex1],
						 genparticles_py[muonIndex1],
						 genparticles_pz[muonIndex1],
						 genparticles_e[muonIndex1]);
	 
	 unsigned int muonIndex2 = genMuons.at(1);
	 TLorentzVector muonLV2; muonLV2.SetXYZT(genparticles_px[muonIndex2],
						 genparticles_py[muonIndex2],
						 genparticles_pz[muonIndex2],
						 genparticles_e[muonIndex2]);
	 
	 
	 
	 unsigned int muonIndexLeading = muonIndex1;
	 if (muonLV2.Pt()>muonLV1.Pt()) {
	   muonIndexLeading = muonIndex2;
	 }
	 
	 for (unsigned int iPion=0; iPion<genPions.size(); ++iPion) {
	   unsigned int pionIndex = genPions.at(iPion);
	   TLorentzVector pionLV; pionLV.SetXYZT(genparticles_px[pionIndex],
						 genparticles_py[pionIndex],
						 genparticles_pz[pionIndex],
						 genparticles_e[pionIndex]);
	   pionPtH->Fill(pionLV.Pt(),weight);
	   float dRMuonPion = deltaR(muonLV1.Eta(),muonLV1.Phi(),
				     pionLV.Eta(),pionLV.Phi());
	   float dRMuon2Pion = deltaR(muonLV2.Eta(),muonLV2.Phi(),
				      pionLV.Eta(),pionLV.Phi());

	   if (dRMuon2Pion<dRMuonPion) dRMuonPion = dRMuon2Pion;
	   deltaRMuonPionH->Fill(dRMuonPion,weight);
	 }
	 for (unsigned int iMuon=0; iMuon<genMuons.size(); ++iMuon) {
	   unsigned int muonIndex = genMuons.at(iMuon);
	   TLorentzVector muLV; muLV.SetXYZT(genparticles_px[muonIndex],
					     genparticles_py[muonIndex],
					     genparticles_pz[muonIndex],
					     genparticles_e[muonIndex]);
	   
	   int qmuon = 1;
	   if (genparticles_pdgid[muonIndex]==13) qmuon = -1;
	   unsigned int nparticlesIp = 0;
	   unsigned int nparticles = 0;
	   unsigned int ntracks = 0;
	   unsigned int ntracksIp = 0;

	   for (unsigned int iP=0; iP<genparticles_count; ++iP) {
	     int pdgPart = TMath::Abs(genparticles_pdgid[iP]);
	     float v0[3];
	     float ppart[3];
	     float v[3];
	     if (genparticles_status[iP]==1&&(pdgPart==211||pdgPart==11)) {
	       TLorentzVector partLV;
	       partLV.SetXYZT(genparticles_px[iP],
			      genparticles_py[iP],
			      genparticles_pz[iP],
			      genparticles_e[iP]);
	       
	       v0[0] = vgen_x;
	       v0[1] = vgen_y;
	       v0[2] = vgen_z;
	       
	       ppart[0] = partLV.Px();
	       ppart[1] = partLV.Py();
	       ppart[2] = partLV.Pz();
	       
	       v[0] = genparticles_vx[iP];
	       v[1] = genparticles_vy[iP];
	       v[2] = genparticles_vz[iP];
	       
	       float ptPart = partLV.Pt();
	       float etaPart = TMath::Abs(partLV.Eta());
	       
	       float dxyPart = 0;
	       float dzPart = 0;
	       float ipPart = impactParameter(v0, v, ppart, dxyPart, dzPart);
	       float dRpart = deltaR(muLV.Eta(),muLV.Phi(),
				     partLV.Eta(),partLV.Phi());

	       if (ptPart>ptTrkLooseCut && 
		   etaPart<etaTrkCut && 
		   dRpart<dRIsoMuon &&dRpart>0.001) {
		 nparticles++;
		 //		   if (dxyPart<dxyTrkLooseCut&&dzPart<dzTrkLooseCut)
		 if (ipPart<dzTrkLooseCut)
		   nparticlesIp++;
	       }
	     }
	   }

	   for (unsigned int iT=0; iT<track_count; ++iT) {
	     float chargeTrk = track_charge[iT];
	     if (chargeTrk>0.5) {
	       if (track_pt[iT]>ptTrkLooseCut&&TMath::Abs(track_eta[iT])<etaTrkCut) {
		 float dRtrack = deltaR(track_eta[iT],track_phi[iT],
					muLV.Eta(),muLV.Phi());
		 if (dRtrack<dRIsoMuon&&dRtrack>0.001) {
		   ntracks++;
		   if (abs(track_dxy[iT])<dxyTrkLooseCut&&TMath::Abs(track_dz[iT])<dzTrkLooseCut)
		     ntracksIp++;
		 }
	       }
	     }
	   }	   
	   //	   TLorentzVector mupionLV = muLV + pionLV;
	   if (muonIndex==muonIndexLeading) {
	     genmu1_Pt = muLV.Pt();
	     genmu1_P = muLV.P();
	     genmu1_Eta = muLV.Eta();
	     genmu1_Q = qmuon;
	     genmu1_ntrk = ntracks; 
	     genmu1_ntrkIp = ntracksIp;
	     genmu1_npart = nparticles; 
	     genmu1_npartIp = nparticlesIp;
	     /*	     gentrk1_Pt = pionLV.Pt();
	     	     gentrk1_P = pionLV.P();
	     	     gentrk1_Eta = pionLV.Eta();
	     	     gentrk1_Q = -1.;
	     	     genmutrk1_DR = dRMuonPion;
	     	     genmutrk1_Pt = mupionLV.Pt();
	     	     genmutrk1_P = mupionLV.P();
	     	     genmutrk1_Eta = mupionLV.Eta();
	     	     genmutrk1_M = mupionLV.M();
	     */
	   }
	   else {
	     genmu2_Pt = muLV.Pt();
	     genmu2_P = muLV.P();
	     genmu2_Eta = muLV.Eta();
	     genmu2_Q = qmuon;
	     genmu2_npart = nparticles;
	     genmu2_npartIp = nparticlesIp;
	     genmu2_ntrk = ntracks;
	     genmu2_ntrkIp = ntracksIp;
	     /*
	       gentrk2_Pt = pionLV.Pt();
	       gentrk2_P = pionLV.P();
	       gentrk2_Eta = pionLV.Eta();
	       gentrk2_Q = -1.;
	       genmutrk2_DR = dRMuonPion;
	       genmutrk2_Pt = mupionLV.Pt();
	       genmutrk2_P = mupionLV.P();
	       genmutrk2_Eta = mupionLV.Eta();
	       genmutrk2_M = mupionLV.M();
	     */
	   }
	   muonPtH->Fill(muLV.Pt(),weight);
	   TLorentzVector dimuonsLV = muonLV1 + muonLV2;
	   dimuons_DR = deltaR(muonLV1.Eta(),muonLV1.Phi(),
			       muonLV2.Eta(),muonLV2.Phi());
	   dimuons_P = dimuonsLV.P();
	   dimuons_Pt = dimuonsLV.Pt();
	   npu = numtruepileupinteractions;
	   dimuonsTree->Fill();
	 }
       }
     }

     // ****************************************************
     // *********** Higgs pT reweighting *******************
     // ****************************************************

     if (HiggsFound) {
       TLorentzVector higgsLV; higgsLV.SetXYZT(genparticles_px[higgsIndex],
					       genparticles_py[higgsIndex],
					       genparticles_pz[higgsIndex],
					       genparticles_e[higgsIndex]);
       higgsPt = higgsLV.Pt();
       higgsTree->Fill();
       if (applyHiggsPtWeight) {
	 double HiggsPtForWeighting = higgsPt;
	 if (higgsPt>500) HiggsPtForWeighting = 499;
	 double higgsPtWeight = 1;
	 if (isVH) {
	   if (isWplus)
	     higgsPtWeight = higgsPt_WPlusH->GetBinContent(higgsPt_WPlusH->FindBin(HiggsPtForWeighting));
	   if (isWminus)
	     higgsPtWeight = higgsPt_WMinusH->GetBinContent(higgsPt_WMinusH->FindBin(HiggsPtForWeighting));
	   if (isZ)
	     higgsPtWeight = higgsPt_ZH->GetBinContent(higgsPt_ZH->FindBin(HiggsPtForWeighting));
	 }
	 else {
	   higgsPtWeight = higgsPtH->GetBinContent(higgsPtH->FindBin(HiggsPtForWeighting));
	 }
	 weight *= higgsPtWeight;
	 //	 std::cout << "Higgs : " << genparticles_pdgid[higgsIndex] 
	 //		   << "    pT = " << higgsLV.Pt()
	 //		   << "    eta = " << higgsLV.Eta() << " weight = " << higgsPtWeight << std::endl;

       }
     }

     counter_InputEventsH->Fill(1.0,weight);

     if (isData) {
	if (applyGoodRunSelection) {
	  bool lumi = false;
	  int n=event_run;
	  int lum = event_luminosityblock;
	  
	  std::string num = std::to_string(n);
	  std::string lnum = std::to_string(lum);
	  for(const auto& a : periods)
	    {
	      if ( num.c_str() ==  a.name ) {
		//	      std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
		//std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
		
		for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
		  
		  //   cout<<b->lower<<"  "<<b->bigger<<endl;
		  if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
		}
		auto last = std::prev(a.ranges.end());
		// std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
		if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
	      }
	    }
	  if (!lumi) continue;
	}
     }

     // MET filters
     bool passedMETFilters = passedFilters(flags,metfilters);
     if (!passedMETFilters) {
       //       std::cout << "MET filters not passed : run = " << event_run << "    event = " << event_nr << std::endl;
       continue;
     }

     bool vetoEvent = false;
     if (applyHEM) {
       bool vetoElectron = false;
       for (unsigned int ie=0; ie<electron_count; ++ie) {
	 bool electronMvaId = electron_mva_wp80_noIso_Fall17_v2[ie];        
	 if (!electronMvaId) continue;
	 if (fabs(electron_dxy[ie]) >= dxyElectronCut) continue;
	 if (fabs(electron_dz[ie]) >= dzElectronCut) continue;
	 if (!electron_pass_conversion[ie]) continue;
	 if (electron_nmissinginnerhits[ie] > 1) continue;
	 float eA = getEffectiveArea(electron_superclusterEta[ie]);
	 float chargedIso = electron_r03_sumChargedHadronPt[ie]; 
	 float neutralIso = electron_r03_sumNeutralHadronEt[ie] + electron_r03_sumPhotonEt[ie]
	   - eA*rho;
	 if (neutralIso<0.0) neutralIso = 0;
	 float relIso = (chargedIso+neutralIso)/electron_pt[ie];
	 if (relIso>isoElectronCut) continue;
	 bool vetoObj = electron_pt[ie]>ptEleHEM;
	 vetoObj = vetoObj && electron_eta[ie]>etaMinEleHEM;
	 vetoObj = vetoObj && electron_eta[ie]<etaMaxEleHEM;
	 vetoObj = vetoObj && electron_phi[ie]>phiMinEleHEM;
	 vetoObj = vetoObj && electron_phi[ie]<phiMaxEleHEM;
	 if (vetoObj) {
	   //	   std::cout << "HEM : Electron Pt = " << electron_pt[ie]
	   //		     << "   Eta = " << electron_eta[ie]
	   //		     << "   Phi = " << electron_phi[ie] << std::endl;
	   vetoElectron = true;
	   break;
	 }
       }

       bool vetoJet = false; 
       for (unsigned int jet=0; jet<pfjet_count; ++jet) {

	 bool jetId = tightJetID(pfjet_e[jet],
				 pfjet_eta[jet],
				 pfjet_neutralhadronicenergy[jet],
				 pfjet_neutralemenergy[jet],
				 pfjet_muonenergy[jet],
				 pfjet_chargedhadronicenergy[jet],
				 pfjet_chargedemenergy[jet],
				 pfjet_neutralmulti[jet],
				 pfjet_chargedmulti[jet],
				 era);
	   
	 if (!jetId) continue;
	 bool vetoObj = pfjet_pt[jet]>ptJetHEM;
	 vetoObj = vetoObj && pfjet_eta[jet]>etaMinJetHEM;
	 vetoObj = vetoObj && pfjet_eta[jet]<etaMaxJetHEM;
	 vetoObj = vetoObj && pfjet_phi[jet]>phiMinJetHEM;
	 vetoObj = vetoObj && pfjet_phi[jet]<phiMaxJetHEM;
	 float dphi = dPhiFrom2P(pfjet_px[jet],pfjet_py[jet],metx,mety);
	 vetoObj = vetoObj && dphi<dphiHEM;
	 if (vetoObj) {
	   //	   std::cout << "HEM : Jet Pt = " << pfjet_pt[jet]
	   //		     << "   Eta = " << pfjet_eta[jet]
	   //		     << "   Phi = " << pfjet_phi[jet] << std::endl;
	   vetoJet = true;
	   break;
	 }
       }

       vetoEvent = vetoJet || vetoElectron;
       if (vetoEvent) {
	 if (isData) {
	   //	   std::cout << "HEM : event is rejected -> run = " << event_run << std::endl;
	   //	   std::cout << std::endl;
	   if (event_run>=runHEM) continue;
	 }  
	 else {
	   weight *= weightHEM;
	 }
       }
       
     }

     if (applyPrefire) {
       if (!isData) {
	 weight *= prefiringweight;
	 //	 std::cout << "prefiring : " << prefiringweight
	 //		   << "  + " << prefiringweightup
	 //		   << "  - " << prefiringweightdown << std::endl;
       }
     }

     puWeightH->Fill(puweight,1.0);
     weight *= puweight;

     // checking if dimuon trigger bit is ON (obsolete : do just trigger match)
     /*
     bool isDimuonTrigger = false;
     bool triggerFound = false;
     for (std::map<string,int>::iterator it=hltriggerresults->begin(); it!=hltriggerresults->end(); ++it) {
       TString trigName(it->first);
       if (trigName.Contains(DiMuonTriggerName)) {
	 //	 std::cout << trigName << " : " << it->second << std::endl;
	 if (it->second==1)
	   isDimuonTrigger = true;
	 triggerFound = true;
       }
     }
     if (!triggerFound) 
       std::cout << "HLT path " << DiMuonTriggerName << " is not found" << std::endl;

     if (!isDimuonTrigger) continue;

     unsigned int ntrig = hltriggerresults->size();
     std::cout << "ntrig = " << ntrig << std::endl;
     for (std::map<string,int>::iterator it=hltriggerresults->begin(); it!=hltriggerresults->end(); ++it) 
       std::cout << it->first << "  :  "  << it->second << std::endl;
     std::cout << std::endl;

     unsigned int npres = hltriggerprescales->size();
     std::cout << "npres = " << npres << std::endl;
     for (std::map<string,int>::iterator it=hltriggerprescales->begin(); it!=hltriggerprescales->end(); ++it) 
       std::cout << it->first << "  :  "  << it->second << std::endl;
     std::cout << std::endl;
     */

     // ******************
     // applying btag veto
     // ******************
    
     // std::cout << "Ok 1" << std::endl;

     float weight_btag = 1.0; 
     float weight_btag_up = 1.0;
     float weight_btag_down = 1.0;

     float weight_btag_corr_up = 1.0;
     float weight_btag_corr_down = 1.0;

     float weight_btag_uncorr_up = 1.0;
     float weight_btag_uncorr_down = 1.0;

     float weight_mistag_up = 1.0;
     float weight_mistag_down = 1.0;

     float weight_mistag_corr_up = 1.0;
     float weight_mistag_corr_down = 1.0;

     float weight_mistag_uncorr_up = 1.0;
     float weight_mistag_uncorr_down = 1.0;

     if (ApplyBTagVeto) {
       int nBTagDiscriminant1 = -1;
       int nBTagDiscriminant2 = -1;
       int nBTagDiscriminant3 = -1;
       unsigned int num_btags = 0;
       float Pdata = 1.;
       float Pdata_lf_up = 1.;
       float Pdata_lf_down = 1.;
       float Pdata_hf_up = 1.;
       float Pdata_hf_down = 1.;

       float Pdata_lf_corr_up = 1.;
       float Pdata_lf_corr_down = 1.;
       float Pdata_hf_corr_up = 1.;
       float Pdata_hf_corr_down = 1.;

       float Pdata_lf_uncorr_up = 1.;
       float Pdata_lf_uncorr_down = 1.;
       float Pdata_hf_uncorr_up = 1.;
       float Pdata_hf_uncorr_down = 1.;

       float Pmc = 1.;
       float Pdata_lf = 1.;
       float Pdata_hf = 1.;

       unsigned int nbtags = 0;
       unsigned int nbtags_jesUp = 0;
       unsigned int nbtags_jesDown =0;
       for (unsigned int ibtag=0; ibtag<btagdiscriminators->size(); ibtag++) {
	 TString discr(btagdiscriminators->at(ibtag));
	 if (discr==BTagDiscriminator1) {
	   nBTagDiscriminant1 = ibtag;
	   num_btags++;
	 }
	 if ((BTagAlgorithm == "pfDeepCSVJetTags" || BTagAlgorithm == "pfDeepFlavourJetTags") && discr == BTagDiscriminator2) {
	   nBTagDiscriminant2 = ibtag;
	   num_btags++;
	 }
	 if (BTagAlgorithm == "pfDeepFlavourJetTags" && discr == BTagDiscriminator3) {
	   nBTagDiscriminant3 = ibtag;
	   num_btags++;
	 }
       }
       bool isCorrectBTag = true;
       if (num_btags==0) {
	 std::cout << "No discriminants are found for " << BTagAlgorithm << std::endl;
	 isCorrectBTag = false;
       }
       if (BTagAlgorithm=="pfDeepCSVJetTags" && num_btags!=2) {
	 std::cout << "Numbers of discriminators for " << BTagAlgorithm << " = " << num_btags 
		   << "   should be 2 " << std::endl;
	 isCorrectBTag = false;
       }
       if (BTagAlgorithm=="pfDeepFlavourJetTags" && num_btags!=3) {
	 std::cout << "Numbers of discriminators for " << BTagAlgorithm << " = " << num_btags 
		   << "   should be 3" << std::endl;
	 isCorrectBTag = false;
       }
       t_njets = 0;
       if (isCorrectBTag) {
	 for (unsigned int jet=0; jet<pfjet_count; ++jet) {
	   
	   float absEta = TMath::Abs(pfjet_eta[jet]);
	   if (absEta>bjetEta) continue;

	   float JetPtForBTag = pfjet_pt[jet];
	   float JetEtaForBTag = absEta;
	   float jetEta = pfjet_eta[jet];

	   //
	   //	   if (pfjet_pt[jet]<bjetPt) 
	   //	     std::cout << "jet : pT = " << pfjet_pt[jet]
	   //		       << "  eta = " << pfjet_eta[jet]
	   //		       << "  phi = " << pfjet_phi[jet] << std::endl;
	   //

	   bool jetId = tightJetID(pfjet_e[jet],
				   pfjet_eta[jet],
				   pfjet_neutralhadronicenergy[jet],
				   pfjet_neutralemenergy[jet],
				   pfjet_muonenergy[jet],
				   pfjet_chargedhadronicenergy[jet],
				   pfjet_chargedemenergy[jet],
				   pfjet_neutralmulti[jet],
				   pfjet_chargedmulti[jet],
				   era);
	   
	   if (!jetId) continue;

	   if (JetPtForBTag > MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	   if (JetPtForBTag < MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	   if (JetEtaForBTag > MaxBJetEta) JetEtaForBTag = MaxBJetEta - 0.01;
	   if (JetEtaForBTag < MinBJetEta) JetEtaForBTag = MinBJetEta + 0.01;
	   
	   float btagDiscr = pfjet_btag[jet][nBTagDiscriminant1];
	   if (BTagAlgorithm=="pfDeepFlavourJetTags"||BTagAlgorithm=="pfDeepCSVJetTags")
	     btagDiscr += pfjet_btag[jet][nBTagDiscriminant2];
	   if (BTagAlgorithm=="pfDeepFlavourJetTags")
	     btagDiscr += pfjet_btag[jet][nBTagDiscriminant3];

	   bool tagged = btagDiscr>btagCut;

	   // counting b-tagged jets for both data and MC
           if (tagged){
	     if (isData) {
	       if (pfjet_pt[jet]>bjetPt) {
		 nbtags++;
		 nbtags_jesUp++;
		 nbtags_jesUp++;
	       }
	     }
	     else { 
	       float pfjet_pt_jesUp = (1.0+pfjet_jecUncertainty[jet])*pfjet_pt[jet];
	       float pfjet_pt_jesDown = (1.0-pfjet_jecUncertainty[jet])*pfjet_pt[jet];
	       if (pfjet_pt[jet]>bjetPt) nbtags++; // central estimate
	       if (pfjet_pt_jesUp>bjetPt) nbtags_jesUp++; // jesUp
	       if (pfjet_pt_jesDown>bjetPt) nbtags_jesDown++; // jesDown
	     }
	   }
	   
	   // BTag correction
	   if (!isData && applyBTagSF) {

	     int flavor = TMath::Abs(pfjet_flavour[jet]);

	     float tageff = tagEff_Light->GetBinContent(tagEff_Light->FindBin(JetPtForBTag, JetEtaForBTag));
	     float jet_scalefactor = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG,JetEtaForBTag,JetPtForBTag);

	     float jet_scalefactor_up = reader_Light.eval_auto_bounds("up",BTagEntry::FLAV_UDSG,JetEtaForBTag,JetPtForBTag);
	     float jet_scalefactor_down = reader_Light.eval_auto_bounds("down",BTagEntry::FLAV_UDSG,JetEtaForBTag,JetPtForBTag);

	     float jet_scalefactor_corr_up = reader_Light.eval_auto_bounds("up_correlated",BTagEntry::FLAV_UDSG,JetEtaForBTag,JetPtForBTag);
	     float jet_scalefactor_corr_down = reader_Light.eval_auto_bounds("down_correlated",BTagEntry::FLAV_UDSG,JetEtaForBTag,JetPtForBTag);

	     float jet_scalefactor_uncorr_up = reader_Light.eval_auto_bounds("up_uncorrelated",BTagEntry::FLAV_UDSG,JetEtaForBTag,JetPtForBTag);
	     float jet_scalefactor_uncorr_down = reader_Light.eval_auto_bounds("down_uncorrelated",BTagEntry::FLAV_UDSG,JetEtaForBTag,JetPtForBTag);

	     if (flavor==4) { 
	       tageff = tagEff_C->GetBinContent(tagEff_C->FindBin(JetPtForBTag,JetEtaForBTag));
	       jet_scalefactor = reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C,JetEtaForBTag,JetPtForBTag);

	       jet_scalefactor_up = reader_C.eval_auto_bounds("up",BTagEntry::FLAV_C,JetEtaForBTag,JetPtForBTag);
	       jet_scalefactor_down = reader_C.eval_auto_bounds("down",BTagEntry::FLAV_C,JetEtaForBTag,JetPtForBTag);

	       jet_scalefactor_corr_up = reader_C.eval_auto_bounds("up_correlated",BTagEntry::FLAV_C,JetEtaForBTag,JetPtForBTag);
	       jet_scalefactor_corr_down = reader_C.eval_auto_bounds("down_correlated",BTagEntry::FLAV_C,JetEtaForBTag,JetPtForBTag);

	       jet_scalefactor_uncorr_up = reader_C.eval_auto_bounds("up_uncorrelated",BTagEntry::FLAV_C,JetEtaForBTag,JetPtForBTag);
	       jet_scalefactor_uncorr_down = reader_C.eval_auto_bounds("down_uncorrelated",BTagEntry::FLAV_C,JetEtaForBTag,JetPtForBTag);

	     }

	     if (flavor==5) { 
	       tageff = tagEff_B->GetBinContent(tagEff_B->FindBin(JetPtForBTag,JetEtaForBTag));
	       jet_scalefactor = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B,JetEtaForBTag,JetPtForBTag);

	       jet_scalefactor_up = reader_B.eval_auto_bounds("up",BTagEntry::FLAV_B,JetEtaForBTag,JetPtForBTag);
	       jet_scalefactor_down = reader_B.eval_auto_bounds("down",BTagEntry::FLAV_B,JetEtaForBTag,JetPtForBTag);

	       jet_scalefactor_corr_up = reader_B.eval_auto_bounds("up_correlated",BTagEntry::FLAV_B,JetEtaForBTag,JetPtForBTag);
	       jet_scalefactor_corr_down = reader_B.eval_auto_bounds("down_correlated",BTagEntry::FLAV_B,JetEtaForBTag,JetPtForBTag);

	       jet_scalefactor_uncorr_up = reader_B.eval_auto_bounds("up_uncorrelated",BTagEntry::FLAV_B,JetEtaForBTag,JetPtForBTag);
	       jet_scalefactor_uncorr_down = reader_B.eval_auto_bounds("down_uncorrelated",BTagEntry::FLAV_B,JetEtaForBTag,JetPtForBTag);

	     }

	     if (pfjet_pt[jet]>bjetPt) {
	       if (tagged) {

		 Pmc = Pmc*tageff;
		 Pdata = Pdata*jet_scalefactor*tageff;
		 
		 if (flavor==4||flavor==5) {
		   Pdata_hf      = Pdata_hf*jet_scalefactor*tageff;

		   Pdata_hf_up   = Pdata_hf_up*jet_scalefactor_up*tageff;
		   Pdata_hf_down = Pdata_hf_down*jet_scalefactor_down*tageff;

		   Pdata_hf_corr_up   = Pdata_hf_corr_up*jet_scalefactor_corr_up*tageff;
		   Pdata_hf_corr_down = Pdata_hf_corr_down*jet_scalefactor_corr_down*tageff;

		   Pdata_hf_uncorr_up   = Pdata_hf_uncorr_up*jet_scalefactor_uncorr_up*tageff;
		   Pdata_hf_uncorr_down = Pdata_hf_uncorr_down*jet_scalefactor_uncorr_down*tageff;

		 }
		 else {
		   Pdata_lf = Pdata_lf*jet_scalefactor*tageff;

		   Pdata_lf_up = Pdata_lf_up*jet_scalefactor_up*tageff;
		   Pdata_lf_down = Pdata_lf_down*jet_scalefactor_down*tageff;

		   Pdata_lf_corr_up = Pdata_lf_corr_up*jet_scalefactor_corr_up*tageff;
		   Pdata_lf_corr_down = Pdata_lf_corr_down*jet_scalefactor_corr_down*tageff;

		   Pdata_lf_uncorr_up = Pdata_lf_uncorr_up*jet_scalefactor_uncorr_up*tageff;
		   Pdata_lf_uncorr_down = Pdata_lf_uncorr_down*jet_scalefactor_uncorr_down*tageff;

		 }
	       }
	       else {

		 Pmc = Pmc*(1-tageff);
		 Pdata = Pdata*(1.0-jet_scalefactor*tageff);
		 
		 if (flavor==4 || flavor==5) {
		   Pdata_hf = Pdata_hf*(1-jet_scalefactor*tageff);

		   Pdata_hf_up = Pdata_hf_up*(1.0-jet_scalefactor_up*tageff);
		   Pdata_hf_down = Pdata_hf_down*(1.0-jet_scalefactor_down*tageff);

		   Pdata_hf_corr_up = Pdata_hf_corr_up*(1.0-jet_scalefactor_corr_up*tageff);
		   Pdata_hf_corr_down = Pdata_hf_corr_down*(1.0-jet_scalefactor_corr_down*tageff);

		   Pdata_hf_uncorr_up = Pdata_hf_uncorr_up*(1.0-jet_scalefactor_uncorr_up*tageff);
		   Pdata_hf_uncorr_down = Pdata_hf_uncorr_down*(1.0-jet_scalefactor_uncorr_down*tageff);

		 }
		 else {

		   Pdata_lf = Pdata_lf*(1-jet_scalefactor*tageff);

		   Pdata_lf_up = Pdata_lf_up*(1.0-jet_scalefactor_up*tageff);
		   Pdata_lf_down = Pdata_lf_down*(1.0-jet_scalefactor_down*tageff);

		   Pdata_lf_corr_up = Pdata_lf_corr_up*(1.0-jet_scalefactor_corr_up*tageff);
		   Pdata_lf_corr_down = Pdata_lf_corr_down*(1.0-jet_scalefactor_corr_down*tageff);

		   Pdata_lf_uncorr_up = Pdata_lf_uncorr_up*(1.0-jet_scalefactor_uncorr_up*tageff);
		   Pdata_lf_uncorr_down = Pdata_lf_uncorr_down*(1.0-jet_scalefactor_uncorr_down*tageff);

		 }
	       }
	     }
	   }
	 }
       }

       if (!isData && applyBTagSF) { 
	 weight_btag = Pdata/Pmc;
	 // some protection against low weights ->
	 if (Pdata_hf<1e-4) {Pdata_hf=1e-4; Pdata_hf_up=1e-4; Pdata_hf_down=1e-4;}
	 if (Pdata_lf<1e-4) {Pdata_lf=1e-4; Pdata_lf_up=1e-4; Pdata_lf_down=1e-4;} 
	 weight *= weight_btag;

	 weight_btag_up = Pdata_hf_up/Pdata_hf;
	 weight_btag_down = Pdata_hf_down/Pdata_hf;

	 weight_btag_corr_up = Pdata_hf_corr_up/Pdata_hf;
	 weight_btag_corr_down = Pdata_hf_corr_down/Pdata_hf;

	 weight_btag_uncorr_up = Pdata_hf_uncorr_up/Pdata_hf;
	 weight_btag_uncorr_down = Pdata_hf_uncorr_down/Pdata_hf;

	 weight_mistag_up = Pdata_lf_up/Pdata_lf;
	 weight_mistag_down = Pdata_lf_down/Pdata_lf;

	 weight_mistag_corr_up = Pdata_lf_corr_up/Pdata_lf;
	 weight_mistag_corr_down = Pdata_lf_corr_down/Pdata_lf;

	 weight_mistag_uncorr_up = Pdata_lf_uncorr_up/Pdata_lf;
	 weight_mistag_uncorr_down = Pdata_lf_uncorr_down/Pdata_lf;

	 /*
	 std::cout << "nbtags = " << nbtags << std::endl;
	 std::cout << "Pdata = " << Pdata 
		   << "  Pdata_lf = " << Pdata_lf 
		   << "  Pdata_hf = " << Pdata_hf
		   << "  Pdata_lf*Pdata_hf = " << Pdata_lf*Pdata_hf << std::endl;
	 std::cout << "weight_btag = " << weight_btag << std::endl;
	 std::cout << "weight_btag_corr     down/up = " << weight_btag_corr_down << "/" << weight_btag_corr_up << std::endl;
	 std::cout << "weight_btag_uncorr   down/up = " << weight_btag_uncorr_down << "/" << weight_btag_uncorr_up << std::endl; 
	 std::cout << "weight_mistag_corr   down/up = " << weight_mistag_corr_down << "/" << weight_mistag_corr_up << std::endl;
	 std::cout << "weight_mistag_uncorr down/up = " << weight_mistag_uncorr_down << "/" << weight_mistag_uncorr_up << std::endl;
	 */
       }

       //       if (nbtags!=nbtags_jesDown||nbtags) 
       //	 std::cout << "nbtags = " << nbtags 
       //		   << "   nbtags_jesDown = " << nbtags_jesDown
       //		   << "   nbtags_jesUp = " << nbtags_jesUp 
       //		   << std::endl;
       if (nbtags==0) counter_btagH->Fill(1.,weight);
       if (nbtags_jesUp==0) counter_btag_jesUpH->Fill(1.,weight);
       if (nbtags_jesDown==0) counter_btag_jesDownH->Fill(1.,weight);

       if (nbtags>0) continue;

     }
     //     std::cout << "Ok 2" << std::endl;
     // finding HLT filters in the HLT Filter library
     unsigned int nMu8Leg   = 0;
     unsigned int nMu17Leg  = 0;
     unsigned int nDZFilter = 0;
     unsigned int nSSFilter = 0;
     bool isMu8Leg = false;
     bool isMu17Leg = false;
     bool isDZFilter = false;
     bool isSSFilter = false;

     unsigned int nfilters = hltfilters->size();
     for (unsigned int i=0; i<nfilters; ++i) {
       //       std::cout << hltfilters->at(i) << std::endl;
       TString HLTFilter(hltfilters->at(i));
       if (HLTFilter==MuonHighPtFilterName) {
	 nMu17Leg = i;
	 isMu17Leg = true;
	 //	 std::cout << HLTFilter << ":" << i << std::endl;
       }
       if (HLTFilter==MuonLowPtFilterName1||HLTFilter==MuonLowPtFilterName2) {
	 nMu8Leg = i;
	 isMu8Leg = true;
	 //	 std::cout << HLTFilter << ":" << i << std::endl;
       }
       if (HLTFilter==DiMuonDzFilterName) {
	 nDZFilter = i;
	 isDZFilter = true;
	 //	 std::cout << HLTFilter << ":" << i << std::endl;
       }
       if (HLTFilter==DiMuonSameSignFilterName) {
	 nSSFilter = i;
	 isSSFilter = true;
	 //	 std::cout << HLTFilter << ":" << i <<  std::endl;
       }
     }
     //     std::cout << std::endl;
     
     if (!isMu17Leg) {
       cout << "Filter " << MuonHighPtFilterName << " not found " << endl;
       exit(-1);
     }
     if (!isMu8Leg) {
       cout << "Filters " << MuonLowPtFilterName1 
	    << " or " << MuonLowPtFilterName2
	    << " not found " << endl;
       exit(-1);
     }
     if (!isDZFilter) {
       cout << "Filter " << DiMuonDzFilterName << " not found " << endl;
       exit(-1);
     }
     if (!isSSFilter) {
       cout << "Filter " << DiMuonSameSignFilterName << " not found " << endl;
       exit(-1);
     }


     // *************************
     // selecting generator muons
     // *************************
     for (unsigned int iP=0; iP<genparticles_count; ++iP) {
       if (genparticles_status[iP]==1) {
	 int pdgid = TMath::Abs(genparticles_pdgid[iP]);
	 if (pdgid==13) {
	   gen_muons.push_back(iP);
	 }
       }
     }
     //     std::cout << "+++++++++++++++++++++++++++++++++++" << std::endl;
     //     for(UInt_t i=0;i<muon_count;i++){
     //       std::cout << " " << i << " " << muon_px[i] << " " << muon_py[i] << " " << muon_pz[i] << std::endl;       
     //     }     
     //     std::cout << std::endl;
     // **********************
     // selecting good muons
     // **********************
     double momSF = 1.0;
     vector<unsigned int> muons; muons.clear();
     for(UInt_t i=0;i<muon_count;i++){
       bool muonID = muon_isMedium[i]; // MC 
       // obsolete
       // if (isData) {
       //   if (event_run >= 278820 && muon_isMedium[i]) muonID = true; // Run2016G-H
       //   if (event_run < 278820  && muon_isICHEP[i]) muonID = true; // Run2016B-F
       // }
       if (!muonID) continue;
       if(fabs(muon_dxy[i])>dxyMuonCut) continue;
       if(fabs(muon_dz[i])>dzMuonCut) continue;
       
       muon_pt_uncorr[i] = muon_pt[i];
       muon_px_uncorr[i] = muon_px[i];
       muon_py_uncorr[i] = muon_py[i];
       muon_pz_uncorr[i] = muon_pz[i];

       if (applyRoccoR) {
	 double pt = muon_pt[i];
	 double eta = muon_eta[i];
	 double phi = muon_phi[i];
	 int Q = 1;
	 if (muon_charge[i]<-0.5) Q = -1;
	 if (isData) 
	   momSF = rc.kScaleDT(Q,pt,eta,phi);
	 else {
	   double genPt = pt;
	   double dRmin = 0.2;
	   for (unsigned int igen=0; igen<gen_muons.size(); ++igen) {
	     unsigned int ipart = gen_muons[igen];
	     TLorentzVector partLV; partLV.SetXYZM(genparticles_px[ipart],
						   genparticles_py[ipart],
						   genparticles_pz[ipart],
						   genparticles_e[ipart]
						   );
	     double dR = deltaR(muon_eta[i],muon_phi[i],
				partLV.Eta(),partLV.Phi());
	     if (dR<dRmin) {
	       genPt = partLV.Pt();
	       dRmin = dR;
	     }
	     
	   }
	   momSF = rc.kSpreadMC(Q,pt,eta,phi,genPt);
       
	 }
	 // don't allow for large corretions
	 if (momSF<=0.8) momSF = 0.8;
	 if (momSF>=1.2) momSF = 1.2;
	 
	 //	 float muon_pT = muon_pt[i];

	 
	 muon_px[i] *= momSF;
	 muon_py[i] *= momSF;
	 muon_pz[i] *= momSF;
	 muon_pt[i] *= momSF;
	 
       }
       
       if(muon_pt[i]<ptMuonLowCut) continue;
       if(fabs(muon_eta[i])>etaMuonLowCut) continue;
       //       cout << "SumPt(0p3) = " << muon_r03_sumChargedHadronPt[i] 
       //	    << "   SumPt(0p4) = " << muon_r04_sumChargedHadronPt[i] << endl;
       muons.push_back(i);
     }
     
     nGoodMuonsH->Fill(float(muons.size()),weight);
      
     if (muons.size()<2) continue; // quit event if number of good muons < 2
     counter_MuonSizeGTE2H->Fill(1.,weight);

     // *************************
     // selection of dimuon pairs 
     // *************************

     float maxPtSum = -1;
     int iLeading = -1;
     int iTrailing = -1;
     for (unsigned int i1=0; i1<muons.size()-1; ++i1) {
       int index1 = muons.at(i1);
       for (unsigned int i2=i1+1; i2<muons.size(); ++i2) {
	 int index2 = muons.at(i2);
	 float ptSum = muon_pt[index1] + muon_pt[index2];
	 float charge = muon_charge[index1] *  muon_charge[index2];
	 bool chargeSelection = charge<0;
	 if (sameSign)
	   chargeSelection = charge>0;
	 if (!chargeSelection) continue;
	 float dRmuons = deltaR(muon_eta[index1],muon_phi[index1],
				muon_eta[index2],muon_phi[index2]);
	 
	 if (dRmuons<dRMuonsCut) continue;
	 bool mu1MatchMu17 = false;
	 bool mu1MatchMu8  = false;
	 bool mu1MatchDz   = false;
	 bool mu1MatchSS   = false;
	 for (unsigned int iT=0; iT<trigobject_count; ++iT) {
	   float dRtrig = deltaR(muon_eta[index1],muon_phi[index1],
				 trigobject_eta[iT],trigobject_phi[iT]);
	   if (dRtrig>DRTrigMatch) continue;
	   if (trigobject_filters[iT][nMu17Leg])
	     mu1MatchMu17 = true;
	   if (trigobject_filters[iT][nMu8Leg])
	     mu1MatchMu8 = true;
	   if (trigobject_filters[iT][nDZFilter])
	     mu1MatchDz = true;
	   if (trigobject_filters[iT][nSSFilter])
	     mu1MatchSS = true;
	 }
	 bool mu2MatchMu17 = false;
	 bool mu2MatchMu8  = false;
	 bool mu2MatchDz   = false;
	 bool mu2MatchSS   = false;
	 for (unsigned int iT=0; iT<trigobject_count; ++iT) {
	   float dRtrig = deltaR(muon_eta[index2],muon_phi[index2],
				 trigobject_eta[iT],trigobject_phi[iT]);
	   if (dRtrig>DRTrigMatch) continue;
	   if (trigobject_filters[iT][nMu17Leg])
	     mu2MatchMu17 = true;
	   if (trigobject_filters[iT][nMu8Leg])
	     mu2MatchMu8 = true;
	   if (trigobject_filters[iT][nDZFilter])
	     mu2MatchDz = true;
	   if (trigobject_filters[iT][nSSFilter])
	     mu2MatchSS = true;
	 }

	 bool mu1PtHigh = muon_pt[index1]>ptMuonHighCut && fabs(muon_eta[index1])<etaMuonHighCut;
	 bool mu1PtLow  = muon_pt[index1]>ptMuonLowCut && fabs(muon_eta[index1])<etaMuonLowCut;
	 bool mu2PtHigh = muon_pt[index2]>ptMuonHighCut && fabs(muon_eta[index2])<etaMuonHighCut;
	 bool mu2PtLow  = muon_pt[index2]>ptMuonLowCut && fabs(muon_eta[index2])<etaMuonLowCut;
	 
	 // trigger condition
	 bool isTriggerMatched = true;
	 if (applyTriggerMatch) {
	   isTriggerMatched = (mu1MatchMu17&&mu2MatchMu8&&mu1PtHigh&&mu2PtLow)||(mu1MatchMu8&&mu2MatchMu17&&mu1PtLow&&mu2PtHigh);
	   if (isData) {
	     isTriggerMatched = isTriggerMatched && mu1MatchSS && mu2MatchSS;
	     if (era==2016) {
	       if (event_run<=274442||event_run>=280919) // when dZ filter is present
		 isTriggerMatched = isTriggerMatched && mu1MatchDz && mu2MatchDz;
	     }
	     if (era==2017) // DZ filter 
	       isTriggerMatched = isTriggerMatched && mu1MatchDz && mu2MatchDz;
	   }
	 }
	 else {
	   isTriggerMatched = (mu1PtHigh&&mu2PtLow) || (mu1PtLow&&mu2PtHigh);
	 }
	 if (!isTriggerMatched) continue;
	 if (ptSum>maxPtSum) { // choose the mair with maximum Sum(pT)
	   maxPtSum = ptSum;
	   if (muon_pt[index1]>muon_pt[index2]) {
	     iLeading = index1;
	     iTrailing = index2;
	   }
	   else {
	     iLeading = index2;
	     iTrailing = index1;
	   }
	 }
       }
     }

     if (iLeading<0) continue;
     if (iTrailing<0) continue;

     double triggerWeight = 1;
     double idLeadingWeight = 1;
     double idTrailingWeight = 1;
     if (!isData) { // trigger efficiency here
       double ptLeading = muon_pt[iLeading];
       double etaLeading = muon_eta[iLeading];
       double ptTrailing = muon_pt[iTrailing];
       double etaTrailing = muon_eta[iTrailing];

       double isoLeading = muon_r03_sumChargedHadronPt[iLeading]/ptLeading;
       double isoTrailing = muon_r03_sumChargedHadronPt[iTrailing]/ptTrailing;

       if (isoLeading<=0) isoLeading = 0.01;
       if (isoLeading>=muonIsoMax) isoLeading = muonIsoMax - 0.01;

       if (isoTrailing<=0) isoTrailing = 0.01;
       if (isoTrailing>=muonIsoMax) isoTrailing = muonIsoMax - 0.01;

       double iso17DataLeading  = 1.0;
       double iso17DataTrailing = 1.0;
       double iso17MCLeading    = 1.0;
       double iso17MCTrailing   = 1.0;

       double iso8DataLeading  = 1.0;
       double iso8DataTrailing = 1.0;
       double iso8MCLeading    = 1.0;
       double iso8MCTrailing   = 1.0;

       if (era==2017) {
	 iso17DataLeading  = Mu17Iso_data->GetBinContent(Mu17Iso_data->FindBin(isoLeading));
	 iso17DataTrailing = Mu17Iso_data->GetBinContent(Mu17Iso_data->FindBin(isoTrailing));
	 iso17MCLeading    = Mu17Iso_mc->GetBinContent(Mu17Iso_mc->FindBin(isoLeading));
	 iso17MCTrailing   = Mu17Iso_mc->GetBinContent(Mu17Iso_mc->FindBin(isoTrailing));

	 iso8DataLeading  = Mu8Iso_data->GetBinContent(Mu8Iso_data->FindBin(isoLeading));
	 iso8DataTrailing = Mu8Iso_data->GetBinContent(Mu8Iso_data->FindBin(isoTrailing));
	 iso8MCLeading    = Mu8Iso_mc->GetBinContent(Mu8Iso_mc->FindBin(isoLeading));
	 iso8MCTrailing   = Mu8Iso_mc->GetBinContent(Mu8Iso_mc->FindBin(isoTrailing));
       }

       if (ptLeading<10.0) ptLeading = 10.01;
       if (ptTrailing<10.0) ptTrailing = 10.01;
       
       //       std::cout << "before weighting " << std::endl;

       correctionWS->var("m_pt")->setVal(ptLeading);
       correctionWS->var("m_eta")->setVal(etaLeading);
       double idLeadingW  = correctionWS->function("m_id_ic_ratio")->getVal();
       double trkLeadingW = correctionWS->function("m_trk_ratio")->getVal();
       idLeadingWeight = idLeadingW * trkLeadingW;
       correctionWS->var("m_pt")->setVal(ptTrailing);
       correctionWS->var("m_eta")->setVal(etaTrailing);
       double idTrailingW = correctionWS->function("m_id_ic_ratio")->getVal();
       double trkTrailingW = correctionWS->function("m_trk_ratio")->getVal();
       idTrailingWeight = idTrailingW * trkTrailingW;

       //       std::cout << "after weighting" << std::endl;

       double effMu17dataTrailing = SF_muon17->get_EfficiencyData(ptTrailing,etaTrailing); 
       double effMu8dataTrailing = SF_muon8->get_EfficiencyData(ptTrailing,etaTrailing); 
       double effMu17dataLeading = SF_muon17->get_EfficiencyData(ptLeading,etaLeading); 
       double effMu8dataLeading = SF_muon8->get_EfficiencyData(ptLeading,etaLeading); 

       double effMu17MCTrailing = SF_muon17->get_EfficiencyMC(ptTrailing,etaTrailing); 
       double effMu8MCTrailing = SF_muon8->get_EfficiencyMC(ptTrailing,etaTrailing); 
       double effMu17MCLeading = SF_muon17->get_EfficiencyMC(ptLeading,etaLeading); 
       double effMu8MCLeading = SF_muon8->get_EfficiencyMC(ptLeading,etaLeading); 

       if (era==2017) {

	 effMu17dataTrailing *= iso17DataTrailing;
	 effMu8dataTrailing  *= iso8DataTrailing;

	 effMu17dataLeading  *= iso17DataLeading;
	 effMu8dataLeading   *= iso8DataLeading;

	 effMu17MCTrailing *= iso17MCTrailing;
	 effMu8MCTrailing  *= iso8MCTrailing;

	 effMu17MCLeading  *= iso17MCLeading;
	 effMu8MCLeading   *= iso8MCLeading;

       }

       /*
       if (era==2017) {
	 std::cout << "iso(leading)  = " << isoLeading << std::endl;
	 std::cout << "iso(trailing) = " << isoTrailing << std::endl;
	 std::cout << "effIso17(leading) : data = " <<  iso17DataLeading 
		   << "  mc = " << iso17MCLeading << std::endl;
	 std::cout << "effIso8(leading) : data = " <<  iso8DataLeading 
		   << "  mc = " << iso8MCLeading << std::endl;
	 std::cout << "effIso17(trailing) : data = " <<  iso17DataTrailing 
		   << "  mc = " << iso17MCTrailing << std::endl;
	 std::cout << "effIso8(trailing) : data = " <<  iso8DataTrailing
		   << "  mc = " << iso8MCTrailing << std::endl;
       }
       */

       double trigWeightData = effMu17dataLeading*effMu8dataTrailing + effMu17dataTrailing*effMu8dataLeading - effMu17dataLeading*effMu17dataTrailing;
       double trigWeightMC = effMu17MCLeading*effMu8MCTrailing + effMu17MCTrailing*effMu8MCLeading - effMu17MCLeading*effMu17MCTrailing;

       if (applyTriggerMatch) {
	 if (trigWeightMC>0)
	   triggerWeight = trigWeightData/trigWeightMC;
       }
       else
	 triggerWeight = trigWeightData;

       triggerWeight *= effDzSS;

       /*       
       std::cout << "pT(leading) = " << ptLeading
		 << "   eta(leading) = " << etaLeading
		 << "   pT(trailing) = " << ptTrailing
		 << "   eta(trailing) = " << etaTrailing << std::endl;
       std::cout << "IdW(leading) = " << idLeadingW
		 << "   TrkW(leading) = " << trkLeadingW 
		 << "   IdW(trailing) = " << idTrailingW
		 << "   TrkW(trailing) = " << trkTrailingW << std::endl;
       std::cout << "Trigger weight = " << triggerWeight << std::endl;
       */

     }
     triggerWeightH->Fill(triggerWeight,1.0);
     weight *= triggerWeight;
     weight *= idLeadingWeight;
     weight *= idTrailingWeight;

     //     std::cout << "Ok 6" << std::endl;
       

     // dimuon selection passed 
     TLorentzVector LeadingMuon4; LeadingMuon4.SetXYZM(muon_px[iLeading],
						       muon_py[iLeading],
						       muon_pz[iLeading],
						       MuMass);
     
     TLorentzVector TrailingMuon4; TrailingMuon4.SetXYZM(muon_px[iTrailing],
							 muon_py[iTrailing],
							 muon_pz[iTrailing],
							 MuMass);

     // uncorrected 4-vectors (needed for track selection)
     TLorentzVector LeadingMuon4_uncorr; LeadingMuon4_uncorr.SetXYZM(muon_px_uncorr[iLeading],
								     muon_py_uncorr[iLeading],
								     muon_pz_uncorr[iLeading],
								     MuMass);
     
     TLorentzVector TrailingMuon4_uncorr; TrailingMuon4_uncorr.SetXYZM(muon_px_uncorr[iTrailing],
								       muon_py_uncorr[iTrailing],
								       muon_pz_uncorr[iTrailing],
								       MuMass);
     
     
     TLorentzVector diMuon4 = LeadingMuon4 + TrailingMuon4;
     
     float dimuonMass = diMuon4.M();

     // trigger weight (obsolete)
     //     if (!isData) {
     //       float effMu17Leading  = MuLegEfficiency(LeadingMuon4.Pt(),LeadingMuon4.Eta(),17.0,2.4);
     //       float effMu8Leading   = MuLegEfficiency(LeadingMuon4.Pt(),LeadingMuon4.Eta(), 8.0,2.4);
     //       float effMu17Trailing = MuLegEfficiency(TrailingMuon4.Pt(),TrailingMuon4.Eta(),17.0,2.4);
     //       float effMu8Trailing  = MuLegEfficiency(TrailingMuon4.Pt(),TrailingMuon4.Eta(), 8.0,2.4);
     //       float trigWeight = 0.935*(effMu17Leading*effMu8Trailing+effMu8Leading*effMu17Trailing-effMu17Leading*effMu17Trailing);
     //       std::cout << "trig weight = " << trigWeight << std::endl;
     //       weight *= trigWeight; 
     //     }

     dimuonMassH->Fill(dimuonMass,weight);
     ptLeadingMuH->Fill(muon_pt[iLeading],weight);
     ptTrailingMuH->Fill(muon_pt[iTrailing],weight);
     etaLeadingMuH->Fill(muon_eta[iLeading],weight);
     etaTrailingMuH->Fill(muon_eta[iTrailing],weight);
     counter_MuonKinematicsH->Fill(1.0,weight); 

     float iso1 = muon_r03_sumChargedHadronPt[iLeading]/muon_pt[iLeading];
     float iso2 = muon_r03_sumChargedHadronPt[iTrailing]/muon_pt[iTrailing];

     if (iso1<0.4&&iso2<0.4) {
       dimuonMassIsoH->Fill(dimuonMass,weight);
       ptLeadingIsoMuH->Fill(muon_pt[iLeading],weight);
       ptTrailingIsoMuH->Fill(muon_pt[iLeading],weight);
     }

     if (iso1>0.4&&iso2>0.4) {
       dimuonMassNonIsoH->Fill(dimuonMass,weight);
       ptLeadingNonIsoMuH->Fill(muon_pt[iLeading],weight);
       ptTrailingNonIsoMuH->Fill(muon_pt[iLeading],weight);
     }


     TLorentzVector Met4; Met4.SetXYZM(metx,mety,0,0);


     // counting tracks around each muon
     std::vector<unsigned int> trkLeadingMu; trkLeadingMu.clear(); // all tracks
     std::vector<unsigned int> trkTrailingMu; trkTrailingMu.clear(); // all tracks
     std::vector<unsigned int> trkSigLeadingMu; trkSigLeadingMu.clear(); // signal tracks
     std::vector<unsigned int> trkSigTrailingMu; trkSigTrailingMu.clear(); // signal tracks
     std::vector<double> trkSigDRLeadingMu; trkSigDRLeadingMu.clear();
     std::vector<double> trkSigDRTrailingMu; trkSigDRTrailingMu.clear();
     unsigned int hardestTrkLeading = 0; // index of hardest track around leading mu
     unsigned int hardestTrkTrailing = 0; // index of hardest track around trailing mu
     unsigned int softestTrkLeading = 0; // index of softest track around leading mu
     unsigned int softestTrkTrailing = 0; // index of softest track around trailing mu
     

     float ptHardestLeading = 0;
     float ptHardestTrailing = 0;
     float ptSoftestLeading = 1e+10;
     float ptSoftestTrailing = 1e+10;
     std::vector<unsigned int> Soft_trkLeadingMu; Soft_trkLeadingMu.clear(); // Lead Muon tracks control region
     std::vector<unsigned int> Soft_trkTrailingMu; Soft_trkTrailingMu.clear(); // Trail Muon tracks control region

     
     for (unsigned int iTrk=0; iTrk<track_count; ++iTrk) {
       if (fabs(track_charge[iTrk])<0.1) continue; // make sure we are not taking neutral stuff
       if (fabs(track_dxy[iTrk])>dxyTrkLooseCut) continue;
       if (fabs(track_dz[iTrk])>dzTrkLooseCut) continue;
       if (fabs(track_eta[iTrk])>etaTrkCut) continue;
       if (fabs(track_pt[iTrk])<ptTrkLooseCut) continue;
       
       TLorentzVector trk4; trk4.SetXYZM(track_px[iTrk],
					 track_py[iTrk],
					 track_pz[iTrk],
					 track_mass[iTrk]);
       
       TLorentzVector leadingMuDiff = LeadingMuon4_uncorr - trk4;
       if (leadingMuDiff.P()>0.1) { // track is not leading muon
	 float drTrkMu = deltaR(muon_eta[iLeading],muon_phi[iLeading],
				track_eta[iTrk],   track_phi[iTrk]);
	 float qTrkLeadingMu = track_charge[iTrk]*muon_charge[iLeading];
	 if (drTrkMu<dRIsoMuon){
	   trkLeadingMu.push_back(iTrk);
	   if (track_pt[iTrk]>ptTrkLooseCut && track_pt[iTrk]< ptTrkCut)
	     Soft_trkLeadingMu.push_back(iTrk);
	 }
	 if (drTrkMu<dRIsoMuon && qTrkLeadingMu<0 && fabs(track_dxy[iTrk])<dxyTrkCut && fabs(track_dz[iTrk])<dzTrkCut && track_pt[iTrk]>ptTrkCut) {
	   trkSigLeadingMu.push_back(iTrk);
	   trkSigDRLeadingMu.push_back(drTrkMu);
	   if (track_pt[iTrk]>ptHardestLeading) {
	     ptHardestLeading = track_pt[iTrk];
	     hardestTrkLeading = iTrk;
	   }
	   if (track_pt[iTrk]<ptSoftestLeading) {
	     ptSoftestLeading = track_pt[iTrk];
	     softestTrkLeading = iTrk;
	   }
	 }
       }
       
       TLorentzVector trailingMuDiff = TrailingMuon4_uncorr - trk4;
       if (trailingMuDiff.P()>0.1) { // track is not trailing muon
	 float drTrkMu = deltaR(muon_eta[iTrailing],muon_phi[iTrailing],
				track_eta[iTrk],track_phi[iTrk]);
	 float qTrkTrailingMu = track_charge[iTrk]*muon_charge[iTrailing];         
	 if (drTrkMu<dRIsoMuon){
	   trkTrailingMu.push_back(iTrk);
	   if (track_pt[iTrk] > ptTrkLooseCut && track_pt[iTrk]< ptTrkCut)
	     Soft_trkTrailingMu.push_back(iTrk);
	 }
	 if (drTrkMu<dRIsoMuon && qTrkTrailingMu<0 && fabs(track_dxy[iTrk])<dxyTrkCut && fabs(track_dz[iTrk])<dzTrkCut && track_pt[iTrk]>ptTrkCut) {
	   trkSigTrailingMu.push_back(iTrk);
	   trkSigDRTrailingMu.push_back(drTrkMu);
	   if (track_pt[iTrk]>ptHardestTrailing) {
	     ptHardestTrailing = track_pt[iTrk];
	     hardestTrkTrailing = iTrk;
	   }
	   if (track_pt[iTrk]<ptSoftestTrailing) {
	     ptSoftestTrailing = track_pt[iTrk];
	     softestTrkTrailing = iTrk;
	   }
	   
	 }
       }
       
     }
     
     nTracksLeadingMuH->Fill(float(trkLeadingMu.size()),weight);
     nTracksTrailingMuH->Fill(float(trkTrailingMu.size()),weight);
     nSigTracksLeadingMuH->Fill(float(trkSigLeadingMu.size()),weight);
     nSigTracksTrailingMuH->Fill(float(trkSigTrailingMu.size()),weight);
     nSoftTracksLeadingMuH->Fill(float(Soft_trkLeadingMu.size()),weight);
     nSoftTracksTrailingMuH->Fill(float(Soft_trkTrailingMu.size()),weight);

     float chargedIsoLeading = muon_chargedHadIso[iLeading];
     for (unsigned int iTrk=0; iTrk<trkSigLeadingMu.size(); ++iTrk) {
       unsigned int indexTrack = trkSigLeadingMu.at(iTrk);       
       if (trkSigDRLeadingMu.at(iTrk)<0.4) {
	 TLorentzVector trk4; trk4.SetXYZM(track_px[indexTrack],
					   track_py[indexTrack],
					   track_pz[indexTrack],
					   track_mass[indexTrack]);
	 chargedIsoLeading -= trk4.Pt();
       }
     }
     if (chargedIsoLeading<0.0) chargedIsoLeading = 0;
     float neutralIsoLeading = muon_neutralHadIso[iLeading] + muon_photonIso[iLeading] - 0.5*muon_puIso[iLeading];
     if (neutralIsoLeading<0) neutralIsoLeading = 0;
     float isoLeading = (chargedIsoLeading+neutralIsoLeading)/muon_pt[iLeading];
     
     float chargedIsoTrailing = muon_chargedHadIso[iTrailing];
     for (unsigned int iTrk=0; iTrk<trkSigTrailingMu.size(); ++iTrk) {
       unsigned int indexTrack = trkSigTrailingMu.at(iTrk);       
       if (trkSigDRTrailingMu.at(iTrk)<0.4) {
	 TLorentzVector trk4; trk4.SetXYZM(track_px[indexTrack],
					   track_py[indexTrack],
					   track_pz[indexTrack],
					   track_mass[indexTrack]);
	 chargedIsoTrailing -= trk4.Pt();
       }
     }
     if (chargedIsoTrailing<0.0) chargedIsoTrailing = 0;
     float neutralIsoTrailing = muon_neutralHadIso[iTrailing] + muon_photonIso[iTrailing] - 0.5*muon_puIso[iTrailing];
     if (neutralIsoTrailing<0) neutralIsoTrailing = 0;
     float isoTrailing = (chargedIsoTrailing+neutralIsoTrailing)/muon_pt[iTrailing];

     // ********************
     // Check isolations -->
     // ********************    
     /*
     if (trkSigLeadingMu.size()==1) {
       unsigned int iTrkLeading  = trkSigLeadingMu.at(0);
       double ptRatioLeading = track_pt[iTrkLeading]/muon_pt[iLeading];
       double dRTrkMuonLeading = deltaR(muon_eta[iLeading],muon_phi[iLeading],
					track_eta[iTrkLeading],track_phi[iTrkLeading]);
       if (ptRatioLeading>0.5) {
	 printf("Leading ---->\n");
	 printf(" ptRatio = %5.2f   dR(trk,mu) = %5.2f\n",ptRatioLeading,dRTrkMuonLeading);
	 printf(" muon  : pT  = %5.2f  eta = %5.2f  phi = %5.2f  q = %2i\n",
		muon_pt[iLeading],muon_eta[iLeading],muon_phi[iLeading],int(muon_charge[iLeading]));
	 printf(" track : pT  = %5.2f  eta = %5.2f  phi = %5.2f  q = %2i  ID = %4i\n",
		track_pt[iTrkLeading],
		track_eta[iTrkLeading],
		track_phi[iTrkLeading],
		int(track_charge[iTrkLeading]),
		track_ID[iTrkLeading]);
	 printf(" iso   : had = %5.2f\n",muon_chargedHadIso[iLeading]);
	 printf(" ptHad : 0p3 = %5.2f  0p4 = %5.2f\n",
		muon_r03_sumChargedHadronPt[iLeading],
		muon_r04_sumChargedHadronPt[iLeading]);
	 printf(" ptAll : 0p3 = %5.2f  0p4 = %5.2f\n",
		muon_r03_sumChargedParticlePt[iLeading],
		muon_r04_sumChargedParticlePt[iLeading]);
	 for (unsigned int im=0; im<muon_count; ++im) {
	   printf(" muon %2i -> pT = %5.2f  eta = %5.2f  phi = %5.2f iso = %5.2f\n",
		  im,muon_pt[im],muon_eta[im],muon_phi[im],muon_r03_sumChargedHadronPt[im]);
	 }
	 unsigned int trigobj = 0;
	 for (unsigned int it=0; it<trigobject_count; ++it) {
	   bool trigObj_Mu17 = trigobject_filters[it][nMu17Leg];
	   bool trigObj_Mu8  = trigobject_filters[it][nMu8Leg];
	   bool trigObj_SS   = trigobject_filters[it][nSSFilter];
	   bool trigObj_DZ   = trigobject_filters[it][nDZFilter];
	   bool validObject  = (trigObj_Mu17||trigObj_Mu8)&&trigObj_SS&&trigObj_DZ;
	   if (validObject) {
	     printf(" trig %2i -> pT = %5.2f  eta = %5.2f  phi = %5.2f\n",
		    trigobj,trigobject_pt[it],trigobject_eta[it],trigobject_phi[it]);
	     trigobj++;
	   }
	 }
	 std::cout << std::endl;
       }     
     }
    
     if (trkSigTrailingMu.size()==1) {
       unsigned int iTrkTrailing = trkSigTrailingMu.at(0); 
       double ptRatioTrailing = track_pt[iTrkTrailing]/muon_pt[iTrailing];
       double dRTrkMuonTrailing = deltaR(muon_eta[iTrailing],muon_phi[iTrailing],
					track_eta[iTrkTrailing],track_phi[iTrkTrailing]);
       if (ptRatioTrailing>0.5) {
	 printf("Trailing ---->\n");
	 printf(" ptRatio = %5.2f   dR(trk,mu) = %5.2f\n",ptRatioTrailing,dRTrkMuonTrailing);
	 printf(" muon  : pT  = %5.2f  eta = %5.2f  phi = %5.2f  q = %2i\n",
		muon_pt[iTrailing],muon_eta[iTrailing],muon_phi[iTrailing],int(muon_charge[iTrailing]));
	 printf(" track : pT  = %5.2f  eta = %5.2f  phi = %5.2f  q = %2i  ID = %4i\n",
		track_pt[iTrkTrailing],
		track_eta[iTrkTrailing],
		track_phi[iTrkTrailing],
		int(track_charge[iTrkTrailing]),
		track_ID[iTrkTrailing]);
	 printf(" iso   : had = %5.2f\n",muon_chargedHadIso[iTrailing]);
	 printf(" ptHad : 0p3 = %5.2f  0p4 = %5.2f\n",
		muon_r03_sumChargedHadronPt[iTrailing],
		muon_r04_sumChargedHadronPt[iTrailing]);
	 printf(" ptAll : 0p3 = %5.2f  0p4 = %5.2f\n",
		muon_r03_sumChargedParticlePt[iTrailing],
		muon_r04_sumChargedParticlePt[iTrailing]);
	 for (unsigned int im=0; im<muon_count; ++im) {
	   printf(" muon %2i -> pT = %5.2f  eta = %5.2f  phi = %5.2f iso = %5.2f\n",
		  im,muon_pt[im],muon_eta[im],muon_phi[im],muon_r03_sumChargedHadronPt[im]);
	 }
	 unsigned int trigobj = 0;
	 for (unsigned int it=0; it<trigobject_count; ++it) {
	   bool trigObj_Mu17 = trigobject_filters[it][nMu17Leg];
	   bool trigObj_Mu8  = trigobject_filters[it][nMu8Leg];
	   bool trigObj_SS   = trigobject_filters[it][nSSFilter];
	   bool trigObj_DZ   = trigobject_filters[it][nDZFilter];
	   bool validObject  = (trigObj_Mu17||trigObj_Mu8)&&trigObj_SS&&trigObj_DZ;
	   if (validObject) {
	     printf(" trig %2i -> pT = %5.2f  eta = %5.2f  phi = %5.2f\n",
		    trigobj,trigobject_pt[it],trigobject_eta[it],trigobject_phi[it]);
	     trigobj++;
	   }
	 }
	 std::cout << std::endl;
       }     
     }
     */
     float isoTrkLeading = muon_r03_sumChargedParticlePt[iLeading]/muon_pt[iLeading];
     float isoHadLeading = muon_r03_sumChargedHadronPt[iLeading]/muon_pt[iLeading];

     float isoTrkTrailing = muon_r03_sumChargedParticlePt[iTrailing]/muon_pt[iTrailing];
     float isoHadTrailing = muon_r03_sumChargedHadronPt[iTrailing]/muon_pt[iTrailing];
     
     IsoTrkLeadingMuH->Fill(isoTrkLeading,weight);
     IsoTrkTrailingMuH->Fill(isoTrkTrailing,weight);
     IsoTrkMuH->Fill(isoTrkLeading,weight);
     IsoTrkMuH->Fill(isoTrkTrailing,weight);

     IsoHadLeadingMuH->Fill(isoHadLeading,weight);
     IsoHadTrailingMuH->Fill(isoHadTrailing,weight);
     IsoHadMuH->Fill(isoHadLeading,weight);
     IsoHadMuH->Fill(isoHadTrailing,weight);

     if (trkSigLeadingMu.size()==1) {
       unsigned int iTrk = trkSigLeadingMu.at(0);
       int id = TMath::Abs(track_ID[iTrk]);
       float pttrk = track_pt[iTrk];
       float ptratio = pttrk/muon_pt[iLeading];
       ptRatioTrkLeadingMuH->Fill(ptratio,weight);
       ptRatioTrkMuH->Fill(ptratio,weight);
       if (id!=13) { 
	 ptRatioHadLeadingMuH->Fill(ptratio,weight);
	 ptRatioHadMuH->Fill(ptratio,weight);	 
       }
     }

     if (trkSigTrailingMu.size()==1) {
       unsigned int iTrk = trkSigTrailingMu.at(0);
       int id = TMath::Abs(track_ID[iTrk]);
       float pttrk = track_pt[iTrk];
       float ptratio = pttrk/muon_pt[iTrailing];
       ptRatioTrkTrailingMuH->Fill(ptratio,weight);
       ptRatioTrkMuH->Fill(ptratio,weight);
       if (id!=13) { 
	 ptRatioHadTrailingMuH->Fill(ptratio,weight);
	 ptRatioHadMuH->Fill(ptratio,weight);	 
       }
     }

     // definition of signal muon+track
     bool signalLeadingMu  = trkLeadingMu.size()==1 && trkSigLeadingMu.size()==1;
     bool signalTrailingMu = trkTrailingMu.size()==1 && trkSigTrailingMu.size()==1;

     double weightTrkLeading = 1;
     double weightTrkTrailing = 1;
     if (trkSigLeadingMu.size()>0&&!isData) {
       unsigned int iTrkLeading = trkSigLeadingMu.at(0);
       int absPdgId = TMath::Abs(track_ID[iTrkLeading]);
       if (absPdgId==11) {
	 correctionWS->var("e_pt")->setVal(track_pt[iTrkLeading]);
	 correctionWS->var("e_eta")->setVal(track_eta[iTrkLeading]);
	 weightTrkLeading = correctionWS->function("e_trk_ratio")->getVal();
       }
       else {
	 correctionWS->var("m_pt")->setVal(track_pt[iTrkLeading]);
         correctionWS->var("m_eta")->setVal(track_eta[iTrkLeading]);
         weightTrkLeading = correctionWS->function("m_trk_ratio")->getVal();
       }
       //       std::cout << "track around leading : pT = " << track_pt[iTrkLeading]
       //		 << "  eta = " << track_eta[iTrkLeading] 
       //		 << "  Id = " << track_ID[iTrkLeading]
       //		 << "  weight(trk) = " << weightTrkLeading << std::endl;
	
     }
     weight *= weightTrkLeading;

     if (trkSigTrailingMu.size()>0&&!isData) {
       unsigned int iTrkTrailing = trkSigTrailingMu.at(0);
       int absPdgId = TMath::Abs(track_ID[iTrkTrailing]);
       if (absPdgId==11) {
         correctionWS->var("e_pt")->setVal(track_pt[iTrkTrailing]);
         correctionWS->var("e_eta")->setVal(track_eta[iTrkTrailing]);
         weightTrkTrailing = correctionWS->function("e_trk_ratio")->getVal();
       }
       else {
         correctionWS->var("m_pt")->setVal(track_pt[iTrkTrailing]);
         correctionWS->var("m_eta")->setVal(track_eta[iTrkTrailing]);
         weightTrkTrailing = correctionWS->function("m_trk_ratio")->getVal();
       }
       //       std::cout << "track around trailing : pT = " << track_pt[iTrkTrailing]
       //                 << "  eta = " << track_eta[iTrkTrailing]
       //                 << "  Id = " << track_ID[iTrkTrailing]
       //                 << "  weight(trk) = " << weightTrkTrailing << std::endl;
     }
     //     std::cout << std::endl;
     weight *= weightTrkTrailing;


     if (trkLeadingMu.size()==1) {

       unsigned int iTrkLeading = trkLeadingMu.at(0);

       ptTrackLeadingMuH->Fill(track_pt[iTrkLeading],weight);
       etaTrackLeadingMuH->Fill(track_eta[iTrkLeading],weight);
       dxyTrackLeadingMuH->Fill(track_dxy[iTrkLeading],weight);
       dzTrackLeadingMuH->Fill(track_dz[iTrkLeading],weight);

       ptTrackN1H->Fill(track_pt[iTrkLeading],weight);
       etaTrackN1H->Fill(track_eta[iTrkLeading],weight);
       dxyTrackN1H->Fill(track_dxy[iTrkLeading],weight);
       dzTrackN1H->Fill(track_dz[iTrkLeading],weight);
       
     }

     if (trkTrailingMu.size()==1) {

       unsigned int iTrkTrailing = trkTrailingMu.at(0);

       ptTrackTrailingMuH->Fill(track_pt[iTrkTrailing],weight);
       etaTrackTrailingMuH->Fill(track_eta[iTrkTrailing],weight);
       dxyTrackTrailingMuH->Fill(track_dxy[iTrkTrailing],weight);
       dzTrackTrailingMuH->Fill(track_dz[iTrkTrailing],weight);

       ptTrackN1H->Fill(track_pt[iTrkTrailing],weight);
       etaTrackN1H->Fill(track_eta[iTrkTrailing],weight);
       dxyTrackN1H->Fill(track_dxy[iTrkTrailing],weight);
       dzTrackN1H->Fill(track_dz[iTrkTrailing],weight);
       
     }

     if (signalLeadingMu) {
       unsigned int iTrkLeading = trkLeadingMu.at(0);

       ptTrackH->Fill(track_pt[iTrkLeading],weight);
       etaTrackH->Fill(track_eta[iTrkLeading],weight);
       dxyTrackH->Fill(track_dxy[iTrkLeading],weight);
       dzTrackH->Fill(track_dz[iTrkLeading],weight);

       counter_nMuTrackSigH->Fill(1.0,weight);                 

     }

     if (signalTrailingMu) {
       unsigned int iTrkTrailing = trkTrailingMu.at(0);
       
       ptTrackH->Fill(track_pt[iTrkTrailing],weight);
       etaTrackH->Fill(track_eta[iTrkTrailing],weight);
       dxyTrackH->Fill(track_dxy[iTrkTrailing],weight);
       dzTrackH->Fill(track_dz[iTrkTrailing],weight);

       counter_nMuTrackSigH->Fill(1.0,weight);                 

     }
     
     // sideband N23 and N45
     bool isN23leading  = (trkLeadingMu.size()==1&&trkSigLeadingMu.size()==1) && (trkTrailingMu.size()==2||trkTrailingMu.size()==3);
     bool isN45leading  = (trkLeadingMu.size()==1&&trkSigLeadingMu.size()==1) && (trkTrailingMu.size()==4||trkTrailingMu.size()==5);
     bool isN23trailing = (trkTrailingMu.size()==1&&trkSigTrailingMu.size()==1) && (trkLeadingMu.size()==2||trkLeadingMu.size()==3); 
     bool isN45trailing = (trkTrailingMu.size()==1&&trkSigTrailingMu.size()==1) && (trkLeadingMu.size()==4||trkLeadingMu.size()==5); 

     
     // sidebands Ntrk23
     bool isNtrk23leading  = trkSigLeadingMu.size()>0 && (trkTrailingMu.size()==2||trkTrailingMu.size()==3);
     bool isNtrk23trailing = trkSigTrailingMu.size()>0 && (trkLeadingMu.size()==2||trkLeadingMu.size()==3);
     
     // sidebands Ntrk1
     bool isNtrk1leading  = trkSigLeadingMu.size()>0 && trkTrailingMu.size()==1;
     bool isNtrk1trailing = trkSigTrailingMu.size()>0 && trkLeadingMu.size()==1;
     
     // signal region
     bool signalRegion = (trkLeadingMu.size()==1&&trkSigLeadingMu.size()==1) && (trkSigTrailingMu.size()==1&&trkTrailingMu.size()==1);
     
     // sidebands N23 (see definition of this sideband in HIG-14-019)
     if (isN23leading) {
       int iTrk = trkSigLeadingMu[0];
       TLorentzVector Track4; Track4.SetXYZM(track_px[iTrk],
					     track_py[iTrk],
					     track_pz[iTrk],
					     track_mass[iTrk]);
       TLorentzVector TrackPlusMuon4 = LeadingMuon4 + Track4;
       float deltaRMuTrk = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),
				  Track4.Eta(),Track4.Phi());
       float mass = TrackPlusMuon4.M();
       InvMassN23H->Fill(mass,weight);
     }


     if (isN45leading) {
       int iTrk = trkSigLeadingMu[0];
       TLorentzVector Track4; Track4.SetXYZM(track_px[iTrk],
					     track_py[iTrk],
					     track_pz[iTrk],
					     track_mass[iTrk]);
       TLorentzVector TrackPlusMuon4 = LeadingMuon4 + Track4;
       float deltaRMuTrk = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),
				  Track4.Eta(),Track4.Phi());
       float mass = TrackPlusMuon4.M();
       InvMassN45H->Fill(mass,weight);
     }
     
     if (isN23trailing) {
       int iTrk = trkSigTrailingMu[0];
       TLorentzVector Track4; Track4.SetXYZM(track_px[iTrk],
					     track_py[iTrk],
					     track_pz[iTrk],
					     track_mass[iTrk]);
       TLorentzVector TrackPlusMuon4 = TrailingMuon4 + Track4;
       float deltaRMuTrk = deltaR(TrailingMuon4.Eta(),TrailingMuon4.Phi(),
				  Track4.Eta(),Track4.Phi());
       float mass = TrackPlusMuon4.M();
       InvMassN23H->Fill(mass,weight);
     }      
     
     if (isN45trailing) {
       int iTrk = trkSigTrailingMu[0];
       TLorentzVector Track4; Track4.SetXYZM(track_px[iTrk],
					     track_py[iTrk],
					     track_pz[iTrk],
					     track_mass[iTrk]);
       TLorentzVector TrackPlusMuon4 = TrailingMuon4 + Track4;
       float deltaRMuTrk = deltaR(TrailingMuon4.Eta(),TrailingMuon4.Phi(),
				  Track4.Eta(),Track4.Phi());
       float mass = TrackPlusMuon4.M();
       InvMassN45H->Fill(mass,weight);
     }      
     

     // sidebands Ntrk23 (see definition of this sideband in HIG-14-019)
     if (isNtrk23leading) {
       TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(track_px[hardestTrkLeading],
							   track_py[hardestTrkLeading],
							   track_pz[hardestTrkLeading],
							   track_mass[hardestTrkLeading]);
       TLorentzVector HardestTrackPlusMuon4 = LeadingMuon4 + HardestTrack4;
       float mass = HardestTrackPlusMuon4.M();
       InvMassHardestNtrk23leadingH->Fill(mass,weight);
       InvMassHardestNtrk23H->Fill(mass,weight);
       TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(track_px[softestTrkLeading],
							   track_py[softestTrkLeading],
							   track_pz[softestTrkLeading],
							   track_mass[softestTrkLeading]);
       TLorentzVector SoftestTrackPlusMuon4 = LeadingMuon4 + SoftestTrack4;
       mass = SoftestTrackPlusMuon4.M();
       InvMassSoftestNtrk23leadingH->Fill(mass,weight);
       InvMassSoftestNtrk23H->Fill(mass,weight);
     }      
     if (isNtrk23trailing) {
       TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(track_px[hardestTrkTrailing],
							   track_py[hardestTrkTrailing],
							   track_pz[hardestTrkTrailing],
							   track_mass[hardestTrkTrailing]);
       TLorentzVector HardestTrackPlusMuon4 = TrailingMuon4 + HardestTrack4;
       float mass = HardestTrackPlusMuon4.M();
       InvMassHardestNtrk23trailingH->Fill(mass,weight);
       InvMassHardestNtrk23H->Fill(mass,weight);
       TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(track_px[softestTrkTrailing],
							   track_py[softestTrkTrailing],
							   track_pz[softestTrkTrailing],
							   track_mass[softestTrkTrailing]);
       TLorentzVector SoftestTrackPlusMuon4 = TrailingMuon4 + SoftestTrack4;
       mass = SoftestTrackPlusMuon4.M();
       InvMassSoftestNtrk23trailingH->Fill(mass,weight);
       InvMassSoftestNtrk23H->Fill(mass,weight);
     }      
      
     // sidebands Ntrk1 (see definition of this sideband in HIG-14-019)
     if (isNtrk1leading) {
       TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(track_px[hardestTrkLeading],
							   track_py[hardestTrkLeading],
							   track_pz[hardestTrkLeading],
							   track_mass[hardestTrkLeading]);
       TLorentzVector HardestTrackPlusMuon4 = LeadingMuon4 + HardestTrack4;
       float mass = HardestTrackPlusMuon4.M();
       InvMassHardestNtrk1leadingH->Fill(mass,weight);
       InvMassHardestNtrk1H->Fill(mass,weight);
       TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(track_px[softestTrkLeading],
							   track_py[softestTrkLeading],
							   track_pz[softestTrkLeading],
							   track_mass[softestTrkLeading]);
       TLorentzVector SoftestTrackPlusMuon4 = LeadingMuon4 + SoftestTrack4;
       mass = SoftestTrackPlusMuon4.M();
       InvMassSoftestNtrk1leadingH->Fill(mass,weight);
       InvMassSoftestNtrk1H->Fill(mass,weight);
     }      
     if (isNtrk1trailing) {
       TLorentzVector HardestTrack4; HardestTrack4.SetXYZM(track_px[hardestTrkTrailing],
							   track_py[hardestTrkTrailing],
							   track_pz[hardestTrkTrailing],
							   track_mass[hardestTrkTrailing]);
       TLorentzVector HardestTrackPlusMuon4 = TrailingMuon4 + HardestTrack4;
       float mass = HardestTrackPlusMuon4.M();
       InvMassHardestNtrk1trailingH->Fill(mass,weight);
       InvMassHardestNtrk1H->Fill(mass,weight);
       TLorentzVector SoftestTrack4; SoftestTrack4.SetXYZM(track_px[softestTrkTrailing],
							   track_py[softestTrkTrailing],
							   track_pz[softestTrkTrailing],
							   track_mass[softestTrkTrailing]);
       TLorentzVector SoftestTrackPlusMuon4 = TrailingMuon4 + SoftestTrack4;
       mass = SoftestTrackPlusMuon4.M();
       InvMassSoftestNtrk1trailingH->Fill(mass,weight);
       InvMassSoftestNtrk1H->Fill(mass,weight);
     }      

     // **************
     // signal region
     // *************
     if (signalRegion) {

       //       if (vetoEvent) {
       //       	 std::cout << "Event is rejected by HEM " << std::endl; 
       //       }
       //       else {
       //       	 std::cout << "Event is selected... " << std::endl;
       //       }

       counter_FinalEventsH->Fill(1.,weight);

       counter_btagCorrUp->Fill(1.,weight*weight_btag_corr_up);
       counter_btagCorrDown->Fill(1.,weight*weight_btag_corr_down);
       counter_mistagCorrUp->Fill(1.,weight*weight_mistag_corr_up);
       counter_mistagCorrDown->Fill(1.,weight*weight_mistag_corr_down);

       counter_btagUncorrUp->Fill(1.,weight*weight_btag_uncorr_up);
       counter_btagUncorrDown->Fill(1.,weight*weight_btag_uncorr_down);
       counter_mistagUncorrUp->Fill(1.,weight*weight_mistag_uncorr_up);
       counter_mistagUncorrDown->Fill(1.,weight*weight_mistag_uncorr_down);

       float prefireUp = 1.0;
       float prefireDown = 1.0;
       if (prefiringweight>0.01) {
	 prefireUp = prefiringweightup/prefiringweight;
	 prefireDown = prefiringweightdown/prefiringweight;
       }

       counter_prefireUp->Fill(1.,weight*prefireUp);
       counter_prefireDown->Fill(1.,weight*prefireDown);

       int iTrkLeading = trkSigLeadingMu[0];
       TLorentzVector TrackLeading4; TrackLeading4.SetXYZM(track_px[iTrkLeading],
							   track_py[iTrkLeading],
							   track_pz[iTrkLeading],
							   track_mass[iTrkLeading]);
       TLorentzVector MuonTrackLeading4 = LeadingMuon4 + TrackLeading4;
       float massTrkMuLeading = MuonTrackLeading4.M();
       
       int iTrkTrailing = trkSigTrailingMu[0];
       TLorentzVector TrackTrailing4; TrackTrailing4.SetXYZM(track_px[iTrkTrailing],
							     track_py[iTrkTrailing],
							     track_pz[iTrkTrailing],
							     track_mass[iTrkTrailing]);
       TLorentzVector MuonTrackTrailing4 = TrailingMuon4 + TrackTrailing4;
       float massTrkMuTrailing = MuonTrackTrailing4.M();
       
       TLorentzVector Visible4 = MuonTrackTrailing4 + MuonTrackLeading4;
       TLorentzVector Total4 = Visible4 + Met4;

       float dRLeadingMuTrk = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),
				     TrackLeading4.Eta(),TrackLeading4.Phi());


       float dRTrailingMuTrk = deltaR(TrailingMuon4.Eta(),TrailingMuon4.Phi(),
				      TrackTrailing4.Eta(),TrackTrailing4.Phi());

       // *************************
       // applying track id/iso SF
       // *************************

       double weight_trk_leading_iso_Up = 1;
       double weight_trk_leading_iso_Down = 1;
       double weight_trk_trailing_iso_Up = 1;
       double weight_trk_trailing_iso_Down = 1;
       double weight_trk_iso_Up = 1;
       double weight_trk_iso_Down = 1;

       if (!isData) {

         double SFIsoLeading = trackIsoSF;
         double SFIsoTrailing = trackIsoSF;

         double derivativeLeading[2];
         derivativeLeading[0] = 1;
         derivativeLeading[1] = TrackLeading4.Pt();

	 double derivativeTrailing[2];
	 derivativeTrailing[0] = 1;
	 derivativeTrailing[1] = TrackTrailing4.Pt();	 

         double ParError[2];
         ParError[0] = InterceptErr;
         ParError[1] = SlopeErr;

         double CorrMatrix[2][2] = {{1.0, Correlation}, {Correlation, 1.0}};

	 double sumdeltaLeading = 0; double sumdeltaTrailing = 0;

         for(int npar=0;npar<2;npar++) {
             for(int mpar=0;mpar<2;mpar++) {
                  sumdeltaLeading = sumdeltaLeading + derivativeLeading[npar]*derivativeLeading[mpar]*ParError[npar]*ParError[mpar]*CorrMatrix[npar][mpar];
                  sumdeltaTrailing = sumdeltaTrailing + derivativeTrailing[npar]*derivativeTrailing[mpar]*ParError[npar]*ParError[mpar]*CorrMatrix[npar][mpar];
             }
         }

         double DeltaSFIsoLeading = TMath::Sqrt(sumdeltaLeading);
         double DeltaSFIsoTrailing = TMath::Sqrt(sumdeltaTrailing);

         weight_trk_leading_iso_Up = SFIsoLeading + DeltaSFIsoLeading;
         weight_trk_leading_iso_Down = SFIsoLeading - DeltaSFIsoLeading;
         weight_trk_trailing_iso_Up = SFIsoTrailing + DeltaSFIsoTrailing;
         weight_trk_trailing_iso_Down = SFIsoTrailing - DeltaSFIsoTrailing;

         weight_trk_iso_Up = weight_trk_leading_iso_Up * weight_trk_trailing_iso_Up;
         weight_trk_iso_Down = weight_trk_leading_iso_Down * weight_trk_trailing_iso_Down;

         weight *= SFIsoLeading*SFIsoTrailing;

       } 

       MetSelH->Fill(Met4.Pt(),weight);
       mTtotSelH->Fill(Total4.M(),weight);
       dimuonMassSelH->Fill(diMuon4.M(),weight);
       invMass2Mu2TrkSelH->Fill(Visible4.M(),weight);

       ptLeadingMuSelH->Fill(LeadingMuon4.Pt(),weight);
       ptTrailingMuSelH->Fill(TrailingMuon4.Pt(),weight);
       ptLeadingMuTrkSelH->Fill(TrackLeading4.Pt(),weight);
       ptTrailingMuTrkSelH->Fill(TrackTrailing4.Pt(),weight);

       ptMuSelH->Fill(LeadingMuon4.Pt(),weight);
       ptMuSelH->Fill(TrailingMuon4.Pt(),weight);
       ptTrkSelH->Fill(TrackLeading4.Pt(),weight);
       ptTrkSelH->Fill(TrackTrailing4.Pt(),weight);

       etaLeadingMuSelH->Fill(LeadingMuon4.Eta(),weight);
       etaTrailingMuSelH->Fill(TrailingMuon4.Eta(),weight);
       etaLeadingMuTrkSelH->Fill(TrackLeading4.Eta(),weight);
       etaTrailingMuTrkSelH->Fill(TrackTrailing4.Eta(),weight);

       etaMuSelH->Fill(LeadingMuon4.Eta(),weight);
       etaMuSelH->Fill(TrailingMuon4.Eta(),weight);
       
       etaTrkSelH->Fill(TrackLeading4.Eta(),weight);
       etaTrkSelH->Fill(TrackTrailing4.Eta(),weight);

       dRMuTrkSelH->Fill(dRLeadingMuTrk,weight);
       dRMuTrkSelH->Fill(dRTrailingMuTrk,weight);

       float masshigh = massTrkMuLeading;
       float masslow = massTrkMuTrailing;
       
       if (masshigh<masslow) {
	 masshigh = massTrkMuTrailing;
	 masslow = massTrkMuLeading;
       }

       InvMassLeadingH->Fill(massTrkMuLeading,weight);
       InvMassTrailingH->Fill(massTrkMuTrailing,weight);
       InvMassH->Fill(massTrkMuLeading,weight);
       InvMassH->Fill(massTrkMuTrailing,weight);
       InvMass2DH->Fill(masslow,masshigh, weight);

       // btag variations for heavy flavours : up ->
       InvMassH_btagUp->Fill(massTrkMuLeading,weight*weight_btag_up);
       InvMassH_btagUp->Fill(massTrkMuTrailing,weight*weight_btag_up);
       InvMass2DH_btagUp->Fill(masslow,masshigh,weight*weight_btag_up);
       // btag variations for heavy flavours : down ->
       InvMassH_btagDown->Fill(massTrkMuLeading,weight*weight_btag_down);
       InvMassH_btagDown->Fill(massTrkMuTrailing,weight*weight_btag_down);
       InvMass2DH_btagDown->Fill(masslow,masshigh,weight*weight_btag_down);

       // mistag variations for light flavours : up ->
       InvMassH_mistagUp->Fill(massTrkMuLeading,weight*weight_mistag_up);
       InvMassH_mistagUp->Fill(massTrkMuTrailing,weight*weight_mistag_up);
       InvMass2DH_mistagUp->Fill(masslow,masshigh,weight*weight_mistag_up);
       // mistag variations for light flavours : down ->
       InvMassH_mistagDown->Fill(massTrkMuLeading,weight*weight_mistag_down);
       InvMassH_mistagDown->Fill(massTrkMuTrailing,weight*weight_mistag_down);
       InvMass2DH_mistagDown->Fill(masslow,masshigh,weight*weight_mistag_down);

       // trk iso sf variations : up ->
       InvMassH_trkIsoUp->Fill(massTrkMuLeading,weight*weight_trk_leading_iso_Up);
       InvMassH_trkIsoUp->Fill(massTrkMuTrailing,weight*weight_trk_trailing_iso_Up);
       InvMass2DH_trkIsoUp->Fill(massTrkMuLeading, massTrkMuTrailing,weight*weight_trk_iso_Up);
       // trk iso sf variations : down ->
       InvMassH_trkIsoDown->Fill(massTrkMuLeading,weight*weight_trk_leading_iso_Down);
       InvMassH_trkIsoDown->Fill(massTrkMuTrailing,weight*weight_trk_trailing_iso_Down);
       InvMass2DH_trkIsoDown->Fill(masslow, masshigh,weight*weight_trk_iso_Down);

       // prefiring variations : up ->
       InvMassH_prefireUp->Fill(massTrkMuLeading,weight*prefireUp);
       InvMassH_prefireUp->Fill(massTrkMuTrailing,weight*prefireUp);
       InvMass2DH_prefireUp->Fill(masslow, masshigh,weight*prefireUp);
       // prefiring variations : down ->
       InvMassH_prefireDown->Fill(massTrkMuLeading,weight*prefireDown);
       InvMassH_prefireDown->Fill(massTrkMuTrailing,weight*prefireDown);
       InvMass2DH_prefireDown->Fill(masslow, masshigh,weight*prefireDown);
       
     }

     // ******************************************************************
     // Control region to study mass correlations (RegionA in HIG-14-019)
     // *****************************************************************
     
     bool bkgdLeadingMu = 
       (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==1 && trkLeadingMu.size()==2) ||
       (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==2 && trkLeadingMu.size()==3);
     
     bool bkgdTrailingMu = 
       (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==1 && trkTrailingMu.size()==2) ||
       (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==2 && trkTrailingMu.size()==3);

     bool bkgdXLeadingMu = 
       (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==1 && trkLeadingMu.size()==2) ||
       (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==2 && trkLeadingMu.size()==3) ||
       (trkSigLeadingMu.size()==1 && Soft_trkLeadingMu.size()==3 && trkLeadingMu.size()==4);
     
     bool bkgdXTrailingMu = 
       (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==1 && trkTrailingMu.size()==2) ||
       (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==2 && trkTrailingMu.size()==3) ||
       (trkSigTrailingMu.size()==1 && Soft_trkTrailingMu.size()==3 && trkTrailingMu.size()==4);

     bool ControlAll = (signalLeadingMu&&bkgdTrailingMu) || 
       (signalTrailingMu&&bkgdLeadingMu) || (bkgdLeadingMu&&bkgdTrailingMu);
     
     // * Now we can use this boolean to select bkg sideband
     // * where correlation coefficients are computed
     // * It is sufficient to use only one boolean - ControlAll
     if(ControlAll){
       // leading muon and associated track
       int iTrkLeading = trkSigLeadingMu[0];
       TLorentzVector TrackLeading4; TrackLeading4.SetXYZM(track_px[iTrkLeading],
							   track_py[iTrkLeading],
							   track_pz[iTrkLeading],
							   track_mass[iTrkLeading]);
       TLorentzVector TrackPlusLeadingMuon4 = LeadingMuon4 + TrackLeading4;
       float massLeadingMuonTrk = TrackPlusLeadingMuon4.M();
       
       // trailing muon and associated track
       int iTrkTrailing = trkSigTrailingMu[0];
       TLorentzVector TrackTrailing4; TrackTrailing4.SetXYZM(track_px[iTrkTrailing],
							     track_py[iTrkTrailing],
							     track_pz[iTrkTrailing],
							     track_mass[iTrkTrailing]);
       TLorentzVector TrackPlusTrailingMuon4 = TrailingMuon4 + TrackTrailing4;
       float massTrailingMuonTrk = TrackPlusTrailingMuon4.M();
       
       float masshigh = massLeadingMuonTrk;
       float masslow = massTrailingMuonTrk;
       
       if (masshigh<masslow) {
	 masshigh = massTrailingMuonTrk;
	 masslow = massLeadingMuonTrk;
       }
       
       // filling histograms
       InvMassTrackPlusMuon1D_ControlH->Fill(massLeadingMuonTrk,weight);
       InvMassTrackPlusMuon1D_ControlH->Fill(massTrailingMuonTrk,weight);
       InvMassTrackPlusMuon2D_ControlH->Fill(masslow, masshigh, weight);

       counter_ControlEventsH->Fill(1.0,weight);

     }

     // ********** ControlX ****************************
     bool ControlAllX = bkgdXLeadingMu&&bkgdXTrailingMu;
     if(ControlAllX){
       // leading muon and associated track
       int iTrkLeading = trkSigLeadingMu[0];
       TLorentzVector TrackLeading4; TrackLeading4.SetXYZM(track_px[iTrkLeading],
							   track_py[iTrkLeading],
							   track_pz[iTrkLeading],
							   track_mass[iTrkLeading]);
       TLorentzVector TrackPlusLeadingMuon4 = LeadingMuon4 + TrackLeading4;
       float massLeadingMuonTrk = TrackPlusLeadingMuon4.M();
       
       // trailing muon and associated track
       int iTrkTrailing = trkSigTrailingMu[0];
       TLorentzVector TrackTrailing4; TrackTrailing4.SetXYZM(track_px[iTrkTrailing],
							     track_py[iTrkTrailing],
							     track_pz[iTrkTrailing],
							     track_mass[iTrkTrailing]);
       TLorentzVector TrackPlusTrailingMuon4 = TrailingMuon4 + TrackTrailing4;
       float massTrailingMuonTrk = TrackPlusTrailingMuon4.M();
       
       float masshigh = massLeadingMuonTrk;
       float masslow = massTrailingMuonTrk;
       
       if (masshigh<masslow) {
	 masshigh = massTrailingMuonTrk;
	 masslow = massLeadingMuonTrk;
       }
       
       // filling histograms
       InvMassTrackPlusMuon1D_ControlXH->Fill(massLeadingMuonTrk,weight);
       InvMassTrackPlusMuon1D_ControlXH->Fill(massTrailingMuonTrk,weight);
       InvMassTrackPlusMuon2D_ControlXH->Fill(masslow, masshigh, weight);

       counter_ControlXEventsH->Fill(1.0,weight);

     }
     
     // ********* ControlY *********************
     bool ControlAllY = (signalLeadingMu&&bkgdXTrailingMu) ||
       (signalTrailingMu&&bkgdXLeadingMu) || (bkgdXLeadingMu&&bkgdXTrailingMu);
     if(ControlAllY){
       // leading muon and associated track
       int iTrkLeading = trkSigLeadingMu[0];
       TLorentzVector TrackLeading4; TrackLeading4.SetXYZM(track_px[iTrkLeading],
							   track_py[iTrkLeading],
							   track_pz[iTrkLeading],
							   track_mass[iTrkLeading]);
       TLorentzVector TrackPlusLeadingMuon4 = LeadingMuon4 + TrackLeading4;
       float massLeadingMuonTrk = TrackPlusLeadingMuon4.M();
       
       // trailing muon and associated track
       int iTrkTrailing = trkSigTrailingMu[0];
       TLorentzVector TrackTrailing4; TrackTrailing4.SetXYZM(track_px[iTrkTrailing],
							     track_py[iTrkTrailing],
							     track_pz[iTrkTrailing],
							     track_mass[iTrkTrailing]);
       TLorentzVector TrackPlusTrailingMuon4 = TrailingMuon4 + TrackTrailing4;
       float massTrailingMuonTrk = TrackPlusTrailingMuon4.M();

       float masshigh = massLeadingMuonTrk;
       float masslow = massTrailingMuonTrk;
       
       if (masshigh<masslow) {
	 masshigh = massTrailingMuonTrk;
	 masslow = massLeadingMuonTrk;
       }
       
       // filling histograms
       InvMassTrackPlusMuon1D_ControlYH->Fill(massLeadingMuonTrk,weight);
       InvMassTrackPlusMuon1D_ControlYH->Fill(massTrailingMuonTrk,weight);
       InvMassTrackPlusMuon2D_ControlYH->Fill(masslow, masshigh, weight);

       counter_ControlYEventsH->Fill(1.0,weight);

       // filling tuple if m1,m2 > 4 GeV
       if (massLeadingMuonTrk>4.0&&massTrailingMuonTrk>4.0) {
	 t_event = event_nr;
	 t_run = event_run;
	 t_nmu = 0;
	 for (auto muon : muons) {
	   if (t_nmu<100) {
	     t_mupt[t_nmu] = muon_pt[muon]; 
	     t_mueta[t_nmu] = muon_eta[muon];
	     t_muphi[t_nmu] = muon_phi[muon];
	     t_muq[t_nmu] = muon_charge[muon];
	     t_mudxy[t_nmu] = muon_dxy[muon];
	     t_mudz[t_nmu] = muon_dz[muon];
	     float chargedIso = muon_chargedHadIso[muon]; 
	     float neutralIso = muon_neutralHadIso[muon] + muon_photonIso[muon] - 0.5*muon_puIso[muon];
	     if (neutralIso<0) neutralIso = 0;
	     t_muiso[t_nmu] = (chargedIso+neutralIso)/muon_pt[muon];
	     int imu = muon; 
	     if (imu==iLeading) t_mu1_index = t_nmu;
	     if (imu==iTrailing) t_mu2_index = t_nmu;
	     t_nmu++;
	   }	 
	 }
	 
	 t_mu1_mutrkmass = massLeadingMuonTrk;
	 t_mu1_trkpt = TrackLeading4.Pt();
	 t_mu1_trketa = TrackLeading4.Eta();
	 t_mu1_trkphi = TrackLeading4.Phi();
	 t_mu1_nsoft = Soft_trkLeadingMu.size();
	 
	 t_mu2_mutrkmass = massTrailingMuonTrk;
	 t_mu2_trkpt = TrackTrailing4.Pt();
	 t_mu2_trketa = TrackTrailing4.Eta();
	 t_mu2_trkphi = TrackTrailing4.Phi();
	 t_mu2_nsoft = Soft_trkTrailingMu.size();

	 t_njets = 0;
	 t_njetspt20 = 0;       
	 
	 for (unsigned int jet=0; jet<pfjet_count; ++jet) {
	   
	   float absJetEta = TMath::Abs(pfjet_eta[jet]);
	   float jetPt = pfjet_pt[jet];
	   
	   bool jetId = tightJetID(pfjet_e[jet],
				   pfjet_eta[jet],
				   pfjet_neutralhadronicenergy[jet],
				   pfjet_neutralemenergy[jet],
				   pfjet_muonenergy[jet],
				   pfjet_chargedhadronicenergy[jet],
				   pfjet_chargedemenergy[jet],
				   pfjet_neutralmulti[jet],
				   pfjet_chargedmulti[jet],
				   era);
	   if (!jetId) continue;
	   float dR_leading = deltaR(LeadingMuon4.Eta(),LeadingMuon4.Phi(),
				     pfjet_eta[jet],pfjet_phi[jet]);
	   if (dR_leading<0.4) continue;
	   float dR_trailing = deltaR(TrailingMuon4.Eta(),TrailingMuon4.Phi(),
				      pfjet_eta[jet],pfjet_phi[jet]);
	   if (dR_trailing<0.4) continue;
	   if (jetPt>30.&&absJetEta<4.7) t_njets++;
	   if (jetPt>20.&&absJetEta<2.4) t_njetspt20++;
	 }
	 
	 t_metfilters = passedFilters(flags,metfilters);
	 t_badmufilter = passedFilters(flags,badmufilter);
	 t_met = met;
	 t_metphi = metphi;
	 
	 tuple->Fill();
       }

     }
     
   } // icand loop
   
   delete tree_;
   file_->Close();
   delete file_;
   
  }// filelist loop
  
  file->cd("");
  file->Write();
  file->Close();
  
  //delete file;
  }// int main loop 
