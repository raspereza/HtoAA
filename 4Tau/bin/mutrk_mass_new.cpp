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
#include "TRandom.h"
#include "TRandom3.h"
#include "HtoAA/Utilities/interface/json.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "HtoAA/Utilities/interface/PileUp.h"
#include "HtoAA/Utilities/interface/QCDModel.h"
#include "HtoAA/Utilities/interface/functions.h"

using namespace std;

int main(int argc, char * argv[]) {
  
  if (argc<2) {
    std::cout << "Usage of the program : mutrk_mass [config_file] [file_list]" << std::endl;
    exit(1);
  }


  // **** configuration
  Config cfg(argv[1]);

  const bool isData = cfg.get<bool>("IsData");
  const bool debug  = cfg.get<bool>("Debug");

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

  const string MuonHighPtTriggerFile = cfg.get<string>("MuonHighPtTriggerEff");
  const string MuonLowPtTriggerFile = cfg.get<string>("MuonLowPtTriggerEff");


  // ********** end of configuration *******************

  std::ifstream fileList(argv[2]);

  // event info
  ULong64_t event_nr;
  unsigned int event_run;
  unsigned int event_luminosityblock;

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
  bool muon_isPF[1000];
  bool muon_isGlobal[1000];
  bool muon_isTracker[1000];
  bool muon_isTight[1000];
  bool muon_isLoose[1000];
  bool muon_isMedium[1000];

  UInt_t genparticles_count;
  Float_t genparticles_e[1000];
  Float_t genparticles_px[1000];
  Float_t genparticles_py[1000];
  Float_t genparticles_pz[1000];
  Int_t genparticles_pdgid[1000];
  Int_t genparticles_status[1000];
  UInt_t genparticles_info[1000];

  UInt_t genjets_count;
  Float_t genjets_e[100];
  Float_t genjets_px[100];
  Float_t genjets_py[100];
  Float_t genjets_pz[100];
  Float_t genjets_pt[100];
  Float_t genjets_eta[100];
  Float_t genjets_phi[100];
  Int_t genjets_pdgid[100];
  Int_t genjets_status[100];

  float genweight;

  float metx;
  float mety;
  float met;
  float metphi;

  // Jets
  unsigned int pfjet_count;
  float pfjet_e[200];
  float pfjet_px[200];
  float pfjet_py[200];
  float pfjet_pz[200];
  float pfjet_pt[200];
  float pfjet_eta[200];
  float pfjet_phi[200];
  int pfjet_flavour[200];
  
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

  std::string rootFileName(argv[2]);
  
  std::string chainName("makeroottree/AC1B");
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
  TH1D * nGoodIsoMuonsH = new TH1D("nGoodIsoMuonsH","",11,-0.5,10.5);
  TH1D * nGoodLooseIsoMuonsH = new TH1D("nGoodLooseIsoMuonsH","",11,-0.5,10.5);
  TH1D * nGoodSbMuonsH = new TH1D("nGoodSbMuonsH","",11,-0.5,10.5);

  TH1D * puWeightH = new TH1D("puWeightH","",250,0,5);
  TH1D * triggerWeightH = new TH1D("triggerWeightH","",100,0,2);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,0.,2.);

  // histograms a selection

  TH1D * ptMuH = new TH1D("ptMuH","",400,0,400);
  TH1D * etaMuH = new TH1D("etaMuH","",48,-2.4,2.4);
  TH1D * dxyMuH = new TH1D("dxyMuH","",200,-0.5,0.5);
  TH1D * dzMuH = new TH1D("dzMuH","",200,-1,1);

  TH1D * nTracksMuH = new TH1D("nTracksMuH","",21,-0.5,20.5);
  TH1D * nSignalTracksMuH = new TH1D("nSignalTracksMuH","",21,-0.5,20.5);
  TH1D * nSoftTracksMuH = new TH1D("nSoftTracksMuH","",21,-0.5,20.5);
  
  // isolated muon-tracks
  TH1D * ptIsoTrackH = new TH1D("ptIsoTrackH","",100,0,100);
  TH1D * etaIsoTrackH = new TH1D("etaIsoTrackH","",48,-2.4,2.4);
  TH1D * dxyIsoTrackH = new TH1D("dxyIsoTrackH","",200,-0.5,0.5);
  TH1D * dzIsoTrackH = new TH1D("dzIsoTrackH","",200,-1,1);

  TH1D * ptIsoMuH = new TH1D("ptIsoMuH","",400,0,400);
  TH1D * etaIsoMuH = new TH1D("etaIsoMuH","",48,-2.4,2.4);
  TH1D * dxyIsoMuH = new TH1D("dxyIsoMuH","",200,-0.5,0.5);
  TH1D * dzIsoMuH = new TH1D("dzIsoMuH","",200,-1,1);

  TH1D * deltaRMuTrkIsoH = new TH1D("deltaRMuTrkIsoH","",100,0,1.0);
  
  // muon-tracks from Bkgd CR "LooseIso"
  TH1D * ptLooseIsoTrackH = new TH1D("ptLooseIsoTrackH","",100,0,100);
  TH1D * etaLooseIsoTrackH = new TH1D("etaLooseIsoTrackH","",48,-2.4,2.4);
  TH1D * dxyLooseIsoTrackH = new TH1D("dxyLooseIsoTrackH","",200,-0.5,0.5);
  TH1D * dzLooseIsoTrackH = new TH1D("dzLooseIsoTrackH","",200,-1,1);

  TH1D * ptLooseIsoMuH = new TH1D("ptLooseIsoMuH","",400,0,400);
  TH1D * etaLooseIsoMuH = new TH1D("etaLooseIsoMuH","",48,-2.4,2.4);
  TH1D * dxyLooseIsoMuH = new TH1D("dxyLooseIsoMuH","",200,-0.5,0.5);
  TH1D * dzLooseIsoMuH = new TH1D("dzLooseIsoMuH","",200,-1,1);

  TH1D * deltaRMuTrkLooseIsoH = new TH1D("deltaRMuTrkLooseIsoH","",100,0,1.0);

  // muon-tracks from Bkgd CR "LooseIso" and signal region "Iso"
  // adopted name Sb
  TH1D * ptSbTrackH = new TH1D("ptSbTrackH","",100,0,100);
  TH1D * etaSbTrackH = new TH1D("etaSbTrackH","",48,-2.4,2.4);
  TH1D * dxySbTrackH = new TH1D("dxySbTrackH","",200,-0.5,0.5);
  TH1D * dzSbTrackH = new TH1D("dzSbTrackH","",200,-1,1);

  TH1D * ptSbMuH = new TH1D("ptSbMuH","",400,0,400);
  TH1D * etaSbMuH = new TH1D("etaSbMuH","",48,-2.4,2.4);
  TH1D * dxySbMuH = new TH1D("dxySbMuH","",200,-0.5,0.5);
  TH1D * dzSbMuH = new TH1D("dzSbMuH","",200,-1,1);

  TH1D * deltaRMuTrkSbH = new TH1D("deltaRMuTrkSbH","",100,0,1.0);

  int nPartonMomBins = 3;
  float partonMomBins[4] = {0,50,100,150};

  int nMuonMomBins = 3;
  float muonMomBins[4] = {0,30,50,150};

  TString partonFlavor[4] = {"uds","g","c","b"};
  TString muonPartonNetCharge[2] = {"opposite","same"};
  TString partonMomRange[3] = {"Lt50","50to100","Gt100"};
  TString muonMomRange[3] = {"Lt30","30to50","Gt50"};
  TString RegionNames[3] = {"Iso","LooseIso","Sb"};
  TString MuTypes[2] = {"HighMu","LowMu"};

  TH1D * PartonMomBinsH = new TH1D("PartonMomBinsH","",nPartonMomBins,partonMomBins);
  TH1D * MuonPartonNetChargeH = new TH1D("MuonPartonNetChargeH","",2,0,2);
  TH1D * PartonFlavorH = new TH1D("PartonFlavorH","",4,0,4);
  for (int iB=0; iB<nPartonMomBins; ++iB)
    PartonMomBinsH->GetXaxis()->SetBinLabel(iB+1,partonMomRange[iB]);
  for (int iB=0; iB<4; ++iB)
    PartonFlavorH->GetXaxis()->SetBinLabel(iB+1,partonFlavor[iB]);
  for (int iB=0; iB<2; ++iB)
    MuonPartonNetChargeH->GetXaxis()->SetBinLabel(iB+1,muonPartonNetCharge[iB]);

  // Signal region
  TH1D * InvMassHighMuIsoH = new TH1D("InvMassHighMuIsoH","",200,0.,20.);
  TH1D * ModelInvMassHighMuIsoH = new TH1D("ModelInvMassHighMuIsoH","",200,0.,20.);
  TH1D * InvMassHighMuIsoFlavorChargeH[4][2];
  TH1D * ModelInvMassHighMuIsoFlavorChargeH[4][2];
  TH1D * PartonMomHighMuIsoFlavorChargeH[4][2];
  TH1D * InvMassHighMuIsoFlavorChargeMomH[4][2][10];

  TH1D * InvMassLowMuIsoH = new TH1D("InvMassLowMuIsoH","",200,0.,20.);
  TH1D * ModelInvMassLowMuIsoH = new TH1D("ModelInvMassLowMuIsoH","",200,0.,20.);
  TH1D * InvMassLowMuIsoFlavorChargeH[4][2];
  TH1D * ModelInvMassLowMuIsoFlavorChargeH[4][2];
  TH1D * PartonMomLowMuIsoFlavorChargeH[4][2];
  TH1D * InvMassLowMuIsoFlavorChargeMomH[4][2][10];

  TH1D * InvMassIsoH = new TH1D("InvMassIsoH","",200,0.,20.);
  TH1D * ModelInvMassIsoH = new TH1D("ModelInvMassIsoH","",200,0.,20.);

  TH1D * InvMassDimuonIsoH = new TH1D("InvMassDimuonIsoH","",200,0.,20.);
  TH1D * ModelInvMassDimuonIsoH = new TH1D("ModelInvMassDimuonIsoH","",200,0.,20.);

  // LooseIso region
  TH1D * InvMassHighMuLooseIsoH = new TH1D("InvMassHighMuLooseIsoH","",200,0.,20.);
  TH1D * ModelInvMassHighMuLooseIsoH = new TH1D("ModelInvMassHighMuLooseIsoH","",200,0.,20.);
  TH1D * InvMassHighMuLooseIsoAllH = new TH1D("InvMassHighMuLooseIsoAllH","",200,0.,20.);
  TH1D * InvMassHighMuLooseIsoFlavorChargeH[4][2];
  TH1D * ModelInvMassHighMuLooseIsoFlavorChargeH[4][2];
  TH1D * PartonMomHighMuLooseIsoFlavorChargeH[4][2];
  TH1D * InvMassHighMuLooseIsoFlavorChargeMomH[4][2][10];

  TH1D * InvMassLowMuLooseIsoH = new TH1D("InvMassLowMuLooseIsoH","",200,0.,20.);
  TH1D * ModelInvMassLowMuLooseIsoH = new TH1D("ModelInvMassLowMuLooseIsoH","",200,0.,20.);
  TH1D * InvMassLowMuLooseIsoAllH = new TH1D("InvMassLowMuLooseIsoAllH","",200,0.,20.);
  TH1D * InvMassLowMuLooseIsoFlavorChargeH[4][2];
  TH1D * ModelInvMassLowMuLooseIsoFlavorChargeH[4][2];
  TH1D * PartonMomLowMuLooseIsoFlavorChargeH[4][2];
  TH1D * InvMassLowMuLooseIsoFlavorChargeMomH[4][2][10];

  TH1D * InvMassLooseIsoH = new TH1D("InvMassLooseIsoH","",200,0.,20.);
  TH1D * ModelInvMassLooseIsoH = new TH1D("ModelInvMassLooseIsoH","",200,0.,20.);

  TH1D * InvMassDimuonLooseIsoH = new TH1D("InvMassDimuonLooseIsoH","",200,0.,20.);
  TH1D * ModelInvMassDimuonLooseIsoH = new TH1D("ModelInvMassDimuonLooseIsoH","",200,0.,20.);

  // Sb region (Signal + Background region)
  TH1D * InvMassHighMuSbH = new TH1D("InvMassHighMuSbH","",200,0.,20.);
  TH1D * ModelInvMassHighMuSbH = new TH1D("ModelInvMassHighMuSbH","",200,0.,20.);
  TH1D * InvMassHighMuSbAllH = new TH1D("InvMassHighMuSbAllH","",200,0.,20.);
  TH1D * InvMassHighMuSbFlavorChargeH[4][2];
  TH1D * ModelInvMassHighMuSbFlavorChargeH[4][2];
  TH1D * PartonMomHighMuSbFlavorChargeH[4][2];
  TH1D * InvMassHighMuSbFlavorChargeMomH[4][2][10];

  TH1D * InvMassLowMuSbH = new TH1D("InvMassLowMuSbH","",200,0.,20.);
  TH1D * ModelInvMassLowMuSbH = new TH1D("ModelInvMassLowMuSbH","",200,0.,20.);
  TH1D * InvMassLowMuSbAllH = new TH1D("InvMassLowMuSbAllH","",200,0.,20.);
  TH1D * InvMassLowMuSbFlavorChargeH[4][2];
  TH1D * ModelInvMassLowMuSbFlavorChargeH[4][2];
  TH1D * PartonMomLowMuSbFlavorChargeH[4][2];
  TH1D * InvMassLowMuSbFlavorChargeMomH[4][2][10];

  // momentum
  TH1D * PartonMomFlavorH[4];
  TH1D * PartonMomHighMuFlavorChargeH[4][2];
  TH1D * PartonMomLowMuFlavorChargeH[4][2];

  // Unmatched muons
  TH1D * MuonMomHighMuUnmatchedH = new TH1D("MuonMomHighMuUnmatchedH","",500,0.,500.);
  TH1D * MuonMomLowMuUnmatchedH = new TH1D("MuonMomLowMuUnmatchedH","",500,0.,500.);

  TH1D * MuonMomHighMuIsoUnmatchedH = new TH1D("MuonMomHighMuIsoUnmatchedH","",500,0.,500.);
  TH1D * MuonMomHighMuLooseIsoUnmatchedH = new TH1D("MuonMomHighMuLooseIsoUnmatchedH","",500,0.,500.);
  TH1D * MuonMomHighMuSbUnmatchedH = new TH1D("MuonMomHighMuSbUnmatchedH","",500,0.,500.);
  TH1D * InvMassHighMuIsoMomUnmatchedH[10];
  TH1D * InvMassHighMuLooseIsoMomUnmatchedH[10];
  TH1D * InvMassHighMuSbMomUnmatchedH[10];

  TH1D * MuonMomLowMuIsoUnmatchedH = new TH1D("MuonMomLowMuIsoUnmatchedH","",500,0.,500.);
  TH1D * MuonMomLowMuLooseIsoUnmatchedH = new TH1D("MuonMomLowMuLooseIsoUnmatchedH","",500,0.,500.);
  TH1D * MuonMomLowMuSbUnmatchedH = new TH1D("MuonMomLowMuSbUnmatchedH","",500,0.,500.);
  TH1D * InvMassLowMuIsoMomUnmatchedH[10];
  TH1D * InvMassLowMuLooseIsoMomUnmatchedH[10];
  TH1D * InvMassLowMuSbMomUnmatchedH[10];  


  // other distributions 
  TH1D * PartonMultiplicityH = new TH1D("PartonMultiplicityH","",20,-0.5,19.5);
  TH1D * PFJetMultiplicityH  = new TH1D("PFJetMultiplicityH","",20,-0.5,19.5);
  TH1D * GenJetMultiplicityH = new TH1D("GenJetMultiplicityH","",20,-0.5,19.5);

  TH1D * PartonMultiplicityMuH = new TH1D("PartonMultiplicityMuH","",20,-0.5,19.5);
  TH1D * PFJetMultiplicityMuH  = new TH1D("PFJetMultiplicityMuH","",20,-0.5,19.5);
  TH1D * GenJetMultiplicityMuH = new TH1D("GenJetMultiplicityMuH","",20,-0.5,19.5);
  TH1D * deltaRPartonMuH = new TH1D("deltaRPartonMuH","",100,0.,1.);

  TH1D * PartonMultiplicityMuIsoH = new TH1D("PartonMultiplicityMuIsoH","",20,-0.5,19.5);
  TH1D * PFJetMultiplicityMuIsoH  = new TH1D("PFJetMultiplicityMuIsoH","",20,-0.5,19.5);
  TH1D * GenJetMultiplicityMuIsoH = new TH1D("GenJetMultiplicityMuIsoH","",20,-0.5,19.5);
  TH1D * deltaRPartonMuIsoH = new TH1D("deltaRPartonMuIsoH","",100,0.,1.);

  TH1D * PartonMultiplicityMuLooseIsoH = new TH1D("PartonMultiplicityMuLooseIsoH","",20,-0.5,19.5);
  TH1D * PFJetMultiplicityMuLooseIsoH  = new TH1D("PFJetMultiplicityMuLooseIsoH","",20,-0.5,19.5);
  TH1D * GenJetMultiplicityMuLooseIsoH = new TH1D("GenJetMultiplicityMuLooseIsoH","",20,-0.5,19.5);
  TH1D * deltaRPartonMuLooseIsoH = new TH1D("deltaRPartonMuLooseIsoH","",100,0.,1.);

  TH1D * PartonMultiplicityMuSbH = new TH1D("PartonMultiplicityMuSbH","",20,-0.5,19.5);
  TH1D * PFJetMultiplicityMuSbH  = new TH1D("PFJetMultiplicityMuSbH","",20,-0.5,19.5);
  TH1D * GenJetMultiplicityMuSbH = new TH1D("GenJetMultiplicityMuSbH","",20,-0.5,19.5);
  TH1D * deltaRPartonMuSbH = new TH1D("deltaRPartonMuSbH","",100,0.,1.);
  
  TH1D * InvMassH = new TH1D("InvMassH","",200,0.,20.);
  TH2D * InvMass2DH = new TH2D("InvMass2DH","",200,0.,20.,200,0.,20.);

  // Correlation Plots
  TH1D * InvMass_ControlXH = new TH1D("InvMass_ControlXH","",200,0.,20.); 
  TH2D * InvMass2D_ControlXH = new TH2D("InvMass2D_ControlXH","",200,0.,20.,200,0.,20.);
   
  TH1D * InvMass_ControlYH = new TH1D("InvMass_ControlYH","",200,0.,20.); 
  TH2D * InvMass2D_ControlYH = new TH2D("InvMass2D_ControlYH","",200,0.,20.,200,0.,20.);

  for (int iM=0; iM<nMuonMomBins; ++iM) {
    InvMassLowMuIsoMomUnmatchedH[iM] = 
      new TH1D("InvMassLowMuIso_"+muonMomRange[iM]+"_UnmatchedH","",200,0.,20.);
    InvMassLowMuLooseIsoMomUnmatchedH[iM] = 
      new TH1D("InvMassLowMuLooseIso_"+muonMomRange[iM]+"_UnmatchedH","",200,0.,20.);
    InvMassLowMuSbMomUnmatchedH[iM] = 
      new TH1D("InvMassLowMuSb_"+muonMomRange[iM]+"_UnmatchedH","",200,0.,20.);
    InvMassHighMuIsoMomUnmatchedH[iM] = 
      new TH1D("InvMassHighMuIso_"+muonMomRange[iM]+"_UnmatchedH","",200,0.,20.);
    InvMassHighMuLooseIsoMomUnmatchedH[iM] = 
      new TH1D("InvMassHighMuLooseIso_"+muonMomRange[iM]+"_UnmatchedH","",200,0.,20.);
    InvMassHighMuSbMomUnmatchedH[iM] = 
      new TH1D("InvMassHighMuSb_"+muonMomRange[iM]+"_UnmatchedH","",200,0.,20.);
  }

  for (int iF=0; iF<4; ++iF) {

    PartonMomFlavorH[iF] = new TH1D("partonMomFlavor_"+partonFlavor[iF],"",500,0.,5000.);
    

    for (int iQ=0; iQ<2; ++iQ) {

      PartonMomHighMuFlavorChargeH[iF][iQ] = new TH1D("partonMomHighMu_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      PartonMomLowMuFlavorChargeH[iF][iQ] = new TH1D("partonMomLowMu_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      
      InvMassHighMuIsoFlavorChargeH[iF][iQ] = new TH1D("InvMassHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      InvMassLowMuIsoFlavorChargeH[iF][iQ] = new TH1D("InvMassLowMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      ModelInvMassHighMuIsoFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      ModelInvMassLowMuIsoFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassLowMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      PartonMomHighMuIsoFlavorChargeH[iF][iQ] = new TH1D("PartonMomHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      PartonMomLowMuIsoFlavorChargeH[iF][iQ] = new TH1D("PartonMomLowMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);

      InvMassHighMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("InvMassHighMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      InvMassLowMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("InvMassLowMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      ModelInvMassHighMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassHighMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      ModelInvMassLowMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassLowMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      PartonMomHighMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("PartonMomHighMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      PartonMomLowMuLooseIsoFlavorChargeH[iF][iQ] = new TH1D("PartonMomLowMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);

      InvMassHighMuSbFlavorChargeH[iF][iQ] = new TH1D("InvMassHighMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      InvMassLowMuSbFlavorChargeH[iF][iQ] = new TH1D("InvMassLowMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      ModelInvMassHighMuSbFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassHighMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      ModelInvMassLowMuSbFlavorChargeH[iF][iQ] = new TH1D("ModelInvMassLowMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",200,0.,20.);
      PartonMomHighMuSbFlavorChargeH[iF][iQ] = new TH1D("PartonMomHighMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);
      PartonMomLowMuSbFlavorChargeH[iF][iQ] = new TH1D("PartonMomLowMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ],"",500,0.,5000.);

      for (int iM=0; iM<nPartonMomBins; ++iM) {

	InvMassLowMuIsoFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassLowMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",200,0.,20.);
	InvMassLowMuLooseIsoFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassLowMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",200,0.,20.);
	InvMassLowMuSbFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassLowMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",200,0.,20.);

	InvMassHighMuIsoFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassHighMuIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",200,0.,20.);
	InvMassHighMuLooseIsoFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassHighMuLooseIso_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",200,0.,20.);
	InvMassHighMuSbFlavorChargeMomH[iF][iQ][iM] = new TH1D("InvMassHighMuSb_"+partonFlavor[iF]+"_"+muonPartonNetCharge[iQ]+"_"+partonMomRange[iM],"",200,0.,20.);

      }      
      
    }
  }
  // for muons selected with same-sign requirement
  TH1D * partonMuSS[2][4][2][10];
  TH1D * partonMuSSPass[2][4][2][10]; // 3 bins : SR, SB, SR+SB
  // unmatched
  TH1D * unmatchedMuSS[2][10];
  TH1D * unmatchedMuSSPass[2][10];
 
  for (unsigned int mu=0; mu<2; ++mu) {
    for (int iMom=0; iMom<nMuonMomBins; ++iMom) {
      TString histName = "unmatchedMuSS_" + MuTypes[mu] + muonMomRange[iMom];
      unmatchedMuSS[mu][iMom] = new TH1D(histName,"",1,0.,1.);
      unmatchedMuSSPass[mu][iMom] = new TH1D(histName+"_passed","",3,-1.,2.);
    }
    for (int iMom=0; iMom<nPartonMomBins; ++iMom) {
      for (unsigned int iF=0; iF<4; ++iF) {
	for (unsigned int iQ=0; iQ<2; ++iQ) {
	  TString histName = "partonMuSS_" + MuTypes[mu] + partonFlavor[iF] + muonPartonNetCharge[iQ] + muonMomRange[iMom];
	  
	  partonMuSS[mu][iF][iQ][iMom] = new TH1D(histName,"",1,0.,1.);
	  partonMuSSPass[mu][iF][iQ][iMom] = new TH1D(histName+"_passed","",3,-1.,2.);
	}
      }
    }
  }

  string cmsswBase = (getenv ("CMSSW_BASE"));

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

  std::string qcdFileName = cfg.get<string>("qcdModelFileName");
  TString QCDFileName(qcdFileName);

  // QCD Model
  TString fileNameQCDModel = TString(cmsswBase)+TString("/src/HtoAA/data/"+QCDFileName);
  QCDModel * qcdModel = new QCDModel(fileNameQCDModel);

  // PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/HtoAA/data/PileUpDistrib/"+PileUpDataFile,"read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/HtoAA/data/PileUpDistrib/"+PileUpMCFile, "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  // Trigger efficiencies

  ScaleFactor * SF_muonHighPt = new ScaleFactor();
  SF_muonHighPt->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonHighPtTriggerFile));
  ScaleFactor * SF_muonLowPt = new ScaleFactor();
  SF_muonLowPt->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonLowPtTriggerFile));

  TString filen;
  int iFiles = 0;
  int events = 0;
  while (fileList >> filen) {
   iFiles++;
   cout << "file " << iFiles << " : " << filen << endl;
   
   TFile * file_ = TFile::Open(TString(filen));
   if (file_==NULL) continue;

   TTree * tree_ = (TTree*)file_->Get(TString(chainName));
   
   if (tree_==NULL) continue;

   tree_->SetMaxVirtualSize(3000000);
   // event info
   tree_->SetBranchAddress("event_nr", &event_nr);
   tree_->SetBranchAddress("event_run", &event_run);
   tree_->SetBranchAddress("event_luminosityblock", &event_luminosityblock);

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
   tree_->SetBranchAddress("muon_chargedHadIso", muon_chargedHadIso);
   tree_->SetBranchAddress("muon_neutralHadIso", muon_neutralHadIso);
   tree_->SetBranchAddress("muon_photonIso", muon_photonIso);
   tree_->SetBranchAddress("muon_puIso", muon_puIso);
   tree_->SetBranchAddress("muon_isMedium", muon_isMedium);

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

   // jets
   tree_->SetBranchAddress("pfjet_count", &pfjet_count);
   tree_->SetBranchAddress("pfjet_e", pfjet_e);
   tree_->SetBranchAddress("pfjet_px", pfjet_px);
   tree_->SetBranchAddress("pfjet_py", pfjet_py);
   tree_->SetBranchAddress("pfjet_pz", pfjet_pz);
   tree_->SetBranchAddress("pfjet_pt", pfjet_pt);
   tree_->SetBranchAddress("pfjet_eta", pfjet_eta);
   tree_->SetBranchAddress("pfjet_phi", pfjet_phi);
   tree_->SetBranchAddress("pfjet_flavour", pfjet_flavour);

   // genjets
   tree_->SetBranchAddress("genjets_count",&genjets_count);
   tree_->SetBranchAddress("genjets_e",genjets_e);
   tree_->SetBranchAddress("genjets_px",genjets_px);
   tree_->SetBranchAddress("genjets_py",genjets_py);
   tree_->SetBranchAddress("genjets_pz",genjets_pz);
   tree_->SetBranchAddress("genjets_pt",genjets_pt);
   tree_->SetBranchAddress("genjets_eta",genjets_eta);
   tree_->SetBranchAddress("genjets_phi",genjets_phi);
   tree_->SetBranchAddress("genjets_pdgid",genjets_pdgid);
   tree_->SetBranchAddress("genjets_status",genjets_status);

   // Additional trigger objects
   tree_->SetBranchAddress("run_hltfilters",&hltfilters);
   //   tree_->SetBranchAddress("run_btagdiscriminators", &run_btagdiscriminators);
   tree_->SetBranchAddress("hltriggerresults",&hltriggerresults);
   tree_->SetBranchAddress("hltriggerprescales",&hltriggerprescales);

   tree_->SetBranchAddress("numtruepileupinteractions",&numtruepileupinteractions);

   if (!isData) {
     tree_->SetBranchAddress("genweight",&genweight);
     tree_->SetBranchAddress("genparticles_count", &genparticles_count);
     tree_->SetBranchAddress("genparticles_e", genparticles_e);
     tree_->SetBranchAddress("genparticles_px", genparticles_px);
     tree_->SetBranchAddress("genparticles_py", genparticles_py);
     tree_->SetBranchAddress("genparticles_pz", genparticles_pz);
     tree_->SetBranchAddress("genparticles_pdgid", genparticles_pdgid);
     tree_->SetBranchAddress("genparticles_status", genparticles_status);
     tree_->SetBranchAddress("genparticles_info", genparticles_info);
   }   

   int numberOfCandidates = tree_->GetEntries();

   std::cout << "number of events = " << numberOfCandidates << std::endl;
   
   TRandom3 rand;

   for (int iCand=0; iCand<numberOfCandidates; iCand++) {
     
     tree_->GetEntry(iCand);

     events++;
     if (events%10000==0) cout << "   processed events : " << events << endl;

     float weight = 1;
     if (!isData) {
       weight *= genweight;
     }


     histWeightsH->Fill(1.0,weight);

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

     float puweight = 1;
     if (!isData) {
       puweight = float(PUofficial->get_PUweight(double(numtruepileupinteractions)));
       //       std::cout << "n(true interactions) = " << numtruepileupinteractions << "   :  PU weight = " << puweight << std::endl; 
     }
     puWeightH->Fill(puweight,1.0);
     weight *= puweight;

     if (debug) {
       std::cout << std::endl;
       std::cout << "+++++++++++++++++++++++++++++" << std::endl;
       std::cout << std::endl;
       std::cout << "genjets : " << genjets_count << std::endl;
     }
     int partons = 0;
     int pfjets  = 0;
     int genjets = 0;
     if (!isData) {
       for (unsigned int igen=0; igen<genjets_count; ++igen) {
	 int pdgId = genjets_pdgid[igen];
	 TLorentzVector genjetLV; genjetLV.SetXYZT(genjets_px[igen],
						   genjets_py[igen],
						   genjets_pz[igen],
						   genjets_e[igen]);
       
	 if (genjetLV.Pt()>10&&fabs(genjetLV.Eta())<3.0) {
	   genjets++;
	   if (debug)
	     printf("  flavor = %3i   pT = %7.2f   eta = %5.2f   phi = %5.2f   status = %3i\n",
		    genjets_pdgid[igen],genjetLV.Pt(),genjetLV.Eta(),genjetLV.Phi(),genjets_status[igen]);
	 }
       }
     }
     if (debug) {
       std::cout << std::endl;
       std::cout << "jets -> " << std::endl;
     }
     std::vector<int> partonPdgId; partonPdgId.clear();
     std::vector<TLorentzVector> partonLV; partonLV.clear();
     for (unsigned ijet=0; ijet<pfjet_count; ++ijet) {
       if (pfjet_pt[ijet]>10&&fabs(pfjet_eta[ijet])<3.0) {
	 pfjets++;
	 //	 if (pfjet_flavour[ijet]==0)
	 //	   cout << "++++++++++++++ " << pfjet_flavour[ijet] << endl;
	 if (pfjet_flavour[ijet]!=0&&!isData) { 
	   partons++;
	   int absPdgId = TMath::Abs(pfjet_flavour[ijet]);
	   int iflav = 0;
	   if (absPdgId==21) iflav = 1;
	   if (absPdgId==4)  iflav = 2;
	   if (absPdgId==5)  iflav = 3;
	   TLorentzVector jetLV; jetLV.SetXYZT(pfjet_px[ijet],
					       pfjet_py[ijet],
					       pfjet_pz[ijet],
					       pfjet_e[ijet]);
	   TLorentzVector partLV = jetLV;
	   float dRMin = 0.5;
	   for (unsigned int igen=0; igen<genjets_count; ++igen) {
	     TLorentzVector genjetLV; genjetLV.SetXYZT(genjets_px[igen],
						       genjets_py[igen],
						       genjets_pz[igen],
						       genjets_e[igen]);
	     float dRJets = deltaR(jetLV.Eta(),jetLV.Phi(),
				   genjetLV.Eta(),genjetLV.Phi());
	     if (dRJets<dRMin) {
	       dRMin = dRJets;
	       partLV = genjetLV;
	     }	     
	   }
	   partonPdgId.push_back(pfjet_flavour[ijet]);
	   partonLV.push_back(partLV);
	   PartonMomFlavorH[iflav]->Fill(partLV.P(),weight);
	 }
	 if (debug)
	   printf("  flavor = %3i   pT = %7.2f   eta = %5.2f   phi = %5.2f\n",
		  pfjet_flavour[ijet],pfjet_pt[ijet],pfjet_eta[ijet],pfjet_phi[ijet]);
       
       }
     }
    
     // filling histograms 
     PartonMultiplicityH->Fill(float(partons),weight);
     PFJetMultiplicityH->Fill(float(pfjets),weight);
     GenJetMultiplicityH->Fill(float(genjets),weight);

     // ********************
     // selecting good muons
     // ********************
     vector<unsigned int> muons; muons.clear();
     vector<int> muons_flavour; muons_flavour.clear();
     vector<int> muons_pdgid; muons_pdgid.clear();
     vector<int> muons_mom; muons_mom.clear();
     vector<int> muons_net; muons_net.clear();
     vector<int> muons_type; muons_type.clear();
     vector<int> muons_region; muons_region.clear();
     vector<float> muons_mutrkmass; muons_mutrkmass.clear();

     for(UInt_t i=0;i<muon_count;i++){
       bool muonID = muon_isMedium[i]; // MC 
       if (!muonID) continue;
       if(fabs(muon_dxy[i])>dxyMuonCut) continue;
       if(fabs(muon_dz[i])>dzMuonCut) continue;
       if(muon_pt[i]<ptMuonLowCut) continue;
       if(fabs(muon_eta[i])>etaMuonLowCut) continue;
       //      cout << "muon pt = " << muon_pt[i] << endl;
       muons.push_back(i);
     }
     
     nGoodMuonsH->Fill(float(muons.size()),weight);
      
     if (muons.size()<1) continue; // quit event if number of good muons < 1

     int nIsoMuons = 0;
     int nLooseIsoMuons = 0;
     int nSbMuons = 0;

     if (debug) {
       std::cout << std::endl;
       std::cout << "muons -> " << std::endl;
     }

     for (unsigned int imu=0; imu<muons.size(); ++imu) {

       unsigned int index = muons.at(imu);

       if (debug) {
	 std::cout << " muon index : " << index << std::endl;
       }

       bool muHighPassed = muon_pt[index]>ptMuonHighCut && 
	 fabs(muon_eta[index])<etaMuonHighCut;
       bool muLowPassed  = !muHighPassed;

       // Muon
       TLorentzVector Muon4; Muon4.SetXYZM(muon_px[index],
					   muon_py[index],
					   muon_pz[index],
					   MuMass);
       // determine flavour of jet
       float dRmin = 0.5;
       int flavour = -1;
       float qnet  = 0;
       int pdgId = 0;
       bool matchedParton = false;
       TLorentzVector matchedPartonLV;
       for (unsigned int ip=0; ip<partonPdgId.size(); ++ip) {
	 TLorentzVector partLV = partonLV.at(ip);
	 float drJetMuon = deltaR(muon_eta[index],muon_phi[index],
				  partLV.Eta(),partLV.Phi());
	 if (drJetMuon<dRmin) {
	   dRmin = drJetMuon;
	   int absFlav = TMath::Abs(partonPdgId.at(ip));
	   pdgId = partonPdgId.at(ip);
	   flavour = 0;
	   if (absFlav==21) 
	     flavour = 1;
	   if (absFlav==4)
	     flavour = 2;
	   if (absFlav==5)
	     flavour = 3;
	   qnet = float(muon_charge[index])*float(pdgId);
	   matchedParton = true;
	   matchedPartonLV = partLV;
	 }
       }

       if (debug) {
	 std::cout << " jet flavor : " << flavour << std::endl;
       }

       int net = 0;
       if (qnet>0.0) net = 1;
       if (flavour==1) net = 0;       

       // counting tracks around each muon
       std::vector<unsigned int> trkMu; trkMu.clear(); // all tracks
       std::vector<unsigned int> trkSignalMu; trkSignalMu.clear(); // signal tracks
       std::vector<unsigned int> trkSoftMu; trkSoftMu.clear(); // tracks control region

       TLorentzVector trackLV; trackLV.SetXYZM(0,0,0,PionMass);
     
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
	 
	 TLorentzVector MuDiff = Muon4 - trk4;
	 if (MuDiff.P()>0.1) { // track is not leading muon
	   float drTrkMu = deltaR(muon_eta[index],muon_phi[index],
				  track_eta[iTrk],   track_phi[iTrk]);
	   float qTrkMu = track_charge[iTrk]*muon_charge[index];
	   if (drTrkMu<dRIsoMuon){
	     trkMu.push_back(iTrk);
	     if (track_pt[iTrk]>ptTrkLooseCut && track_pt[iTrk]< ptTrkCut)
	       trkSoftMu.push_back(iTrk);
	   }
	   if (drTrkMu<dRIsoMuon && qTrkMu<0 
	       && fabs(track_dxy[iTrk])<dxyTrkCut 
	       && fabs(track_dz[iTrk])<dzTrkCut 
	       && track_pt[iTrk]>ptTrkCut) {
	     trkSignalMu.push_back(iTrk);
	     trackLV = trk4;
	   }
	 }
       }

       

       TLorentzVector muonTrkLV = trackLV + Muon4;
       float muonTrkMass = muonTrkLV.M();
       float deltaRMuonTrk = 1.1;
       unsigned int indexTrk = 0;
       if (trkSignalMu.size()==1) { 
	 indexTrk = trkSignalMu.at(0);
	 deltaRMuonTrk = deltaR(Muon4.Eta(),Muon4.Phi(),
				trackLV.Eta(),trackLV.Phi());
       }

       if (debug) {
	 printf("  pT = %6.2f   eta = %5.2f   phi = %5.2f   jetFlavor = %3i  nTrk = %2i  nSigTrk = %2i\n",
		Muon4.Pt(),Muon4.Eta(),Muon4.Phi(),pdgId,int(trkMu.size()),int(trkSignalMu.size()));
       }

       bool sigMu = trkSignalMu.size()==1 && trkMu.size()==1;

       bool bkgdMu = 
	 (trkSignalMu.size()==1 && trkSoftMu.size()==1 && trkMu.size()==2) ||
	 (trkSignalMu.size()==1 && trkSoftMu.size()==2 && trkMu.size()==3) ||
	 (trkSignalMu.size()==1 && trkSoftMu.size()==3 && trkMu.size()==4);

       bool sbMu = sigMu || bkgdMu;
       
       PartonMultiplicityMuH->Fill(float(partons),weight);
       PFJetMultiplicityMuH->Fill(float(pfjets),weight);
       GenJetMultiplicityMuH->Fill(float(genjets),weight);

       ptMuH->Fill(muon_pt[index],weight);
       etaMuH->Fill(muon_eta[index],weight);
       dxyMuH->Fill(muon_dxy[index],weight);
       dzMuH->Fill(muon_dz[index],weight);
       
       nTracksMuH->Fill(float(trkMu.size()),weight);
       nSoftTracksMuH->Fill(float(trkSoftMu.size()),weight);
       nSignalTracksMuH->Fill(float(trkSoftMu.size()),weight);

       int momBin = 0;
       int partonMomBin = 0;
       int muonMomBin = 0;
       if (matchedParton) {
	 partonMomBin = binNumber(TMath::Min(float(matchedPartonLV.P()),float(partonMomBins[nPartonMomBins]-0.1)),nPartonMomBins,partonMomBins);
	 momBin = partonMomBin;
	 deltaRPartonMuH->Fill(dRmin,weight);
	 if (muHighPassed) PartonMomHighMuFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weight);
	 if (muLowPassed) PartonMomLowMuFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weight);
	 if (muHighPassed) {
	   for (int iM=0; iM<20; ++iM) {
	     double mass = double(iM) + double(0.5);
	     int muType = 0;
	     int ireg = 0;
	     double pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
	     ModelInvMassHighMuIsoH->Fill(mass,weight*pdf);
	     ModelInvMassHighMuIsoFlavorChargeH[flavour][net]->Fill(mass,weight*pdf);
	     ireg = 1;
             pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
             ModelInvMassHighMuLooseIsoH->Fill(mass,weight*pdf);
             ModelInvMassHighMuLooseIsoFlavorChargeH[flavour][net]->Fill(mass,weight*pdf);
	     ireg = 2;
             pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
             ModelInvMassHighMuSbH->Fill(mass,weight*pdf);
             ModelInvMassHighMuSbFlavorChargeH[flavour][net]->Fill(mass,weight*pdf);
	   }
	 }
	 if (muLowPassed) {
	   for (int iM=0; iM<20; ++iM) {
	     double mass = double(iM) + double(0.5);
	     int muType = 1;
	     int ireg = 0;
	     double pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
	     ModelInvMassLowMuIsoH->Fill(mass,weight*pdf);
	     ModelInvMassLowMuIsoFlavorChargeH[flavour][net]->Fill(mass,weight*pdf);
	     ireg = 1;
             pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
             ModelInvMassLowMuLooseIsoH->Fill(mass,weight*pdf);
             ModelInvMassLowMuLooseIsoFlavorChargeH[flavour][net]->Fill(mass,weight*pdf);
	     ireg = 2;
             pdf = qcdModel->getMuMassPdf(partonMomBin,muType,ireg,flavour,net,mass);
             ModelInvMassLowMuSbH->Fill(mass,weight*pdf);
             ModelInvMassLowMuSbFlavorChargeH[flavour][net]->Fill(mass,weight*pdf);
	   }
	 }
       }
       else {
	 muonMomBin = binNumber(TMath::Min(float(Muon4.P()),float(muonMomBins[nMuonMomBins]-0.1)),nMuonMomBins,muonMomBins);
	 momBin=muonMomBin;
	 if (muHighPassed) MuonMomHighMuUnmatchedH->Fill(Muon4.P(),weight);
	 if (muLowPassed) MuonMomLowMuUnmatchedH->Fill(Muon4.P(),weight);
       }

       // 
       // save muon info
       //
       int muonType = 0;
       if (muLowPassed) muonType = 1;
       muons_mom.push_back(momBin);
       muons_flavour.push_back(flavour);
       muons_net.push_back(net);
       muons_type.push_back(muonType);
       int muon_region = -1;
       if (sigMu) muon_region = 0;
       if (bkgdMu) muon_region = 1;
       muons_region.push_back(muon_region);
       muons_mutrkmass.push_back(muonTrkMass);
       muons_pdgid.push_back(pdgId);

       //       cout << "flav = " << flavour << "  momBin = " << momBin << "  muType = " << muonType << "  net = " << net << std::endl;

       
       //       cout << "after saving" << endl;


       if (sigMu) {
	 
	 nIsoMuons++;
	 ptIsoTrackH->Fill(trackLV.Pt(),weight);
	 etaIsoTrackH->Fill(trackLV.Eta(),weight);
	 dxyIsoTrackH->Fill(track_dxy[indexTrk],weight);
	 dzIsoTrackH->Fill(track_dz[indexTrk],weight);
	 ptIsoMuH->Fill(Muon4.Pt(),weight);
	 etaIsoMuH->Fill(Muon4.Eta(),weight);
	 dxyIsoMuH->Fill(muon_dxy[index],weight);
	 dzIsoMuH->Fill(muon_dz[index],weight);
	 PartonMultiplicityMuIsoH->Fill(float(partons),weight);
	 PFJetMultiplicityMuIsoH->Fill(float(pfjets),weight);
	 GenJetMultiplicityMuIsoH->Fill(float(genjets),weight);
	 deltaRMuTrkIsoH->Fill(deltaRMuonTrk,weight);


	 if (matchedParton) {
	   deltaRPartonMuIsoH->Fill(dRmin,weight);
	   if (muHighPassed) {
	     InvMassHighMuIsoFlavorChargeH[flavour][net]->Fill(muonTrkMass,weight);
	     PartonMomHighMuIsoFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weight);
	     InvMassHighMuIsoFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weight);
	   }
	   if (muLowPassed) { 
	     InvMassLowMuIsoFlavorChargeH[flavour][net]->Fill(muonTrkMass,weight);
             PartonMomLowMuIsoFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weight);
             InvMassLowMuIsoFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weight);
	   }
	 }
	 else {
	   if (muHighPassed) {
	     MuonMomHighMuIsoUnmatchedH->Fill(Muon4.P(),weight);
	     InvMassHighMuIsoMomUnmatchedH[muonMomBin]->Fill(muonTrkMass,weight);
	   }
	   if (muLowPassed) {
	     MuonMomLowMuIsoUnmatchedH->Fill(Muon4.P(),weight);
	     InvMassLowMuIsoMomUnmatchedH[muonMomBin]->Fill(muonTrkMass,weight);
	   }

	 }
	 
       }

       //       std::cout << "Ok3  " << muonMomBin << "  " << matchedParton << std::endl;

       if (bkgdMu) {

	 nLooseIsoMuons++;
	 ptLooseIsoTrackH->Fill(trackLV.Pt(),weight);
         etaLooseIsoTrackH->Fill(trackLV.Eta(),weight);
         dxyLooseIsoTrackH->Fill(track_dxy[indexTrk],weight);
         dzLooseIsoTrackH->Fill(track_dz[indexTrk],weight);
         ptLooseIsoMuH->Fill(Muon4.Pt(),weight);
         etaLooseIsoMuH->Fill(Muon4.Eta(),weight);
         dxyLooseIsoMuH->Fill(muon_dxy[index],weight);
         dzLooseIsoMuH->Fill(muon_dz[index],weight);
         PartonMultiplicityMuLooseIsoH->Fill(float(partons),weight);
         PFJetMultiplicityMuLooseIsoH->Fill(float(pfjets),weight);
         GenJetMultiplicityMuLooseIsoH->Fill(float(genjets),weight);
	 deltaRMuTrkLooseIsoH->Fill(deltaRMuonTrk,weight);


	 if (matchedParton) {
           deltaRPartonMuLooseIsoH->Fill(dRmin,weight);
	   if (muHighPassed) {
             InvMassHighMuLooseIsoFlavorChargeH[flavour][net]->Fill(muonTrkMass,weight);
             PartonMomHighMuLooseIsoFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weight);
             InvMassHighMuLooseIsoFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weight);
           }
           if (muLowPassed) {
             InvMassLowMuLooseIsoFlavorChargeH[flavour][net]->Fill(muonTrkMass,weight);
             PartonMomLowMuLooseIsoFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weight);
             InvMassLowMuLooseIsoFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weight);
           }

         }
	 else {
	   if (muHighPassed) {
	     MuonMomHighMuLooseIsoUnmatchedH->Fill(Muon4.P(),weight);
	     InvMassHighMuLooseIsoMomUnmatchedH[muonMomBin]->Fill(muonTrkMass,weight);
	   }
	   if (muLowPassed) {
	     MuonMomLowMuLooseIsoUnmatchedH->Fill(Muon4.P(),weight);
	     InvMassLowMuLooseIsoMomUnmatchedH[muonMomBin]->Fill(muonTrkMass,weight);
	   }

	 }
       }


       if (sbMu) {

	 nSbMuons++;
	 ptSbTrackH->Fill(trackLV.Pt(),weight);
         etaSbTrackH->Fill(trackLV.Eta(),weight);
         dxySbTrackH->Fill(track_dxy[indexTrk],weight);
         dzSbTrackH->Fill(track_dz[indexTrk],weight);
         ptSbMuH->Fill(Muon4.Pt(),weight);
         etaSbMuH->Fill(Muon4.Eta(),weight);
         dxySbMuH->Fill(muon_dxy[index],weight);
         dzSbMuH->Fill(muon_dz[index],weight);
         PartonMultiplicityMuSbH->Fill(float(partons),weight);
         PFJetMultiplicityMuSbH->Fill(float(pfjets),weight);
         GenJetMultiplicityMuSbH->Fill(float(genjets),weight);
	 deltaRMuTrkSbH->Fill(deltaRMuonTrk,weight);

 	 if (matchedParton) {
           deltaRPartonMuSbH->Fill(dRmin,weight);
           if (muHighPassed) {
             InvMassHighMuSbFlavorChargeH[flavour][net]->Fill(muonTrkMass,weight);
             PartonMomHighMuSbFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weight);
             InvMassHighMuSbFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weight);
           }
           if (muLowPassed) {
             InvMassLowMuSbFlavorChargeH[flavour][net]->Fill(muonTrkMass,weight);
             PartonMomLowMuSbFlavorChargeH[flavour][net]->Fill(matchedPartonLV.P(),weight);
             InvMassLowMuSbFlavorChargeMomH[flavour][net][partonMomBin]->Fill(muonTrkMass,weight);
           }
         }
	 else {
	   if (muHighPassed) {
	     MuonMomHighMuSbUnmatchedH->Fill(Muon4.P(),weight);
	     InvMassHighMuSbMomUnmatchedH[muonMomBin]->Fill(muonTrkMass,weight);
	   }
	   if (muLowPassed) {
	     MuonMomLowMuSbUnmatchedH->Fill(Muon4.P(),weight);
	     InvMassLowMuSbMomUnmatchedH[muonMomBin]->Fill(muonTrkMass,weight);
	   }
	   
	 }
       }

     }

     nGoodIsoMuonsH->Fill(float(nIsoMuons),weight);
     nGoodLooseIsoMuonsH->Fill(float(nLooseIsoMuons),weight);
     nGoodSbMuonsH->Fill(float(nSbMuons),weight);

     for (unsigned int im=0; im<muons.size();++im) {


       int muonType = muons_type[im];
       int net = muons_net[im];
       int flavour = muons_flavour[im];
       int momBin = muons_mom[im];
       int region = muons_region[im];
       float mutrk_mass = muons_mutrkmass[im];
       int pdgid = muons_pdgid[im];

       if (region==0) {
	 InvMassIsoH->Fill(mutrk_mass,weight);
	 if (muonType==0) InvMassHighMuIsoH->Fill(mutrk_mass,weight);
	 if (muonType==1) InvMassLowMuIsoH->Fill(mutrk_mass,weight);
       }

       if (region==1) {
	 InvMassLooseIsoH->Fill(mutrk_mass,weight);
	 if (muonType==0) InvMassHighMuLooseIsoH->Fill(mutrk_mass,weight);
	 if (muonType==1) InvMassLowMuLooseIsoH->Fill(mutrk_mass,weight);
       }
       
       for (int iM=0; iM<20; ++iM) {
	 double mass = double(iM) + double(0.5);
	 int ireg = 0;
	 double pdf = qcdModel->getMuMassPdf(momBin,muonType,ireg,flavour,net,mass);
	 ModelInvMassIsoH->Fill(mass,weight*pdf);
	 ireg = 1;
	 pdf = qcdModel->getMuMassPdf(momBin,muonType,ireg,flavour,net,mass);
	 ModelInvMassLooseIsoH->Fill(mass,weight*pdf);
       }

     }

     // find muons
     if (muons.size()<2) continue;
     //     cout << muons.size() << endl;
     unsigned int mu1 = 0;
     unsigned int mu2 = 0;
     bool foundMuons = false;
     double ptmax = 0;
     for (unsigned int i1=0; i1<muons.size()-1; ++i1) {
       for (unsigned int i2=i1+1; i2<muons.size(); ++i2) {
	 unsigned int index1 = muons[i1];
	 unsigned int index2 = muons[i2];

	 //	 cout << index1 << " " << index2 << endl;

	 bool os = muon_charge[index1]*muon_charge[index2]<0;
	 if (os) continue;
	 bool high1 = muon_pt[index1]>ptMuonHighCut && fabs(muon_eta[index1])<etaMuonHighCut;
	 bool high2 = muon_pt[index2]>ptMuonHighCut && fabs(muon_eta[index2])<etaMuonHighCut;
	 bool passed = high1 || high2;
	 if (!passed) continue;
	 double dR = deltaR(muon_eta[index1],muon_phi[index1],
			    muon_eta[index2],muon_phi[index2]);
	 if (dR<dRMuonsCut) continue;
	 double ptsum = muon_pt[index1]+muon_pt[index2];
	 if (ptsum>ptmax) {
	   foundMuons = true;
	   mu1 = i1;
	   mu2 = i2;
	 }
       }
     }
     
     if (foundMuons) {


       int muonType1 = muons_type[mu1];
       int net1 = muons_net[mu1];
       int flavour1 = muons_flavour[mu1];
       int mom1 = muons_mom[mu1];
       int region1 = muons_region[mu1];
       float mutrk_mass1 = muons_mutrkmass[mu1];
       int pdgid1 = muons_pdgid[mu1];

       int muonType2 = muons_type[mu2];
       int net2 = muons_net[mu2];
       int flavour2 = muons_flavour[mu2];
       int mom2 = muons_mom[mu2];
       int region2 = muons_region[mu2];
       float mutrk_mass2 = muons_mutrkmass[mu2];
       int pdgid2 = muons_pdgid[mu2];

       if (flavour1>=0) {
	 partonMuSS[muonType1][flavour1][net1][mom1]->Fill(0.5,weight);
	 partonMuSSPass[muonType1][flavour1][net1][mom1]->Fill(0.5+double(region1),weight);
       }
       else {
	 unmatchedMuSS[muonType1][mom1]->Fill(0.5,weight);
	 unmatchedMuSSPass[muonType1][mom1]->Fill(0.5+double(region1),weight);	 
       }

       if (flavour2>=0) {
	 partonMuSS[muonType2][flavour2][net2][mom2]->Fill(0.5,weight);
	 partonMuSSPass[muonType2][flavour2][net2][mom2]->Fill(0.5+double(region2),weight);
       }
       else {
	 unmatchedMuSS[muonType2][mom2]->Fill(0.5,weight);
	 unmatchedMuSSPass[muonType2][mom2]->Fill(0.5+double(region2),weight);	 
       }

       unsigned int index1 = muons[mu1];
       unsigned int index2 = muons[mu2];

       if (region1==0)
	 InvMassDimuonIsoH->Fill(mutrk_mass1,weight);
       if (region1==1)
	 InvMassDimuonLooseIsoH->Fill(mutrk_mass1,weight);

       if (region2==0)
	 InvMassDimuonIsoH->Fill(mutrk_mass2,weight);
       if (region2==1)
	 InvMassDimuonLooseIsoH->Fill(mutrk_mass2,weight);


       if (region1==0||region2==0) {
	 cout << "muons found " << endl;
	 cout << "mu1 = " << mu1 
	      << "  pdgid = " << pdgid1
	      << "  flavor = " << flavour1
	      << "  charge = " << muon_charge[index1]
	      << "  net = " << net1 	   
	      << "  pt  = " << muon_pt[index1] 
	      << "  type = " << muonType1 
	      << "  reg = " << region1 << endl;
	 cout << "mu2 = " << mu2 
	      << "  pdgid = " << pdgid2 
	      << "  flavour = " << flavour2
	      << "  charge = " << muon_charge[index2]
	      << "  net = " << net2 
	      << "  pt  = " << muon_pt[index2] 
	      << "  type = " << muonType2 
	      << "  reg = " << region2 << endl;
       }

       for (unsigned int imass=0; imass<20; ++imass) {
	 double mass = 0.5 * double(imass);
	 int ireg = 0;
	 double pdf1 = qcdModel->getSSMuMassPdf(mom1,muonType1,ireg,flavour1,net1,mass);
	 double pdf2 = qcdModel->getSSMuMassPdf(mom2,muonType2,ireg,flavour2,net2,mass);
	 ModelInvMassDimuonIsoH->Fill(mass,weight*pdf1);
	 ModelInvMassDimuonIsoH->Fill(mass,weight*pdf2);
	 ireg = 1;
	 pdf1 = qcdModel->getSSMuMassPdf(mom1,muonType1,ireg,flavour1,net1,mass);
	 pdf2 = qcdModel->getSSMuMassPdf(mom2,muonType2,ireg,flavour2,net2,mass);
	 ModelInvMassDimuonLooseIsoH->Fill(mass,weight*pdf1);
	 ModelInvMassDimuonLooseIsoH->Fill(mass,weight*pdf2);
       }

       // SR and ControlX
       vector<TH1D*> array_1d = {InvMassH,InvMass_ControlXH};
       vector<TH2D*> array_2d = {InvMass2DH,InvMass2D_ControlXH};
       for (int ireg=0; ireg<2; ++ireg) {
	 double probMu1 = qcdModel->getProbSSIsoMu(mom1,muonType1,ireg,flavour1,net1);
	 double probMu2 = qcdModel->getProbSSIsoMu(mom2,muonType2,ireg,flavour2,net2);
	 double probDoubleMu = probMu1*probMu2;
	 

	 for (unsigned int imass=0; imass<20; ++imass) {
	   double mass = 0.5 * double(imass);
	   double pdf1 = qcdModel->getMassPdf(mom1,muonType1,ireg,flavour1,net1,mass);
	   array_1d[ireg]->Fill(mass,weight*probDoubleMu*pdf1);
	   double pdf2 = qcdModel->getMassPdf(mom2,muonType2,ireg,flavour2,net2,mass);
	   array_1d[ireg]->Fill(mass,weight*probDoubleMu*pdf2);	 
	   for (unsigned int imass2=0; imass2<20; ++imass2) {
	     double mass2 = 0.5 + double(imass2);
	     pdf2 = qcdModel->getMassPdf(mom2,muonType2,ireg,flavour2,net2,mass2);
	     array_2d[ireg]->Fill(mass,mass2,weight*probDoubleMu*pdf1*pdf2);
	   }
	 }
       }

       // ControlY
       for (int ireg1=0; ireg1<2; ++ireg1) {
	 for (int ireg2=0; ireg2<2; ++ireg2) {
	   if (ireg1==0&&ireg2==0) continue;
	   double probMu1 = qcdModel->getProbSSIsoMu(mom1,muonType1,ireg1,flavour1,net1);
	   double probMu2 = qcdModel->getProbSSIsoMu(mom2,muonType2,ireg2,flavour2,net2);
	   double probDoubleMu = probMu1*probMu2;
	   for (unsigned int imass=0; imass<20; ++imass) {
	     double mass = 0.5 * double(imass);
	     double pdf1 = qcdModel->getMassPdf(mom1,muonType1,ireg1,flavour1,net1,mass);
	     InvMass_ControlYH->Fill(mass,weight*probDoubleMu*pdf1);
	     double pdf2 = qcdModel->getMassPdf(mom2,muonType2,ireg2,flavour2,net2,mass);
	     InvMass_ControlYH->Fill(mass,weight*probDoubleMu*pdf2);
	     for (unsigned int imass2=0; imass2<20; ++imass2) {
	       double mass2 = 0.5 + double(imass2);
	       pdf2 = qcdModel->getMassPdf(mom2,muonType2,ireg2,flavour2,net2,mass2);
	       InvMass2D_ControlYH->Fill(mass,mass2,weight*probDoubleMu*pdf1*pdf2);
	     }
	   }
	 }
       }

     }

   }
   delete tree_;
   file_->Close();
   delete file_;
   
  }// filelist loop
  
  file->cd("");
  file->Write();
  file->Close();
  
  //delete file;
}// int main loop 

 

