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
#include "HtoAA/Utilities/src/Config.cc"
#include "TRandom.h"
#include "TRandom3.h"
#include "HtoAA/Utilities/interface/json.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "HtoAA/Utilities/interface/PileUp.h"
//#include "HTT-utilities/RecoilCorrections_KIT/interface/RecoilCorrector.h"
//#include "HTT-utilities/RecoilCorrections_KIT/interface/MEtSys.h"
#include "TSystem.h"
#include "HtoAA/Utilities/interface/functions.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

//using namespace kit;
using namespace std;

int main(int argc, char * argv[]) {
  
  if (argc<2) {
    std::cout << "Usage of the program : Hto4TausAnalysis [file_list]" << std::endl;
    std::cout << "file_list : file list of RooT files to be processed. To run on Data the string has to include the string \"Data\"." << std::endl;
    exit(1);
  }

  string cmsswBase = (getenv ("CMSSW_BASE"));

  // **** configuration
  Config cfg(argv[1]);

  const string jsonFile = cfg.get<string>("jsonFile");
  
  // kinematic cuts on muons
  const float ptMuonCut      = cfg.get<float>("ptMuonLowCut");
  const float ptMuonLooseCut = cfg.get<float>("ptMuonLooseCut");
  const float etaMuonCut  = cfg.get<float>("etaMuonHighCut");
  const float dxyMuonCut  = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut   = cfg.get<float>("dzMuonCut");
  const float isoMuonCut  = cfg.get<float>("isoMuonCut");
  const float isoMuonLooseCut  = cfg.get<float>("isoMuonLooseCut");
  const float dRdimuon = cfg.get<float>("dRdimuon");

  // topological cuts
  const float dPhiMuonTrkCut   = cfg.get<float>("dPhiMuonTrkCut");

  // track selection
  const float dRIsoMuon       = cfg.get<float>("dRIsoMuon");

  const float ptTrkCut        = cfg.get<float>("ptTrkCut");
  const float ptTrkLooseCut   = cfg.get<float>("ptTrkLooseCut");
  const float etaTrkCut       = cfg.get<float>("etaTrkCut");
  const float dxyTrkLooseCut  = cfg.get<float>("dxyTrkLooseCut");
  const float dxyTrkCut       = cfg.get<float>("dxyTrkCut");
  const float dzTrkLooseCut   = cfg.get<float>("dzTrkLooseCut");
  const float dzTrkCut        = cfg.get<float>("dzTrkCut");

  // jets
  const float jetEtaCut = cfg.get<float>("jetEtaCut"); 
  const float jetPtCut  = cfg.get<float>("jetPtCut");
  const float deltaRJetLeptonCut = cfg.get<float>("deltaRJetLeptonCut");

  // mT cuts
  const float mTLowCut = cfg.get<float>("mTLowCut");
  const float mTHighCut = cfg.get<float>("mTHighCut");

  // trigger
  const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");
  const vector<string> singleMuonTriggerNames = cfg.get<vector<string>>("SingleMuonTriggerNames");
  const vector<string> singleMuonFilterNames = cfg.get<vector<string>>("SingleMuonFilterNames");

  std::cout << std::endl;
  vector<TString> SingleMuonTriggerNames;
  for (auto singleMuonTriggerName : singleMuonTriggerNames) {
    std::cout << "HLT path " << singleMuonTriggerName << " is included " << std::endl;
    SingleMuonTriggerNames.push_back(TString(singleMuonTriggerName));
  }
  std::cout << std::endl;
  vector<TString> SingleMuonFilterNames; 
  for (auto singleMuonFilterName : singleMuonFilterNames) {
    std::cout << "HLT filter " << singleMuonFilterName << " is included " << std::endl;
    SingleMuonFilterNames.push_back(TString(singleMuonFilterName));
  }
  std::cout << std::endl;

  // zptweight file
  // const string ZptweightFile = cfg.get<string>("ZptweightFile");

  // trigger matching
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch"); 

  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");

  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);

  const string MuonTriggerEffFile = cfg.get<string>("MuonTriggerEffFile");
  const string workspaceFileName  = cfg.get<string>("WorkspaceFileName");
  const string zptWorkspaceFileName = cfg.get<string>("ZPtWorkspaceFileName");
  TString WorkspaceFileName = TString(cmsswBase)+"/src/HtoAA/data/"+TString(workspaceFileName);
  TString ZptWorkspaceFileName = TString(cmsswBase)+"/src/HtoAA/data/"+TString(zptWorkspaceFileName);

  const bool puppiMET = cfg.get<bool>("PuppiMET");

  // sys uncertainties
  //  const float jetES = cfg.get<float>("JetES"); 
  //  const float unclusteredES = cfg.get<float>("UnclusteredES"); 
  //  const float zPtReweighting = cfg.get<float>("ZPtReweighting"); 
  //  const float topPtReweighting = cfg.get<float>("TopPtReweighting"); 
  //
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
  UInt_t muon_nMuonStations[100];
  UInt_t muon_nMuonHits[100];
  UInt_t muon_nPixelHits[100];
  UInt_t muon_nTrackerHits[100];
  float muon_px[100];
  float muon_py[100];
  float muon_pz[100];
  float muon_pt[100];
  float muon_eta[100];
  float muon_phi[100];
  float muon_pterror[100];
  float muon_chi2[100];
  float muon_ndof[100];
  float muon_charge[100];
  float muon_dxy[100];
  float muon_dxyerr[100];
  float muon_dz[100];
  float muon_dzerr[100];
  float muon_chargedHadIso[100];
  float muon_neutralHadIso[100];
  float muon_photonIso[100];
  float muon_puIso[100];
  bool muon_isPF[100];
  bool muon_isGlobal[100];
  bool muon_isTracker[100];
  bool muon_isTight[100];
  bool muon_isLoose[100];
  bool muon_isMedium[100];

  UInt_t genparticles_count;
  Float_t genparticles_e[1000];
  Float_t genparticles_px[1000];
  Float_t genparticles_py[1000];
  Float_t genparticles_pz[1000];
  Int_t genparticles_pdgid[1000];
  Int_t genparticles_status[1000];
  UInt_t genparticles_info[1000];
  Int_t genparticles_fromHardProcess[1000];
  Int_t genparticles_fromHardProcessBeforeFSR[1000];
  Int_t genparticles_isDecayedLeptonHadron[1000];
  Int_t genparticles_isDirectHadronDecayProduct[1000];
  Int_t genparticles_isDirectHardProcessTauDecayProduct[1000];
  Int_t genparticles_isDirectPromptTauDecayProduct[1000];
  Int_t genparticles_isDirectTauDecayProduct[1000];
  Int_t genparticles_isFirstCopy[1000];
  Int_t genparticles_isHardProcess[1000];
  Int_t genparticles_isHardProcessTauDecayProduct[1000];
  Int_t genparticles_isLastCopy[1000];
  Int_t genparticles_isLastCopyBeforeFSR[1000];
  Int_t genparticles_isPrompt[1000];
  Int_t genparticles_isPromptTauDecayProduct[1000];
  Int_t genparticles_isTauDecayProduct[1000];

  UInt_t          gentau_count;
  Float_t         gentau_visible_e[100]; 
  Float_t         gentau_visible_px[100];
  Float_t         gentau_visible_py[100];
  Float_t         gentau_visible_pz[100];
  Int_t           gentau_decayMode[100]; 
  Int_t           gentau_status[100]; 
  Int_t           gentau_isPrompt[100]; 
  Int_t           gentau_fromHardProcess[100];

  float genweight;

  UInt_t genparticles_noutgoing;

  float metx;
  float mety;
  float met;
  float metphi;
  float metx_JetEnDown;
  float mety_JetEnDown;
  float metx_JetEnUp;
  float mety_JetEnUp;
  float metx_UnclusteredEnDown;
  float mety_UnclusteredEnDown;
  float metx_UnclusteredEnUp;
  float mety_UnclusteredEnUp;
  
   // Trigger
  unsigned int trigobject_count;
  float trigobject_px[1000];
  float trigobject_py[1000];
  float trigobject_pz[1000];
  float trigobject_pt[1000];
  float  trigobject_eta[1000];
  float trigobject_phi[1000];
  bool trigobject_filters[1000][200];

  // pat jets                                                                                                                                                                                                      
  UInt_t pfjet_count;
  Float_t pfjet_e[1000];
  Float_t pfjet_px[1000];
  Float_t pfjet_py[1000];
  Float_t pfjet_pz[1000];
  Float_t pfjet_pt[1000];
  Float_t pfjet_eta[1000];
  Float_t pfjet_phi[1000];
  Float_t pfjet_energycorr[1000];

  Float_t pfjet_neutralhadronicenergy[1000];
  Float_t pfjet_chargedhadronicenergy[1000];
  Float_t pfjet_neutralemenergy[1000];
  Float_t pfjet_chargedemenergy[1000];
  Float_t pfjet_muonenergy[1000];
  Float_t pfjet_chargedmuonenergy[1000];
  UInt_t pfjet_chargedmulti[1000];
  UInt_t pfjet_neutralmulti[1000];
  UInt_t pfjet_chargedhadronmulti[1000];

  float numtruepileupinteractions;

  //unsigned int iLeadingPosTrig = 0;
  //vector<bool> trigobject_filter; trigobject_filter.clear();

  std::map<std::string, int> * hltriggerresults = new std::map<std::string, int>() ;
  std::map<std::string, int> * hltriggerprescales = new std::map<std::string, int>() ;
  std::vector<std::string>   * hltfilters = new std::vector<std::string>();

  std::string rootFileName(argv[2]);
  
  std::string chainName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");
  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;

  bool isData = TStrName.Contains("SingleMuon") || TStrName.Contains("DoubleMuon");
  bool isTOP  = TStrName.Contains("TTTo");
  bool isZTT  = TStrName.Contains("DYJetsToTT");
  bool isDY   = TStrName.Contains("DYJets");

  bool applyGoodRunSelection = isData;

  if (isData) {
    std::cout << "=========================" << std::endl;
    std::cout << "=== Running on Data =====" << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << std::endl;
  }
  else {
    std::cout << "=====================" << std::endl;
    std::cout << "=== Running on MC ===" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << std::endl; 
    if (isDY) {
      std::cout << "=== Running on DYJets =====" << std::endl;
      std::cout << "activated Z pt reweighting " << std::endl;
      std::cout << std::endl;
      if (isZTT) {
	std::cout<< "== Selecting Z->tautau ===" << std::endl;
	std::cout<<std::endl;
      }
    }
    if (isTOP) {
      std::cout<< "===== Running on TTbar =====" << std::endl;
      std::cout<< "activated Top pt reweighting" << std::endl;
      std::cout<< std::endl;
    }

  }

  TString FullName = TStrName;      
  
  TFile * file = new TFile(FullName+TString(".root"),"recreate");

  file->cd("");
  
  // weights
  TH1D * puWeightH = new TH1D("puWeightH","",250,0,5);
  TH1D * zptWeightH = new TH1D("zptWeightH","",250,0,5);
  TH1D * topWeightH = new TH1D("topWeightH","",250,0,5);
  TH1D * muTrigWeightH = new TH1D("muTrigWeightH","",100,0,2);
  TH1D * muIdWeightH = new TH1D("muIdWeightH","",100,0,2);
  TH1D * histWeightsH = new TH1D("histWeightsH","",1,0.,2.);

  // histograms not after final selection
  TH1D * nCloseTrksH         = new TH1D("nCloseTrksH","",20,0,20);
  TH1D * nCloseTrksH_2p5to5  = new TH1D("nCloseTrksH_2p5to5","",20,0,20);
  TH1D * nCloseTrksH_5to10 = new TH1D("nCloseTrksH_5to10","",20,0,20);
  TH1D * nCloseTrksH_10to15  = new TH1D("nCloseTrksH_10to15","",20,0,20);
  TH1D * nCloseTrksH_15to20  = new TH1D("nCloseTrksH_15to20","",20,0,20);
  TH1D * nCloseTrksH_20toInf = new TH1D("nCloseTrksH_20toInf","",20,0,20);

  // histograms after final inclusive selection in the signal region
  // isolated muon, OS muon-track pair, low mT region
  TH1D * mTH         = new TH1D("mTH","",500,0,500);
  TH1D * mTH_2p5to5  = new TH1D("mTH_2p5to5","",500,0,500);
  TH1D * mTH_5to10   = new TH1D("mTH_5to10","",500,0,500);
  TH1D * mTH_10to15  = new TH1D("mTH_10to15","",500,0,500);
  TH1D * mTH_15to20  = new TH1D("mTH_15to20","",500,0,500);
  TH1D * mTH_20toInf = new TH1D("mTH_20toInf","",500,0,500);

  TH1D * ptMuonH         = new TH1D("ptMuonH","",100,0,500);
  TH1D * etaMuonH = new TH1D("etaMuonH","",25,-2.5,2.5);
  TH1D * ptMuonH_2p5to5  = new TH1D("ptMuonH_2p5to5","",100,0,500);
  TH1D * ptMuonH_5to10 = new TH1D("ptMuonH_5to10","",100,0,500);
  TH1D * ptMuonH_10to15  = new TH1D("ptMuonH_10to15","",100,0,500);
  TH1D * ptMuonH_15to20  = new TH1D("ptMuonH_15to20","",100,0,500);
  TH1D * ptMuonH_20toInf = new TH1D("ptMuonH_20toInf","",100,0,500);

  TH1D * ptTrackH = new TH1D("ptTrackH","",400,0,400);
  TH1D * ptTrackLowH = new TH1D("ptTrackLowH","",100,0,50);
  TH1D * etaTrackH = new TH1D("etaTrackH","",25,-2.5,2.5);
  TH1D * dzTrackH = new TH1D("dzTrackH","",100,0,1.);
  TH1D * dxyTrackH = new TH1D("dxyTrackH","",100,0,1.0);

  TH1D * metH         = new TH1D("metH","",500,0,500);
  TH1D * metH_2p5to5  = new TH1D("metH_2p5to5","",500,0,500);
  TH1D * metH_5to10 = new TH1D("metH_5to10","",500,0,500);
  TH1D * metH_10to15  = new TH1D("metH_10to15","",500,0,500);
  TH1D * metH_15to20  = new TH1D("metH_15to20","",500,0,500);
  TH1D * metH_20toInf = new TH1D("metH_20toInf","",500,0,500);

  TH1D * dPhiMuonTrkH = new TH1D("dPhiMuonTrkH","",32,0,3.15);
  TH1D * dPhiMuonTrkH_2p5to5  = new TH1D("dPhiMuonTrkH_2p5to5","",32,0,3.15);
  TH1D * dPhiMuonTrkH_5to10 = new TH1D("dPhiMuonTrkH_5to10","",32,0,3.15);
  TH1D * dPhiMuonTrkH_10to15  = new TH1D("dPhiMuonTrkH_10to15","",32,0,3.15);
  TH1D * dPhiMuonTrkH_15to20  = new TH1D("dPhiMuonTrkH_15to20","",32,0,3.15);
  TH1D * dPhiMuonTrkH_20toInf = new TH1D("dPhiMuonTrkH_20toInf","",32,0,3.15);

  TH1D * dPhiMuonMetH = new TH1D("dPhiMuonMetH","",32,0,3.15);
  TH1D * dPhiMuonMetH_2p5to5  = new TH1D("dPhiMuonMetH_2p5to5","",32,0,3.15);
  TH1D * dPhiMuonMetH_5to10 = new TH1D("dPhiMuonMetH_5to10","",32,0,3.15);
  TH1D * dPhiMuonMetH_10to15  = new TH1D("dPhiMuonMetH_10to15","",32,0,3.15);
  TH1D * dPhiMuonMetH_15to20  = new TH1D("dPhiMuonMetH_15to20","",32,0,3.15);
  TH1D * dPhiMuonMetH_20toInf = new TH1D("dPhiMuonMetH_20toInf","",32,0,3.15);

  TH1D * dPhiTrkMetH = new TH1D("dPhiTrkMetH","",32,0,3.15);
  TH1D * dPhiTrkMetH_2p5to5  = new TH1D("dPhiTrkMetH_2p5to5","",32,0,3.15);
  TH1D * dPhiTrkMetH_5to10 = new TH1D("dPhiTrkMetH_5to10","",32,0,3.15);
  TH1D * dPhiTrkMetH_10to15  = new TH1D("dPhiTrkMetH_10to15","",32,0,3.15);
  TH1D * dPhiTrkMetH_15to20  = new TH1D("dPhiTrkMetH_15to20","",32,0,3.15);
  TH1D * dPhiTrkMetH_20toInf = new TH1D("dPhiTrkMetH_20toInf","",32,0,3.15);

  TH1D * dPhiMuonJetH = new TH1D("dPhiMuonJetH","",32,0,3.15);
  TH1D * dPhiMuonJetH_2p5to5  = new TH1D("dPhiMuonJetH_2p5to5","",32,0,3.15);
  TH1D * dPhiMuonJetH_5to10 = new TH1D("dPhiMuonJetH_5to10","",32,0,3.15);
  TH1D * dPhiMuonJetH_10to15  = new TH1D("dPhiMuonJetH_10to15","",32,0,3.15);
  TH1D * dPhiMuonJetH_15to20  = new TH1D("dPhiMuonJetH_15to20","",32,0,3.15);
  TH1D * dPhiMuonJetH_20toInf = new TH1D("dPhiMuonJetH_20toInf","",32,0,3.15);

  TH1D * dzetaH         = new TH1D("dzetaH","",300,-200,100);
  TH1D * dzetaH_2p5to5  = new TH1D("dzetaH_2p5to5","",300,-200,100);
  TH1D * dzetaH_5to10   = new TH1D("dzetaH_5to10","",300,-200,100);
  TH1D * dzetaH_10to15  = new TH1D("dzetaH_10to15","",300,-200,100);
  TH1D * dzetaH_15to20  = new TH1D("dzetaH_15to20","",300,-200,100);
  TH1D * dzetaH_20toInf = new TH1D("dzetaH_20toInf","",300,-200,100);

  TH1D * nJetH = new TH1D("nJetH","",10,0,10);

  // ------------------------------------------------------
  // Z->mumu control region
  TH1D * massMuMuH    = new TH1D("massMuMuH","",60,60.,120.);

  // ********************************************
  // * mu-trk invariant mass in various regions *
  // ********************************************
  std::vector<TString> muon_isolations = {"iso","antiiso"};
  std::vector<TString> mt_regions = {"lowMT","highMT"};
  std::vector<TString> sign_regions = {"OS","SS"};
  std::vector<TString> pt_bins {"2p5to5","5to10","10to15","15to20","20toInf"};
  std::map<TString,TH1D*> mass_histos;
  TString baseName("invMassMuTrk");
  for (auto mt_region : mt_regions) {
    for (auto muon_isolation : muon_isolations) {
      for (auto sign_region : sign_regions) {
	TString histName = baseName + "_" + mt_region + "_" + muon_isolation + "_" + sign_region;
	mass_histos[histName] = new TH1D(histName,"",500,0,500);
	for (auto pt_bin : pt_bins) {
	  histName = baseName + "_" + pt_bin + "_" + mt_region +"_" + muon_isolation + "_" + sign_region;
	  mass_histos[histName] =new TH1D(histName,"",500,0,500);
	}
      }
    }
  }

  // ------------------------------
  // ---------- tuple -------------
  // ------------------------------

  bool os_tree; // charge(mu)*charge

  float mt_tree; // mT
  float mtJetEnUp_tree; // mT (JetEnUp)
  float mtJetEnDown_tree; // mT (JetEnDown)
  float mtUnclEnUp_tree; // mT (UnclEnUp)
  float mtUnclEnDown_tree; // mT (UnclEnDown)

  float met_tree; // met
  float metJetEnUp_tree; // met (JetEnUp)
  float metJetEnDown_tree; // met (JetEnDown)
  float metUnclEnUp_tree; // met (UnclEnUp)
  float metUnclEnDown_tree; // met (UnclEnDown)

  float dzeta_tree; // dzeta
  float dzetaJetEnUp_tree; // dzeta (JetEnUp)
  float dzetaJetEnDown_tree; // dzeta (JetEnDown)
  float dzetaUnclEnUp_tree; // dzeta (UnclEnUp)
  float dzetaUnclEnDown_tree; // dzeta (UnclEnDown)

  unsigned int njets_tree; // njets

  float massmutrk_tree; // mass(mu,trk)
  float dphimutrk_tree; // dphi(mu,trk)

  float muiso_tree; // muon isolation
  float mupt_tree; // muon pT
  float mueta_tree; // muon eta

  float trkpt_tree; // trk pT
  float trketa_tree; // trk eta
  int trktyp_tree; // gen matching of trk 

  float weight_tree; // total weight
  float zptweight_tree; // zpt weight
  float topptweight_tree; // top weight
  float puweight_tree; // puweight
  float genweight_tree; // genweight
  float muweight_tree; // mu id/iso weight

  unsigned int npart_tree; // number of partons

  bool cone0p15_tree; // isolation in random cone : dR(trk,cone)=0.15
  bool cone0p30_tree; // isolation in random cone : dR(trk,cone)=0.30
  bool cone0p45_tree; // isolation in random cone : dR(trk,cone)=0.45

  // systematic variations 
  std::vector<TString> uncNames = {"JetEnUp","JetEnDown","UnclEnUp","UnclEnDown"};  

  // nutrk tree used for creation of datacards
  TTree * mutrkTree = new TTree("mutrkTree","Muon Track Tree");

  mutrkTree->Branch("os",&os_tree,"os/O");

  mutrkTree->Branch("mt",&mt_tree,"mt/F"); 
  mutrkTree->Branch("mtJetEnUp",&mtJetEnUp_tree,"mtJetEnUp/F");
  mutrkTree->Branch("mtJetEnDown",&mtJetEnDown_tree,"mtJetEnDown/F");
  mutrkTree->Branch("mtUnclEnUp",&mtUnclEnUp_tree,"mtUnclEnUp/F");
  mutrkTree->Branch("mtUnclEnDown",&mtUnclEnDown_tree,"mtUnclEnDown/F");

  mutrkTree->Branch("met",&met_tree,"met/F");
  mutrkTree->Branch("metJetEnUp",&metJetEnUp_tree,"metJetEnUp/F");
  mutrkTree->Branch("metJetEnDown",&metJetEnDown_tree,"metJetEnDown/F");
  mutrkTree->Branch("metUnclEnUp",&metUnclEnUp_tree,"metUnclEnUp/F");
  mutrkTree->Branch("metUnclEnDown",&metUnclEnDown_tree,"metUnclEnDown/F");

  mutrkTree->Branch("dzeta",&dzeta_tree,"dzeta/F");
  mutrkTree->Branch("dzetaJetEnUp",&dzetaJetEnUp_tree,"dzetaJetEnUp/F");
  mutrkTree->Branch("dzetaJetEnDown",&dzetaJetEnDown_tree,"dzetaJetEnDown/F");
  mutrkTree->Branch("dzetaUnclEnUp",&dzetaUnclEnUp_tree,"dzetaUnclEnUp/F");
  mutrkTree->Branch("dzetaUnclEnDown",&dzetaUnclEnDown_tree,"dzetaUnclEnDown/F");

  mutrkTree->Branch("njets",&njets_tree,"njets/i");
  mutrkTree->Branch("massmutrk",&massmutrk_tree,"massmutrk/F");
  mutrkTree->Branch("dphimutrk",&dphimutrk_tree,"dphimutrk/F");
  mutrkTree->Branch("trktyp",&trktyp_tree,"trktyp/I");
  mutrkTree->Branch("trkpt",&trkpt_tree,"trkpt/F");
  mutrkTree->Branch("trketa",&trketa_tree,"trketa/F");
  mutrkTree->Branch("muiso",&muiso_tree,"muiso/F");
  mutrkTree->Branch("mupt",&mupt_tree,"mupt/F");
  mutrkTree->Branch("mueta",&mueta_tree,"mueta/F");
  mutrkTree->Branch("weight",&weight_tree,"weight/F");
  mutrkTree->Branch("zptweight",&zptweight_tree,"zptweight/F");
  mutrkTree->Branch("topptweight",&topptweight_tree,"topptweight/F");
  mutrkTree->Branch("genweight",&genweight_tree,"genweight/F");
  mutrkTree->Branch("puweight",&puweight_tree,"puweight/F");
  mutrkTree->Branch("muweight",&muweight_tree,"muweight/F");
  mutrkTree->Branch("npart",&npart_tree,"npart/i");

  mutrkTree->Branch("cone0p15",&cone0p15_tree,"cone0p15/O");
  mutrkTree->Branch("cone0p30",&cone0p30_tree,"cone0p30/O");
  mutrkTree->Branch("cone0p45",&cone0p45_tree,"cone0p45/O");

  // Load CrystalBallEfficiency class
  TString pathToCrystalLib = (TString) cmsswBase + "/src/HTT-utilities/CorrectionsWorkspace/CrystalBallEfficiency_cxx.so";
  int openSuccessful = gSystem->Load( pathToCrystalLib );
  if (openSuccessful !=0 ) {
    cout<<pathToCrystalLib<<" not found. Please create this file by running \"root -l -q CrystalBallEfficiency.cxx++\" in src/HTT-utilities/CorrectionsWorkspace/. "<<endl;
    exit( -1 );
  }

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

  // Trigger efficiencies
  //  ScaleFactor * SF_muonId = new ScaleFactor();
  //  SF_muonId->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdEffFile));
  TFile * correctionWorkSpaceFile = new TFile(WorkspaceFileName);
  RooWorkspace *correctionWS = (RooWorkspace*)correctionWorkSpaceFile->Get("w");

  // ZPt corrections
  TFile * zptWorkSpaceFile = new TFile(ZptWorkspaceFileName);
  RooWorkspace *zptCorr = (RooWorkspace*)zptWorkSpaceFile->Get("w");

  ScaleFactor * SF_muonTrigger = new ScaleFactor();
  SF_muonTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTriggerEffFile));

  // Zpt reweighting for LO DY samples from correction WS                    
  //  TString correctionsWorkspaceFileName = TString(cmsswBase)+"/src/"+ZptweightFile;
  //  TFile * correctionWorkSpaceFile = new TFile(correctionsWorkspaceFileName);
  //  RooWorkspace *correctionWS = (RooWorkspace*)correctionWorkSpaceFile->Get("w");

  TRandom3 *rand = new TRandom3();

  TString filen;
  int iFiles = 0;
  int events = 0;

  while (fileList >> filen) {
   iFiles++;
   cout << "file " << iFiles << " : " << filen << endl;
   
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

   TTree * tree_ = (TTree*)file_->Get(TString(chainName));
   TTree * inittree_ = (TTree*)file_->Get(TString(initNtupleName));
   if (tree_==NULL) continue;
   if (inittree_==NULL) continue;
   
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
   if (puppiMET) {
     tree_->SetBranchAddress("puppimet_ex", &metx);
     tree_->SetBranchAddress("puppimet_ey", &mety);
     tree_->SetBranchAddress("puppimet_pt", &met);
     tree_->SetBranchAddress("puppimet_phi",&metphi);
     tree_->SetBranchAddress("puppimet_ex_JetEnDown",&metx_JetEnDown);
     tree_->SetBranchAddress("puppimet_ey_JetEnDown",&mety_JetEnDown);
     tree_->SetBranchAddress("puppimet_ex_JetEnUp",&metx_JetEnUp);
     tree_->SetBranchAddress("puppimet_ey_JetEnUp",&mety_JetEnUp);
     tree_->SetBranchAddress("puppimet_ex_UnclusteredEnDown",&metx_UnclusteredEnDown);
     tree_->SetBranchAddress("puppimet_ey_UnclusteredEnDown",&mety_UnclusteredEnDown);
     tree_->SetBranchAddress("puppimet_ex_UnclusteredEnUp",&metx_UnclusteredEnUp);
     tree_->SetBranchAddress("puppimet_ey_UnclusteredEnUp",&mety_UnclusteredEnUp);
   }
   else {
     tree_->SetBranchAddress("pfmetcorr_ex", &metx);
     tree_->SetBranchAddress("pfmetcorr_ey", &mety);
     tree_->SetBranchAddress("pfmetcorr_pt", &met);
     tree_->SetBranchAddress("pfmetcorr_phi",&metphi);
     tree_->SetBranchAddress("pfmetcorr_ex_JetEnDown",&metx_JetEnDown);
     tree_->SetBranchAddress("pfmetcorr_ey_JetEnDown",&mety_JetEnDown);
     tree_->SetBranchAddress("pfmetcorr_ex_JetEnUp",&metx_JetEnUp);
     tree_->SetBranchAddress("pfmetcorr_ey_JetEnUp",&mety_JetEnUp);
     tree_->SetBranchAddress("pfmetcorr_ex_UnclusteredEnDown",&metx_UnclusteredEnDown);
     tree_->SetBranchAddress("pfmetcorr_ey_UnclusteredEnDown",&mety_UnclusteredEnDown);
     tree_->SetBranchAddress("pfmetcorr_ex_UnclusteredEnUp",&metx_UnclusteredEnUp);
     tree_->SetBranchAddress("pfmetcorr_ey_UnclusteredEnUp",&mety_UnclusteredEnUp);
   }

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
   //   tree_->SetBranchAddress("run_btagdiscriminators", &run_btagdiscriminators);
   tree_->SetBranchAddress("hltriggerresults",&hltriggerresults);
   tree_->SetBranchAddress("hltriggerprescales",&hltriggerprescales);

   tree_->SetBranchAddress("numtruepileupinteractions",&numtruepileupinteractions);

   if (!isData) {
     inittree_->SetBranchAddress("genweight",&genweight);
     tree_->SetBranchAddress("genweight",&genweight);
     tree_->SetBranchAddress("genparticles_count", &genparticles_count);
     tree_->SetBranchAddress("genparticles_e", genparticles_e);
     tree_->SetBranchAddress("genparticles_px", genparticles_px);
     tree_->SetBranchAddress("genparticles_py", genparticles_py);
     tree_->SetBranchAddress("genparticles_pz", genparticles_pz);
     tree_->SetBranchAddress("genparticles_pdgid", genparticles_pdgid);
     tree_->SetBranchAddress("genparticles_status", genparticles_status);
     tree_->SetBranchAddress("genparticles_info", genparticles_info);
     tree_->SetBranchAddress("genparticles_fromHardProcess", genparticles_fromHardProcess);
     tree_->SetBranchAddress("genparticles_isDirectHardProcessTauDecayProduct", genparticles_isDirectHardProcessTauDecayProduct);
     tree_->SetBranchAddress("gentau_count", &gentau_count);
     tree_->SetBranchAddress("gentau_visible_e", gentau_visible_e);
     tree_->SetBranchAddress("gentau_visible_px", gentau_visible_px);
     tree_->SetBranchAddress("gentau_visible_py", gentau_visible_py);
     tree_->SetBranchAddress("gentau_visible_pz", gentau_visible_pz);
     tree_->SetBranchAddress("gentau_status", gentau_status);
     tree_->SetBranchAddress("gentau_fromHardProcess", gentau_fromHardProcess);
     tree_->SetBranchAddress("gentau_isPrompt", gentau_isPrompt);
     tree_->SetBranchAddress("gentau_decayMode", gentau_decayMode);
     tree_->SetBranchAddress("genparticles_noutgoing", &genparticles_noutgoing);
   }   

   tree_->SetBranchAddress("pfjet_count",&pfjet_count);
   tree_->SetBranchAddress("pfjet_e",pfjet_e);
   tree_->SetBranchAddress("pfjet_px",pfjet_px);
   tree_->SetBranchAddress("pfjet_py",pfjet_py);
   tree_->SetBranchAddress("pfjet_pz",pfjet_pz);
   tree_->SetBranchAddress("pfjet_pt",pfjet_pt);
   tree_->SetBranchAddress("pfjet_eta",pfjet_eta);
   tree_->SetBranchAddress("pfjet_phi",pfjet_phi);
   tree_->SetBranchAddress("pfjet_energycorr",pfjet_energycorr);
   tree_->SetBranchAddress("pfjet_neutralhadronicenergy",pfjet_neutralhadronicenergy);
   tree_->SetBranchAddress("pfjet_chargedhadronicenergy",pfjet_chargedhadronicenergy);
   tree_->SetBranchAddress("pfjet_neutralemenergy",pfjet_neutralemenergy);
   tree_->SetBranchAddress("pfjet_chargedemenergy",pfjet_chargedemenergy);
   tree_->SetBranchAddress("pfjet_muonenergy",pfjet_muonenergy);
   tree_->SetBranchAddress("pfjet_chargedmulti",pfjet_chargedmulti);
   tree_->SetBranchAddress("pfjet_neutralmulti",pfjet_neutralmulti);
   tree_->SetBranchAddress("pfjet_chargedhadronmulti",pfjet_chargedhadronmulti);

   int numberOfCandidates = tree_->GetEntries();

   std::cout << "number of events = " << numberOfCandidates << std::endl;

   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   //%%%%%%%%%%%%  Event Loop %%%%%%%%%%%%%%%%
   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for (int iCand=0; iCand<inittree_->GetEntries(); iCand++) {
     
     inittree_->GetEntry(iCand);

     float generator_weight = 1.0;
     if (!isData) {
       if (genweight<0.0) 
	 generator_weight = -1.0;
     }
     histWeightsH->Fill(1.0,generator_weight);
     
   }

   delete inittree_;

   bool isFirst = true;

   for (int iCand=0; iCand<numberOfCandidates; iCand++) {
     
     tree_->GetEntry(iCand);

     events++;
     if (events%10000==0) cout << "   processed events : " << events << endl;

     float weight = 1.0;
     float generator_weight = 1.0;
     if (!isData) {
       if (genweight<0.)
	 generator_weight = -1.0;
       weight *= generator_weight;
       //       std::cout << "genweight = " << genweight << std::endl;
     }

     // ---------- Systematic uncertainties -----------
     
     map<TString, TLorentzVector> metSys;

     float met_JetEnUp = TMath::Sqrt(metx_JetEnUp*metx_JetEnUp+
				     mety_JetEnUp*mety_JetEnUp);
     metSys["JetEnUp"] = TLorentzVector(metx_JetEnUp,
					mety_JetEnUp,
					0.,
					met_JetEnUp); 


     float met_JetEnDown = TMath::Sqrt(metx_JetEnDown*metx_JetEnDown+
				       mety_JetEnDown*mety_JetEnDown);
     metSys["JetEnDown"] = TLorentzVector(metx_JetEnDown,
					  mety_JetEnDown,
					  0.,
					  met_JetEnDown);


     float met_UnclEnUp = TMath::Sqrt(metx_UnclusteredEnUp*metx_UnclusteredEnUp+
					mety_UnclusteredEnUp*mety_UnclusteredEnUp);
     metSys["UnclEnUp"] = TLorentzVector(metx_UnclusteredEnUp,
					 mety_UnclusteredEnUp,
					 0.,
					   met_UnclEnUp);
     

     float met_UnclEnDown = TMath::Sqrt(metx_UnclusteredEnDown*metx_UnclusteredEnDown+
					mety_UnclusteredEnDown*mety_UnclusteredEnDown);
     metSys["UnclEnDown"] = TLorentzVector(metx_UnclusteredEnDown,
					   mety_UnclusteredEnDown,
					   0.,
					   met_UnclEnDown);
     // -----------------------------------------------------------------------------------------------


     float topPt = -1;
     float antitopPt = -1;
     TLorentzVector genBosonLV;genBosonLV.SetXYZT(0,0,0,0);
     TLorentzVector genVisBosonLV; genVisBosonLV.SetXYZT(0,0,0,0);

     if (!isData) {
       //       std::cout << "Generated particles = " << genparticles_count << std::endl;
       for (unsigned int igen=0; igen<genparticles_count; ++igen) {
	 TLorentzVector genLV; genLV.SetXYZT(genparticles_px[igen],
					     genparticles_py[igen],
					     genparticles_pz[igen],
					     genparticles_e[igen]);
	 
	 if (genparticles_pdgid[igen]==6){
	   topPt = TMath::Sqrt(genparticles_px[igen]*genparticles_px[igen]+
			       genparticles_py[igen]*genparticles_py[igen]);
	 }
	 if (genparticles_pdgid[igen]==-6){
	   antitopPt = TMath::Sqrt(genparticles_px[igen]*genparticles_px[igen]+
				   genparticles_py[igen]*genparticles_py[igen]);
	 }
	 bool fromHardProcessFinalState = genparticles_fromHardProcess[igen]&&genparticles_status[igen]==1;
	 bool isMuon = false;
	 bool isElectron = false;
	 bool isNeutrino = false;
	 bool isDirectHardProcessTauDecayProduct = genparticles_isDirectHardProcessTauDecayProduct[igen];

	 if (abs(genparticles_pdgid[igen])==11) isElectron = true;
	 if (abs(genparticles_pdgid[igen])==13) isMuon = true;
	 if (abs(genparticles_pdgid[igen])==12||
	     abs(genparticles_pdgid[igen])==14||
	     abs(genparticles_pdgid[igen])==16) isNeutrino = true;

	 bool isBoson = (fromHardProcessFinalState && (isMuon || isElectron || isNeutrino)) || isDirectHardProcessTauDecayProduct;
	 bool isVisibleBoson = (fromHardProcessFinalState && (isMuon || isElectron)) || (isDirectHardProcessTauDecayProduct && !isNeutrino);

	 if (isBoson)        genBosonLV += genLV;
	 if (isVisibleBoson) genVisBosonLV += genLV;
       }
     }
     
     
     float zptmassweight = 1;     
     float topptweight = 1;
     if (!isData) {
       
       if (isDY) { // applying Z pt mass weights
	 float bosonMass = genBosonLV.M();
	 float bosonPt = genBosonLV.Pt();
	 if (bosonMass>50.0) {
	   float bosonMassX = bosonMass;
	   float bosonPtX = bosonPt;
	   if (bosonMassX>1000.) bosonMassX = 1000.;
	   if (bosonPtX<1.)      bosonPtX = 1.;
	   if (bosonPtX>1000.)   bosonPtX = 1000.;
	   zptCorr->var("z_gen_pt")->setVal(bosonPtX);
	   zptCorr->var("z_gen_mass")->setVal(bosonMassX);
	   zptmassweight = zptCorr->function("zptmass_weight_nom")->getVal();
	   //	   std::cout << "bosonMass = " << bosonMass << "  bosonPt = " << bosonPt << " " << "  zptweight = " << zptmassweight << std::endl;
	 } 
	 weight *= zptmassweight; 
       }
     
       if (isTOP) { // applying top pT weights
	 if (topPt>0&&antitopPt>0) topptweight = topPtWeight(topPt,antitopPt,true);
	 weight *= topptweight;
       }
     }       
     zptWeightH->Fill(zptmassweight,1.0);
     topWeightH->Fill(topptweight,1.0);
     
     if (isData) {
       bool lumi = false;
       int n=event_run;
       int lum = event_luminosityblock;
       
       std::string num = std::to_string(n);
       std::string lnum = std::to_string(lum);
       for(const auto& a : periods)
	 {
	   if ( num.c_str() ==  a.name ) {
	     //	     std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
	     //	     std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	     
	     for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
	       
	       //	       cout<<b->lower<<"  "<<b->bigger<<endl;
	       if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
	     }
	     auto last = std::prev(a.ranges.end());
	     // std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	     if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
	   }
	 }
       if (!lumi) { 
	 //	 std::cout << "run lumi " << n << " " << lum << "rejected " << std::endl;
	 continue;
       }
     }
   

     float puweight = 1;
     if (!isData) {
       puweight = float(PUofficial->get_PUweight(double(numtruepileupinteractions)));
       //       std::cout << "n(true interactions) = " << numtruepileupinteractions << "   :  PU weight = " << puweight << std::endl; 
     }
     puWeightH->Fill(puweight,1.0);
     weight *= puweight;

     // finding HLT filters in the HLT Filter library
     std::map<TString, unsigned int> singleMuonFilters; 
     std::map<TString, bool> isSingleMuonFilterFound;
     for (auto SingleMuonFilterName : SingleMuonFilterNames) {
       singleMuonFilters[SingleMuonFilterName] = 0;
       isSingleMuonFilterFound[SingleMuonFilterName] = false;
     }

     unsigned int nfilters = hltfilters->size();
     for (unsigned int i=0; i<nfilters; ++i) {
       //       std::cout << hltfilters->at(i) << std::endl;
       TString HLTFilter(hltfilters->at(i));
       for (auto SingleMuonFilterName : SingleMuonFilterNames) {
	 if (HLTFilter==SingleMuonFilterName) {
	   singleMuonFilters[SingleMuonFilterName] = i;
	   isSingleMuonFilterFound[SingleMuonFilterName] = true;
	   if (isFirst)
	     std::cout << SingleMuonFilterName << " has index " << i << std::endl;
	 }
       }
     }

     bool allFiltersNotFound = true;
     for (auto SingleMuonFilterName : SingleMuonFilterNames) {
       bool isSingleMuonFilter = isSingleMuonFilterFound[SingleMuonFilterName];
       if (isSingleMuonFilter) {
	 allFiltersNotFound = false;
       }
       else {
	 if (isFirst) 
	   std::cout << SingleMuonFilterName << "not found in event " << std::endl;
       }
     }

     if (allFiltersNotFound) {
       std::cout << "Non of specified filters is found... quitting" << std::endl;
       exit(-1);
     }
     
     isFirst = false;

     // Check if it is ztt or zll (for dy sample only)
     if(isDY){
       bool isDYtoTT = false;
       int nPromptElectrons = 0;
       int nPromptMuons     = 0;
       for (unsigned int igen=0; igen<genparticles_count; ++igen) {
	 if(abs(genparticles_pdgid[igen])==11&&genparticles_fromHardProcess[igen]&&genparticles_status[igen]==1){
	   nPromptElectrons += 1;
	 }
	 
	 if(abs(genparticles_pdgid[igen])==13&&genparticles_fromHardProcess[igen]&&genparticles_status[igen]==1){
	   nPromptMuons += 1;
	 }
       }

       if(nPromptMuons != 2 && nPromptElectrons != 2) isDYtoTT = true;
       if(isZTT && !isDYtoTT) continue;
       if(!isZTT && isDYtoTT) continue;
     }

     // ********************
     // selecting good muons
     // ********************
     vector<unsigned int> muons; muons.clear();
     vector<unsigned int> looseMuons; looseMuons.clear();
     float maxMuonPt = -10;
     float muonIsolation = -10;
     unsigned int indexMu = 0;
     for(UInt_t i=0;i<muon_count;i++){
       bool muonID = muon_isMedium[i]; // MC 
       if(!muonID) continue;
       if(fabs(muon_dxy[i])>dxyMuonCut) continue;
       if(fabs(muon_dz[i])>dzMuonCut) continue;
       if(muon_pt[i]<ptMuonLooseCut) continue;
       if(fabs(muon_eta[i])>etaMuonCut) continue;
       float neutralHadIsoMu = muon_neutralHadIso[i];
       float photonIsoMu = muon_photonIso[i];
       float chargedHadIsoMu = muon_chargedHadIso[i];
       float puIsoMu = muon_puIso[i];
       float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
       neutralIsoMu = TMath::Max(float(0),neutralIsoMu);
       float absIsoMu = chargedHadIsoMu + neutralIsoMu;
       float relIsoMu = absIsoMu/muon_pt[i];
       if (relIsoMu<isoMuonCut)	{ 
	 muons.push_back(i);
       }
       if (relIsoMu<0.4) {
	 if (muon_pt[i]>maxMuonPt) {
	   maxMuonPt = muon_pt[i];
	   muonIsolation = relIsoMu;
	   indexMu = i;
	 }
       }
       if (relIsoMu<isoMuonLooseCut) looseMuons.push_back(i);
     }
  
     // *****************************************************
     // *************** Z->mumu selection *******************
     // *****************************************************
     if (muons.size()>1) {

       float maxPairPT = -1;
       unsigned int indexMu1 = 0;
       unsigned int indexMu2 = 0;
       for (unsigned int i1 = 0; i1<muons.size()-1; ++i1) {
	 unsigned int index1 = muons.at(i1);
	 bool muMatched1 = false;
	 for (unsigned int iT=0; iT<trigobject_count; ++iT) {
	   float dRtrig = deltaR(muon_eta[index1],muon_phi[index1],
				 trigobject_eta[iT],trigobject_phi[iT]);
	   if (dRtrig>DRTrigMatch) continue;
	   for (auto SingleMuonFilterName : SingleMuonFilterNames) {
	     bool isSingleMuonFilter = isSingleMuonFilterFound[SingleMuonFilterName];
	     if (!isSingleMuonFilter) continue;
	     unsigned int singleMuonFilter = singleMuonFilters[SingleMuonFilterName];
	     if (trigobject_filters[iT][singleMuonFilter]) { 
	       muMatched1 = true;
	     }
	   }
	   if (muMatched1) break;
	 }

	 for (unsigned int i2=i1+1; i2<muons.size(); ++i2) {
	   unsigned int index2 = muons.at(i2);
	   bool muMatched2 = false;
	   for (unsigned int iT=0; iT<trigobject_count; ++iT) {
	     float dRtrig = deltaR(muon_eta[index2],muon_phi[index2],
				   trigobject_eta[iT],trigobject_phi[iT]);
	     if (dRtrig>DRTrigMatch) continue;
	     for (auto SingleMuonFilterName : SingleMuonFilterNames) {
	       bool isSingleMuonFilter = isSingleMuonFilterFound[SingleMuonFilterName];
	       if (!isSingleMuonFilter) continue;
	       unsigned int singleMuonFilter = singleMuonFilters[SingleMuonFilterName];
	       if (trigobject_filters[iT][singleMuonFilter]) {
		 muMatched2 = true;
	       }
	     }
	     if (muMatched2) break;
	   }
	   
	   bool chargePair = (muon_charge[index1]*muon_charge[index2]) < 0.0;
	   bool triggerPair = (muMatched1 && muon_pt[index1]>ptMuonCut) || 
	     (muMatched2 && muon_pt[index2]>ptMuonCut);
	   float dRpair = deltaR(muon_eta[index1],muon_phi[index1],
				 muon_eta[index2],muon_phi[index2]);
	   bool deltaRPair = dRpair > dRdimuon;
	   bool pairSelected = deltaRPair && triggerPair && chargePair;

	   if (pairSelected) {
	     float sumpT = muon_pt[index1]+muon_pt[index2];
	     if (sumpT>maxPairPT) {
	       maxPairPT = sumpT;
	       if (muon_pt[index1]>muon_pt[index2]) {
		 indexMu1 = index1;
		 indexMu2 = index2;
	       }
	       else {
		 indexMu1 = index2;
		 indexMu2 = index1;
	       }
	     }
	   }

	 } // loop over second muon
       } // loop over first muon       

       if (maxPairPT>0) {
	 float mumu_weight = 1.0;
	 if (!isData) {
	   correctionWS->var("m_pt")->setVal(muon_pt[indexMu1]);
	   correctionWS->var("m_eta")->setVal(muon_eta[indexMu1]);
	   float weight_ID1 = correctionWS->function("m_idiso_ic_ratio")->getVal();
	   weight_ID1 *= correctionWS->function("m_trk_ratio")->getVal();
	   float trig_eff_Data1 = correctionWS->function("m_trg_ic_data")->getVal();
	   float trig_eff_MC1   = correctionWS->function("m_trg_ic_mc")->getVal();
	   
	   correctionWS->var("m_pt")->setVal(muon_pt[indexMu2]);
	   correctionWS->var("m_eta")->setVal(muon_eta[indexMu2]);
	   float weight_ID2 = correctionWS->function("m_idiso_ic_ratio")->getVal();
	   weight_ID2 *= correctionWS->function("m_trk_ratio")->getVal();
	   float trig_eff_Data2 = correctionWS->function("m_trg_ic_data")->getVal();
	   float trig_eff_MC2   = correctionWS->function("m_trg_ic_mc")->getVal();

	   //	   std::cout << "eta(mu1) = " << muon_eta[indexMu1]
	   //		     << "  eta(mu2) = " << muon_eta[indexMu2] << std::endl;
	   //	   float trig_eff_Data1 = SF_muonTrigger->get_EfficiencyData(muon_pt[indexMu1],muon_eta[indexMu1]);
	   //	   float trig_eff_Data2 = SF_muonTrigger->get_EfficiencyData(muon_pt[indexMu2],muon_eta[indexMu2]);
	   //	   float trig_eff_MC1 = SF_muonTrigger->get_EfficiencyMC(muon_pt[indexMu1],muon_eta[indexMu1]);
	   //	   float trig_eff_MC2 = SF_muonTrigger->get_EfficiencyMC(muon_pt[indexMu2],muon_eta[indexMu2]);
	   float trig_eff_Data = 1.0 - (1.0-trig_eff_Data1)*(1.0-trig_eff_Data2);
	   float trig_eff_MC = 1.0 - (1.0-trig_eff_MC1)*(1.0-trig_eff_MC2);
	   float trig_SF = 1.0;
	   if (trig_eff_MC>0.01) trig_SF = trig_eff_Data/trig_eff_MC;
	   mumu_weight = weight_ID1 * weight_ID2 * trig_SF;
	   //	   std::cout << "mu1 weight = " << weight_ID1 << std::endl;
	   //	   std::cout << "mu2 weight = " << weight_ID2 << std::endl;
	   //	   std::cout << "trg weight = " << trig_SF << std::endl;
	   //	   std::cout << "mumu_weight = " << mumu_weight << std::endl; 
	 }
	 TLorentzVector mu1LV; mu1LV.SetXYZM(muon_px[indexMu1],muon_py[indexMu1],muon_pz[indexMu1],MuMass);
	 TLorentzVector mu2LV; mu2LV.SetXYZM(muon_px[indexMu2],muon_py[indexMu2],muon_pz[indexMu2],MuMass);
	 float mumuMass = (mu1LV+mu2LV).M();
	 //	 std::cout << "mumuMass = " << mumuMass << std::endl;
	 if (isData) {
	   massMuMuH->Fill(mumuMass,1.);
	 }
	 else {
	   //	   std::cout << "weight = " << weight*mumu_weight << std::endl;
	   massMuMuH->Fill(mumuMass,weight*mumu_weight);
	 }
	 /*
	 float neutralHadIsoMu = muon_neutralHadIso[indexMu1];
	 float photonIsoMu = muon_photonIso[indexMu1];
	 float chargedHadIsoMu = muon_chargedHadIso[indexMu1];
	 float puIsoMu = muon_puIso[indexMu1];
	 float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	 neutralIsoMu = TMath::Max(float(0),neutralIsoMu);
	 float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	 float relIsoMu = absIsoMu/muon_pt[indexMu2];
	 std::cout << "1 : pt = " << muon_pt[indexMu1]
		   << "  eta = " << muon_eta[indexMu1]
		   << "  phi = " << muon_phi[indexMu1]
		   << "  Q = " << muon_charge[indexMu1]
		   << "  ID = " << muon_isMedium[indexMu1]
		   << "  dxy = " << muon_dxy[indexMu1]
		   << "  dz  = " << muon_dz[indexMu1]
		   << "  iso = " << relIsoMu << std::endl;

	 neutralHadIsoMu = muon_neutralHadIso[indexMu2];
	 photonIsoMu = muon_photonIso[indexMu2];
	 chargedHadIsoMu = muon_chargedHadIso[indexMu2];
	 puIsoMu = muon_puIso[indexMu2];
	 neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	 neutralIsoMu = TMath::Max(float(0),neutralIsoMu);
	 absIsoMu = chargedHadIsoMu + neutralIsoMu;
	 relIsoMu = absIsoMu/muon_pt[indexMu2];
	 std::cout << "1 : pt = " << muon_pt[indexMu2]
		   << "  eta = " << muon_eta[indexMu2]
		   << "  phi = " << muon_phi[indexMu2]
		   << "  Q = " << muon_charge[indexMu2]
		   << "  ID = " << muon_isMedium[indexMu2]
		   << "  dxy = " << muon_dxy[indexMu2]
		   << "  dz  = " << muon_dz[indexMu2]
		   << "  iso = " << relIsoMu << std::endl;
	 std::cout << "isData = " << isData << std::endl;
	 std::cout << "mass = " << mumuMass << std::endl;
	 std::cout << std::endl;
	 */
       }
     }
    
     // *******************************************
     // ********* muon+track selection ************
     // *******************************************
     if (looseMuons.size()>1) continue; // second muon veto
     if (maxMuonPt<0) continue;
     
     TLorentzVector muonLV; muonLV.SetXYZM(muon_px[indexMu],
					   muon_py[indexMu],
					   muon_pz[indexMu],
					   MuMass);

     if (muonLV.Pt()<ptMuonCut) continue;

     bool muMatchedTrigger = false;
     for (unsigned int iT=0; iT<trigobject_count; ++iT) {
       float dRtrig = deltaR(muon_eta[indexMu],muon_phi[indexMu],
			     trigobject_eta[iT],trigobject_phi[iT]);
       if (dRtrig>DRTrigMatch) continue;
       for (auto SingleMuonFilterName : SingleMuonFilterNames) {
	 bool isSingleMuonFilter = isSingleMuonFilterFound[SingleMuonFilterName];
	 if (!isSingleMuonFilter) continue;
	 unsigned int singleMuonFilter = singleMuonFilters[SingleMuonFilterName];
	 if (trigobject_filters[iT][singleMuonFilter]) {
	   muMatchedTrigger = true;
	 }
       }
       if (muMatchedTrigger) break;
     }
     if (applyTriggerMatch && !muMatchedTrigger) continue;

     // ****************************
     // selecting candidate tracks *
     // ****************************
     bool trackFound = false;
     double ptTrkMax = -1;
     TLorentzVector trackLV;
     unsigned int indexTrack = 0;
     for (unsigned int iTrk=0; iTrk<track_count; ++iTrk) {
       if (fabs(track_charge[iTrk])<0.1) continue; // make sure we are not taking neutral stuff
       if (!track_highPurity[iTrk]) continue;
       if (fabs(track_dxy[iTrk])>dxyTrkCut) continue;
       if (fabs(track_dz[iTrk])>dzTrkCut) continue;
       if (fabs(track_eta[iTrk])>etaTrkCut) continue;
       if (fabs(track_pt[iTrk])<ptTrkCut) continue;
       if (abs(track_ID[iTrk])==13) continue; // not muons

       TLorentzVector trk4; trk4.SetXYZM(track_px[iTrk],track_py[iTrk],track_pz[iTrk],track_mass[iTrk]);
       
       double dPhi = TVector2::Phi_mpi_pi(muon_phi[indexMu] - track_phi[iTrk]);
       double dEta = muon_eta[indexMu] - track_eta[iTrk];
       double dR   = sqrt(dPhi*dPhi + dEta*dEta);
       if(dR < 0.1) continue;

       bool netCharge = double(muon_charge[indexMu])*double(track_charge[iTrk]) < -0.5;

       float dphiTrkMu = dPhiFrom2P(muonLV.Px(),muonLV.Py(),trk4.Px(),trk4.Py());
       if (dphiTrkMu<dPhiMuonTrkCut) continue;

       // now check isolation
       bool isolation = true;
       for (unsigned int jTrk=0; jTrk<track_count; ++jTrk)  {
	 if (iTrk==jTrk) continue;
	 if (fabs(track_charge[jTrk])<0.1) continue; // make sure we are not taking neutral stuff
	 if (fabs(track_dxy[jTrk])>dxyTrkLooseCut) continue;
	 if (fabs(track_dz[jTrk])>dzTrkLooseCut) continue;
	 if (fabs(track_eta[jTrk])>etaTrkCut) continue;
	 if (fabs(track_pt[jTrk])<ptTrkLooseCut) continue;
	 float dRtracks = deltaR(trk4.Eta(),trk4.Phi(),track_eta[jTrk],track_phi[jTrk]);
	 if (dRtracks<dRIsoMuon) {
	   isolation = false;
	   break;
	 }
       }

       if(!isolation) continue;

       if (trk4.Pt()>ptTrkMax) {
	 trackFound = true;
	 indexTrack = iTrk;
	 ptTrkMax = trk4.Pt();
	 trackLV = trk4;
       }
     }

     if (!trackFound) continue;

     rand->SetSeed((int)((trackLV.Eta()+2.4)*100000));;
     double rannum = rand->Rndm();
     //     std::cout << "rannum = " << rannum << std::endl;
     double phi = 2.0*TMath::Pi()*rannum;

     double etaCone0p15 = trackLV.Eta() + 0.15*TMath::Cos(phi);
     double phiCone0p15 = trackLV.Phi() + 0.15*TMath::Sin(phi);

     double etaCone0p30 = trackLV.Eta() + 0.30*TMath::Cos(phi);
     double phiCone0p30 = trackLV.Phi() + 0.30*TMath::Sin(phi);

     double etaCone0p45 = trackLV.Eta() + 0.45*TMath::Cos(phi);
     double phiCone0p45 = trackLV.Phi() + 0.45*TMath::Sin(phi);

     double isolation0p15 = true;
     double isolation0p30 = true;
     double isolation0p45 = true;

     for (unsigned int jTrk=0; jTrk<track_count; ++jTrk)  {
       if (indexTrack==jTrk) continue;
       if (fabs(track_charge[jTrk])<0.1) continue; // make sure we are not taking neutral stuff
       if (fabs(track_dxy[jTrk])>dxyTrkLooseCut) continue;
       if (fabs(track_dz[jTrk])>dzTrkLooseCut) continue;
       if (fabs(track_eta[jTrk])>etaTrkCut) continue;
       if (fabs(track_pt[jTrk])<ptTrkLooseCut) continue;
       float dRtrk0p15 = deltaR(etaCone0p15,phiCone0p15,track_eta[jTrk],track_phi[jTrk]);
       float dRtrk0p30 = deltaR(etaCone0p30,phiCone0p30,track_eta[jTrk],track_phi[jTrk]);
       float dRtrk0p45 = deltaR(etaCone0p45,phiCone0p45,track_eta[jTrk],track_phi[jTrk]);
       if (dRtrk0p15<dRIsoMuon) {
	 isolation0p15 = false;
       }
       if (dRtrk0p30<dRIsoMuon) {
	 isolation0p30 = false;
       }
       if (dRtrk0p45<dRIsoMuon) {
	 isolation0p45 = false;
       }
     }


     double invMassMuTrk = (muonLV + trackLV).M();

     // jets 
     unsigned int njets = 0;
     double leadingJetPt = 0;
     TLorentzVector jetLV; 
     for (unsigned int ijet=0; ijet<pfjet_count; ++ijet) {
       float absJetEta = fabs(pfjet_eta[ijet]);
       float jetEta = pfjet_eta[ijet];
       float jetPt = pfjet_pt[ijet];
       
       if (absJetEta>jetEtaCut) continue;
       if (jetPt<jetPtCut) continue;
       bool jetId = looseId(pfjet_e[ijet],
			    pfjet_eta[ijet],
			    pfjet_energycorr[ijet],
			    pfjet_chargedhadronicenergy[ijet],
			    pfjet_neutralhadronicenergy[ijet],
			    pfjet_neutralemenergy[ijet],
			    pfjet_chargedemenergy[ijet],
			    pfjet_muonenergy[ijet],
			    pfjet_chargedmulti[ijet],
			    pfjet_neutralmulti[ijet]);
       if (!jetId) continue;

       float dRMuJet = deltaR(pfjet_eta[ijet],pfjet_phi[ijet],muonLV.Eta(),muonLV.Phi());
       if (dRMuJet<deltaRJetLeptonCut) continue;
       
       float dRTrkJet = deltaR(pfjet_eta[ijet],pfjet_phi[ijet],trackLV.Eta(),trackLV.Phi());
       if (dRTrkJet<deltaRJetLeptonCut) continue;

       if(pfjet_pt[ijet] > leadingJetPt){
	 jetLV.SetPxPyPzE(pfjet_px[ijet],pfjet_py[ijet],pfjet_pz[ijet],pfjet_e[ijet]);
	 leadingJetPt = pfjet_pt[ijet];
       }
       

       njets++;
     }
    
     TLorentzVector Met4; Met4.SetXYZM(metx,mety,0,0);

     float dPhiMuTrk = dPhiFrom2P(muonLV.Px(),muonLV.Py(),trackLV.Px(),trackLV.Py());
     float dPhiMuJet = dPhiFrom2P(muonLV.Px(),muonLV.Py(),jetLV.Px(),jetLV.Py());
     float dPhiMuMet = dPhiFrom2P(muonLV.Px(),muonLV.Py(),metx,mety);
     float dPhiTrkMet = dPhiFrom2P(trackLV.Px(),trackLV.Py(),metx,mety);
     float mT = TMath::Sqrt(2*met*muonLV.Pt()*(1-TMath::Cos(dPhiMuMet)));

     // systematics
     std::map<TString,float> mTsys;
     for (auto unc : uncNames) {
       TLorentzVector metsys = metSys[unc];
       float dphi = dPhiFrom2P(muonLV.Px(),muonLV.Py(),metsys.Px(),metsys.Py());
       mTsys[unc] = TMath::Sqrt(2*metsys.Pt()*muonLV.Pt()*(1-TMath::Cos(dphi)));
     }
     
     // compute dzeta
     float dzeta = 0;
     float pzetamiss = 0;

     // bisector of electron and muon transverse momenta
     float trackUnitX = trackLV.Px()/trackLV.Pt();
     float trackUnitY = trackLV.Py()/trackLV.Pt();
     float muonUnitX = muonLV.Px()/muonLV.Pt();
     float muonUnitY = muonLV.Py()/muonLV.Pt();
     float zetaX = trackUnitX + muonUnitX;
     float zetaY = trackUnitY + muonUnitY;
     float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);
     zetaX = zetaX/normZeta;
     zetaY = zetaY/normZeta;
     
     float vectorVisX = muonLV.Px() + trackLV.Px();
     float vectorVisY = muonLV.Py() + trackLV.Py();
     float pzetavis = vectorVisX*zetaX + vectorVisY*zetaY;
     
     computeDzeta(metx,mety,zetaX,zetaY,pzetavis,pzetamiss,dzeta);

     map<TString,float> dzetaSys;
     for (auto unc : uncNames) {
       float dzeta_x = 0;
       float pzetamiss_x = 0;
       computeDzeta(metx,mety,zetaX,zetaY,pzetavis,pzetamiss_x,dzeta_x);
       dzetaSys[unc] = dzeta_x;
     }

     // muon id/iso/trigger weights here
     float muTrigWeight = 1;
     float muIdWeight   = 1;
     if (!isData) { 

       correctionWS->var("m_pt")->setVal(muonLV.Pt());
       correctionWS->var("m_eta")->setVal(muonLV.Eta());
       
       muIdWeight = correctionWS->function("m_idiso_ic_ratio")->getVal();
       muIdWeight *= correctionWS->function("m_trk_ratio")->getVal();

       //       std::cout << "eta(muon) = " << muonLV.Eta() << std::endl;
       //       float trigWeightData = SF_muonTrigger->get_EfficiencyData(muonLV.Pt(),muonLV.Eta());
       //       float trigWeightMC = SF_muonTrigger->get_EfficiencyMC(muonLV.Pt(),muonLV.Eta());

       if (applyTriggerMatch) 
	 muTrigWeight = correctionWS->function("m_trg_ic_ratio")->getVal();
       else
	 muTrigWeight = correctionWS->function("m_trg_ic_data")->getVal();

       /*
	 std::cout << "pT(muon) = " << muonLV.Pt()
	 << "   eta(muon) = " << muonLV.Eta() << std::endl;
	 std::cout << "Trigger weight = " << muTrigWeight << std::endl;
	 std::cout << "Id/Iso  weight = " << muIdWeight << std::endl;
       */

     }

     muTrigWeightH->Fill(muTrigWeight,1.0);
     weight *= muTrigWeight;
     muIdWeightH->Fill(muIdWeight,1.0);
     weight *= muIdWeight;

     os_tree = double(muon_charge[indexMu])*double(track_charge[indexTrack]) < 0.5;

     mt_tree = mT;
     mtJetEnUp_tree = mTsys["JetEnUp"];
     mtJetEnDown_tree = mTsys["JetEnDown"];
     mtUnclEnUp_tree = mTsys["UnclEnUp"];
     mtUnclEnDown_tree = mTsys["UnclEnDown"];

     met_tree = met;
     metJetEnUp_tree = metSys["JetEnUp"].Pt();
     metJetEnDown_tree = metSys["JetEnDown"].Pt();
     metUnclEnUp_tree = metSys["UnclEnUp"].Pt();
     metUnclEnDown_tree = metSys["UnclEnDown"].Pt();

     dzeta_tree = dzeta;
     dzetaJetEnUp_tree = dzetaSys["JetEnUp"];
     dzetaJetEnDown_tree = dzetaSys["JetEnDown"];
     dzetaUnclEnUp_tree = dzetaSys["UnclEnUp"];
     dzetaUnclEnDown_tree = dzetaSys["UnclEnDown"];

     njets_tree = njets;
     dphimutrk_tree = dPhiMuTrk;
     massmutrk_tree = invMassMuTrk;
     trkpt_tree = trackLV.Pt();
     trketa_tree = trackLV.Eta();
     mupt_tree = muonLV.Pt();
     mueta_tree = muonLV.Eta();
     muiso_tree = muonIsolation;
     weight_tree = weight;
     muweight_tree = muTrigWeight*muIdWeight;
     genweight_tree = generator_weight;
     puweight_tree = puweight;
     zptweight_tree = zptmassweight;
     topptweight_tree = topptweight;
     npart_tree = genparticles_noutgoing;

     cone0p15_tree = isolation0p15;
     cone0p30_tree = isolation0p30;
     cone0p45_tree = isolation0p45;

     //     std::cout << std::endl;
     // type of track : 0 - jet, 11 - muon, 13 - electron, 15 - tau
     trktyp_tree = 0;
     if (!isData) {
       float minDR = 0.3;
       for (unsigned int igen=0; igen<genparticles_count; ++igen) {
	 //	 std::cout << genparticles_pdgid[igen] << std::endl;
	 int pdgid = TMath::Abs(genparticles_pdgid[igen]);
	 TLorentzVector genLV; genLV.SetXYZT(genparticles_px[igen],
					     genparticles_py[igen],
					     genparticles_pz[igen],
					     genparticles_e[igen]); 
	 if (pdgid==11||pdgid==13) {
	   if (genparticles_fromHardProcess[igen] && genparticles_status[igen]==1) {
	     float dR = deltaR(trackLV.Eta(),trackLV.Phi(),
			       genLV.Eta(),genLV.Phi());
	     if (dR<minDR) {
	       minDR = dR;
	       trktyp_tree = pdgid;
	     }
	   }
	 }
       }
       
       for (unsigned int igen=0; igen<gentau_count; ++igen) {
	 //	 std::cout << "gentau : " << igen 
	 //	   << "   isPrompt = " << gentau_isPrompt[igen] << std::endl;
	 TLorentzVector genLV; genLV.SetXYZT(gentau_visible_px[igen],
					     gentau_visible_py[igen],
					     gentau_visible_pz[igen],
					     gentau_visible_e[igen]); 
	 if (gentau_isPrompt[igen]==1) {
	   float dR = deltaR(trackLV.Eta(),trackLV.Phi(),
			     genLV.Eta(),genLV.Phi());

	   //	   std::cout << "dR : " << dR << std::endl;
	   if (dR<minDR) {
	     minDR = dR;
	     trktyp_tree = 15;
	   }
	 }
       }
     }
     
     mutrkTree->Fill();

     bool lowMT = mT < mTLowCut;
     bool highMT = mT > mTHighCut;
     bool muIso = muiso_tree < 0.15;

     // signal region
     if (muIso&&os_tree) {

       mTH->Fill(mT,weight);
       
       if(trackLV.Pt()<5.){
	 mTH_2p5to5->Fill(mT,weight);
       }
       else if(trackLV.Pt()<10.){
	 mTH_5to10->Fill(mT,weight);
       }
       else if(trackLV.Pt()<15.){
	 mTH_10to15->Fill(mT,weight);
       }
       else if(trackLV.Pt()<20.){
	 mTH_15to20->Fill(mT,weight);
       }
       else if(trackLV.Pt()>20.){
	 mTH_20toInf->Fill(mT,weight);
       }

       // low mT region
       if (lowMT) {

	 ptMuonH->Fill(muonLV.Pt(),weight);
	 etaMuonH->Fill(muonLV.Eta(),weight);
	 ptTrackH->Fill(trackLV.Pt(),weight);
	 ptTrackLowH->Fill(trackLV.Pt(),weight);
	 etaTrackH->Fill(trackLV.Eta(),weight);
	 metH->Fill(met,weight);
	 dzetaH->Fill(dzeta,weight);
	 dPhiMuonTrkH->Fill(dPhiMuTrk,weight);
	 dPhiMuonMetH->Fill(dPhiMuMet,weight);
	 if(njets>0) dPhiMuonJetH->Fill(dPhiMuJet,weight);
	 dPhiTrkMetH->Fill(dPhiTrkMet,weight);
	 nJetH->Fill(njets,weight);
	 dzTrackH->Fill(track_dz[indexTrack],weight);
	 dxyTrackH->Fill(track_dxy[indexTrack],weight);
	 
	 if(trackLV.Pt()<5.){
	   ptMuonH_2p5to5->Fill(muonLV.Pt(),weight);
	   metH_2p5to5->Fill(met,weight);
	   dzetaH_2p5to5->Fill(dzeta,weight);
	   dPhiMuonTrkH_2p5to5->Fill(dPhiMuTrk,weight);
	   dPhiMuonMetH_2p5to5->Fill(dPhiMuMet,weight);
	   if(njets>0) dPhiMuonJetH_2p5to5->Fill(dPhiMuJet,weight);
	   dPhiTrkMetH_2p5to5->Fill(dPhiTrkMet,weight);
	 }
	 else if(trackLV.Pt()<10.){	 
	   ptMuonH_5to10->Fill(muonLV.Pt(),weight);
	   metH_5to10->Fill(met,weight);
	   dzetaH_5to10->Fill(dzeta,weight);
	   dPhiMuonTrkH_5to10->Fill(dPhiMuTrk,weight);
	   dPhiMuonMetH_5to10->Fill(dPhiMuMet,weight);
	   if(njets>0) dPhiMuonJetH_5to10->Fill(dPhiMuJet,weight);
	   dPhiTrkMetH_5to10->Fill(dPhiTrkMet,weight);
	 }
	 else if(trackLV.Pt()<15.){
	   ptMuonH_10to15->Fill(muonLV.Pt(),weight);
	   metH_10to15->Fill(met,weight);
	   dzetaH_10to15->Fill(dzeta,weight);
	   dPhiMuonTrkH_10to15->Fill(dPhiMuTrk,weight);
	   dPhiMuonMetH_10to15->Fill(dPhiMuMet,weight);
	   if(njets>0) dPhiMuonJetH_10to15->Fill(dPhiMuJet,weight);
	   dPhiTrkMetH_10to15->Fill(dPhiTrkMet,weight);
	 }
	 else if(trackLV.Pt()<20.){
	   ptMuonH_15to20->Fill(muonLV.Pt(),weight);
	   metH_15to20->Fill(met,weight);
	   dzetaH_15to20->Fill(dzeta,weight);
	   dPhiMuonTrkH_15to20->Fill(dPhiMuTrk,weight);
	   dPhiMuonMetH_15to20->Fill(dPhiMuMet,weight);
	   if(njets>0) dPhiMuonJetH_15to20->Fill(dPhiMuJet,weight);
	   dPhiTrkMetH_15to20->Fill(dPhiTrkMet,weight);
	 }
	 else if(trackLV.Pt()>20.){
	   ptMuonH_20toInf->Fill(muonLV.Pt(),weight);
	   metH_20toInf->Fill(met,weight);
	   dzetaH_20toInf->Fill(dzeta,weight);
	   dPhiMuonTrkH_20toInf->Fill(dPhiMuTrk,weight);
	   dPhiMuonMetH_20toInf->Fill(dPhiMuMet,weight);
	   if(njets>0) dPhiMuonJetH_20toInf->Fill(dPhiMuJet,weight);
	   dPhiTrkMetH_20toInf->Fill(dPhiTrkMet,weight);
	 }
       }
     }

     TString mt_region("");
     if (lowMT) mt_region = "lowMT";
     else if (highMT) mt_region = "highMT";
     else mt_region = "";

     TString muon_isolation("");
     if (muiso_tree<0.15) muon_isolation = "iso";
     else if (muiso_tree>0.2) muon_isolation = "antiiso";

     TString pt_bin("");
     if (trackLV.Pt()<5.) pt_bin = "2p5to5";
     else if (trackLV.Pt()<10.) pt_bin = "5to10";
     else if (trackLV.Pt()<15.) pt_bin = "10to15";
     else if (trackLV.Pt()<20.) pt_bin = "15to20";
     else if (trackLV.Pt()>20.) pt_bin = "20toInf";

     TString sign_region("");
     if (os_tree) sign_region = "OS";
     else sign_region = "SS";

     if ((mt_region!="")&&
	 (muon_isolation!="")&&
	 (sign_region!="")&&
	 (pt_bin!="")) {
       TString histName = baseName + "_" + mt_region + "_" + muon_isolation + "_" + sign_region;
       mass_histos[histName]->Fill(invMassMuTrk,weight);
       histName = baseName + "_" + pt_bin + "_" + mt_region + "_" + muon_isolation + "_" + sign_region;
       mass_histos[histName]->Fill(invMassMuTrk,weight);
     }

   } // loop over tree entries    
   
   delete tree_;
   file_->Close();
   delete file_;
   
  } // filelist loop
  
  file->cd("");
  file->Write();
  file->Close();
  
  //delete file;
  
}// int main  

 
