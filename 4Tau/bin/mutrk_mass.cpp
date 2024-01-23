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

#include "HtoAA/Utilities/interface/Jets.h"
#include "HtoAA/Utilities/interface/QCDModelDefinitions.h"
#include "HtoAA/Utilities/interface/QCDModel.h"
#include "HtoAA/Utilities/interface/functions.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "CondFormats/BTauObjects/interface/BTagEntry.h"

using namespace std;

int main(int argc, char * argv[]) {

  if (argc<2) {
    std::cout << "Usage of the program : mutrk_mass [config_file] [file_list]" << std::endl;
    exit(1);
  }

  // **** configuration
  Config cfg(argv[1]);

  const int  debug  = cfg.get<int>("Debug");
  const int  era = cfg.get<int>("Era");

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

  // trigger
  const string dimuonTriggerName = cfg.get<string>("DiMuonTriggerName");
  const string muonHighPtFilterName = cfg.get<string>("MuonHighPtFilterName");
  const string muonLowPtFilterName1 = cfg.get<string>("MuonLowPtFilterName1");
  const string muonLowPtFilterName2 = cfg.get<string>("MuonLowPtFilterName2");
  const string dimuonDzFilterName = cfg.get<string>("DimuonDzFilterName");
  const string dimuonSameSignFilterName = cfg.get<string>("DimuonSameSignFilterName");

  // trigger info
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

  // 
  const bool ApplyBTagVeto = cfg.get<bool>("ApplyBTagVeto");
  const string bTagAlgorithm = cfg.get<string>("BTagAlgorithm");
  const string bTagDiscriminator1 = cfg.get<string>("BTagDiscriminator1");
  const string bTagDiscriminator2 = cfg.get<string>("BTagDiscriminator2");
  const string bTagDiscriminator3 = cfg.get<string>("BTagDiscriminator3");
  const float btagCut = cfg.get<float>("BTagCut");
  const float bjetEta = cfg.get<float>("BJetEta");
  const float bjetPt = cfg.get<float>("BJetPt");
  
  // rebinning of pdf
  vector<double> bins = cfg.get<vector<double>>("bins");

  TString BTagAlgorithm(bTagAlgorithm);
  TString BTagDiscriminator1(bTagDiscriminator1);
  TString BTagDiscriminator2(bTagDiscriminator2);
  TString BTagDiscriminator3(bTagDiscriminator3);

  const bool applyBTagSF = cfg.get<bool>("ApplyBTagSF");

  // BTag SF file
  const string BtagSfFile = cfg.get<string>("BtagSfFile");
  BTagCalibration calib = BTagCalibration(bTagAlgorithm, BtagSfFile);
  BTagCalibrationReader reader_B = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central",{"up","down"});
  BTagCalibrationReader reader_C = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central",{"up","down"});
  BTagCalibrationReader reader_Light = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central",{"up","down"});
  reader_B.load(calib, BTagEntry::FLAV_B, "comb");
  reader_C.load(calib, BTagEntry::FLAV_C, "comb");
  reader_Light.load(calib, BTagEntry::FLAV_UDSG, "incl");
  
  // BTAG efficiency for various flavours ->
  TString fileBtagEff = (TString)cfg.get<string>("BtagMCeffFile");
  TFile *fileTagging  = new TFile(fileBtagEff);
  TH2F  *tagEff_B     = (TH2F*)fileTagging->Get("btag_eff_b");
  TH2F  *tagEff_C     = (TH2F*)fileTagging->Get("btag_eff_c");
  TH2F  *tagEff_Light = (TH2F*)fileTagging->Get("btag_eff_oth");
  TRandom3 *rand = new TRandom3();

  float MaxBJetPt = 1000.;
  float MinBJetPt = 20.;
  float MaxBJetEta = 2.4;
  float MinBJetEta = 0.0;

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
  float pfjet_neutralhadronicenergy[200];
  float pfjet_chargedhadronicenergy[200];
  float pfjet_neutralemenergy[200];
  float pfjet_chargedemenergy[200];
  float pfjet_muonenergy[200];
  float pfjet_chargedmuonenergy[200];
  UInt_t pfjet_chargedmulti[200];
  UInt_t pfjet_neutralmulti[200];
  UInt_t pfjet_chargedhadronmulti[200];
  Int_t  pfjet_flavour[200];
  float pfjet_btag[200][10];
  float pfjet_jecUncertainty[200];
  
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

  std::string rootFileName(argv[2]);
  
  TString chainName("makeroottree/AC1B");
  TString initChainName("initroottree/AC1B");
  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;
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

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++ Definitions are in HtoAA/Utilities/interface/functions.h ++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TH1D * PartonMomBinsH = new TH1D("PartonMomBinsH","",nPartMom,partonMomBins);
  TH1D * MuonMomBinsH = new TH1D("MuonMomBinsH","",nMuMom,muonMomBins);
  TH1D * MuonPartonNetChargeH = new TH1D("MuonPartonNetChargeH","",nNetQ,0.,float(nNetQ));
  TH1D * PartonFlavorH = new TH1D("PartonFlavorH","",nFlav,0.,float(nFlav));
  TH1D * RegionsH = new TH1D("Regions","",nReg,0.,float(nReg));

  for (unsigned int iB=0; iB<nMuMom; ++iB)
    MuonMomBinsH->GetXaxis()->SetBinLabel(iB+1,muonMomRange[iB]);
  for (unsigned int iB=0; iB<nPartMom; ++iB)
    PartonMomBinsH->GetXaxis()->SetBinLabel(iB+1,partonMomRange[iB]);
  for (unsigned int iB=0; iB<nFlav; ++iB)
    PartonFlavorH->GetXaxis()->SetBinLabel(iB+1,partonFlavor[iB]);
  for (unsigned int iB=0; iB<nNetQ; ++iB)
    MuonPartonNetChargeH->GetXaxis()->SetBinLabel(iB+1,muonPartonNetCharge[iB]);
  for (unsigned int iB=0; iB<nReg; ++iB)
    RegionsH->GetXaxis()->SetBinLabel(iB+1,Regions[iB]);


  double binWidth = (xMax-xMin)/double(nBins);
  unsigned int nbins = nBins;
  double xbins[200]; 
  if (bins.size()>1) {
    nbins = bins.size() - 1;
    for (unsigned int ib=0; ib<bins.size(); ++ib) {
      xbins[ib] = bins[ib];
    }
  }
  else {
    nbins = nBins;    
    for (unsigned int ib=0; ib<=nbins; ++ib) {
      xbins[ib] = xMin + double(ib)*binWidth; 
    }
  }



  TH1D * massHist = new TH1D("massHist","",nbins,xbins);
  std::cout << std::endl;
  std::cout << "Number of bins = " << nbins << std::endl;
  std::cout << "bin boundaries : " << std::endl;
  for (unsigned int ib=1; ib<=nbins; ++ib) {
    double lower = massHist->GetXaxis()->GetBinLowEdge(ib);
    double upper = massHist->GetXaxis()->GetBinLowEdge(ib+1);
    printf("[%4.1f,%4.1f]\n",lower,upper);
  }
  std::cout << std::endl;

  // testing model in inclusive muon sample
  // single muon passing selection
  TH1D * InvMassIsoH = new TH1D("InvMassIsoH","",nbins,xbins);
  TH1D * ModelInvMassIsoH = new TH1D("ModelInvMassIsoH","",nbins,xbins);

  TH1D * InvMassLooseIsoH = new TH1D("InvMassLooseIsoH","",nbins,xbins);
  TH1D * ModelInvMassLooseIsoH = new TH1D("ModelInvMassLooseIsoH","",nbins,xbins);

  // testing model in sample of SS muons
  TH1D * InvMassDimuonIsoH = new TH1D("InvMassDimuonIsoH","",nbins,xbins);
  TH1D * HybridModelInvMassDimuonIsoH = new TH1D("HybridModelInvMassDimuonIsoH","",nbins,xbins);
  TH1D * InclusiveModelInvMassDimuonIsoH = new TH1D("InclusiveModelInvMassDimuonIsoH","",nbins,xbins);
  TH1D * ClosureInvMassDimuonIsoH = new TH1D("ClosureInvMassDimuonIsoH","",nbins,xbins);

  TH1D * InvMassDimuonLooseIsoH = new TH1D("InvMassDimuonLooseIsoH","",nbins,xbins);
  TH1D * HybridModelInvMassDimuonLooseIsoH = new TH1D("HybridModelInvMassDimuonLooseIsoH","",nbins,xbins);
  TH1D * InclusiveModelInvMassDimuonLooseIsoH = new TH1D("InclusiveModelInvMassDimuonLooseIsoH","",nbins,xbins);
  TH1D * ClosureInvMassDimuonLooseIsoH = new TH1D("ClosureInvMassDimuonLooseIsoH","",nbins,xbins);


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

  // Signal region 
  TH1D * InvMassH = new TH1D("InvMassH","",nbins,xbins);
  TH2D * InvMass2DH = new TH2D("InvMass2DH","",nbins,xbins,nbins,xbins);

  TH1D * ClosureInvMassH = new TH1D("ClosureInvMassH","",nbins,xbins);
  TH2D * ClosureInvMass2DH = new TH2D("ClosureInvMass2DH","",nbins,xbins,nbins,xbins);

  TH1D * HybridModelInvMassH = new TH1D("HybridModelInvMassH","",nbins,xbins);
  TH2D * HybridModelInvMass2DH = new TH2D("HybridModelInvMass2DH","",nbins,xbins,nbins,xbins);

  // Background sideband (Y)   
  TH1D * InvMass_ControlYH = new TH1D("InvMass_ControlYH","",nbins,xbins); 
  TH2D * InvMass2D_ControlYH = new TH2D("InvMass2D_ControlYH","",nbins,xbins,nbins,xbins);

  TH1D * ClosureInvMass_ControlYH = new TH1D("ClosureInvMass_ControlYH","",nbins,xbins); 
  TH2D * ClosureInvMass2D_ControlYH = new TH2D("ClosureInvMass2D_ControlYH","",nbins,xbins,nbins,xbins);

  TH1D * HybridModelInvMass_ControlYH = new TH1D("HybridModelInvMass_ControlYH","",nbins,xbins); 
  TH2D * HybridModelInvMass2D_ControlYH = new TH2D("HybridModelInvMass2D_ControlYH","",nbins,xbins,nbins,xbins);

  // Background sideband (X)
  TH1D * InvMass_ControlXH = new TH1D("InvMass_ControlXH","",nbins,xbins); 
  TH2D * InvMass2D_ControlXH = new TH2D("InvMass2D_ControlXH","",nbins,xbins,nbins,xbins);

  TH1D * ClosureInvMass_ControlXH = new TH1D("ClosureInvMass_ControlXH","",nbins,xbins); 
  TH2D * ClosureInvMass2D_ControlXH = new TH2D("ClosureInvMass2D_ControlXH","",nbins,xbins,nbins,xbins);

  TH1D * HybridModelInvMass_ControlXH = new TH1D("HybridModelInvMass_ControlXH","",nbins,xbins); 
  TH2D * HybridModelInvMass2D_ControlXH = new TH2D("HybridModelInvMass2D_ControlXH","",nbins,xbins,nbins,xbins);



  // for muons selected with same-sign requirement
  // probability to pass selection as a function of
  // flavour, net charged (mu,parton), parton P, muon pT
  TH1D * partonMu_SS = new TH1D("partonMu_SS","",1,0.,1.);  
  TH1D * partonMuSelectedIso_SS = new TH1D("partonMuSelectedIso_SS","",1,0.,1.); 
  TH1D * partonMuModelledIso_SS = new TH1D("partonMuModelledIso_SS","",1,0.,1.);
  TH1D * partonMuSelectedLooseIso_SS = new TH1D("partonMuSelectedLooseIso_SS","",1,0.,1.); 
  TH1D * partonMuModelledLooseIso_SS = new TH1D("partonMuModelledLooseIso_SS","",1,0.,1.);


  // inclusive muon region
  TH1D * partonMu = new TH1D("partonMu","",1,0.,1.); 
  TH1D * partonMuSelectedIso = new TH1D("partonMuSelectedIso","",1,0.,1.); 
  TH1D * partonMuModelledIso = new TH1D("partonMuModelledIso","",1,0.,1.);
  TH1D * partonMuSelectedLooseIso = new TH1D("partonMuSelectedLooseIso","",1,0.,1.); 
  TH1D * partonMuModelledLooseIso = new TH1D("partonMuModelledLooseIso","",1,0.,1.);

  TH1D * partonMuProbe_SS[nFlav][nNetQ][nPartMom][nMuMom];
  TH1D * partonMuPass_SS[nFlav][nNetQ][nPartMom][nMuMom][nReg]; 
  TH1D * MassMuTrk_SS[nFlav][nNetQ][nPartMom][nMuMom][nReg];

  TH1D * partonMuProbe[nFlav][nNetQ][nPartMom][nMuMom];
  TH1D * partonMuPass[nFlav][nNetQ][nPartMom][nMuMom][nReg]; 
  TH1D * MassMuTrk[nFlav][nNetQ][nPartMom][nMuMom][nReg];
 
  for (unsigned int iF=0; iF<nFlav; ++iF) {
    for (unsigned int iQ=0; iQ<nNetQ; ++iQ) {
      for (unsigned int iMom=0; iMom<nPartMom; ++iMom) {
	for (unsigned int mu=0; mu<nMuMom; ++mu) {
	    TString name = 
	      partonFlavor[iF] + "_" + 
	      muonPartonNetCharge[iQ] + "_" +
	      partonMomRange[iMom] + "_" +
	      muonMomRange[mu];

	    TString histName = "partonMuProbe_" + name + "_SS";
	    partonMuProbe_SS[iF][iQ][iMom][mu] = new TH1D(histName,"",1,0.,1.);

	    histName = "partonMuProbe_" + name;
	    partonMuProbe[iF][iQ][iMom][mu] = new TH1D(histName,"",1,0.,1.);


	  for (unsigned int iR=0; iR<nReg; ++iR) {
	    name = 
	      partonFlavor[iF] + "_" + 
	      muonPartonNetCharge[iQ] + "_" +
	      partonMomRange[iMom] + "_" +
	      muonMomRange[mu] + "_" +
	      Regions[iR] ;
	  
	    // mass distributions
	    histName = "MuTrkMass_" + name + "_SS";
	    MassMuTrk_SS[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",nBins,xMin,xMax);
      
	    histName = "MuTrkMass_" + name;
	    MassMuTrk[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",nBins,xMin,xMax);
	    
	    // probability to pass selection
	    histName = "partonMuPass_" + name + "_SS";
	    partonMuPass_SS[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",1,0.,1.);
	    
	    histName = "partonMuPass_" + name;
	    partonMuPass[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",1,0.,1.);
	  }
	}
      }
    }
  }

  string cmsswBase = (getenv ("CMSSW_BASE"));

  // Run-lumi selector
  std::string qcdFileName = cfg.get<string>("qcdModelFileName");
  TString QCDFileName(qcdFileName);

  // QCD Model
  TString fileNameQCDModel = TString(cmsswBase)+TString("/src/HtoAA/data/")+QCDFileName;
  QCDModel * qcdModel = NULL;
  bool applyQCDModel = QCDFileName.Contains(".root");
  if (applyQCDModel) 
    qcdModel = new QCDModel(fileNameQCDModel,bins);

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

   TTree * _inittree = (TTree*)file_->Get(initChainName);
   if (_inittree!=NULL) {
     Float_t Genweight;
     _inittree->SetBranchAddress("genweight",&Genweight);
     Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
     std::cout << "Number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
     for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
       _inittree->GetEntry(iEntry);
       if (Genweight>0.0)
	 histWeightsH->Fill(1.,1.);
       else
	 histWeightsH->Fill(1.,-1.);
     }
   }  

   TTree * tree_ = (TTree*)file_->Get(chainName);
   
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
   tree_->SetBranchAddress("pfjet_neutralhadronicenergy",pfjet_neutralhadronicenergy);
   tree_->SetBranchAddress("pfjet_chargedhadronicenergy",pfjet_chargedhadronicenergy);
   tree_->SetBranchAddress("pfjet_neutralemenergy",pfjet_neutralemenergy);
   tree_->SetBranchAddress("pfjet_chargedemenergy",pfjet_chargedemenergy);
   tree_->SetBranchAddress("pfjet_muonenergy",pfjet_muonenergy);
   tree_->SetBranchAddress("pfjet_chargedmuonenergy",pfjet_chargedmuonenergy);
   tree_->SetBranchAddress("pfjet_chargedmulti",pfjet_chargedmulti);
   tree_->SetBranchAddress("pfjet_neutralmulti",pfjet_neutralmulti);
   tree_->SetBranchAddress("pfjet_btag",pfjet_btag);
   tree_->SetBranchAddress("pfjet_jecUncertainty",pfjet_jecUncertainty);


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

   // Additional objects
   tree_->SetBranchAddress("run_hltfilters",&hltfilters);
   tree_->SetBranchAddress("run_btagdiscriminators", &btagdiscriminators);
   tree_->SetBranchAddress("hltriggerresults",&hltriggerresults);
   tree_->SetBranchAddress("hltriggerprescales",&hltriggerprescales);

   tree_->SetBranchAddress("numtruepileupinteractions",&numtruepileupinteractions);

   tree_->SetBranchAddress("genweight",&genweight);
   tree_->SetBranchAddress("genparticles_count", &genparticles_count);
   tree_->SetBranchAddress("genparticles_e", genparticles_e);
   tree_->SetBranchAddress("genparticles_px", genparticles_px);
   tree_->SetBranchAddress("genparticles_py", genparticles_py);
   tree_->SetBranchAddress("genparticles_pz", genparticles_pz);
   tree_->SetBranchAddress("genparticles_pdgid", genparticles_pdgid);
   tree_->SetBranchAddress("genparticles_status", genparticles_status);
   tree_->SetBranchAddress("genparticles_info", genparticles_info);

   int numberOfCandidates = tree_->GetEntries();

   std::cout << "number of events = " << numberOfCandidates << std::endl;
   
   for (int iCand=0; iCand<numberOfCandidates; iCand++) {
     
     tree_->GetEntry(iCand);

     events++;
     if (events%10000==0) cout << "   processed events : " << events << endl;
     
     float weight = 1.;
     if (genweight<0.)
       weight = -1.;
          
     float puweight = float(PUofficial->get_PUweight(double(numtruepileupinteractions)));
     puWeightH->Fill(puweight,1.0);
     weight *= puweight;
     
     if (debug>2) {
       std::cout << std::endl;
       std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
       std::cout << std::endl;
       std::cout << "run = " << event_run << "  number = " << event_nr << std::endl;
       std::cout << std::endl;
       std::cout << "B-tagging " << std::endl;
     }
     
     if (ApplyBTagVeto) {

       int nBTagDiscriminant1 = -1;
       int nBTagDiscriminant2 = -1;
       int nBTagDiscriminant3 = -1;

       unsigned int num_btags = 0;
       unsigned int nbtags = 0;
       float Pdata = 1.;
       float Pmc = 1.;
       
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
       
       if (num_btags==0) {
	 std::cout << "No discriminants are found for " << BTagAlgorithm << std::endl;
	 exit(-1);
       }
       if (BTagAlgorithm=="pfDeepCSVJetTags" && num_btags!=2) {
	 std::cout << "Numbers of discriminators for " << BTagAlgorithm << " = " << num_btags 
		   << "   should be 2 " << std::endl;
	 exit(-1);
       }
       if (BTagAlgorithm=="pfDeepFlavourJetTags" && num_btags!=3) {
	 std::cout << "Numbers of discriminators for " << BTagAlgorithm << " = " << num_btags 
		   << "   should be 3" << std::endl;
	 exit(-1);
       }
       
       for (unsigned int jet=0; jet<pfjet_count; ++jet) {
	   
	 float absEta = TMath::Abs(pfjet_eta[jet]);
	 float JetPtForBTag = pfjet_pt[jet];
	 float JetEtaForBTag = absEta;
	 float jetEta = pfjet_eta[jet];
	 
	 if (JetEtaForBTag>bjetEta) continue;
	 if (JetPtForBTag<bjetPt) continue;

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
	 if (tagged) nbtags++;
	 
	 if (debug>2) 
	   printf("jet pT = %6.2f  eta = %6.2f  discr = %5.3f  : tagged = %1i\n",
		  pfjet_pt[jet],pfjet_eta[jet],btagDiscr,tagged);

	 // BTag scale factors
	 if (applyBTagSF) {
	   
	   int flavor = TMath::Abs(pfjet_flavour[jet]);
	   float tageff = tagEff_Light->GetBinContent(tagEff_Light->FindBin(JetPtForBTag, JetEtaForBTag));
	   float jet_scalefactor = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG,JetEtaForBTag,JetPtForBTag);
	   if (flavor==4) { 
	     tageff = tagEff_C->GetBinContent(tagEff_C->FindBin(JetPtForBTag,JetEtaForBTag));
	     jet_scalefactor = reader_C.eval_auto_bounds("central",BTagEntry::FLAV_C,JetEtaForBTag,JetPtForBTag);
	   }
	   if (flavor==5) { 
	     tageff = tagEff_B->GetBinContent(tagEff_B->FindBin(JetPtForBTag,JetEtaForBTag));
	     jet_scalefactor = reader_B.eval_auto_bounds("central",BTagEntry::FLAV_B,JetEtaForBTag,JetPtForBTag);
	   }

	   if (tagged) {  
	     Pmc = Pmc*tageff;
	     Pdata = Pdata*jet_scalefactor*tageff;	     
	   }
	   else {
	     Pmc = Pmc*(1-tageff);
	     Pdata = Pdata*(1.0-jet_scalefactor*tageff);
	   }
	 }	 
       }
       if (debug>2) {
	 std::cout << std::endl;
	 std::cout << "Number of tagged jets = " << nbtags << std::endl;
       }
       
       if (nbtags>0) continue;
       float weight_btag = Pdata/Pmc;
       weight *= weight_btag;

     }
     // 
     // end of btag veto
     //

     std::vector<int> partonPdgId; partonPdgId.clear();
     std::vector<TLorentzVector> partonLV; partonLV.clear();
     unsigned int partons = 0;
     unsigned int pfjets = 0;

     if (debug>2) {
       std::cout << std::endl;
       std::cout << "jets -> " << std::endl;
     }

     for (unsigned ijet=0; ijet<pfjet_count; ++ijet) {
       if (fabs(pfjet_eta[ijet])<5.0) {
	 pfjets++;
	 int absPdgId = TMath::Abs(pfjet_flavour[ijet]);
	 if (pfjet_flavour[ijet]!=0) 
	   partons++;

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
	 if (debug>2)
	   printf("%2i  flavor = %3i   pT = %7.2f   eta = %5.2f\n",
		  pfjets,pfjet_flavour[ijet],pfjet_pt[ijet],pfjet_eta[ijet]);
       } 
     }
    
     // filling histograms 
     PartonMultiplicityH->Fill(float(partons),weight);
     PFJetMultiplicityH->Fill(float(pfjets),weight);

     // ********************
     // selecting good muons
     // ********************
     vector<unsigned int> muons; muons.clear();
     vector<unsigned int> muons_flavour; muons_flavour.clear();
     vector<int> muons_pdgid; muons_pdgid.clear();
     vector<unsigned int> muons_partmom; muons_partmom.clear();
     vector<unsigned int> muons_muonmom; muons_muonmom.clear();
     vector<unsigned int> muons_net; muons_net.clear();
     vector<int> muons_region; muons_region.clear();
     vector<float> muons_mutrkmass; muons_mutrkmass.clear();

     for(UInt_t i=0;i<muon_count;i++){
       bool muonID = muon_isMedium[i]; // MC 
       if (!muonID) continue;
       if(fabs(muon_dxy[i])>dxyMuonCut) continue;
       if(fabs(muon_dz[i])>dzMuonCut) continue;
       if(muon_pt[i]<ptMuonLowCut) continue;
       if(fabs(muon_eta[i])>etaMuonLowCut) continue;
       muons.push_back(i);
     }
     
     nGoodMuonsH->Fill(float(muons.size()),weight);
      
     if (muons.size()<1) continue; // quit event if number of good muons < 1

     int nIsoMuons = 0;
     int nLooseIsoMuons = 0;
     int nSbMuons = 0;

     if (debug>1) {
       std::cout << std::endl;
       std::cout << "loop over muons -> " << std::endl;
     }

     for (unsigned int imu=0; imu<muons.size(); ++imu) {

       unsigned int index = muons.at(imu);


       // Muon
       TLorentzVector Muon4; Muon4.SetXYZM(muon_px[index],
					   muon_py[index],
					   muon_pz[index],
					   MuMass);
       // determine flavour of jet
       float dRmin = 0.5;
       unsigned int flavour = 0;
       float qnet  = 0;
       int pdgId = 0;
       TLorentzVector matchedPartonLV = Muon4;
       for (unsigned int ip=0; ip<partonPdgId.size(); ++ip) {
	 TLorentzVector partLV = partonLV.at(ip);
	 float drJetMuon = deltaR(muon_eta[index],muon_phi[index],
				  partLV.Eta(),partLV.Phi());
	 if (drJetMuon<dRmin) {
	   dRmin = drJetMuon;
	   int absFlav = TMath::Abs(partonPdgId.at(ip));
	   pdgId = partonPdgId.at(ip);	   
	   flavour = 0;
	   /*
	   if (absFlav==21)
	     flavour = 1;
	   else if (absFlav>=1&&absFlav<=3)
	     flavour = 2;
	   else if (absFlav==4)
	     flavour = 3;
	   else if (absFlav==5)
	     flavour = 4;
	   */
	   if (absFlav>=1&&absFlav<=3)
	     flavour = 1;
	   else if (absFlav==4)
	     flavour = 2;
	   else if (absFlav==5)
	     flavour = 3;
	   qnet = float(muon_charge[index])*float(pdgId);

	   //	   if (flavour==0||flavour==1) qnet = -1.0;
	   if (flavour==0) qnet = -1.0;

	   matchedPartonLV = partLV;
	 }
       }

       unsigned int net = 0;
       if (qnet>0.0) net = 1;

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
				  track_eta[iTrk],track_phi[iTrk]);
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

       bool sigMu = trkSignalMu.size()==1 && trkMu.size()==1;

       bool bkgdMu = 
	 (trkSignalMu.size()==1 && trkSoftMu.size()==1 && trkMu.size()==2) ||
	 (trkSignalMu.size()==1 && trkSoftMu.size()==2 && trkMu.size()==3) ||
	 (trkSignalMu.size()==1 && trkSoftMu.size()==3 && trkMu.size()==4);

       bool sbMu = sigMu || bkgdMu;
       
       PartonMultiplicityMuH->Fill(float(partons),weight);
       PFJetMultiplicityMuH->Fill(float(pfjets),weight);

       ptMuH->Fill(muon_pt[index],weight);
       etaMuH->Fill(muon_eta[index],weight);
       dxyMuH->Fill(muon_dxy[index],weight);
       dzMuH->Fill(muon_dz[index],weight);
       
       nTracksMuH->Fill(float(trkMu.size()),weight);
       nSoftTracksMuH->Fill(float(trkSoftMu.size()),weight);
       nSignalTracksMuH->Fill(float(trkSoftMu.size()),weight);
       deltaRPartonMuH->Fill(dRmin,weight);

       // bin of parton pT
       float partonP = matchedPartonLV.Pt();
       if (partonP>partonMomBins[nPartMom]) partonP = partonMomBins[nPartMom] - 1.0;
       unsigned int partonMomBin = PartonMomBinsH->FindBin(partonP) - 1;
       //       std::cout << "mom(parton) = " << partonP << "  bin = " << partonMomBin << std::endl;

       // bin of muon pT
       float muonP = Muon4.Pt();
       if (muonP>muonMomBins[nMuMom]) muonP = muonMomBins[nMuMom] - 1.0;
       unsigned int muonMomBin = MuonMomBinsH->FindBin(muonP) - 1;
       //       std::cout << "mom(muon) = " << muonP << "  bin = " << muonMomBin << std::endl;

       // --------------------
       // -- save muon info --
       // --------------------
       muons_partmom.push_back(partonMomBin);
       muons_flavour.push_back(flavour);
       muons_muonmom.push_back(muonMomBin);
       muons_net.push_back(net);
       muons_mutrkmass.push_back(muonTrkMass);
       muons_pdgid.push_back(pdgId);
       int muon_region = -1; // doesn't pass selection
       if (sigMu) muon_region = 1; // signal region
       if (bkgdMu) muon_region = 0; // background region
       muons_region.push_back(muon_region);

       if (debug>1) {
	 printf("%2i -> pT(mu) = %3.0f : eta(mu) = %5.2f : Pt = %3.0f : flavour = %1i : pdgId = %2i : netQ = %1i\n",
		imu,Muon4.Pt(),Muon4.Eta(),matchedPartonLV.Pt(),flavour,pdgId,net); 
	 printf("       region : %2i : muonMomBin = %1i : partonMomBin = %1i\n",muon_region,muonMomBin,partonMomBin);
	 if (muon_region>=0) 
	   printf("       mass(mu,trk) = %5.2F\n",muonTrkMass);
       }

       partonMu->Fill(0.5,weight);
       partonMuProbe[flavour][net][partonMomBin][muonMomBin]->Fill(0.5,weight);

       if (sigMu) {

	 partonMuSelectedIso->Fill(0.5,weight);

	 InvMassIsoH->Fill(muonTrkMass,weight);
	 MassMuTrk[flavour][net][partonMomBin][muonMomBin][1]->Fill(muonTrkMass,weight);
	 partonMuPass[flavour][net][partonMomBin][muonMomBin][1]->Fill(0.5,weight);

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
	 deltaRMuTrkIsoH->Fill(deltaRMuonTrk,weight);
       }

       if (bkgdMu) {

	 partonMuSelectedLooseIso->Fill(0.5,weight);

	 InvMassLooseIsoH->Fill(muonTrkMass,weight);
	 MassMuTrk[flavour][net][partonMomBin][muonMomBin][0]->Fill(muonTrkMass,weight);
	 partonMuPass[flavour][net][partonMomBin][muonMomBin][0]->Fill(0.5,weight);

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
	 deltaRMuTrkLooseIsoH->Fill(deltaRMuonTrk,weight);
	 deltaRPartonMuLooseIsoH->Fill(dRmin,weight);
       }

       // applying QCD Model to the inclusive selected sample of muons
       if (applyQCDModel) {	 
	 //	 std::cout << "model is applied" << std::endl;
	 int ireg = 1; // Iso
	 double prob = qcdModel->getProb(partonMomBin,muonMomBin,ireg,flavour,net,true);
	 partonMuModelledIso->Fill(0.5,weight*prob);
	 for (unsigned int imass=1; imass<=nbins; ++imass) {
	   double mass = massHist->GetBinCenter(imass);
	   double pdf = qcdModel->getMassPdf(partonMomBin,muonMomBin,ireg,flavour,net,imass,true);
	   ModelInvMassIsoH->Fill(mass,weight*pdf*prob);
	 }

	 ireg = 0;
	 prob = qcdModel->getProb(partonMomBin,muonMomBin,ireg,flavour,net,true);
	 partonMuModelledLooseIso->Fill(0.5,weight*prob);
	 for (unsigned int imass=1; imass<=nbins; ++imass) {
	   double mass = massHist->GetBinCenter(imass);
	   double pdf = qcdModel->getMassPdf(partonMomBin,muonMomBin,ireg,flavour,net,imass,true);
	   ModelInvMassLooseIsoH->Fill(mass,weight*pdf*prob);
	 }
       }

     }

     nGoodIsoMuonsH->Fill(float(nIsoMuons),weight);
     nGoodLooseIsoMuonsH->Fill(float(nLooseIsoMuons),weight);

     // *****************************************
     // **** finding pair of same sign muons ****
     // *****************************************
     if (muons.size()<2) continue;
     unsigned int mu1 = 0;
     unsigned int mu2 = 0;
     double ptmax = -1;
     for (unsigned int i1=0; i1<muons.size()-1; ++i1) {
       for (unsigned int i2=i1+1; i2<muons.size(); ++i2) {

	 unsigned int index1 = muons[i1];
	 unsigned int index2 = muons[i2];

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
	   ptmax = ptsum;
	   mu1 = i1;
	   mu2 = i2;
	   if (muon_pt[index2]>muon_pt[index1]) {
	     mu1 = i2;
	     mu2 = i1;
	   }
	 }
       }
     }

     if (ptmax>10.) {

       unsigned int index1 = muons[mu1];
       unsigned int index2 = muons[mu2];

       int muonmom1 = muons_muonmom[mu1];
       int net1 = muons_net[mu1];
       int flavour1 = muons_flavour[mu1];
       int partmom1 = muons_partmom[mu1];
       int region1 = muons_region[mu1];
       float mutrk_mass1 = muons_mutrkmass[mu1];
       int pdgid1 = muons_pdgid[mu1];

       int muonmom2 = muons_muonmom[mu2];
       int net2 = muons_net[mu2];
       int flavour2 = muons_flavour[mu2];
       int partmom2 = muons_partmom[mu2];
       int region2 = muons_region[mu2];
       float mutrk_mass2 = muons_mutrkmass[mu2];
       int pdgid2 = muons_pdgid[mu2];

       if (debug>0) {
	 cout << endl;
	 cout << "same sign uons found ---> " << endl;
	 cout << "muon 1 : "
	      << "  pdgid = " << pdgid1
	      << "  flavor = " << flavour1
	      << "  charge = " << muon_charge[index1]
	      << "  net = " << net1 	   
	      << "  muonmom = " << muonmom1
	      << "  partmom = " << partmom1 
	      << "  reg = " << region1;
	 if (region1>=0) cout << "  mass(mu,trk) = " << mutrk_mass1;
	 cout << endl;
	 cout << "muon 2 : "
	      << "  pdgid = " << pdgid2 
	      << "  flavour = " << flavour2
	      << "  charge = " << muon_charge[index2]
	      << "  net = " << net2 
	      << "  muonmom  = " << muonmom2 
	      << "  partmom = " << partmom2 
	      << "  reg = " << region2;
	 if (region2>=0) cout << "  mass(mu,trk) = " << mutrk_mass2;
	 cout << endl;
	 cout << endl;
       }

       partonMu_SS->Fill(0.5,2.0*weight);
       partonMuProbe_SS[flavour1][net1][partmom1][muonmom1]->Fill(0.5,weight);
       partonMuProbe_SS[flavour2][net2][partmom2][muonmom2]->Fill(0.5,weight);

       if (region1==1)
	 partonMuSelectedIso_SS->Fill(0.5,weight);
       if (region2==1)
	 partonMuSelectedIso_SS->Fill(0.5,weight);

       if (region1==0)
	 partonMuSelectedLooseIso_SS->Fill(0.5,weight);
       if (region2==0)
	 partonMuSelectedLooseIso_SS->Fill(0.5,weight);

       if (applyQCDModel) {

	 int ireg = 1;
	 double prob1 = qcdModel->getProb(partmom1,muonmom1,ireg,flavour1,net1,false);
	 double prob2 = qcdModel->getProb(partmom2,muonmom2,ireg,flavour2,net2,false);
	 partonMuModelledIso_SS->Fill(0.5,prob1*weight);
	 partonMuModelledIso_SS->Fill(0.5,prob2*weight);

	 ireg = 0;
	 prob1 = qcdModel->getProb(partmom1,muonmom1,ireg,flavour1,net1,false);
	 prob2 = qcdModel->getProb(partmom2,muonmom2,ireg,flavour2,net2,false);
	 partonMuModelledLooseIso_SS->Fill(0.5,prob1*weight);
	 partonMuModelledLooseIso_SS->Fill(0.5,prob2*weight);

       }

       if (region1>=0) {
	 partonMuPass_SS[flavour1][net1][partmom1][muonmom1][region1]->Fill(0.5,weight);
	 MassMuTrk_SS[flavour1][net1][partmom1][muonmom1][region1]->Fill(mutrk_mass1,weight);
       } 
       if (region2>=0) {
	 partonMuPass_SS[flavour2][net2][partmom2][muonmom2][region2]->Fill(0.5,weight);
	 MassMuTrk_SS[flavour2][net2][partmom2][muonmom2][region2]->Fill(mutrk_mass2,weight);
       }

       if (region1==0) {
	 InvMassDimuonLooseIsoH->Fill(mutrk_mass1,weight);	 
       }
       if (region1==1)
	 InvMassDimuonIsoH->Fill(mutrk_mass1,weight);

       if (region2==0)
	 InvMassDimuonLooseIsoH->Fill(mutrk_mass2,weight);
       if (region2==1)
	 InvMassDimuonIsoH->Fill(mutrk_mass2,weight);

       if (applyQCDModel) {

	 // ---------------
	 // QCD model test
	 // ---------------
	 for (unsigned int imass=1; imass<=nbins; ++imass) {
	   double mass = massHist->GetBinCenter(imass);
	   // closure 
	   // Iso region
	   int ireg = 1;
	   double prob1 = qcdModel->getProb(partmom1,muonmom1,ireg,flavour1,net1,false);
	   double prob2 = qcdModel->getProb(partmom2,muonmom2,ireg,flavour2,net2,false);
	   double pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg,flavour1,net1,imass,false);
	   double pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass,false);
	   ClosureInvMassDimuonIsoH->Fill(mass,weight*pdf1*prob1);
	   ClosureInvMassDimuonIsoH->Fill(mass,weight*pdf2*prob2);
	   // LooseIso region
	   ireg = 0;
	   prob1 = qcdModel->getProb(partmom1,muonmom1,ireg,flavour1,net1,false);
	   pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg,flavour1,net1,imass,false);
	   prob2 = qcdModel->getProb(partmom2,muonmom2,ireg,flavour2,net2,false);
	   pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass,false);
	   ClosureInvMassDimuonLooseIsoH->Fill(mass,weight*pdf1*prob1);
	   ClosureInvMassDimuonLooseIsoH->Fill(mass,weight*pdf2*prob2);

	   // hybrid model
	   // Iso region
	   ireg = 1;
	   prob1 = qcdModel->getProb(partmom1,muonmom1,ireg,flavour1,net1,false);
	   pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg,flavour1,net1,imass,true);
	   prob2 = qcdModel->getProb(partmom2,muonmom2,ireg,flavour2,net2,false);
	   pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass,true);
	   HybridModelInvMassDimuonIsoH->Fill(mass,weight*pdf1*prob1);
	   HybridModelInvMassDimuonIsoH->Fill(mass,weight*pdf2*prob2);
	   // LooseIso region
	   ireg = 0;
	   prob1 = qcdModel->getProb(partmom1,muonmom1,ireg,flavour1,net1,false);
	   pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg,flavour1,net1,imass,true);
	   prob2 = qcdModel->getProb(partmom2,muonmom2,ireg,flavour2,net2,false);
	   pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass,true);
	   HybridModelInvMassDimuonLooseIsoH->Fill(mass,weight*pdf1*prob1);
	   HybridModelInvMassDimuonLooseIsoH->Fill(mass,weight*pdf2*prob2);

	   // inclusive model
	   // Iso region
	   ireg = 1;
	   prob1 = qcdModel->getProb(partmom1,muonmom1,ireg,flavour1,net1,true);
	   pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg,flavour1,net1,imass,true);
	   prob2 = qcdModel->getProb(partmom2,muonmom2,ireg,flavour2,net2,true);
	   pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass,true);
	   InclusiveModelInvMassDimuonIsoH->Fill(mass,weight*pdf1*prob1);
	   InclusiveModelInvMassDimuonIsoH->Fill(mass,weight*pdf2*prob2);
	   // LooseIso region
	   ireg = 0;
	   prob1 = qcdModel->getProb(partmom1,muonmom1,ireg,flavour1,net1,true);
	   pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg,flavour1,net1,imass,true);
	   prob2 = qcdModel->getProb(partmom2,muonmom2,ireg,flavour2,net2,true);
	   pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass,true);
	   InclusiveModelInvMassDimuonLooseIsoH->Fill(mass,weight*pdf1*prob1);
	   InclusiveModelInvMassDimuonLooseIsoH->Fill(mass,weight*pdf2*prob2);

	 }
	 
	 // SR : Iso
	 // closure
	 int ireg = 1; // signal region!
	 double prob1 = qcdModel->getProb(partmom1,muonmom1,ireg,flavour1,net1,false);
	 double prob2 = qcdModel->getProb(partmom2,muonmom2,ireg,flavour2,net2,false);
	 for (unsigned int imass=1; imass<=nbins; ++imass) {
	   double mass = massHist->GetBinCenter(imass);
	   double pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg,flavour1,net1,imass,false);
	   double pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass,false);
	   ClosureInvMassH->Fill(mass,weight*prob1*prob2*pdf1);
	   ClosureInvMassH->Fill(mass,weight*prob1*prob2*pdf2);	 
	   for (unsigned int imass2=1; imass2<=nbins; ++imass2) {
	     double mass2 = massHist->GetBinCenter(imass2);;
	     pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass2,false);
	     ClosureInvMass2DH->Fill(mass,mass2,weight*prob1*prob2*pdf1*pdf2);
	   }
	 }

	 // SR : Iso
	 // hybrid model
	 for (unsigned int imass=1; imass<=nbins; ++imass) {
           double mass = massHist->GetBinCenter(imass);
	   double pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg,flavour1,net1,imass,true);
	   double pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass,true);
	   HybridModelInvMassH->Fill(mass,weight*prob1*prob2*pdf1);
	   HybridModelInvMassH->Fill(mass,weight*prob1*prob2*pdf2);	 
	   for (unsigned int imass2=1; imass2<=nbins; ++imass2) {
	     double mass2 = massHist->GetBinCenter(imass2);
	     pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass2,true);
	     HybridModelInvMass2DH->Fill(mass,mass2,weight*prob1*prob2*pdf1*pdf2);
	   }
	 }
	 
	 // SR : Iso
	 // inclusive model
	 prob1 = qcdModel->getProb(partmom1,muonmom1,ireg,flavour1,net1,true);
	 prob2 = qcdModel->getProb(partmom2,muonmom2,ireg,flavour2,net2,true);
	 for (unsigned int imass=1; imass<=nbins; ++imass) {
           double mass = massHist->GetBinCenter(imass);
	   double pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg,flavour1,net1,imass,true);
	   double pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass,true);
	   InvMassH->Fill(mass,weight*prob1*prob2*pdf1);
	   InvMassH->Fill(mass,weight*prob1*prob2*pdf2);	 
	   for (unsigned int imass2=1; imass2<=nbins; ++imass2) {
	     double mass2 = massHist->GetBinCenter(imass2);
	     pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg,flavour2,net2,imass2,true);
	     InvMass2DH->Fill(mass,mass2,weight*prob1*prob2*pdf1*pdf2);
	   }
	 }

	 // CR : LooseIso (Y and X)
	 // closure
	 for (int ireg1=0; ireg1<2; ++ireg1) {
	   for (int ireg2=0; ireg2<2; ++ireg2) {
	     if (ireg1==1&&ireg2==1) continue; // signal region (skip)
	     double prob1 = qcdModel->getProb(partmom1,muonmom1,ireg1,flavour1,net1,false);
	     double prob2 = qcdModel->getProb(partmom2,muonmom2,ireg2,flavour2,net2,false);
	     for (unsigned int imass=1; imass<=nbins; ++imass) {
	       double mass = massHist->GetBinCenter(imass);
	       double pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg1,flavour1,net1,imass,false);
	       double pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg2,flavour2,net2,imass,false);
	       ClosureInvMass_ControlYH->Fill(mass,weight*prob1*prob2*pdf1);
	       ClosureInvMass_ControlYH->Fill(mass,weight*prob1*prob2*pdf2);
	       if (ireg1==0&&ireg2==0) {
		 ClosureInvMass_ControlXH->Fill(mass,weight*prob1*prob2*pdf1);
		 ClosureInvMass_ControlXH->Fill(mass,weight*prob1*prob2*pdf2);
	       }
	       for (unsigned int imass2=1; imass2<=nbins; ++imass2) {
		 double mass2 = massHist->GetBinCenter(imass2);
		 pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg2,flavour2,net2,imass2,false);
		 ClosureInvMass2D_ControlYH->Fill(mass,mass2,weight*prob1*prob2*pdf1*pdf2);
		 if (ireg1==0&&ireg2==0)
		   ClosureInvMass2D_ControlXH->Fill(mass,mass2,weight*prob1*prob2*pdf1*pdf2);
	       }
	     }
	   }
	 }

	 // CR : LooseIso (regions Y and X)
	 // hybrid model
	 for (int ireg1=0; ireg1<2; ++ireg1) {
	   for (int ireg2=0; ireg2<2; ++ireg2) {
	     if (ireg1==1&&ireg2==1) continue; // signal region (skip)
	     double prob1 = qcdModel->getProb(partmom1,muonmom1,ireg1,flavour1,net1,false);
	     double prob2 = qcdModel->getProb(partmom2,muonmom2,ireg2,flavour2,net2,false);
	     for (unsigned int imass=1; imass<=nbins; ++imass) {
	       double mass = massHist->GetBinCenter(imass);
	       double pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg1,flavour1,net1,imass,true);
	       double pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg2,flavour2,net2,imass,true);
	       HybridModelInvMass_ControlYH->Fill(mass,weight*prob1*prob2*pdf1);
	       HybridModelInvMass_ControlYH->Fill(mass,weight*prob1*prob2*pdf2);
	       if (ireg1==0&&ireg2==0) {
		 HybridModelInvMass_ControlXH->Fill(mass,weight*prob1*prob2*pdf1);
		 HybridModelInvMass_ControlXH->Fill(mass,weight*prob1*prob2*pdf2);
	       }
	       for (unsigned int imass2=1; imass2<=nbins; ++imass2) {
		 double mass2 = massHist->GetBinCenter(imass2);
		 pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg2,flavour2,net2,imass2,true);
		 HybridModelInvMass2D_ControlYH->Fill(mass,mass2,weight*prob1*prob2*pdf1*pdf2);
		 if (ireg1==0&&ireg2==0)
		   HybridModelInvMass2D_ControlXH->Fill(mass,mass2,weight*prob1*prob2*pdf1*pdf2);
	       }
	     }
	   }
	 }

	 // CR : LooseIso (regions Y and X)
	 // inclusive model
	 for (int ireg1=0; ireg1<2; ++ireg1) {
	   for (int ireg2=0; ireg2<2; ++ireg2) {
	     if (ireg1==1&&ireg2==1) continue; // signal region (not considered)
	     double prob1 = qcdModel->getProb(partmom1,muonmom1,ireg1,flavour1,net1,true);
	     double prob2 = qcdModel->getProb(partmom2,muonmom2,ireg2,flavour2,net2,true);
	     for (unsigned int imass=1; imass<=nbins; ++imass) {
	       double mass = massHist->GetBinCenter(imass);
	       double pdf1 = qcdModel->getMassPdf(partmom1,muonmom1,ireg1,flavour1,net1,imass,true);
	       double pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg2,flavour2,net2,imass,true);
	       InvMass_ControlYH->Fill(mass,weight*prob1*prob2*pdf1);
	       InvMass_ControlYH->Fill(mass,weight*prob1*prob2*pdf2);
	       if (ireg1==0&&ireg2==0) {
		 InvMass_ControlXH->Fill(mass,weight*prob1*prob2*pdf1);
		 InvMass_ControlXH->Fill(mass,weight*prob1*prob2*pdf2);
	       }
	       for (unsigned int imass2=1; imass2<=nbins; ++imass2) {
		 double mass2 = massHist->GetBinCenter(imass2);
		 pdf2 = qcdModel->getMassPdf(partmom2,muonmom2,ireg2,flavour2,net2,imass2,true);
		 InvMass2D_ControlYH->Fill(mass,mass2,weight*prob1*prob2*pdf1*pdf2);
		 if (ireg1==0&&ireg2==0) {
		   InvMass2D_ControlXH->Fill(mass,mass2,weight*prob1*prob2*pdf1*pdf2);
		 }
	       }
	     }
	   }
	 }



       }
     }

   } // event loop
   delete tree_;
   file_->Close();
   delete file_;
   
  } // filelist loop
  
  file->cd("");
  file->Write();
  file->Close();
  
  //delete file;
}// int main loop 

 

