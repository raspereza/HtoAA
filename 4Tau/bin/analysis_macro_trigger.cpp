#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "HtoAA/Utilities/interface/Config.h"
#include "HtoAA/Utilities/interface/json.h"
#include "HtoAA/Utilities/interface/functions.h"
#include "HtoAA/Utilities/interface/PileUp.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  // Data
  const bool isData = cfg.get<bool>("IsData");
  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection"); 

  // kinematic cuts on muons
  const float ptMuonCut      = cfg.get<float>("ptMuonCut");
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float ptMuonTagCut   = cfg.get<float>("ptMuonTagCut");
  const float etaMuonLowCut  = cfg.get<float>("etaMuonLowCut");
  const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
  const float etaMuonTagCut  = cfg.get<float>("etaMuonTagCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonCut     = cfg.get<float>("isoMuonCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");
  const bool applyMuonIso    = cfg.get<bool>("ApplyMuonIso");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dZleptonsCut   = cfg.get<float>("dZleptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");
  const bool oppositeSign    = cfg.get<bool>("OppositeSign");
  int chargeTagMuon          = cfg.get<int>("chargeTagMuon"); // -1/0/+1 (neg,all,pos)

  // triggers
  const bool applyTrigger = cfg.get<bool>("ApplyTrigger");

  const string hltDoubleMu           = cfg.get<string>("HLTDoubleMu");
  const string hltDoubleMuDZ         = cfg.get<string>("HLTDoubleMuDZ");
  const string hltDoubleMuSameSign   = cfg.get<string>("HLTDoubleMuSameSign");
  const string hltDoubleMuSameSignDZ = cfg.get<string>("HLTDoubleMuSameSignDZ");

  const string hltIsoSingleMu   = cfg.get<string>("HLTIsoSingleMu");  

  // HLT filters
  const string hltHighPtLeg   = cfg.get<string>("HLTHighPtLeg");
  const string hltLowPtLeg1   = cfg.get<string>("HLTLowPtLeg1");
  const string hltLowPtLeg2   = cfg.get<string>("HLTLowPtLeg2");
  const string dzFilter       = cfg.get<string>("dzFilter"); 
  const string sameSignFilter = cfg.get<string>("sameSignFilter");

  const string hltIsoSingleMuFilter = cfg.get<string>("HLTIsoSingleMuFilter");

  const string jsonFile = cfg.get<string>("jsonFile");
  const string puDataFile = cfg.get<string>("PileUpDataFile");
  const string puMCFile = cfg.get<string>("PileUpMCFile");

  float dxyTrkLooseCut = cfg.get<float>("dxyTrkLooseCut");
  float dzTrkLooseCut  = cfg.get<float>("dzTrkLooseCut");
  float ptTrkLooseCut  = cfg.get<float>("ptTrkLooseCut");

  float dxyTrkCut = cfg.get<float>("dxyTrkCut");
  float dzTrkCut  = cfg.get<float>("dzTrkCut");
  float ptTrkCut  = cfg.get<float>("ptTrkCut");

  float etaTrkCut = cfg.get<float>("etaTrkCut"); 

  TString PUDataFile(puDataFile);
  TString PUMCFile(puMCFile);

  // convesrsion from string to TString
  TString HLTDoubleMu(hltDoubleMu);
  TString HLTDoubleMuDZ(hltDoubleMuDZ);
  TString HLTDoubleMuSameSign(hltDoubleMuSameSign);
  TString HLTDoubleMuSameSignDZ(hltDoubleMuSameSignDZ);

  TString HLTIsoSingleMu(hltIsoSingleMu);

  TString HLTHighPtLeg(hltHighPtLeg);
  TString HLTLowPtLeg1(hltLowPtLeg1);
  TString HLTLowPtLeg2(hltLowPtLeg2);

  TString DZFilter(dzFilter);
  TString SameSignFilter(sameSignFilter);

  TString HLTIsoSingleMuFilter(hltIsoSingleMuFilter);
  
  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");


  // **** end of configuration

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

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  int presHLTDoubleMu = 0;
  int presHLTDoubleMuDZ = 0;
  int presHLTDoubleMuSameSign = 0;
  int presHLTDoubleMuSameSignDZ = 0;
  int presHLTIsoSingleMu = 0;

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
  Float_t genparticles_vx[1000];
  Float_t genparticles_vy[1000];
  Float_t genparticles_vz[1000];
  Int_t genparticles_pdgid[1000];
  Int_t genparticles_status[1000];
  UInt_t genparticles_info[1000];

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

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1F * histWeightsH = new TH1F("histWeightsH","",1,0.,2.);
  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);
  TH1F * weightsH = new TH1F("weightsH","",1,-0.5,0.5);

  // prescales
  TH1F * prescaleHLTDoubleMuH = new TH1F("prescaleHLTDoubleMuH","",21,-0.5,20.5);
  TH1F * prescaleHLTDoubleMuDZH = new TH1F("prescaleHLTDoubleMuDZH","",21,-0.5,20.5);
  TH1F * prescaleHLTDoubleMuSameSignH = new TH1F("prescaleHLTDoubleMuSameSignH","",21,-0.5,20.5);
  TH1F * prescaleHLTDoubleMuSameSignDZH = new TH1F("prescaleHLTDoubleMuSameSignDZH","",21,-0.5,20.5);
  TH1F * prescaleHLTIsoSingleMuH = new TH1F("prescaleHLTIsoSingleMuH","",21,-0.5,20.5);

  // J/Psi ->
  TH1F * JPsiMassDZFilterPassH =  new TH1F("JPsiMassDZFilterPassH","",200,2,4);
  TH1F * JPsiMassDZFilterFailH =  new TH1F("JPsiMassDZFilterFailH","",200,2,4);

  TH1F * JPsiMassDZFilterDz0to1PassH =  new TH1F("JPsiMassDZFilterDz0to1PassH","",200,2,4);
  TH1F * JPsiMassDZFilterDz0to1FailH =  new TH1F("JPsiMassDZFilterDz0to1FailH","",200,2,4);

  TH1F * JPsiMassDZFilterDz1to2PassH =  new TH1F("JPsiMassDZFilterDz1to2PassH","",200,2,4);
  TH1F * JPsiMassDZFilterDz1to2FailH =  new TH1F("JPsiMassDZFilterDz1to2FailH","",200,2,4);

  TH1F * JPsiMassDZFilterDzGt2PassH =  new TH1F("JPsiMassDZFilterDzGt2PassH","",200,2,4);
  TH1F * JPsiMassDZFilterDzGt2FailH =  new TH1F("JPsiMassDZFilterDzGt2FailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterPassH =  new TH1F("JPsiMassSameSignFilterPassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterFailH =  new TH1F("JPsiMassSameSignFilterFailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterDR0to0p15PassH =  new TH1F("JPsiMassSameSignFilterDR0to0p15PassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterDR0to0p15FailH =  new TH1F("JPsiMassSameSignFilterDR0to0p15FailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterDRGt0p15PassH =  new TH1F("JPsiMassSameSignFilterDRGt0p15PassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterDRGt0p15FailH =  new TH1F("JPsiMassSameSignFilterDRGt0p15FailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterDirIsoPassH =  new TH1F("JPsiMassSameSignFilterDirIsoPassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterDirIsoFailH =  new TH1F("JPsiMassSameSignFilterDirIsoFailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterInvIsoPassH =  new TH1F("JPsiMassSameSignFilterInvIsoPassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterInvIsoFailH =  new TH1F("JPsiMassSameSignFilterInvIsoFailH","",200,2,4);

  TH1F * muIsoLeadJPsiHLTDoubleMuH   = new TH1F("muIsoLeadJPsiHLTDoubleMuH","",200,0,2);
  TH1F * muIsoLeadJPsiHLTDoubleMuDZH = new TH1F("muIsoLeadJPsiHLTDoubleMuDZH","",200,0,2);

  TH1F * muIsoTrailJPsiHLTDoubleMuH   = new TH1F("muIsoTrailJPsiHLTDoubleMuH","",200,0,2);
  TH1F * muIsoTrailJPsiHLTDoubleMuDZH = new TH1F("muIsoTrailJPsiHLTDoubleMuDZH","",200,0,2);

  TH1F * dRmumuJPsiHLTDoubleMuH   = new TH1F("dRmumuJPsiHLTDoubleMuH","",200,0,2);
  TH1F * dRmumuJPsiHLTDoubleMuDZH = new TH1F("dRmumuJPsiHLTDoubleMuDZH","",200,0,2);

  TH1F * dZmumuJPsiHLTDoubleMuH   = new TH1F("dZmumuJPsiHLTDoubleMuH","",100,0,1);
  TH1F * dZmumuJPsiHLTDoubleMuDZH = new TH1F("dZmumuJPsiHLTDoubleMuDZH","",100,0,1);

  // Z ->
  TH1F * ZMassDZFilterPassH =  new TH1F("ZMassDZFilterPassH","",60,60,120);
  TH1F * ZMassDZFilterFailH =  new TH1F("ZMassDZFilterFailH","",60,60,120);

  TH1F * ZMassDZFilterDz0to1PassH =  new TH1F("ZMassDZFilterDz0to1PassH","",60,60,120);
  TH1F * ZMassDZFilterDz0to1FailH =  new TH1F("ZMassDZFilterDz0to1FailH","",60,60,120);

  TH1F * ZMassDZFilterDz1to2PassH =  new TH1F("ZMassDZFilterDz1to2PassH","",60,60,120);
  TH1F * ZMassDZFilterDz1to2FailH =  new TH1F("ZMassDZFilterDz1to2FailH","",60,60,120);

  TH1F * ZMassDZFilterDzGt2PassH =  new TH1F("ZMassDZFilterDzGt2PassH","",60,60,120);
  TH1F * ZMassDZFilterDzGt2FailH =  new TH1F("ZMassDZFilterDzGt2FailH","",60,60,120);

  TH1F * ZMassSameSignFilterPassH =  new TH1F("ZMassSameSignFilterPassH","",60,60,120);
  TH1F * ZMassSameSignFilterFailH =  new TH1F("ZMassSameSignFilterFailH","",60,60,120);

  int nEtaBins = 4;
  float etaBins[5] = {-0.001, 0.9, 1.2, 2.1, 2.4};
  /*
  int nPtBins = 18;
  float ptBins[19] = { 5,  7,  9, 11, 13, 
		      15, 17, 19, 21, 23, 
		      25, 27, 30, 40, 50,
  		      60, 80, 100, 200};
  */

  int nPtBins = 11;
  float ptBins[12] = { 6, 10, 14, 19,  24, 30, 
		      40, 50, 60, 80, 100, 200};

  int nPtBins45 = 9;
  float ptBins45[10] = {21,26,31,36,41,46,51,56,70,1000};

  int nDzBins = 6;
  float dzBins[7] = {-0.001,0.05,0.1,0.15,0.2,0.3,0.5};


  TH1F * ZMassHighPtLegPassH = new TH1F("ZMassHighPtLegPassH","",60,60,120);
  TH1F * ZMassHighPtLegFailH = new TH1F("ZMassHighPtLegFailH","",60,60,120);
  TH1F * ZMassLowPtLegPassH = new TH1F("ZMassLowPtLegPassH","",60,60,120);
  TH1F * ZMassLowPtLegFailH = new TH1F("ZMassLowPtLegFailH","",60,60,120);


  TString EtaBins[4] = {"EtaLt0p9","Eta0p9to1p2","Eta1p2to2p1","EtaGt2p1"};
  /*
  TString PtBins[18] = {"Pt5to7","Pt7to9","Pt9to11",
			"Pt11to13","Pt13to15","Pt15to17","Pt17to19",
			"Pt19to21","Pt21to23","Pt23to25","Pt25to27",
			"Pt27to30","Pt30to40","Pt40to50","Pt50to60",
			"Pt60to80","Pt80to100","PtGt100"};
  */
  TString PtBins[11] = {"Pt6to10" ,"Pt10to14", "Pt14to19","Pt19to24",
			"Pt24to30","Pt30to40", "Pt40to50","Pt50to60",
			"Pt60to80","Pt80to100","PtGt100"};

  TString DzBins[6] = {"Dz0to0p5","Dz0p5to1p0","Dz1p0to1p5",
		       "Dz1p5to2p0","Dz2p0to3p0","DzGt3p0"};
  TString PtBins45[9] = {"Pt21to26","Pt26to31","Pt31to36",
			 "Pt36to41","Pt41to46","Pt46to51",
			 "Pt51to56","Pt56to70","Pt70toInf"};

  TH1F * EtaBinsH = new TH1F("EtaBinsH","",nEtaBins,etaBins);
  TH1F * PtBinsH  = new TH1F("PtBinsH","",nPtBins,ptBins);
  TH1F * DzBinsH  = new TH1F("DzBinsH","",nDzBins,dzBins); 
  TH1F * PtBins45H = new TH1F("PtBins45H","",nPtBins45,ptBins45);

  for (int iB=0; iB<nEtaBins; ++iB)
    EtaBinsH->GetXaxis()->SetBinLabel(iB+1,EtaBins[iB]);

  for (int iB=0; iB<nPtBins; ++iB)
    PtBinsH->GetXaxis()->SetBinLabel(iB+1,PtBins[iB]);

  for (int iB=0; iB<nPtBins45; ++iB)
    PtBins45H->GetXaxis()->SetBinLabel(iB+1,PtBins45[iB]);

  for (int iB=0; iB<nDzBins; ++iB)
    DzBinsH->GetXaxis()->SetBinLabel(iB+1,DzBins[iB]);

  // (Pt,Eta)

  TH1F * ZMassHighPtLegPtEtaPassH[4][18];
  TH1F * ZMassHighPtLegPtEtaFailH[4][18];

  TH1F * ZMassLowPtLegPtEtaPassH[4][18];
  TH1F * ZMassLowPtLegPtEtaFailH[4][18];

  // Eta dependence

  TH1F * ZMassHighPtLegEtaPassH[4];
  TH1F * ZMassHighPtLegEtaFailH[4];

  TH1F * ZMassLowPtLegEtaPassH[4];
  TH1F * ZMassLowPtLegEtaFailH[4];

  // Pt dependence

  TH1F * ZMassHighPtLegPtPassH[18];
  TH1F * ZMassHighPtLegPtFailH[18];

  TH1F * ZMassLowPtLegPtPassH[18];
  TH1F * ZMassLowPtLegPtFailH[18];
  TH1F * JPsiMassDzPassH[6];
  TH1F * JPsiMassDzFailH[6];

  for (int iDz=0; iDz<nDzBins; ++iDz) {
    JPsiMassDzPassH[iDz] = new TH1F("JPsiMass_"+DzBins[iDz]+"_PassH","",200,2,4);
    JPsiMassDzFailH[iDz] = new TH1F("JPsiMass_"+DzBins[iDz]+"_FailH","",200,2,4);
  }

  for (int iEta=0; iEta<nEtaBins; ++iEta) {
    ZMassHighPtLegEtaPassH[iEta] = new TH1F("ZMassHighPtLeg_"+EtaBins[iEta]+"_PassH","",60,60,120);
    ZMassHighPtLegEtaFailH[iEta] = new TH1F("ZMassHighPtLeg_"+EtaBins[iEta]+"_FailH","",60,60,120);
    ZMassLowPtLegEtaPassH[iEta]  = new TH1F("ZMassLowPtLeg_"+EtaBins[iEta]+"_PassH","",60,60,120);
    ZMassLowPtLegEtaFailH[iEta]  = new TH1F("ZMassLowPtLeg_"+EtaBins[iEta]+"_FailH","",60,60,120);
  }

  for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassHighPtLegPtPassH[iPt] = new TH1F("ZMassHighPtLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassHighPtLegPtFailH[iPt] = new TH1F("ZMassHighPtLeg_"+PtBins[iPt]+"_FailH","",60,60,120);
      ZMassLowPtLegPtPassH[iPt]  = new TH1F("ZMassLowPtLeg_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassLowPtLegPtFailH[iPt]  = new TH1F("ZMassLowPtLeg_"+PtBins[iPt]+"_FailH","",60,60,120);
  }


  for (int iEta=0; iEta<nEtaBins; ++iEta) {

    for (int iPt=0; iPt<nPtBins; ++iPt) {
      ZMassHighPtLegPtEtaPassH[iEta][iPt] = new TH1F("ZMassHighPtLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassHighPtLegPtEtaFailH[iEta][iPt] = new TH1F("ZMassHighPtLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_FailH","",60,60,120);
      ZMassLowPtLegPtEtaPassH[iEta][iPt]  = new TH1F("ZMassLowPtLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_PassH","",60,60,120);
      ZMassLowPtLegPtEtaFailH[iEta][iPt]  = new TH1F("ZMassLowPtLeg_"+EtaBins[iEta]+"_"+PtBins[iPt]+"_FailH","",60,60,120);
    }

  }

  unsigned int iRun;
  unsigned int iEvent;
  //  TTree * eventTree = new TTree("eventTree","eventTree");
  //  eventTree->Branch("Run",&iRun,"Run/i");
  //  eventTree->Branch("Event",&iEvent,"Event/i");

  int nFiles = 0;
  int nEvents = 0;

  int selEvents = 0;
  int selEventsHLTDoubleMu = 0;
  int selEventsHLTDoubleMuDZ = 0;
  int selEventsHLTDoubleMuSameSignDZ = 0;

  int selPairs = 0;
  int selPairsHLTDoubleMu = 0;
  int selPairsHLTDoubleMuDZ = 0;
  int selPairsHLTDoubleMuSameSignDZ = 0;

  // PileUp
  PileUp * PUofficial = new PileUp();
  TFile * filePUOfficial_data = new TFile(TString(cmsswBase)+"/src/HtoAA/data/PileUpDistrib/"+PUDataFile,"read");
  TFile * filePUOfficial_MC = new TFile (TString(cmsswBase)+"/src/HtoAA/data/PileUpDistrib/"+PUMCFile, "read");
  TH1D * PUOfficial_data = (TH1D *)filePUOfficial_data->Get("pileup");
  TH1D * PUOfficial_mc = (TH1D *)filePUOfficial_MC->Get("pileup");
  PUofficial->set_h_data(PUOfficial_data);
  PUofficial->set_h_MC(PUOfficial_mc);

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = new TFile(TString(filen),"READ");
    if (file_==NULL) continue;
    if (file_==0) continue;
    if (file_->GetListOfKeys()->GetSize() == 0) continue; 
    if (file_->GetEND() > file_->GetSize()) continue; 

    if (file_->GetSeekKeys()<=file_->GetEND()-file_->GetSize())
      continue;

    if (file_->IsZombie()) {
      cout << "cannot open file " << filen << std::endl;
      continue;
    }
    
    // sum of weights (needed for normalization)
    TTree * _inittree = (TTree*)file_->Get("initroottree/AC1B");
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


    TTree * tree_ = (TTree*)file_->Get("makeroottree/AC1B");
    if (tree_==NULL) continue;    

    TH1D * histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");    
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());    
    std::cout << "      number of input events    = " << NE << std::endl;
    
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    tree_->SetMaxVirtualSize(3000000);
    // event info
    tree_->SetBranchAddress("event_nr", &event_nr);
    tree_->SetBranchAddress("event_run", &event_run);
    tree_->SetBranchAddress("event_luminosityblock", &event_luminosityblock);
    
    // Primary vertex
    tree_->SetBranchAddress("primvertex_x",&primvertex_x);
    tree_->SetBranchAddress("primvertex_y",&primvertex_y);
    tree_->SetBranchAddress("primvertex_z",&primvertex_z);
    
    // Muons
    tree_->SetBranchAddress("muon_count", &muon_count);
    tree_->SetBranchAddress("muon_px", muon_px);
    tree_->SetBranchAddress("muon_py", muon_py);
    tree_->SetBranchAddress("muon_pz", muon_pz);
    tree_->SetBranchAddress("muon_pt", muon_pt);
    tree_->SetBranchAddress("muon_eta", muon_eta);
    tree_->SetBranchAddress("muon_phi", muon_phi);
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
      tree_->SetBranchAddress("genparticles_vx", genparticles_vx);
      tree_->SetBranchAddress("genparticles_vy", genparticles_vy);
      tree_->SetBranchAddress("genparticles_vz", genparticles_vz);
      tree_->SetBranchAddress("genparticles_pdgid", genparticles_pdgid);
      tree_->SetBranchAddress("genparticles_status", genparticles_status);
      tree_->SetBranchAddress("genparticles_info", genparticles_info);
    }   

    int numberOfEntries = tree_->GetEntries();
    
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
    
    for (int iEntry=0; iEntry<numberOfEntries; iEntry++) { 
    
      tree_->GetEntry(iEntry);
      nEvents++;

      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      float weight = 1;

      if (isData) {
	if (applyGoodRunSelection) {
	  bool lumi = false;
	  int n = event_run;
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

      weightsH->Fill(0.0,weight);

      if (!isData) {
	float puWeight =  float(PUofficial->get_PUweight(double(numtruepileupinteractions)));
	weight *= puWeight;
      }
 

      presHLTDoubleMu = 0;
      presHLTDoubleMuDZ = 0;
      presHLTDoubleMuSameSign = 0;
      presHLTDoubleMuSameSignDZ = 0;
      presHLTIsoSingleMu = 0;
      for (std::map<string,int>::iterator it=hltriggerprescales->begin(); it!=hltriggerprescales->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(HLTDoubleMu)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  presHLTDoubleMu = it->second;
	}
	if (trigName.Contains(HLTDoubleMuDZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  presHLTDoubleMuDZ = it->second;
	}
	if (trigName.Contains(HLTDoubleMuSameSign)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  presHLTDoubleMuSameSign = it->second;
	}
	if (trigName.Contains(HLTDoubleMuSameSignDZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  presHLTDoubleMuSameSignDZ = it->second;
	}
	if (trigName.Contains(HLTIsoSingleMu)) {
	  presHLTIsoSingleMu = it->second;
	}
      }
      
      unsigned int nHighPtLeg = 0;
      bool isHighPtLeg = false;
      unsigned int nLowPtLeg = 0;
      bool isLowPtLeg = false;
      unsigned int nDZFilter = 0;
      bool isDZFilter = false;
      unsigned int nSameSignFilter = 0;
      bool isSameSignFilter = false;

      unsigned int nIsoSingleMuFilter = 0;
      bool isIsoSingleMuFilter;

      unsigned int nfilters = hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(hltfilters->at(i));
	if (HLTFilter==HLTHighPtLeg) {
	  nHighPtLeg = i;
	  isHighPtLeg = true;
	}
	if (HLTFilter==HLTLowPtLeg1||HLTFilter==HLTLowPtLeg2) {
	  nLowPtLeg = i;
	  isLowPtLeg = true;
	}
	if (HLTFilter==DZFilter) {
	  nDZFilter = i;
	  isDZFilter = true;
	}
	if (HLTFilter==SameSignFilter) {
	  nSameSignFilter = i;
	  isSameSignFilter = true;
	}
	if (HLTFilter==HLTIsoSingleMuFilter) {
	  nIsoSingleMuFilter = i;
	  isIsoSingleMuFilter = true;
	}
      }
      if (!isHighPtLeg) {
	std::cout << "HLT filter " << HLTHighPtLeg << " not found" << std::endl;
	exit(-1);
      }
      if (!isLowPtLeg) {
	std::cout << "HLT filter " << HLTLowPtLeg1
		  << " or " << HLTLowPtLeg2
		  << " not found" << std::endl;
	exit(-1);
      }
      if (!isDZFilter) {
	std::cout << "HLT filter " << DZFilter << " not found" << std::endl;
	exit(-1);
      }
      if (!isSameSignFilter) {
	std::cout << "HLT filter " << SameSignFilter << " not found" << std::endl;
	exit(-1);
      }
      if (!isIsoSingleMuFilter) {
	std::cout << "HLT filter " << HLTIsoSingleMuFilter << " not found" << std::endl;
        exit(-1);
      }

      prescaleHLTDoubleMuH->Fill(float(presHLTDoubleMu),weight);
      prescaleHLTDoubleMuDZH->Fill(float(presHLTDoubleMuDZ),weight);
      prescaleHLTDoubleMuSameSignH->Fill(float(presHLTDoubleMuSameSign),weight);
      prescaleHLTDoubleMuSameSignDZH->Fill(float(presHLTDoubleMuSameSignDZ),weight);
      prescaleHLTIsoSingleMuH->Fill(float(presHLTIsoSingleMu),weight);
      
      // muon selection
      vector<unsigned int> muons; muons.clear();
      vector<unsigned int> isoMuons; isoMuons.clear();

      for (unsigned int im = 0; im<muon_count; ++im) {
	if (muon_pt[im]<ptMuonCut) continue;
	if (fabs(muon_eta[im])>etaMuonLowCut) continue;
	if (fabs(muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(muon_dz[im])>dzMuonCut) continue;
	if (applyMuonId && !muon_isMedium[im]) continue;
	float neutralHadIsoMu = muon_neutralHadIso[im];
	float photonIsoMu = muon_photonIso[im];
	float chargedHadIsoMu = muon_chargedHadIso[im];
	float puIsoMu = muon_puIso[im];
	
	float neutralIsoMu = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
	neutralIsoMu = TMath::Max(float(0),neutralIsoMu);
	float absIsoMu = chargedHadIsoMu + neutralIsoMu;
	float relIsoMu = absIsoMu/muon_pt[im];
	if (relIsoMu<isoMuonCut) {
	  muons.push_back(im);
	  isoMuons.push_back(im); 
	}
	/*
	float muonEta = muon_eta[im];
	float muonPhi = muon_phi[im];

	unsigned int ntracks = 0;
	unsigned int nsigtracks = 0;
	bool isoTrk = true;
	for (unsigned int iTrk=0; iTrk<track_count; ++iTrk) {
	  if (fabs(track_charge[iTrk])<0.1) continue;
	  if (fabs(track_eta[iTrk])>etaTrkCut) continue;
	  float trackEta = track_eta[iTrk];
	  float trackPhi = track_phi[iTrk];
	  float dR = deltaR(muonEta,muonPhi,trackEta,trackPhi);
	  if (dR>0.01&&dR<0.4) {

	    if (fabs(track_dxy[iTrk])<dxyTrkLooseCut &&
		fabs(track_dz[iTrk])<dzTrkLooseCut &&
		fabs(track_pt[iTrk])>ptTrkLooseCut) ntracks++;

	    if (fabs(track_dxy[iTrk])<dxyTrkCut &&
		fabs(track_dz[iTrk])<dzTrkCut &&
		fabs(track_pt[iTrk])>ptTrkCut) nsigtracks++;	    

	  }
	}
	if (ntracks>1) isoTrk = false;
	if (nsigtracks!=1) isoTrk = false;
	if (isoTrk) isoMuons.push_back(im);
	*/
      }


      bool isPairSelected = false;
      bool isPairSelectedHLTDoubleMu = false;
      bool isPairSelectedHLTDoubleMuDZ = false;
      bool isPairSelectedHLTDoubleMuSameSignDZ = false;

      //      std::cout << "number of tag muons   = " << muons.size() << std::endl;
      //      std::cout << "number of probe muons = " << isoMuons.size() << std::endl;

      // selecting muon pair
      for (unsigned int im1=0; im1<muons.size(); ++im1) {
	int  mu1Index = muons[im1];
	bool mu1MatchHighPt = false;
	bool mu1MatchLowPt  = false;
	bool mu1MatchDz     = false;
	bool mu1MatchSS     = false;
	bool mu1MatchIsoSingleMu = false;
	for (unsigned int iT=0; iT<trigobject_count; ++iT) {
	  float dRtrig = deltaR(muon_eta[mu1Index],muon_phi[mu1Index],
				trigobject_eta[iT],trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (trigobject_filters[iT][nHighPtLeg]) // Muon17 Leg
	    mu1MatchHighPt = true;
	  if (trigobject_filters[iT][nLowPtLeg]) // Muon8 Leg
	    mu1MatchLowPt = true;
	  if (trigobject_filters[iT][nDZFilter]) // DZ filter
	    mu1MatchDz = true;
	  if (trigobject_filters[iT][nSameSignFilter]) // same-sign filter
	    mu1MatchSS = true;
	  if (trigobject_filters[iT][nIsoSingleMuFilter]) // HLT_IsoMu filter
	    mu1MatchIsoSingleMu = true;
	}

	bool mu1HighPt = mu1MatchHighPt && muon_pt[mu1Index]>ptMuonHighCut && fabs(muon_eta[mu1Index])<etaMuonHighCut;
	bool mu1LowPt  = mu1MatchLowPt && muon_pt[mu1Index]>ptMuonLowCut && fabs(muon_eta[mu1Index])<etaMuonLowCut;
	bool mu1IsoSingleMu = mu1MatchIsoSingleMu && muon_pt[mu1Index]>ptMuonTagCut && fabs(muon_eta[mu1Index])<etaMuonTagCut;
	float q1 = muon_charge[mu1Index];
	
	for (unsigned int im2=0; im2<isoMuons.size(); ++im2) {

	  int  mu2Index = isoMuons[im2];
	  if (mu1Index==mu2Index) continue;

	  float q2 = muon_charge[mu2Index];
	  if (oppositeSign && (q1*q2>0)) continue;

	  float dRmumu = deltaR(muon_eta[mu1Index],muon_phi[mu1Index],
				muon_eta[mu2Index],muon_phi[mu2Index]);

	  bool mu2MatchHighPt = false;
	  bool mu2MatchLowPt  = false;
	  bool mu2MatchDz   = false;
	  bool mu2MatchSS   = false;

	  for (unsigned int iT=0; iT<trigobject_count; ++iT) {
	    float dRtrig = deltaR(muon_eta[mu2Index],muon_phi[mu2Index],
				  trigobject_eta[iT],trigobject_phi[iT]);
	    if (dRtrig>DRTrigMatch) continue;
	    if (trigobject_filters[iT][nHighPtLeg]) // Muon17 Leg
              mu2MatchHighPt = true;
	    if (trigobject_filters[iT][nLowPtLeg]) // Muon8 Leg
              mu2MatchLowPt = true;
	    if (trigobject_filters[iT][nDZFilter]) // DZ filter
              mu2MatchDz = true;
	    if (trigobject_filters[iT][nSameSignFilter]) // same-sign filter
              mu2MatchSS = true;
	  }
	  bool mu2HighPt = mu2MatchHighPt && muon_pt[mu2Index]>ptMuonHighCut && fabs(muon_eta[mu2Index])<etaMuonHighCut;
	  bool mu2LowPt  = mu2MatchLowPt && muon_pt[mu2Index]>ptMuonLowCut && fabs(muon_eta[mu2Index])<etaMuonLowCut;

	  bool triggerMatch = (mu1HighPt&&mu2LowPt) || (mu1LowPt&&mu2HighPt);
	  bool triggerMatchDz = triggerMatch && mu1MatchDz && mu2MatchDz;
	  bool triggerMatchSS = triggerMatchDz && mu1MatchSS && mu2MatchSS;


	  float dZ = fabs(muon_dz[mu1Index]-muon_dz[mu2Index]);

	  TLorentzVector mu1lv; mu1lv.SetXYZM(muon_px[mu1Index],
					      muon_py[mu1Index],
					      muon_pz[mu1Index],
					      MuMass);
	  TLorentzVector mu2lv; mu2lv.SetXYZM(muon_px[mu2Index],
                                              muon_py[mu2Index],
                                              muon_pz[mu2Index],
                                              MuMass);

	  float mass = (mu1lv+mu2lv).M();

	  float absIso1 = muon_chargedHadIso[mu1Index];
	  float relIso1 = absIso1/muon_pt[mu1Index];

	  float absIso2 = muon_chargedHadIso[mu2Index];
	  float relIso2 = absIso2/muon_pt[mu2Index];

	  float mu1RelIso = relIso1;

	  if (muon_pt[mu2Index]>muon_pt[mu1Index]) {
	    float temp = relIso1;
	    relIso1 = relIso2;
	    relIso2 = temp;
	  }

	  bool dirIso = (relIso2<isoMuonCut) && (relIso1<isoMuonCut);

	  isPairSelected = true;
	  selPairs++;

	  float DZ = TMath::Min(float(dZ),float(0.5));
	  int DZBin  = binNumber(DZ,nDzBins,dzBins);

	  if (mu1MatchIsoSingleMu && triggerMatch) { // pass HLT_HighPt_LowPt and HLT_SingleMuon
	    isPairSelectedHLTDoubleMu = true;
	    selPairsHLTDoubleMu++;
	    if (mass>3.0&&mass<3.2) {
	      muIsoLeadJPsiHLTDoubleMuH->Fill(relIso1,weight);
	      muIsoTrailJPsiHLTDoubleMuH->Fill(relIso2,weight);
	      dRmumuJPsiHLTDoubleMuH->Fill(dRmumu,weight);
	      dZmumuJPsiHLTDoubleMuH->Fill(dZ,weight);
	    }
	    if (triggerMatchDz) { // pass HLT_HighPt_LowPt_DZ
	      if (dZ<dZleptonsCut) {
		JPsiMassDZFilterPassH->Fill(mass,weight);
		if (dRmumu>dRleptonsCut) ZMassDZFilterPassH->Fill(mass,weight);
	      }
	      JPsiMassDzPassH[DZBin]->Fill(mass,weight);
	      if (dZ<0.1) {
		JPsiMassDZFilterDz0to1PassH->Fill(mass,weight);
		if (dRmumu>dRleptonsCut) ZMassDZFilterDz0to1PassH->Fill(mass,weight);
	      }
	      else if (dZ<0.2) {
		JPsiMassDZFilterDz1to2PassH->Fill(mass,weight); 
                if (dRmumu>dRleptonsCut) ZMassDZFilterDz1to2PassH->Fill(mass,weight);
	      }
	      else {
		JPsiMassDZFilterDzGt2PassH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterDzGt2PassH->Fill(mass,weight);
	      }
	    }
	    else { // fail HLT_HighPt_LowPt_DZ
	      if (dZ<dZleptonsCut) {
		JPsiMassDZFilterFailH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterFailH->Fill(mass,weight);
              }
	      JPsiMassDzFailH[DZBin]->Fill(mass,weight);
	      if (dZ<0.1) {
                JPsiMassDZFilterDz0to1FailH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterDz0to1FailH->Fill(mass,weight);
              }
              else if (dZ<0.2) {
                JPsiMassDZFilterDz1to2FailH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterDz1to2FailH->Fill(mass,weight);
              }
              else {
                JPsiMassDZFilterDzGt2FailH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassDZFilterDzGt2FailH->Fill(mass,weight);
              }
	    }
	  }

	  if (triggerMatchDz && mu1MatchIsoSingleMu) { // pass HLT_HighPt_LowPt_DZ
	    isPairSelectedHLTDoubleMuDZ = true;
            selPairsHLTDoubleMuDZ++;
	    if (mass>3.0&&mass<3.2) {
	      muIsoLeadJPsiHLTDoubleMuDZH->Fill(relIso1,weight);
	      muIsoTrailJPsiHLTDoubleMuDZH->Fill(relIso2,weight);
	      dRmumuJPsiHLTDoubleMuDZH->Fill(dRmumu,weight);
	      dZmumuJPsiHLTDoubleMuDZH->Fill(dZ,weight);
	    }
	    if (triggerMatchSS) { // pass HLT_HighPt_LowPt_SameSign_DZ
              if (dZ<dZleptonsCut) {
                JPsiMassSameSignFilterPassH->Fill(mass,weight);
		if (dirIso)
		  JPsiMassSameSignFilterDirIsoPassH->Fill(mass,weight);
		else 
		  JPsiMassSameSignFilterInvIsoPassH->Fill(mass,weight);
                if (dRmumu>dRleptonsCut) ZMassSameSignFilterPassH->Fill(mass,weight);
		if (dRmumu<0.15)
		  JPsiMassSameSignFilterDR0to0p15PassH->Fill(mass,weight);
		else
		  JPsiMassSameSignFilterDRGt0p15PassH->Fill(mass,weight);
              }
            }
            else { // fail HLT_HighPt_LowPt_SameSign_DZ
              if (dZ<dZleptonsCut) {
                JPsiMassSameSignFilterFailH->Fill(mass,weight);
		if (dirIso)
                  JPsiMassSameSignFilterDirIsoFailH->Fill(mass,weight);
                else
                  JPsiMassSameSignFilterInvIsoFailH->Fill(mass,weight); 
                if (dRmumu>dRleptonsCut) ZMassSameSignFilterFailH->Fill(mass,weight);
		if (dRmumu<0.15)
		  JPsiMassSameSignFilterDR0to0p15FailH->Fill(mass,weight);
		else
		  JPsiMassSameSignFilterDRGt0p15FailH->Fill(mass,weight);
              }
            }
	  }

	  if (triggerMatchSS && mu1MatchIsoSingleMu) {
	    isPairSelectedHLTDoubleMuSameSignDZ = true;
            selPairsHLTDoubleMuSameSignDZ++;
	  }

	  if (mu1IsoSingleMu && mu1RelIso<isoMuonCut && dZ<dZleptonsCut && dRmumu>dRleptonsCut) { // Single muon selection

	    float mu2AbsEta = fabs(muon_eta[mu2Index]);
	    float mu2Pt = TMath::Max(float(5.01),TMath::Min(float(muon_pt[mu2Index]),float(99.9)));
	    int etaBin = binNumber(mu2AbsEta,nEtaBins,etaBins);
	    int ptBin  = binNumber(mu2Pt,nPtBins,ptBins);

	    bool chargeTagPassed = true;
	    if (chargeTagMuon<0 && muon_charge[mu1Index]>0) chargeTagPassed = false;
	    if (chargeTagMuon>0 && muon_charge[mu1Index]<0) chargeTagPassed = false;

	    if (chargeTagPassed) {

	      if (mu2MatchLowPt) { // LowPt Leg
		if (muon_pt[mu2Index]>ptMuonLowCut&&mu2AbsEta<etaMuonLowCut) 
		  ZMassLowPtLegPassH->Fill(mass,weight);
		if (muon_pt[mu2Index]>ptMuonLowCut)
		  ZMassLowPtLegEtaPassH[etaBin]->Fill(mass,weight);
		if (mu2AbsEta<etaMuonLowCut)
		  ZMassLowPtLegPtPassH[ptBin]->Fill(mass,weight);
	      }
	      else {
		if (muon_pt[mu2Index]>ptMuonLowCut&&mu2AbsEta<etaMuonLowCut)
		  ZMassLowPtLegFailH->Fill(mass,weight);
		if (muon_pt[mu2Index]>ptMuonLowCut)
		  ZMassLowPtLegEtaFailH[etaBin]->Fill(mass,weight);
		if (mu2AbsEta<etaMuonLowCut)
		  ZMassLowPtLegPtFailH[ptBin]->Fill(mass,weight);
	      }
	     
	      if (mu2MatchHighPt) { // HighPt Leg
		if (muon_pt[mu2Index]>ptMuonHighCut&&mu2AbsEta<etaMuonHighCut)  
		  ZMassHighPtLegPassH->Fill(mass,weight);
		if (muon_pt[mu2Index]>ptMuonHighCut) 
		  ZMassHighPtLegEtaPassH[etaBin]->Fill(mass,weight);
		if (mu2AbsEta<etaMuonHighCut)
		  ZMassHighPtLegPtPassH[ptBin]->Fill(mass,weight);
	      }
	      else {
		if (muon_pt[mu2Index]>ptMuonHighCut&&mu2AbsEta<etaMuonHighCut)
		  ZMassHighPtLegFailH->Fill(mass,weight);
		if (muon_pt[mu2Index]>ptMuonHighCut) 
		  ZMassHighPtLegEtaFailH[etaBin]->Fill(mass,weight);
		if (mu2AbsEta<etaMuonHighCut)
		  ZMassHighPtLegPtFailH[ptBin]->Fill(mass,weight);
	      }

	      if (mu2AbsEta<etaMuonHighCut) { // bins in (eta,pt)

		if (mu2MatchLowPt)
		  ZMassLowPtLegPtEtaPassH[etaBin][ptBin]->Fill(mass,weight);
		else
		  ZMassLowPtLegPtEtaFailH[etaBin][ptBin]->Fill(mass,weight);

		if (mu2MatchHighPt)
		  ZMassHighPtLegPtEtaPassH[etaBin][ptBin]->Fill(mass,weight);
		else
		  ZMassHighPtLegPtEtaFailH[etaBin][ptBin]->Fill(mass,weight);

	      }

	    }

	  }

	}

      }
    
      if (isPairSelected)
	selEvents++;
      if (isPairSelectedHLTDoubleMu)
	selEventsHLTDoubleMu++;
      if (isPairSelectedHLTDoubleMuDZ)
	selEventsHLTDoubleMuDZ++;
      if (isPairSelectedHLTDoubleMuSameSignDZ)
	selEventsHLTDoubleMuSameSignDZ++;

      
    } // end of file processing (loop over events in one file)
    nFiles++;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << "Total number of selected pairs  = " << selPairs << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events (HLT_Dimuon) = " << selEventsHLTDoubleMu << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Dimuon) = " << selPairsHLTDoubleMu  << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events (HLT_Dimuon_DZ) = " << selEventsHLTDoubleMuDZ << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Dimuon_DZ) = " << selPairsHLTDoubleMuDZ << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events (HLT_Dimuon_SameSign_DZ) = " << selEventsHLTDoubleMuSameSignDZ << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Dimuon_SameSign_DZ) = " << selPairsHLTDoubleMuSameSignDZ << std::endl;
  std::cout << std::endl;
  file->cd("");
  file->Write();
  file->Close();
  delete file;
  
}



