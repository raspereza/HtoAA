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

using namespace std;


int main(int argc, char * argv[]) {
  
  if (argc!=3) {
    std::cout << "Usage of the program : analysis_macro_azh [config] [file_list]" << std::endl;
    std::cout << "file_list : file list of RooT files to be processed. To run on Data the string has to include the string \"Data\"." << std::endl;
    exit(1);
  }


  // **** configuration
  Config cfg(argv[1]);

  string cmsswBase = (getenv ("CMSSW_BASE"));
  const int era = cfg.get<int>("Era");
  // trigger
  const string bTagAlgorithm = cfg.get<string>("BTagAlgorithm");
  const string bTagDiscriminator1 = cfg.get<string>("BTagDiscriminator1");
  const string bTagDiscriminator2 = cfg.get<string>("BTagDiscriminator2");
  const string bTagDiscriminator3 = cfg.get<string>("BTagDiscriminator3");

  const float btagCut = cfg.get<float>("BTagCut");
  const float bjetEtaCut = cfg.get<float>("BJetEtaCut");
  const float bjetPtCut = cfg.get<float>("BJetPtCut");

  const float jetEtaCut = cfg.get<float>("JetEtaCut");
  const float jetPtCut = cfg.get<float>("JetPtCut");

  const float ptLepCut = cfg.get<float>("PtLepCut");
  const float etaLepCut = cfg.get<float>("EtaLepCut");

  const float ptTauMinCut = cfg.get<float>("PtTauMinCut");
  const float ptTauMaxCut = cfg.get<float>("PtTauMaxCut");
  const float etaTauCut = cfg.get<float>("EtaTauCut");

  const float dRJetLep = cfg.get<float>("dRJetLep");

  TString BTagAlgorithm(bTagAlgorithm);
  TString BTagDiscriminator1(bTagDiscriminator1);
  TString BTagDiscriminator2(bTagDiscriminator2);
  TString BTagDiscriminator3(bTagDiscriminator3);

  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");

  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);
  // BTag SF file
  const string BtagSfFile = cfg.get<string>("BtagSfFile");
  BTagCalibration calib = BTagCalibration(bTagAlgorithm, BtagSfFile);
  BTagCalibrationReader reader_B = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central",{"up","down"});
  BTagCalibrationReader reader_C = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central",{"up","down"});
  BTagCalibrationReader reader_Light = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central",{"up","down"}); 

  reader_B.load(calib, BTagEntry::FLAV_B, "comb");
  reader_C.load(calib, BTagEntry::FLAV_C, "comb");
  reader_Light.load(calib, BTagEntry::FLAV_UDSG, "incl");
  
  // BTAG efficiency for various flavours ->
  TString fileBtagEff = (TString)cfg.get<string>("BtagMCeffFile");
  TFile *fileTagging  = new TFile(fileBtagEff);
  TH2F  *tagEff_B     = (TH2F*)fileTagging->Get("btag_eff_b");
  TH2F  *tagEff_C     = (TH2F*)fileTagging->Get("btag_eff_c");
  TH2F  *tagEff_Light = (TH2F*)fileTagging->Get("btag_eff_oth");

  float MaxBJetPt = 1000.;
  float MinBJetPt = 20.;
  float MaxBJetEta = 2.4;
  float MinBJetEta = 0.0;

  // ********** end of configuration *******************

  std::ifstream fileList(argv[2]);
  std::ifstream fileListX(argv[2]);

  Float_t genweight;
  Float_t numtruepileupinteractions;
  UInt_t  genparticles_count;
  Float_t genparticles_e[1000];
  Float_t genparticles_px[1000];
  Float_t genparticles_py[1000];
  Float_t genparticles_pz[1000];
  Int_t   genparticles_pdgid[1000];
  Int_t   genparticles_status[1000];
  UInt_t  genparticles_info[1000];
  Int_t   genparticles_fromHardProcess[1000];

  UInt_t  gentau_count;
  Float_t gentau_charge[100];
  Float_t gentau_visible_e[100]; 
  Float_t gentau_visible_px[100];
  Float_t gentau_visible_py[100];
  Float_t gentau_visible_pz[100];
  Int_t   gentau_decayMode[100]; 
  Int_t   gentau_status[100]; 
  Int_t   gentau_isLastCopy[100];
  Int_t   gentau_fromHardProcess[100];

  UInt_t   pfjet_count;
  Float_t  pfjet_e[200];   //[pfjet_count]
  Float_t  pfjet_pt[200];    //[pfjet_count]
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

  UInt_t   genjets_count;
  Float_t  genjets_pt[200];   //[genjets_count]
  Float_t  genjets_eta[200];   //[genjets_count]
  Float_t  genjets_phi[200];   //[genjets_count]


  std::map<std::string, int> * hltriggerresults = new std::map<std::string, int>() ;
  std::map<std::string, int> * hltriggerprescales = new std::map<std::string, int>() ;
  std::vector<std::string>   * hltfilters = new std::vector<std::string>();
  std::vector<std::string>   * btagdiscriminators = new std::vector<std::string>();

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

  // PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/HtoAA/data/PileUpDistrib/"+PileUpDataFile,"read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/HtoAA/data/PileUpDistrib/"+PileUpMCFile, "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  TFile * file = new TFile(FullName+TString(".root"),"recreate");
  file->cd("");
  
  float weight;
  float mcweight;
  float puweight;
  float btagweight;
  float btagweight_hf_up;
  float btagweight_hf_down;
  float btagweight_lf_up;
  float btagweight_lf_down;

  unsigned int nbtag;
  unsigned int nbtag_jes_up;
  unsigned int nbtag_jes_down;
  unsigned int nbtag_jer_up;
  unsigned int nbtag_jer_down;

  unsigned int njets;
  unsigned int njetspt20;
  float ptjet[100];
  float etajet[100];
  float phijet[100];
  float ejet[100];
  int flavorjet[100];
  bool tagjet[100];
  
  int pdgPosLep;
  float ptPosLep;
  float etaPosLep;
  float phiPosLep;

  int pdgNegLep;
  float ptNegLep;
  float etaNegLep;
  float phiNegLep;

  float ptPosTau;
  float etaPosTau;
  float phiPosTau;
  int dmPosTau;

  float ptNegTau;
  float etaNegTau;
  float phiNegTau;
  int dmNegTau;
  
  TTree * tuple = new TTree("tuple","Tuple");

  tuple->Branch("mcweight",&mcweight,"mcweight/F");
  tuple->Branch("puweight",&puweight,"puweight/F");
  tuple->Branch("weight",&weight,"weight/F");
  tuple->Branch("btagweight",&btagweight,"btagweight/F");
  tuple->Branch("btagweight_hf_up",&btagweight_hf_up,"btagweight_hf_up/F");
  tuple->Branch("btagweight_hf_down",&btagweight_hf_down,"btagweight_hf_down/F");
  tuple->Branch("btagweight_lf_up",&btagweight_lf_up,"btagweight_lf_up/F");
  tuple->Branch("btagweight_lf_down",&btagweight_lf_down,"btagweight_lf_down/F");

  tuple->Branch("pdgPosLep",&pdgPosLep,"pdgPosLep/I");
  tuple->Branch("ptPosLep",&ptPosLep,"ptPosLep/F");
  tuple->Branch("etaPosLep",&etaPosLep,"etaPosLep/F");
  tuple->Branch("phiPosLep",&phiPosLep,"phiPosLep/F");

  tuple->Branch("pdgNegLep",&pdgNegLep,"pdgNegLep/I");
  tuple->Branch("ptNegLep",&ptNegLep,"ptNegLep/F");
  tuple->Branch("etaNegLep",&etaNegLep,"etaNegLep/F");
  tuple->Branch("phiNegLep",&phiNegLep,"phiNegLep/F");

  tuple->Branch("dmPosTau",&dmPosTau,"dmPosTau/I");
  tuple->Branch("ptPosTau",&ptPosTau,"ptPosTau/F");
  tuple->Branch("etaPosTau",&etaPosTau,"etaPosTau/F");
  tuple->Branch("phiPosTau",&phiPosTau,"phiPosTau/F");

  tuple->Branch("dmNegTau",&dmNegTau,"dmNegTau/I");
  tuple->Branch("ptNegTau",&ptNegTau,"ptNegTau/F");
  tuple->Branch("etaNegTau",&etaNegTau,"etaNegTau/F");
  tuple->Branch("phiNegTau",&phiNegTau,"phiNegTau/F");

  tuple->Branch("njets",&njets,"njets/i");
  tuple->Branch("njetspt20",&njetspt20,"njetspt20/i");
  tuple->Branch("nbtag",&nbtag,"nbtag/i");
  tuple->Branch("nbtag_jes_up",&nbtag_jes_up,"nbtag_jes_up/i");
  tuple->Branch("nbtag_jes_down",&nbtag_jes_down,"nbtag_jes_down/i");
  tuple->Branch("nbtag_jer_up",&nbtag_jer_up,"nbtag_jer_up/i");
  tuple->Branch("nbtag_jer_down",&nbtag_jer_down,"nbtag_jer_down/i");

  tuple->Branch("ptjet",ptjet,"ptjet[njets]/F");
  tuple->Branch("etajet",etajet,"etajet[njets]/F");
  tuple->Branch("phijet",phijet,"phijet[njets]/F");
  tuple->Branch("ejet",ejet,"ejet[njets]/F");
  tuple->Branch("flavorjet",flavorjet,"flavorjet[njets]/I");
  tuple->Branch("tagjet",tagjet,"tagjet[njets]/O");

  TString filen;
  int iFiles = 0;
  int events = 0;
  int counterFiles = 0;
  while (fileListX >> filen)
    counterFiles++;

  while (fileList >> filen) {
   iFiles++;
   cout << "file " << iFiles << " out of " << counterFiles << " : " << filen << endl;
   
   TFile * file_ = TFile::Open(TString(filen));
   if (file_==NULL) continue;

   TTree * tree_ = (TTree*)file_->Get(TString(chainName));
   
   if (tree_==NULL) continue;

   tree_->SetMaxVirtualSize(3000000);

   tree_->SetBranchAddress("genweight",&genweight);
   tree_->SetBranchAddress("genparticles_count", &genparticles_count);
   tree_->SetBranchAddress("genparticles_e", genparticles_e);
   tree_->SetBranchAddress("genparticles_px", genparticles_px);
   tree_->SetBranchAddress("genparticles_py", genparticles_py);
   tree_->SetBranchAddress("genparticles_pz", genparticles_pz);
   tree_->SetBranchAddress("genparticles_pdgid", genparticles_pdgid);
   tree_->SetBranchAddress("genparticles_status", genparticles_status);
   tree_->SetBranchAddress("genparticles_fromHardProcess", genparticles_fromHardProcess);
  
   tree_->SetBranchAddress("gentau_count", &gentau_count);
   tree_->SetBranchAddress("gentau_charge", gentau_charge);
   tree_->SetBranchAddress("gentau_visible_e", gentau_visible_e);
   tree_->SetBranchAddress("gentau_visible_px", gentau_visible_px);
   tree_->SetBranchAddress("gentau_visible_py", gentau_visible_py);
   tree_->SetBranchAddress("gentau_visible_pz", gentau_visible_pz);
   tree_->SetBranchAddress("gentau_fromHardProcess", gentau_fromHardProcess);
   tree_->SetBranchAddress("gentau_isLastCopy", gentau_isLastCopy);
   tree_->SetBranchAddress("gentau_decayMode", gentau_decayMode);  

   tree_->SetBranchAddress("pfjet_count",&pfjet_count);
   tree_->SetBranchAddress("pfjet_e",pfjet_e);
   tree_->SetBranchAddress("pfjet_pt",pfjet_pt);
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

   tree_->SetBranchAddress("genjets_count", &genjets_count);
   tree_->SetBranchAddress("genjets_pt", genjets_pt);
   tree_->SetBranchAddress("genjets_eta", genjets_eta);
   tree_->SetBranchAddress("genjets_phi", genjets_phi);

   tree_->SetBranchAddress("numtruepileupinteractions",&numtruepileupinteractions);

   // Additional trigger objects
   tree_->SetBranchAddress("run_hltfilters",&hltfilters);
   tree_->SetBranchAddress("run_btagdiscriminators", &btagdiscriminators);
   tree_->SetBranchAddress("hltriggerresults",&hltriggerresults);
   tree_->SetBranchAddress("hltriggerprescales",&hltriggerprescales);

   int numberOfCandidates = tree_->GetEntries();

   std::cout << "number of events = " << numberOfCandidates << std::endl;
   
   for (int iCand=0; iCand<numberOfCandidates; iCand++) {
     
     tree_->GetEntry(iCand);

     events++;
     if (events%10000==0) cout << "   processed events : " << events << endl;

     mcweight = 1.0;
     if (genweight<0) mcweight = -1.0;
     weight = mcweight;

     puweight = float(PUofficial->get_PUweight(double(numtruepileupinteractions)));
     weight *= puweight;

     ptPosLep = -1.0;
     ptNegLep = -1.0;

     for (unsigned int igen=0; igen<genparticles_count; ++igen) {
       //       std::cout << "genparticle " << igen
       //		 << " status = " << genparticles_status[igen]
       //		 << " fromHardProcess = " << genparticles_fromHardProcess[igen] 
       //		 << " pdgid = " << genparticles_pdgid[igen] 
       //		 << " info " << genparticles_info[igen] << std::endl;
       bool accept = genparticles_status[igen]==1;
       accept = accept && genparticles_fromHardProcess[igen]==1;
       if (!accept) continue;

       if (genparticles_pdgid[igen]==11||genparticles_pdgid[igen]==13) {
	 TLorentzVector genLV;
	 genLV.SetXYZT(genparticles_px[igen],
		       genparticles_py[igen],
		       genparticles_pz[igen],
		       genparticles_e[igen]);
	 if (genLV.Pt()>ptNegLep) {
	   ptNegLep = genLV.Pt();
	   etaNegLep = genLV.Eta();
	   phiNegLep = genLV.Phi();
	   pdgNegLep = genparticles_pdgid[igen];
	 }
       }

       if (genparticles_pdgid[igen]==-11||genparticles_pdgid[igen]==-13) {
	 TLorentzVector genLV;
	 genLV.SetXYZT(genparticles_px[igen],
		       genparticles_py[igen],
		       genparticles_pz[igen],
		       genparticles_e[igen]);
	 if (genLV.Pt()>ptPosLep) {
	   ptPosLep = genLV.Pt();
	   etaPosLep = genLV.Eta();
	   phiPosLep = genLV.Phi();
	   pdgPosLep = genparticles_pdgid[igen];
	 }
       }

     }

     //     std::cout << "ptPosLep = " << ptPosLep << std::endl;
     //     std::cout << "ptNegLep = " << ptNegLep << std::endl;

     if (ptPosLep<ptLepCut) continue;
     if (abs(etaPosLep)>etaLepCut) continue;
 
     if (ptNegLep<ptLepCut) continue;
     if (abs(etaNegLep)>etaLepCut) continue;
 
     bool isZpair = (pdgPosLep==-11 && pdgNegLep==11) || (pdgPosLep==-13 && pdgNegLep==13);
     if (!isZpair) continue;
       
     // ****************************
     // ** selection of tau pair  **
     // ****************************
     ptPosTau = -1.0;
     ptNegTau = -1.0;

     for (unsigned int igen=0; igen<gentau_count; ++igen) {
       //       std::cout << "gentau " << igen 
       //                 << " fromHardProcess = " << gentau_fromHardProcess[igen]	
       //                 << "  lastCopy = " << gentau_isLastCopy[igen] << std::endl;
       bool accept = gentau_fromHardProcess[igen]==1 && gentau_isLastCopy[igen]==1;
       if (!accept) continue;
       TLorentzVector genLV;
       genLV.SetXYZT(gentau_visible_px[igen],
		     gentau_visible_py[igen],
		     gentau_visible_pz[igen],
		     gentau_visible_e[igen]);
       
       float dR = deltaR(etaPosLep,phiPosLep,
			 genLV.Eta(),genLV.Phi());
       if (dR<0.4) continue;
       dR = deltaR(etaNegLep,phiNegLep,
		   genLV.Eta(),genLV.Phi());
       if (dR<0.4) continue;


       if (gentau_charge[igen]<0.) {
	 if (genLV.Pt()>ptNegTau) {
	   ptNegTau = genLV.Pt();
	   etaNegTau = genLV.Eta();
	   phiNegTau = genLV.Phi();
	   dmNegTau = gentau_decayMode[igen];
	 }
       }

       if (gentau_charge[igen]>0.) {
	 if (genLV.Pt()>ptPosTau) {
	   ptPosTau = genLV.Pt();
	   etaPosTau = genLV.Eta();
	   phiPosTau = genLV.Phi();
	   dmPosTau = gentau_decayMode[igen];
	 }
       }

     }

     //     std::cout << "ptPosTau = " << ptPosTau << std::endl;
     //     std::cout << "ptNegTau = " << ptNegTau << std::endl;
     //     std::cout << std::endl;
     

     float ptMaxTau = TMath::Max(ptPosTau,ptNegTau);
     float ptMinTau = TMath::Min(ptPosTau,ptNegTau);
     if (ptMaxTau<ptTauMaxCut) continue;
     if (ptMinTau<ptTauMinCut) continue;

     if (abs(etaPosTau)>etaTauCut) continue; 
     if (abs(etaNegTau)>etaTauCut) continue;
 

     int nBTagDiscriminant1 = -1;
     int nBTagDiscriminant2 = -1;
     int nBTagDiscriminant3 = -1;
     unsigned int num_btags = 0;
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

     njets = 0;
     nbtag = 0;
     nbtag_jes_up = 0;
     nbtag_jes_down = 0;
     nbtag_jer_up = 0;
     nbtag_jer_down = 0;
     njetspt20 = 0;

     btagweight = 1;
     btagweight_lf_up = 1;
     btagweight_lf_down = 1;
     btagweight_hf_up = 1;
     btagweight_hf_down = 1;
     
     for (unsigned int jet=0; jet<pfjet_count; ++jet) {
	 
       float absJetEta = TMath::Abs(pfjet_eta[jet]);
       float jetPt = pfjet_pt[jet];
       if (absJetEta>jetEtaCut) continue;

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

       float dR = deltaR(etaPosLep,phiPosLep,
			 pfjet_eta[jet],pfjet_phi[jet]);
       if (dR<dRJetLep) continue;

       dR = deltaR(etaNegLep,phiNegLep,
		   pfjet_eta[jet],pfjet_phi[jet]);
       if (dR<dRJetLep) continue;

       dR = deltaR(etaPosTau,phiPosTau,
		   pfjet_eta[jet],pfjet_phi[jet]);
       if (dR<dRJetLep) continue;

       dR = deltaR(etaNegTau,phiNegTau,
		   pfjet_eta[jet],pfjet_phi[jet]);
       if (dR<dRJetLep) continue;
       
       float genJetPt = jetPt;
       float dRmin = 0.4;
       for (unsigned int genjet=0; genjet<genjets_count; ++genjet) {
	 float dRjets = deltaR(pfjet_eta[jet],pfjet_phi[jet],
			       genjets_eta[genjet],genjets_phi[genjet]);
	 if (dRjets<dRmin) {
	   dRmin = dRjets;
	   genJetPt = genjets_pt[genjet];
	 }
       }

       bool tagged = false;

       if (absJetEta<bjetEtaCut) {

	 njetspt20++;
	 
	 float btagDiscr = pfjet_btag[jet][nBTagDiscriminant1];
	 if (BTagAlgorithm=="pfDeepFlavourJetTags"||BTagAlgorithm=="pfDeepCSVJetTags")
	   btagDiscr += pfjet_btag[jet][nBTagDiscriminant2];
	 if (BTagAlgorithm=="pfDeepFlavourJetTags")
	   btagDiscr += pfjet_btag[jet][nBTagDiscriminant3];
	 
	 tagged = btagDiscr>btagCut;
	 if (tagged) {
	   if (jetPt>jetPtCut) nbtag++;
	 }

       }
       if (jetPt>jetPtCut) {
	 ptjet[njets] = pfjet_pt[jet];
	 etajet[njets] = pfjet_eta[jet];
	 phijet[njets] = pfjet_phi[jet];
	 ejet[njets] = pfjet_e[jet];
	 flavorjet[njets] = pfjet_flavour[jet];
	 tagjet[njets] = tagged;
	 njets++;
       }
     }
     
     tuple->Fill();
     
   } // icand loop
   
   delete tree_;
   file_->Close();
   delete file_;
   
  }// filelist loop
  
  file->cd("");
  file->Write();
  file->Close();
  
  //delete file;
} // int main loop 

 

