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
//#include "HtoAA/Utilities/src/Config.cc"
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
  
  if (argc<2) {
    std::cout << "Usage of the program : btag_eff [file_list]" << std::endl;
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
  const float bjetEta = cfg.get<float>("BJetEta");
  const float bjetPt = cfg.get<float>("BJetPt");

  TString BTagAlgorithm(bTagAlgorithm);
  TString BTagDiscriminator1(bTagDiscriminator1);
  TString BTagDiscriminator2(bTagDiscriminator2);
  TString BTagDiscriminator3(bTagDiscriminator3);

  const string pileUpDataFile = cfg.get<string>("PileUpDataFileName");
  const string pileUpMCFile = cfg.get<string>("PileUpMCFileName");

  TString PileUpDataFile(pileUpDataFile);
  TString PileUpMCFile(pileUpMCFile);

  // ********** end of configuration *******************

  std::ifstream fileList(argv[2]);
  std::ifstream fileListX(argv[2]);

  Float_t genweight;
  Float_t numtruepileupinteractions;

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
  
  int n_pt_bins = 11;
  float pt_bins[12] ={20,30,40,50,60,80,100,150,200,300,500,1000};
  int n_eta_bins = 4;
  float eta_bins[5] = {0.,0.9,1.2,2.1,2.4};
  
  TFile * file = new TFile(FullName+TString(".root"),"recreate");
  file->cd("");
  TH2D * b_pass = new TH2D("b_pass","b_pass",n_pt_bins,pt_bins,n_eta_bins,eta_bins); 
  TH2D * c_pass = new TH2D("c_pass","c_pass",n_pt_bins,pt_bins,n_eta_bins,eta_bins); 
  TH2D * l_pass = new TH2D("oth_pass","oth_pass",n_pt_bins,pt_bins,n_eta_bins,eta_bins); 

  TH2D * b_all = new TH2D("b_all","b_all",n_pt_bins,pt_bins,n_eta_bins,eta_bins); 
  TH2D * c_all = new TH2D("c_all","c_all",n_pt_bins,pt_bins,n_eta_bins,eta_bins); 
  TH2D * l_all = new TH2D("oth_all","oth_all",n_pt_bins,pt_bins,n_eta_bins,eta_bins); 

  // PU reweighting
  PileUp * PUofficial = new PileUp();
  TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/HtoAA/data/PileUpDistrib/"+PileUpDataFile,"read");
  TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/HtoAA/data/PileUpDistrib/"+PileUpMCFile, "read");
  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  float MaxBJetPt = 1000.;
  float MinBJetPt = 20.;
  float MaxBJetEta = 2.4;
  float MinBJetEta = 0.0;

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

   tree_->SetBranchAddress("genweight",&genweight);
   tree_->SetBranchAddress("numtruepileupinteractions",&numtruepileupinteractions);

   // Additional trigger objects
   tree_->SetBranchAddress("run_hltfilters",&hltfilters);
   tree_->SetBranchAddress("run_btagdiscriminators", &btagdiscriminators);
   tree_->SetBranchAddress("hltriggerresults",&hltriggerresults);
   tree_->SetBranchAddress("hltriggerprescales",&hltriggerprescales);

   int numberOfCandidates = tree_->GetEntries();

   std::cout << "number of events = " << numberOfCandidates << std::endl;
   
   TRandom3 rand;

   for (int iCand=0; iCand<numberOfCandidates; iCand++) {
     
     tree_->GetEntry(iCand);

     events++;
     if (events%10000==0) cout << "   processed events : " << events << endl;

     float weight = genweight;
     float puweight = float(PUofficial->get_PUweight(double(numtruepileupinteractions)));
     weight *= puweight;

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
     if (isCorrectBTag) {
       for (unsigned int jet=0; jet<pfjet_count; ++jet) {
	 
	 float absEta = TMath::Abs(pfjet_eta[jet]);
	 float JetPtForBTag = pfjet_pt[jet];
	 float JetEtaForBTag = absEta;
	 float jetEta = pfjet_eta[jet];
	   
	 if (absEta>bjetEta) continue;
	 if (pfjet_pt[jet]<bjetPt) continue;

	 if (JetPtForBTag > MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	 if (JetPtForBTag < MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	 if (JetEtaForBTag > MaxBJetEta) JetEtaForBTag = MaxBJetEta - 0.01;
	 if (JetEtaForBTag < MinBJetEta) JetEtaForBTag = MinBJetEta + 0.01;

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

	 float btagDiscr = pfjet_btag[jet][nBTagDiscriminant1];
	 if (BTagAlgorithm=="pfDeepFlavourJetTags"||BTagAlgorithm=="pfDeepCSVJetTags")
	   btagDiscr += pfjet_btag[jet][nBTagDiscriminant2];
	 if (BTagAlgorithm=="pfDeepFlavourJetTags")
	   btagDiscr += pfjet_btag[jet][nBTagDiscriminant3];

	 int flavor = TMath::Abs(pfjet_flavour[jet]);
	 bool tagged = btagDiscr>btagCut;

	 if (flavor==5) {
	   b_all->Fill(JetPtForBTag,JetEtaForBTag,weight);
	   if (tagged) b_pass->Fill(JetPtForBTag,JetEtaForBTag,weight);
	 }
	 else if (flavor==4) {
	   c_all->Fill(JetPtForBTag,JetEtaForBTag,weight);
           if (tagged) c_pass->Fill(JetPtForBTag,JetEtaForBTag,weight);
	 }
	 else {
	   l_all->Fill(JetPtForBTag,JetEtaForBTag,weight);
           if (tagged) l_pass->Fill(JetPtForBTag,JetEtaForBTag,weight);
	 }

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
} // int main loop 

 

