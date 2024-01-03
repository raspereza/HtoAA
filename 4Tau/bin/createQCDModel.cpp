#include "HtoAA/Utilities/interface/QCDModelDefinitions.h"
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
#include "TH1.h"
#include "TString.h"
#include "TFile.h"

int main(int argc, char * argv[]) {

  if (argc!=3) {
    std::cout << "Usage of the program : createQCDModel [era] [ptbin]" << std::endl;
    exit(1);
  }

  TString era(argv[1]);
  TString ptbin(argv[2]);
  TString cmsswBase = TString(getenv("CMSSW_BASE"));

  std::vector<TString> ptbins = {
    "pt2p5"
  };
  
  bool ptbinNotFound = true;
  for (auto bin : ptbins) {
    if (ptbin==bin) {
      ptbinNotFound = false;
      break;
    }
  }
  
  if (ptbinNotFound) {
    std::cout << "Uknown ptbin specfied : " << ptbin << std::endl;
    std::cout << "Available options : ";
    for (auto bin : ptbins)
      std::cout << bin << " ";
    std::cout << std::endl;
    exit(-1);
  }

  std::map<TString,double> eraLumi = {
    {"2016_preVFP_UL", 19520},
    {"2016_postVFP_UL",16810},
    {"2017_UL",        41480},
    {"2018_UL",        59830}
  };

  bool eraNotFound = true;
  for (auto lumi : eraLumi) {
    TString ERA = lumi.first;
    if (ERA==era) {
      eraNotFound = false;
      break;
    }
  }

  if (eraNotFound) {
    std::cout << "Uknown era specfied : " << era << std::endl;
    std::cout << "Available options : 2016_preVFP_UL, 2016_postVFP_UL, 2017_UL, 2018_UL" << std::endl;
    exit(-1);
  }

  double lumi = eraLumi[era];

  TString basefolder = "/nfs/dust/cms/user/rasp/Run/Run"+era+"/mutrk/"+ptbin;

  std::map<TString,double> samples = {
    {"QCD_Pt-20To30_MuEnrichedPt5",  558528000*0.00530},
    {"QCD_Pt-30To50_MuEnrichedPt5",  139803000*0.01182},
    {"QCD_Pt-50To80_MuEnrichedPt5",   19222500*0.02276},
    {"QCD_Pt-80To120_MuEnrichedPt5",   2758420*0.03844},
    {"QCD_Pt-120To170_MuEnrichedPt5",   469797*0.05362},
    {"QCD_Pt-170To300_MuEnrichedPt5",   117989*0.07335},
    {"QCD_Pt-300To470_MuEnrichedPt5",   7820.3*0.10196},
    {"QCD_Pt-470To600_MuEnrichedPt5",   645.52*0.12242},
    {"QCD_Pt-600To800_MuEnrichedPt5",  187.109*0.13412},
    {"QCD_Pt-800To1000_MuEnrichedPt5", 32.3486*0.14552}
  };

  TH1D * denomMu = new TH1D("denomMu","",1,0.,1.);
  TH1D * denomMu_SS = new TH1D("denomMu_SS","",1,0.,1.);

  TH1D * probMu[5][2][4][4][2]; 
  TH1D * probMu_SS[5][2][4][4][2]; 
  
  TH1D * Mass[5][2][4][4][2];
  TH1D * Mass_SS[5][2][4][4][2];

  for (unsigned int iF=0; iF<nPartonFlavours; ++iF) {
    for (unsigned int iQ=0; iQ<nMuonPartonNetCharge; ++iQ) {
      for (unsigned int iMom=0; iMom<nPartonMomBins; ++iMom) {
	for (unsigned int mu=0; mu<nMuonMomBins; ++mu) {
	  for (unsigned int iR=0; iR<nRegions; ++iR) {

	    TString name  = 
	      partonFlavor[iF] + "_" + 
	      muonPartonNetCharge[iQ] + "_" +
	      partonMomRange[iMom] + "_" +
	      muonMomRange[mu] + "_" +
	      Regions[iR] ;
	    
	    // mass distributions
	    TString histName = "Mass_" + name + "_SS";
	    Mass_SS[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",nBins,xMin,xMax);
			    
	    histName = "Mass_" + name;
	    Mass[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",nBins,xMin,xMax);
	    
	    // probability to pass selection
	    histName = "probMu_" + name + "_SS";
	    probMu_SS[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",1,0.,1.);
	    
	    histName = "probMu_" + name;
	    probMu[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",1,0.,1.);
	  }
	}
      }
    }
  }

  std::cout << std::endl;
  // run over samples
  for (auto sample : samples) {
    TString samplename = sample.first;
    TString filename = basefolder + "/" + samplename + ".root";
    TFile * file = new TFile(filename);
    if (file->IsZombie()) {
      std::cout << "file " << filename << " does not exist. Quitting" << std::endl;
      exit(-1);
    }
    std::cout << "processing sample " << samplename << std::endl;
    TH1D * histweights = (TH1D*)file->Get("histWeightsH");
    double nevt = histweights->GetSumOfWeights();
    double xsec = sample.second;
    double norm = xsec*lumi/nevt;
 
    TH1D * partonMu = (TH1D*)file->Get("partonMu");
    partonMu->Scale(norm);
    denomMu->Add(denomMu,partonMu);

    TH1D * partonMu_SS = (TH1D*)file->Get("partonMuSS");
    partonMu_SS->Scale(norm);
    denomMu_SS->Add(denomMu_SS,partonMu_SS);

    for (unsigned int iF=0; iF<nPartonFlavours; ++iF) {
      for (unsigned int iQ=0; iQ<nMuonPartonNetCharge; ++iQ) {
	for (unsigned int iMom=0; iMom<nPartonMomBins; ++iMom) {
	  for (unsigned int mu=0; mu<nMuonMomBins; ++mu) {
	    for (unsigned int iR=0; iR<nRegions; ++iR) {
	      TString name  = 
		partonFlavor[iF] + "_" + 
		muonPartonNetCharge[iQ] + "_" +
		partonMomRange[iMom] + "_" +
		muonMomRange[mu] + "_" +
		Regions[iR] ;
	      
	      // mass distributions
	      TString histName = "MuTrkMass_" + name + "_SS";
	      TH1D * MassMuTrk_SS = (TH1D*)file->Get(histName);
	      MassMuTrk_SS->Scale(norm);
	      Mass_SS[iF][iQ][iMom][mu][iR]->Add(Mass_SS[iF][iQ][iMom][mu][iR],MassMuTrk_SS);
      
	      histName = "MuTrkMass_" + name;
	      TH1D * MassMuTrk = (TH1D*)file->Get(histName);
	      MassMuTrk->Scale(norm);
	      Mass[iF][iQ][iMom][mu][iR]->Add(Mass[iF][iQ][iMom][mu][iR],MassMuTrk);
	    
	      // probability to pass selection
	      histName = "partonMu_" + name + "_SS_passed";
	      TH1D * partonMuPass_SS = (TH1D*)file->Get(histName);
	      partonMuPass_SS->Scale(norm);
	      probMu_SS[iF][iQ][iMom][mu][iR]->Add(probMu_SS[iF][iQ][iMom][mu][iR],partonMuPass_SS);
	      
	      histName = "partonMu_" + name + "_passed";
	      TH1D * partonMuPass = (TH1D*)file->Get(histName);
	      partonMuPass->Scale(norm);
	      probMu[iF][iQ][iMom][mu][iR]->Add(probMu[iF][iQ][iMom][mu][iR],partonMuPass);
	      
	    }
	  }
	}
      }
    }

  }

  // denominator
  double scale = 1.0/denomMu->GetSumOfWeights();
  double scale_SS = 0.5/denomMu_SS->GetSumOfWeights(); // filled only once in SS region

  TString outputname = cmsswBase + "/src/HtoAA/data/QCDModel_"+era+"_"+ptbin+".root";
  TFile * outputfile = new TFile(outputname,"recreate");
  outputfile->cd("");

  //  normalizing distributions
  std::cout << std::endl;
  for (unsigned int iF=0; iF<nPartonFlavours; ++iF) {
    for (unsigned int iQ=0; iQ<nMuonPartonNetCharge; ++iQ) {
      for (unsigned int iMom=0; iMom<nPartonMomBins; ++iMom) {
	for (unsigned int mu=0; mu<nMuonMomBins; ++mu) {
	  for (unsigned int iR=0; iR<nRegions; ++iR) {

	    outputfile->cd("");

	    probMu_SS[iF][iQ][iMom][mu][iR]->Scale(scale_SS);
	    probMu_SS[iF][iQ][iMom][mu][iR]->Write(probMu_SS[iF][iQ][iMom][mu][iR]->GetName());
	    double prob_SS = probMu_SS[iF][iQ][iMom][mu][iR]->GetSumOfWeights();

	    probMu[iF][iQ][iMom][mu][iR]->Scale(scale);
	    probMu[iF][iQ][iMom][mu][iR]->Write(probMu[iF][iQ][iMom][mu][iR]->GetName());
	    double prob = probMu[iF][iQ][iMom][mu][iR]->GetSumOfWeights();

	    std::cout << partonFlavor[iF] << ":"
		      << muonPartonNetCharge[iQ] << ":"
		      << partonMomRange[iMom] << ":"
		      << muonMomRange[mu] << ":"
		      << Regions[iR] << "->  prob = " << prob << "  probSS = " << prob_SS << std::endl;
	    
	    double Norm = Mass[iF][iQ][iMom][mu][iR]->GetSumOfWeights();
	    if (Norm>0) Mass[iF][iQ][iMom][mu][iR]->Scale(1./Norm);
	    Mass[iF][iQ][iMom][mu][iR]->Write(Mass[iF][iQ][iMom][mu][iR]->GetName());

	    double Norm_SS = Mass_SS[iF][iQ][iMom][mu][iR]->GetSumOfWeights();
	    if (Norm_SS>0) Mass_SS[iF][iQ][iMom][mu][iR]->Scale(1./Norm_SS);
	    Mass_SS[iF][iQ][iMom][mu][iR]->Write(Mass_SS[iF][iQ][iMom][mu][iR]->GetName());


	  }
	}
      }
    }
  }

  outputfile->Close();
  delete outputfile;
  
}
