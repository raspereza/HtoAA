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
#include "TMath.h"

int main(int argc, char * argv[]) {


  if (argc!=4) {
    std::cout << "Usage of the program : createQCDModel [folder] [era] [bin]" << std::endl;
    exit(1);
  }

  TString folder(argv[1]);
  TString era(argv[2]);
  TString bin(argv[3]);
  TString cmsswBase = TString(getenv("CMSSW_BASE"));

  std::map<TString,double> eraLumi = {
    {"2016_preVFP" ,19520},
    {"2016_postVFP",16810},
    {"2017",        41480},
    {"2018",        59830}
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
    std::cout << "Available options : 2016_preVFP, 2016_postVFP, 2017, 2018" << std::endl;
    exit(-1);
  }

  TString basefolder = folder+"/"+era+"/"+bin;
  double lumi = eraLumi[era];

  std::vector<TString> sampleNames = {
    "QCD_Pt-20To30_MuEnrichedPt5",
    "QCD_Pt-30To50_MuEnrichedPt5",
    "QCD_Pt-50To80_MuEnrichedPt5",
    "QCD_Pt-80To120_MuEnrichedPt5",
    "QCD_Pt-120To170_MuEnrichedPt5",
    "QCD_Pt-170To300_MuEnrichedPt5",
    "QCD_Pt-300To470_MuEnrichedPt5",
    "QCD_Pt-470To600_MuEnrichedPt5",
    "QCD_Pt-600To800_MuEnrichedPt5",
    "QCD_Pt-800To1000_MuEnrichedPt5"
  };

  std::map<TString,double> sample_xsec = {
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

  TH1D * probMu[nFlav][nNetQ][nPartMom][nMuMom][nReg]; 
  TH1D * probMu_SS[nFlav][nNetQ][nPartMom][nMuMom][nReg]; 
  
  TH1D * normMu[nFlav][nNetQ][nPartMom][nMuMom]; 
  TH1D * normMu_SS[nFlav][nNetQ][nPartMom][nMuMom]; 
  
  TH1D * Mass[nFlav][nNetQ][nPartMom][nMuMom][nReg];
  TH1D * Mass_SS[nFlav][nNetQ][nPartMom][nMuMom][nReg];

  for (unsigned int iF=0; iF<nFlav; ++iF) {
    for (unsigned int iQ=0; iQ<nNetQ; ++iQ) {
      for (unsigned int iMom=0; iMom<nPartMom; ++iMom) {
	for (unsigned int mu=0; mu<nMuMom; ++mu) {

	  TString name  = 
	    partonFlavor[iF] + "_" + 
	    muonPartonNetCharge[iQ] + "_" +
	    partonMomRange[iMom] + "_" +
	    muonMomRange[mu];
	  
	  // probability to pass selection
	  TString histName = "normMu_" + name + "_SS";
	  normMu_SS[iF][iQ][iMom][mu] = new TH1D(histName,"",1,0.,1.);
	  
	  histName = "normMu_" + name;
	  normMu[iF][iQ][iMom][mu] = new TH1D(histName,"",1,0.,1.);

	  for (unsigned int iR=0; iR<nReg; ++iR) {

	    name  = 
	      partonFlavor[iF] + "_" + 
	      muonPartonNetCharge[iQ] + "_" +
	      partonMomRange[iMom] + "_" +
	      muonMomRange[mu] + "_" +
	      Regions[iR] ;
	    
	    // mass distributions
	    histName = "Mass_" + name + "_SS";
	    Mass_SS[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",nBins,xMin,xMax);
			    
	    histName = "Mass_" + name;
	    Mass[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",nBins,xMin,xMax);
	    
	    // probability to pass selection
	    histName = "passMu_" + name + "_SS";
	    probMu_SS[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",1,0.,1.);
	    
	    histName = "passMu_" + name;
	    probMu[iF][iQ][iMom][mu][iR] = new TH1D(histName,"",1,0.,1.);

	  }
	}
      }
    }
  }

  std::cout << std::endl;
  // run over samples
  for (auto samplename : sampleNames) {
    TString filename = basefolder + "/" + samplename + ".root";
    TFile * file = new TFile(filename);
    if (file->IsZombie()) {
      std::cout << "file " << filename << " does not exist. Quitting" << std::endl;
      exit(-1);
    }
    TH1D * histweights = (TH1D*)file->Get("histWeightsH");
    double nevt = histweights->GetSumOfWeights();
    std::cout << "processing sample " << samplename << " : " << int(nevt) << std::endl;
    double xsec = sample_xsec[samplename];
    double norm = xsec*lumi/nevt;
 
    TH1D * partonMu = (TH1D*)file->Get("partonMu");
    partonMu->Scale(norm);
    denomMu->Add(denomMu,partonMu);

    TH1D * partonMu_SS = (TH1D*)file->Get("partonMu_SS");
    partonMu_SS->Scale(norm);
    denomMu_SS->Add(denomMu_SS,partonMu_SS);

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

	    TH1D * partonMuProbe_SS = (TH1D*)file->Get(histName);
	    partonMuProbe_SS->Scale(norm);
	    normMu_SS[iF][iQ][iMom][mu]->Add(normMu_SS[iF][iQ][iMom][mu],partonMuProbe_SS);

	    histName = "partonMuProbe_" + name;
	    TH1D * partonMuProbe = (TH1D*)file->Get(histName);	    
	    partonMuProbe->Scale(norm);
	    normMu[iF][iQ][iMom][mu]->Add(normMu[iF][iQ][iMom][mu],partonMuProbe);

	    for (unsigned int iR=0; iR<nReg; ++iR) {
	      name  = 
		partonFlavor[iF] + "_" + 
		muonPartonNetCharge[iQ] + "_" +
		partonMomRange[iMom] + "_" +
		muonMomRange[mu] + "_" +
		Regions[iR] ;
	      
	      // mass distributions
	      histName = "MuTrkMass_" + name + "_SS";
	      TH1D * MassMuTrk_SS = (TH1D*)file->Get(histName);
	      MassMuTrk_SS->Scale(norm);
	      Mass_SS[iF][iQ][iMom][mu][iR]->Add(Mass_SS[iF][iQ][iMom][mu][iR],MassMuTrk_SS);
      
	      histName = "MuTrkMass_" + name;
	      TH1D * MassMuTrk = (TH1D*)file->Get(histName);
	      MassMuTrk->Scale(norm);
	      Mass[iF][iQ][iMom][mu][iR]->Add(Mass[iF][iQ][iMom][mu][iR],MassMuTrk);
	    
	      // probability to pass selection
	      histName = "partonMuPass_" + name + "_SS";
	      TH1D * partonMuPass_SS = (TH1D*)file->Get(histName);
	      partonMuPass_SS->Scale(norm);
	      probMu_SS[iF][iQ][iMom][mu][iR]->Add(probMu_SS[iF][iQ][iMom][mu][iR],partonMuPass_SS);
	      
	      histName = "partonMuPass_" + name;
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
  //  double scale = 1.0/denomMu->GetSumOfWeights();
  //  double scale_SS = 0.5/denomMu_SS->GetSumOfWeights(); // filled only once in SS region

  TString outputname = cmsswBase + "/src/HtoAA/data/QCDModel_"+era+".root";
  TFile * outputfile = new TFile(outputname,"recreate");
  outputfile->cd("");
  denomMu->Write(denomMu->GetName());
  denomMu_SS->Write(denomMu_SS->GetName());

  //  normalizing distributions
  std::cout << std::endl;
  for (unsigned int iF=0; iF<nFlav; ++iF) {
    for (unsigned int iQ=0; iQ<nNetQ; ++iQ) {
      for (unsigned int iMom=0; iMom<nPartMom; ++iMom) {
	for (unsigned int mu=0; mu<nMuMom; ++mu) {
	  outputfile->cd("");
	  normMu_SS[iF][iQ][iMom][mu]->Write(normMu_SS[iF][iQ][iMom][mu]->GetName());
	  normMu[iF][iQ][iMom][mu]->Write(normMu[iF][iQ][iMom][mu]->GetName());

	  for (unsigned int iR=0; iR<nReg; ++iR) {

	    outputfile->cd("");

	    double num = probMu_SS[iF][iQ][iMom][mu][iR]->GetBinContent(1);
	    double den = normMu_SS[iF][iQ][iMom][mu]->GetBinContent(1);
	    double numE = probMu_SS[iF][iQ][iMom][mu][iR]->GetBinError(1);
	    double denE = normMu_SS[iF][iQ][iMom][mu]->GetBinError(1);

	    double denR = 0;
	    double numR = 0;
	    double prob_SS = 0;
	    double probE_SS = 0;
	    double probR_SS = 0;

	    if (num>0) numR = numE/num;
	    if (den>0) {
	      denR = denE/den;
	      prob_SS = num/den;
	      probR_SS = TMath::Sqrt(numR*numR+denR*denR);
	      probE_SS = prob_SS*probR_SS;
	    }
	    probMu_SS[iF][iQ][iMom][mu][iR]->Write(probMu_SS[iF][iQ][iMom][mu][iR]->GetName());

	    num = probMu[iF][iQ][iMom][mu][iR]->GetBinContent(1);
            den = normMu[iF][iQ][iMom][mu]->GetBinContent(1);
	    numE = probMu[iF][iQ][iMom][mu][iR]->GetBinError(1);
            denE = normMu[iF][iQ][iMom][mu]->GetBinError(1);
	    numR = 0;
	    denR = 0;
	    double prob = 0;
	    double probE = 0;
	    double probR = 0;

	    if (num>0) numR = numE/num;
	    if (den>0) {
	      denR = denE/den;
	      prob = num/den;
	      probR = TMath::Sqrt(numR*numR+denR*denR);
	      probE = prob*probR;
	    }
	    probMu[iF][iQ][iMom][mu][iR]->Write(probMu[iF][iQ][iMom][mu][iR]->GetName());

	    std::cout << partonFlavor[iF] << ":"
		      << muonPartonNetCharge[iQ] << ":"
		      << partonMomRange[iMom] << ":"
		      << muonMomRange[mu] << ":"
		      << Regions[iR] << " -> ";
	    printf("%5.3f+/-%5.3f : %5.3f+/-%5.3f\n",prob,probE,prob_SS,probE_SS);	    

	    double Norm = Mass[iF][iQ][iMom][mu][iR]->GetSumOfWeights();
	    if (Norm>0) Mass[iF][iQ][iMom][mu][iR]->Scale(1./Norm);
	    Mass[iF][iQ][iMom][mu][iR]->Write(Mass[iF][iQ][iMom][mu][iR]->GetName());

	    double Norm_SS = Mass_SS[iF][iQ][iMom][mu][iR]->GetSumOfWeights();
	    if (Norm_SS>0) Mass_SS[iF][iQ][iMom][mu][iR]->Scale(1./Norm_SS);
	    Mass_SS[iF][iQ][iMom][mu][iR]->Write(Mass_SS[iF][iQ][iMom][mu][iR]->GetName());


	  }
	}
      }
      std::cout << std::endl;
    }
  }

  outputfile->Close();
  delete outputfile;
  
}
