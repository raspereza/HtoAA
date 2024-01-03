#include "HtoAA/Utilities/interface/QCDModel.h"

QCDModel::QCDModel(TString fileName) {

  file = new TFile(fileName);
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
	    pdfMassSS[iF][iQ][iMom][mu][iR] = (TH1D*)file->Get(histName);
      
	    histName = "Mass_" + name;
	    pdfMass[iF][iQ][iMom][mu][iR] = (TH1D*)file->Get(histName);
	    
	    // probability to pass selection
	    histName = "probMu_" + name + "_SS";
	    probSS[iF][iQ][iMom][mu][iR] = (TH1D*)file->Get(histName);
	    
	    histName = "probMu_" + name;
	    prob[iF][iQ][iMom][mu][iR] = (TH1D*)file->Get(histName);
	  }
	}
      }
    }
  }

}

QCDModel::~QCDModel() {
  delete file;
}

double QCDModel::getProb(int pMom, int muMom, int region, int flav, int qnet, bool inclusive) {

  if (pMom<0) pMom = 0;
  if (pMom>3) pMom = 3;
  if (muMom<0) muMom = 0;
  if (muMom>3) muMom = 3;
  if (region<0) region = 0;
  if (region>1) region = 1;
  if (flav<0) flav = 0;
  if (flav>4) flav = 4;
  if (qnet<0) qnet = 0;
  if (qnet>1) qnet = 1;
  
  double output = 1.0;
  if (inclusive) 
    output = prob[flav][qnet][pMom][muMom][region]->GetBinContent(1);
  else
    output = probSS[flav][qnet][pMom][muMom][region]->GetBinContent(1);

  return output;

}

double QCDModel::getMassPdf(int pMom, int muMom, int region, int flav, int qnet, double mass, bool inclusive) {
  
  if (pMom<0) pMom = 0;
  if (pMom>3) pMom = 3;
  if (muMom<0) muMom = 0;
  if (muMom>3) muMom = 3;
  if (region<0) region = 0;
  if (region>1) region = 1;
  if (flav<0) flav = 0;
  if (flav>4) flav = 4;
  if (qnet<0) qnet = 0;
  if (qnet>1) qnet = 1;

  double output = 1.0;
  if (inclusive) {
    int bin = pdfMass[flav][qnet][pMom][muMom][region]->FindBin(mass);
    output  = pdfMass[flav][qnet][pMom][muMom][region]->GetBinContent(bin);
  }
  else {
    int bin = pdfMassSS[flav][qnet][pMom][muMom][region]->FindBin(mass);
    output  = pdfMassSS[flav][qnet][pMom][muMom][region]->GetBinContent(bin);
  }
  return output;

}

