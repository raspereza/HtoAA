#include "HtoAA/Utilities/interface/QCDModel.h"
#include "HtoAA/Utilities/interface/HtoH.h"

QCDModel::QCDModel(TString fileName, vector<double> bins) {

  double xbins[200];
  int nbins = bins.size()-1;
  for (int ib=0; ib<=nbins; ++ib) {
    xbins[ib] = bins[ib];
    
  }

  file = new TFile(fileName);
  std::cout << fileName << std::endl;
  for (unsigned int iF=0; iF<nFlav; ++iF) {
    for (unsigned int iQ=0; iQ<nNetQ; ++iQ) {
      for (unsigned int iMom=0; iMom<nPartMom; ++iMom) {
	for (unsigned int mu=0; mu<nMuMom; ++mu) {

	  TString name  = 
	    partonFlavor[iF] + "_" + 
	    muonPartonNetCharge[iQ] + "_" +
	    partonMomRange[iMom] + "_" +
	    muonMomRange[mu];

	  // denominators
	  TString histName = "normMu_" + name + "_SS";
	  normSS[iF][iQ][iMom][mu] = (TH1D*)file->Get(histName);
	  
	  histName = "normMu_" + name;
	  norm[iF][iQ][iMom][mu] = (TH1D*)file->Get(histName);

	  for (unsigned int iR=0; iR<nReg; ++iR) {
	    name  = 
	      partonFlavor[iF] + "_" + 
	      muonPartonNetCharge[iQ] + "_" +
	      partonMomRange[iMom] + "_" +
	      muonMomRange[mu] + "_" +
	      Regions[iR] ;

	    // mass distributions
	    histName = "Mass_" + name + "_SS";
	    TH1D * hist = (TH1D*)file->Get(histName);
	    if (hist==NULL) {
	      cout << "histogram " << histName << " is absent" << endl;
	      exit(-1);
	    }
	    pdfMassSS[iF][iQ][iMom][mu][iR] = (TH1D*)TH1DtoTH1D(hist,nbins,xbins,true,"_X");
      
	    histName = "Mass_" + name;
	    hist = (TH1D*)file->Get(histName);
	    if (hist==NULL) {
	      cout << "histogram " << histName << " is absent" << endl;
	      exit(-1);
	    }
	    pdfMass[iF][iQ][iMom][mu][iR] = (TH1D*)TH1DtoTH1D(hist,nbins,xbins,true,"_X");
	    
	    // numerator : pass selection
	    histName = "passMu_" + name + "_SS";
	    hist = (TH1D*)file->Get(histName);
	    if (hist==NULL) {
	      cout << "histogram " << histName << " is absent" << endl;
	      exit(-1);
	    }
	    passSS[iF][iQ][iMom][mu][iR] = (TH1D*)file->Get(histName);
	    
	    histName = "passMu_" + name;
	    hist = (TH1D*)file->Get(histName);
	    if (hist==NULL) {
	      cout << "histogram " << histName << " is absent" << endl;
	      exit(-1);
	    }
	    pass[iF][iQ][iMom][mu][iR] = (TH1D*)file->Get(histName);
	  }
	}
      }
    }
  }

  std::cout << "opened" << std::endl;

}

QCDModel::~QCDModel() {
  delete file;
}

double QCDModel::getProb(unsigned int pMom, 
			 unsigned int muMom, 
			 int region, 
			 unsigned int flav, 
			 unsigned int qnet, 
			 bool inclusive) {

  if (pMom>(nPartMom-1)) pMom = nPartMom - 1;
  if (muMom>(nMuMom-1)) muMom = nMuMom - 1;
  if (region<0) region = 0;
  if (region>1) region = 1;
  if (flav>(nFlav-1)) flav = nFlav - 1;
  if (qnet>1) qnet = 1;
  
  double output = 0.0;
  if (inclusive) {
    double num = pass[flav][qnet][pMom][muMom][region]->GetBinContent(1);
    double den = norm[flav][qnet][pMom][muMom]->GetBinContent(1);
    if (den>0)
      output = num/den;
  }
  else {
    double num = passSS[flav][qnet][pMom][muMom][region]->GetBinContent(1);
    double den = normSS[flav][qnet][pMom][muMom]->GetBinContent(1);
    if (den>0)
      output = num/den;
  }
  return output;

}

double QCDModel::getMassPdf(unsigned int pMom, 
			    unsigned int muMom, 
			    int region, 
			    unsigned int flav, 
			    unsigned int qnet, 
			    unsigned int mass, 
			    bool inclusive) {
  
  if (pMom>(nPartMom-1)) pMom = nPartMom - 1;
  if (muMom>nMuMom-1) muMom = nMuMom - 1;
  if (region<0) region = 0;
  if (region>1) region = 1;
  if (flav>(nFlav-1)) flav = nFlav - 1;
  if (qnet>1) qnet = 1;

  double output = 1.0;
  if (inclusive) {
    output  = pdfMass[flav][qnet][pMom][muMom][region]->GetBinContent(mass);
  }
  else {
    output  = pdfMassSS[flav][qnet][pMom][muMom][region]->GetBinContent(mass);
  }
  return output;

}

