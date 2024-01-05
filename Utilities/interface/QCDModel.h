#ifndef QCDModel_h
#define QCDModel_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "QCDModelDefinitions.h"

using namespace std;

class QCDModel {

 public: 

 QCDModel(TString fileName);
 ~QCDModel();

 // bool incl = True -> probability and mass pdf are taken from inclusive muon sample
 // bool incl = False -> probability and mass pdf are taken from sample of SS muons
 double getProb   (unsigned int pMom, 
		   unsigned int muMom, 
		   int region, 
		   unsigned int flav, 
		   unsigned int qnet, 
		   bool incl); // prob(mu->Iso/LooseIso) 
 double getMassPdf(unsigned int pMom, 
		   unsigned int muMom, 
		   int region, 
		   unsigned int flav, 
		   unsigned int qnet, 
		   double mass, 
		   bool incl); // pdf(m) inclusive

 private:

 TFile * file;

 // 5 flavours 
 // x 2 qnet (net charge of muon parton) 
 // x 4 parton momentum bins 
 // x 4 muon momentum bins 
 // x 2 regions (Iso, LooseIso)
 // Bins are defined in the header file
 // HtoAA/Utilities/interface/QCDModelDefinitions 
 TH1D * pass[nFlav][nNetQ][nPartMom][nMuMom][nReg]; 
 TH1D * passSS[nFlav][nNetQ][nPartMom][nMuMom][nReg];

 TH1D * norm[nFlav][nNetQ][nPartMom][nMuMom]; 
 TH1D * normSS[nFlav][nNetQ][nPartMom][nMuMom]; 

 TH1D * pdfMass[nFlav][nNetQ][nPartMom][nMuMom][nReg]; 
 TH1D * pdfMassSS[nFlav][nNetQ][nPartMom][nMuMom][nReg];


};

#endif
