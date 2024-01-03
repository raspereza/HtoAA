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
 double getProb   (int pMom, int muMom, int region, int flav, int qnet, bool incl); // prob(mu->Iso/LooseIso) 
 double getMassPdf(int pMom, int muMom, int region, int flav, int qnet, double mass, bool incl); // pdf(m) inclusive

 private:

 TFile * file;

 // 5 flavours 
 // x 2 qnet (net charge of muon parton) 
 // x 4 parton momentum bins 
 // x 4 muon momentum bins 
 // x 2 regions (Iso, LooseIso)
 // Bins are defined in the header file
 // HtoAA/Utilities/interface/QCDModelDefinitions 
 TH1D * prob[5][2][4][4][2]; 
 TH1D * probSS[5][2][4][4][2]; 

 TH1D * pdfMass[5][2][4][4][2]; 
 TH1D * pdfMassSS[5][2][4][4][2]; 


};

#endif
