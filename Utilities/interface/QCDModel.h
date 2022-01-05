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

using namespace std;

class QCDModel {

 public: 

 QCDModel(TString fileName);
 ~QCDModel();

 double getProbIsoMu(int iMom, int muType, int ireg, int iflav,int inet); // P(mu->iso mu)
 double getProbSSIsoMu(int iMom, int muType, int ireg, int iflav,int inet); // P(SS mu->iso mu)
 double getMassPdf(int iMom, int muType, int ireg, int iflav, int inet, double mass); // pdf(m)
 double getMuMassPdf(int iMom, int muType, int ireg, int iflav,int inet, double mass); // P(mu->iso mu)*pdf(m)
 double getSSMuMassPdf(int iMom, int muType, int ireg, int iflav,int inet, double mass); // P(SS mu->iso mu)*pdf(m)


 private:

 TFile * file;

 // 3 momentum bins x 2 muon types x 3 regions
 TH2D * ProbIsoMu[3][2][3]; // 2D : [flavor,net-charge]
 TH2D * ProbSSIsoMu[3][2][3]; // 2D : [flavor,net-charge]
 TH1D * ProbIsoMuUnmatched[3][2][3]; // 1D : 1bin
 TH1D * ProbSSIsoMuUnmatched[3][2][3]; // 1D : 1bin
 

 // 3 momentum bins x 2 muon types x 3 regions
 TH3D * pdfMass[3][2][3]; // 3D : [flavor,net-charge,mass]                                              
 TH1D * pdfMassUnmatched[3][2][3]; // 1D : [mass]


};

#endif
