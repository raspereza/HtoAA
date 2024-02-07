#ifndef QCD_MODEL_DEFINITIONS_H
#define QCD_MODEL_DEFINITIONS_H

#include "TString.h"

const int nBins = 100;
const double xMin = 0.;
const double xMax = 20.;

const unsigned int nPartMom = 3;
const float partonMomBins[nPartMom+1] = {0,30,50,1000};

const unsigned int nMuMom = 3;
const float muonMomBins[nMuMom+1] = {0,20,30,1000};

const unsigned int nFlav = 4;
const TString partonFlavor[nFlav] = {"g","uds","c","b"};

const unsigned int nNetQ = 2;
const TString muonPartonNetCharge[nNetQ] = {"opposite","same"};

const TString partonMomRange[nPartMom] = {
  "partPtLt30",
  "partPt30to50",
  "partPtGt50"
};

const TString muonMomRange[nMuMom] = {
  "muonPtLt20",
  "muonPt20to30",
  "muonPtGt30"
};

const unsigned int nReg = 2;
const TString Regions[nReg] = {"LooseIso","Iso"};

#endif
