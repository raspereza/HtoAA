#ifndef QCD_MODEL_DEFINITIONS_H
#define QCD_MODEL_DEFINITIONS_H

#include "TString.h"

const int nBins = 100;
const double xMin = 0.;
const double xMax = 20.;

const unsigned int nPartMom = 4;
const float partonMomBins[nPartMom+1] = {0,30,50,100,1000};

const unsigned int nMuMom = 4;
const float muonMomBins[nMuMom+1] = {0,20,30,50,1000};

const unsigned int nFlav = 5;
const TString partonFlavor[nFlav] = {"undef","g","udsg","c","b"};

const unsigned int nNetQ = 2;
const TString muonPartonNetCharge[nNetQ] = {"opposite","same"};

const TString partonMomRange[nPartMom] = {
  "partPLt30",
  "partP30to50",
  "partP50to100",
  "partPGt100"
};

const TString muonMomRange[nMuMom] = {
  "muonPLt20",
  "muonP20to30",
  "muonP30to50",
  "muonPGt50"
};

const unsigned int nReg = 2;
const TString Regions[nReg] = {"LooseIso","Iso"};

#endif
