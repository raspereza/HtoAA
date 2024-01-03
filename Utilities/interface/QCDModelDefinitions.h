#ifndef QCD_MODEL_DEFINITIONS_H
#define QCD_MODEL_DEFINITIONS_H

#include "TString.h"

const int nBins = 40;
const double xMin = 0.;
const double xMax = 20.;

unsigned int nPartonMomBins = 4;
float partonMomBins[5] = {0,30,50,100,1000};

unsigned int nMuonMomBins = 4;
float muonMomBins[5] = {0,20,30,50,1000};

unsigned int nPartonFlavours = 5;
TString partonFlavor[5] = {"undef","g","uds","c","b"};

unsigned int nMuonPartonNetCharge = 2;
TString muonPartonNetCharge[2] = {"opposite","same"};

TString partonMomRange[4] = {
  "partPtLt30",
  "partPt30to50",
  "partPt50to100",
  "partPtGt100"
};

TString muonMomRange[4] = {
  "muonPtLt20",
  "muonPt20to30",
  "muonPt30to50",
  "muonPtGt50"
};

unsigned int nRegions = 2;
TString Regions[2] = {"LooseIso","Iso"};

#endif
