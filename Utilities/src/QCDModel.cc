#include "HtoAA/Utilities/interface/QCDModel.h"

QCDModel::QCDModel(TString fileName) {

  file = new TFile(fileName);


  TString partonFlavor[4] = {"uds","g","c","b"};
  TString muonPartonNetCharge[2] = {"opposite","same"};
  TString partonMomRange[3] = {"Lt50","50to100","Gt100"};
  TString muonMomRange[3] = {"Lt30","30to50","Gt50"};
  TString muonType[2] = {"HighMu","LowMu"}; 
  TString region[3] = {"Iso","LooseIso","Sb"};

  for (int imom=0; imom<3; ++imom) {
    for (int imu=0; imu<2; ++imu) {
      for (int ireg=0; ireg<3; ++ireg) {
	ProbIsoMu[imom][imu][ireg] = (TH2D*)file->Get("ProbIsoMu_"+partonMomRange[imom]+"_"+muonType[imu]+"_"+region[ireg]);
	ProbSSIsoMu[imom][imu][ireg] = (TH2D*)file->Get("ProbSSIsoMu_"+partonMomRange[imom]+"_"+muonType[imu]+"_"+region[ireg]);
	pdfMass[imom][imu][ireg] = (TH3D*)file->Get("pdfMass_"+partonMomRange[imom]+"_"+muonType[imu]+"_"+region[ireg]);
	ProbIsoMuUnmatched[imom][imu][ireg] = (TH1D*)file->Get("ProbIsoMuUnmatched_"+muonMomRange[imom]+"_"+muonType[imu]+"_"+region[ireg]);
	ProbSSIsoMuUnmatched[imom][imu][ireg] = (TH1D*)file->Get("ProbSSIsoMuUnmatched_"+muonMomRange[imom]+"_"+muonType[imu]+"_"+region[ireg]);
	pdfMassUnmatched[imom][imu][ireg] = (TH1D*)file->Get("pdfMassUnmatched_"+muonMomRange[imom]+"_"+muonType[imu]+"_"+region[ireg]);
	//	std::cout << imom << " " << imu << " " << ireg << " : " << ProbSSIsoMu[imom][imu][ireg] << " " << ProbSSIsoMuUnmatched[imom][imu][ireg] << std::endl;
      }
    }
  }



}

QCDModel::~QCDModel() {

}

double QCDModel::getProbIsoMu(int iMom, int muType, int ireg, int iflav, int inet) {

  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;
  if (ireg<0) ireg = 0;
  if (ireg>2) ireg = 2;
  
  if (iflav==1) inet=0; // gluons

  double output = 1;
  if (iflav<0) {
    output = ProbIsoMuUnmatched[iMom][muType][ireg]->GetBinContent(1);
  }
  else 
    output = ProbIsoMu[iMom][muType][ireg]->GetBinContent(iflav+1,inet+1);

  return output;

}

double QCDModel::getProbSSIsoMu(int iMom, int muType, int ireg, int iflav, int inet) {

  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;
  if (ireg<0) ireg = 0;
  if (ireg>2) ireg = 2;
  
  if (iflav==1) inet=0; // gluons

  double output = 1;
  if (iflav<0) {
    output = ProbSSIsoMuUnmatched[iMom][muType][ireg]->GetBinContent(1);
  }
  else 
    output = ProbSSIsoMu[iMom][muType][ireg]->GetBinContent(iflav+1,inet+1);

  return output;

}

double QCDModel::getMassPdf(int iMom, int muType, int ireg, int iflav, int inet, double mass) {
  
  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;
  if (ireg<0) ireg = 0;
  if (ireg>2) ireg = 2;

  if (iflav==1) inet=0; // gluons 

  int imass = int(mass);

  double output = 1;
  if (iflav<0)
    output = pdfMassUnmatched[iMom][muType][ireg]->GetBinContent(imass+1);
  else
    output = pdfMass[iMom][muType][ireg]->GetBinContent(iflav+1,inet+1,imass+1);
  return output;

}

double QCDModel::getMuMassPdf(int iMom, int muType, int ireg, int iflav,int inet, double mass) {

  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;
  if (ireg<0) ireg = 0;
  if (ireg>2) ireg = 2;

  if (iflav==1) inet=0; // gluons

  double output = getProbIsoMu(iMom,muType,ireg,iflav,inet) *
    getMassPdf(iMom,muType,ireg,iflav,inet,mass);

  return output;

}

double QCDModel::getSSMuMassPdf(int iMom, int muType, int ireg, int iflav,int inet, double mass) {

  if (iMom<0) iMom = 0;
  if (iMom>2) iMom = 2;
  if (muType<0) muType = 0;
  if (muType>1) muType = 1;
  if (ireg<0) ireg = 0;
  if (ireg>2) ireg = 2;

  if (iflav==1) inet=0; // gluons

  double output = getProbSSIsoMu(iMom,muType,ireg,iflav,inet) *
    getMassPdf(iMom,muType,ireg,iflav,inet,mass);

  return output;

}
