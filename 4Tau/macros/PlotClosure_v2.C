#include "HttStylesNew.cc"
#include "HtoH.h"
#include "CMS_lumi.C"

// ClosureInvMassDimuonIso
// InvMassHighMuLooseIso
// InvMassHighMuIso

void PlotClosure_v2(TString histName1 = "partonMuSelectedIso_SS",
		    TString histName2 = "partonMuModelledIso_SS",
		    TString era = "2016_preVFP",
		    TString bin = "bin0p2",
		    bool logY = true) {

  SetStyle();

  std::cout << std::endl;
  std::cout << "Era = " << era << "    bin = " << bin << std::endl;
  std::cout << std::endl;

  TString dir="/nfs/dust/cms/user/rasp/Run/QCDModel/"+era+"/"+bin;

  TString SamplesQCD[10] = {
    "QCD_Pt-20To30_MuEnrichedPt5",   // (1)
    "QCD_Pt-30To50_MuEnrichedPt5",   // (2)
    "QCD_Pt-50To80_MuEnrichedPt5",   // (3)
    "QCD_Pt-80To120_MuEnrichedPt5",  // (4)
    "QCD_Pt-120To170_MuEnrichedPt5", // (5)
    "QCD_Pt-170To300_MuEnrichedPt5", // (6)
    "QCD_Pt-300To470_MuEnrichedPt5", // (7)
    "QCD_Pt-470To600_MuEnrichedPt5", // (8)
    "QCD_Pt-600To800_MuEnrichedPt5", // (9)
    "QCD_Pt-800To1000_MuEnrichedPt5" // (10)
  };

  double xsecQCD[10] = {
    558528000*0.0053,  // (1)
    139803000*0.01182, // (2)
    19222500*0.02276,  // (3)
    2758420*0.03844,   // (4)
    469797*0.05362,    // (5)
    117989*0.07335,    // (6)
    7820.25*0.10196,   // (7)
    645.528*0.12242,   // (8)
    187.109*0.13412,   // (9)
    32.3486*0.14552    // (10)
  };


  TH1D * hist1;
  TH1D * hist2;

  for (int iS=0; iS<10; ++iS) {
    TFile * fileSample = new TFile(dir+"/"+SamplesQCD[iS]+".root");
    TH1D * histWeightsH = (TH1D*)fileSample->Get("histWeightsH");
    double norm = xsecQCD[iS]/histWeightsH->GetSumOfWeights();
    TH1D * hist1Sample = (TH1D*)fileSample->Get(histName1);
    TH1D * hist2Sample = (TH1D*)fileSample->Get(histName2);
    
    if (iS==0) {
      hist1 = hist1Sample;
      hist2 = hist2Sample;
      hist1->Scale(norm);
      hist2->Scale(norm);
    }
    else {
      hist1->Add(hist1,hist1Sample,1.0,norm);
      hist2->Add(hist2,hist2Sample,1.0,norm);
    }
  }

  printf("%s : %3.2f\n",histName1.Data(),hist1->GetSumOfWeights());
  printf("%s : %3.2f\n",histName2.Data(),hist2->GetSumOfWeights());

}
