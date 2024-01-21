// =====================================================
// ========== Computing DY corrections =================
// =====================================================

#include "HtoH.h"

std::map<TString, double> sample_xsec = {
  {"WW_13TeV-pythia8",118.7},
  {"WZ_13TeV-pythia8",27.68},
  {"ZZ_13TeV-pythia8",12.19},
  {"WJetsToLNu",61526.7},
  {"ST_t-channel_top",136.02},
  {"ST_t-channel_antitop",80.95},
  {"ST_tW_top",35.85}, 
  {"ST_tW_antitop",35.85},
  {"TTTo2L2Nu",88.29},
  {"TTToSemiLeptonic",365.35},
  {"TTToHadronic",377.96},
  {"DYJetsToLL_M-10to50",21610.0},
  {"DYJetsToLL_M-50",6077.22},
  {"DYJetsToTT_M-50",6077.22},
  {"DYJetsToLL_M-50_amcatnlo",6077.22},
  {"DYJetsToTT_M-50_amcatnlo",6077.22}
  
};

std::vector<TString > bkgs = {
  "WW_13TeV-pythia8",
  "WZ_13TeV-pythia8",
  "ZZ_13TeV-pythia8",
  "WJetsToLNu",
  "ST_t-channel_top",
  "ST_t-channel_antitop",
  "ST_tW_top",
  "ST_tW_antitop",
  "TTTo2L2Nu",
  "TTToSemiLeptonic",
  "TTToHadronic",
  "DYJetsToTT_M-50_amcatnlo"
};

std::vector<TString > signals = {
  "DYJetsToLL_M-10to50",
  "DYJetsToLL_M-50_amcatnlo"
};

std::map<TString,double> eraLumi = {
  {"2016_preVFP", 19520.},
  {"2016_postVFP",16810.},
  {"2017",        41480.},
  {"2018",        59830.}
};

void ComputeDYcorr(TString era = "2016_postVFP",
		   TString histName = "pTmassMuMuH") {

  std::vector<double> binsPt   = {0.,5.,10.,15.,20.,30.,40.,60.,80.,100.,200.,400.,1000.};
  std::vector<double> binsMass = {50.,60.,70.,80.,85.,90.,95.,100.,110.,120.,150.,200.,1000.};

  gROOT->SetBatch(true);
  
  TString folder = "/nfs/dust/cms/user/rasp/Run/MuTrk/PuppiMET/"+era;

  int nBinsX = binsPt.size()-1;
  int nBinsY = binsMass.size()-1;
  double binsX[100];
  double binsY[100];

  for (int iB=0; iB<=nBinsX; ++iB)
    binsX[iB] = binsPt[iB];

  for (int iB=0; iB<=nBinsY; ++iB)
    binsY[iB] = binsMass[iB];

  double lumi = eraLumi[era];

  TFile * fileData = new TFile(folder+"/Data.root");
  TH2D * histDataOld = (TH2D*)fileData->Get(histName);

  std::cout << std::endl;
  std::cout << "Entries in data : " << histDataOld->GetEntries() << std::endl;
  std::cout << std::endl;


  int nbinsx = histDataOld->GetNbinsX();
  int nbinsy = histDataOld->GetNbinsY();
  double xmin = histDataOld->GetXaxis()->GetBinLowEdge(1);
  double xmax = histDataOld->GetXaxis()->GetBinLowEdge(nbinsx+1);
  double ymin = histDataOld->GetYaxis()->GetBinLowEdge(1);
  double ymax = histDataOld->GetYaxis()->GetBinLowEdge(nbinsy+1);

  std::cout << "Histogram " << histName << " : " 
	    << "nbinsX = " << nbinsx
	    << "  xmin = " << xmin
	    << "  xmax = " << xmax 
	    << "  nbinsY = " << nbinsy
	    << "  ymin = " << ymin
	    << "  ymax = " << ymax 
	    << std::endl;
  std::cout << std::endl;

  TH2D * histData = (TH2D*)TH2DtoTH2D(histDataOld,nBinsX,binsX,nBinsY,binsY,"_Data_DYcorr_"+era);
  std::cout << " yield (1) = " << int(histDataOld->GetSumOfWeights()) << std::endl; 
  std::cout << " yield (2) = " << int(histData->GetSumOfWeights()) << std::endl;; 
  std::cout << std::endl;
  TH2D * histBkg = NULL;
  TH2D * histSig = NULL;
  bool isFirst = true;
  double totYieldBkg = 0;
  double totYieldSig = 0;
  for (auto sample : bkgs) {
    TFile * file = new TFile(folder+"/"+sample+".root");
    TH2D * histOld = (TH2D*)file->Get(histName);
    TH2D * hist = (TH2D*)TH2DtoTH2D(histOld,nBinsX,binsX,nBinsY,binsY,"_"+sample+"_DYcorr_"+era);
    TH1D * eventCount = (TH1D*)file->Get("histWeightsH");
    double nGen = eventCount->GetSumOfWeights();
    double norm = sample_xsec[sample]*lumi/nGen;
    hist->Scale(norm);
    double yield = hist->GetSumOfWeights();
    printf("%30s %8.0f\n",sample.Data(),yield);
    if (isFirst) { 
      histBkg = hist;
      isFirst = false;
    }
    else {
      histBkg->Add(histBkg,hist,1.,1.);
    }
    totYieldBkg += yield;
  }

  isFirst = true;
  for (auto sample : signals) {
    TFile * file = new TFile(folder+"/"+sample+".root");
    TH2D * histOld = (TH2D*)file->Get(histName);
    TH2D * hist = (TH2D*)TH2DtoTH2D(histOld,nBinsX,binsX,nBinsY,binsY,"_"+sample+"_DYcorr_"+era);
    TH1D * eventCount = (TH1D*)file->Get("histWeightsH");
    double nGen = eventCount->GetSumOfWeights();
    double norm = sample_xsec[sample]*lumi/nGen;
    hist->Scale(norm);
    double yield = hist->GetSumOfWeights();
    printf("%30s %8.0f\n",sample.Data(),yield);
    if (isFirst) { 
      histSig = hist;
      isFirst = false;
    }
    else {
      histSig->Add(histSig,hist,1.,1.);
    }
    totYieldSig += yield;
  }

  
  std::cout << std::endl;

  TString name="Total bkg";
  printf("%15s : %8.0f\n",name.Data(),totYieldBkg);

  name = "Total sig";
  printf("%15s : %8.0f\n",name.Data(),totYieldSig);

  double totYield = totYieldBkg + totYieldSig;
  name = "Total";
  printf("%15s : %8.0f\n",name.Data(),totYield);

  name = "Data";
  double totYieldData = histData->GetSumOfWeights();
  printf("%15s : %8.0f\n",name.Data(),totYieldData);
  std::cout << std::endl;
  double ratio = totYieldData/totYield;
  printf("ratio = %5.3f\n",ratio);
  std::cout << std::endl;

  histData->Add(histData,histBkg,1.,-1.);
  
  double normData = 1.0/totYieldData;
  histData->Scale(normData);

  double normZLL = 1.0/totYieldSig;
  histSig->Scale(normZLL);

  TH2D * histRatio = (TH2D*)histSig->Clone("DY_corr");

  for (int iB=1; iB<=nBinsX; ++iB) {
    for (int jB=1; jB<=nBinsY; ++jB) {
      double num = histData->GetBinContent(iB,jB);
      double den = histSig->GetBinContent(iB,jB);
      double ratio = 1.0;
      if (den>0.) {
	ratio = num/den;
      }
      if (ratio>10.) ratio = 10.;
      if (ratio<0.1) ratio = 0.1;
      histRatio->SetBinContent(iB,jB,ratio);
      histRatio->SetBinError(iB,jB,0.);
    }
  }

  TString outfilename = "DYCorr_NLO_"+era+".root";
  TFile * file = new TFile(outfilename,"recreate");
  file->cd("");
  histRatio->Write("DY_NLO");
  file->Close();


}
