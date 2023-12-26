#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <string> 
#include <map>

#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "HtoAA/Utilities/interface/Config.h"
#include "HtoAA/Utilities/src/Config.cc"
#include <iomanip>

using namespace std;

std::map<TString,double> eraLumi = {
  {"2016_preVFP", 19520.},
  {"2016_postVFP",16810.},
  {"2017",        41480.},
  {"2018",        59830.}
};

std::map<TString, TString> LUMI_label = {
  {"2016"     ,"2016, 36.3 fb^{-1}"},
  {"2016_pre" ,"2016, preVFP, 19.5 fb^{-1}"},
  {"2016_post","2016, postVFP, 16.8 fb^{-1}"},
  {"2017"     ,"2017, 41.5 fb^{-1}"},
  {"2018"     ,"2018, 59.8 fb^{-1}"}
};

// =====================================================
// ============ Extracting histograms ==================
// =====================================================
void ExtractHistograms() {

}


// =====================================================
// ==================== Plotting =======================
// =====================================================

void Plot(std::map<TString, TH1D*> histograms,
	  TString xtitle,
	  TString region,
	  TString legend,
	  TString era,
	  TString plotdir) {

  SetStyle();


  TH1D * Data = histograms["Data"];
  int nBins = Data->GetNbinsX();

  for (auto sysSample : sysSamples) {
    TString name = sysSample.first;
    double sys = sysSample.second;
    TH1D * hist = histograms[name];
    for (int iB=0; iB<nBins; ++iB) {
      double x = hist->GetBinContent(iB);
      double e = hist->GetBinError(iB);
      double error = TMath::Sqrt(e*e+x*x*sys*sys);
      hist->SetBinError(iB,error);
    }
  }


  TH1D * ZLL = histograms["ZLL"];
  TH1D * ZTT = histograms["ZTT"];
  TH1D * EWK = histograms["EWK"];
  TH1D * TT  = histograms["TT"];
  TH1D * W   = histograms["WJ"];
  TH1D * QCD = histograms["QCD"];

  std::cout << std::endl;
  std::cout << "Plotting distributions in region " << region << std::endl;
  std::cout << std::endl;
  printf("ZLL   : %7.0f\n",ZLL->GetSumOfWeights());
  printf("EWK   : %7.0f\n",EWK->GetSumOfWeights());
  printf("TT    : %7.0f\n",TT->GetSumOfWeights());
  printf("W     : %7.0f\n",W->GetSumOfWeights());
  printf("QCD   : %7.0f\n",QCD->GetSumOfWeights());
  printf("ZTT   : %7.0f\n",ZTT->GetSumOfWeights());

  TT->Add(TT,EWK,1.,1.);
  QCD->Add(QCD,TT,1.,1.);
  W->Add(W,QCD,1.,1.);
  ZLL->Add(ZLL,W,1.,1.);
  ZTT->Add(ZTT,ZLL,1.,1.);

  TH1D * bkgdErr = (TH1D*)ZTT->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);
  printf("----------------\n");
  printf("Total : %7.0f\n",bkgdErr->GetSumOfWeights());
  printf("Data  : %7.0f\n",histograms["Data"]->GetSumOfWeights());
  std::cout << std::endl;

  
  for (int iB=1; iB<=nBins; ++iB) {
    TT->SetBinError(iB,0);
    EWK->SetBinError(iB,0);
    W->SetBinError(iB,0);
    ZLL->SetBinError(iB,0);
    ZTT->SetBinError(iB,0);
    QCD->SetBinError(iB,0);
  }
  
  InitData(Data);

  InitHist(ZTT,xtitle,"Events",TColor::GetColor("#FFCC66"),1001);
  InitHist(ZLL,"","",TColor::GetColor("#4496C8"),1001);
  InitHist(QCD,"","",TColor::GetColor("#FFCCFF"),1001);
  InitHist(TT,"","",TColor::GetColor("#9999CC"),1001);
  InitHist(W,"","",TColor::GetColor("#c6f74a"),1001);
  InitHist(EWK,"","",TColor::GetColor("#DE5A6A"),1001);

  double ymax = Data->GetMaximum();
  if (ZTT->GetMaximum()>ymax) ymax = ZTT->GetMaximum();
  ZTT->GetYaxis()->SetRangeUser(0.,1.2*ymax);

  TCanvas * canv = MakeCanvas("canv","",600,600);

  ZTT->Draw("h");
  ZLL->Draw("hsame");
  W->Draw("hsame");
  QCD->Draw("hsame");
  TT->Draw("hsame");
  EWK->Draw("hsame");
  bkgdErr->Draw("e2same");
  Data->Draw("e1same");
    
  TLegend * leg = new TLegend(0.6,0.4,0.90,0.75);
  SetLegendStyle(leg);
  leg->SetHeader(legend);
  leg->SetTextSize(0.03);
  leg->AddEntry(Data,"Data","lp");
  leg->AddEntry(ZTT,"Z#rightarrow#tau#tau","f");
  leg->AddEntry(ZLL,"Z#rightarrow#mu#mu","f");
  leg->AddEntry(W,"W+jets","f");
  leg->AddEntry(QCD,"QCD","f");
  leg->AddEntry(TT,"t#bar{t}","f");
  leg->AddEntry(EWK,"electroweak","f");
  leg->Draw();

  lumi_13TeV = LUMI_label[era];
  writeExtraText = true;
  extraText = "Preliminary";
  CMS_lumi(canv,4,33); 

  canv->RedrawAxis();
  canv->Update();
  canv->Print(plotdir+"/"+region+"_"+era+".png");

  std::cout << std::endl;
  delete canv;

}

int main(int argc, char * argv[]) {

  if (argc!=3) {
    std::cout << "Usage : createCardsZtt [config] [era] [ptbin] [mTregion]" << std::endl;
    exit(-1);
  }

  if (fopen(argv[1],"r")==NULL) {
    std::cout << "Configuration file " << argv[1] << " does not exist" << std::endl;
    exit(-1);
  }

  TString era(argv[2]);
  Config cfg(argv[1]);
  const string inputdir = cfg.get<string>("InputFolder");
  const string plotdir = cfg.get<string>("PlotFolder"); 
  const string carddir = cfg.get<string>("DatacardsFolder");
  const float xmin = cfg.get<float>("XMin");
  const float xmax = cfg.get<float>("XMax");
  const int nbins = cfg.get<int>("NBins");
  const string lowMT = cfg.get<string>("CutLowMT");
  const string highMT = cfg.get<string>("CutHighMT");
  const bool simpleCards = cfg.get<bool>("SimpleCards");

  TString CutLowMT = TString("mt<")+TString(lowMT);
  TString CutHighMT = TString("mt>")+TString(highMT);
  TString DatacardsDir(carddir);
  TString PlotDir(plotdir);
  TString dir(inputdir);

  TFile * fileMLFit


  if ()

}
