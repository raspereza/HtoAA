#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"
#include "boost/range/algorithm.hpp"
#include "boost/range/algorithm_ext.hpp"
#include "Plotting.h"
#include "Plotting_Style.h"

void PlotMass2D(
		TString era = "2016",
		bool prefit = false, // prefit (or postfit) distributions
		bool blindData = false, // blind data
		bool drawLeg = true, // draw legend
		bool logY = true // use log scale for
		) {
  

  TH1::SetDefaultSumw2(true);

  std::map<TString, vector<TString> > era_groups = {
    {"2016_postVFP",{"2016_postVFP"}},
    {"2016_preVFP",{"2016_preVFP"}},
    {"2016",{"2016"}},
    {"2017",{"2017"}},
    {"2018",{"2018"}},
    {"Run2",{"2016","2017","2018"}}
  };

  vector<TString> group = era_groups[era];

  std::map<TString,TString> channel_map = 
    {
      {"2016","ch1"},
      {"2017","ch2"},
      {"2018","ch3"}
    };

  std::map<TString,TString> lumi_label = {
    {"2016_preVFP","2016 APV, 19.5 fb^{-1}"},
    {"2016_postVFP","2016 16.8 fb^{-1}"},
    {"2016","2016 36.3 fb^{-1}"},
    {"2017","2017 41.5 fb^{-1}"},
    {"2018","2018 59.8 fb^{-1}"},
    {"Run2","138 fb^{-1}"}
  }; 

  lumi_13TeV = lumi_label[era];

  TString dir("/nfs/dust/cms/user/rasp/Run/HtoAA/stat");

  TH1D * hist5  = new TH1D("sig5","",21,0.,21.);
  TH1D * hist10 = new TH1D("sig10","",21,0.,21.);
  TH1D * hist15 = new TH1D("sig15","",21,0.,21.);

  vector<TString> signals = {"ggh","vbf","vh","tth","mmtt"};
  bool isFirst = true;
  for (auto grp : group) {
    TFile * file5  = new TFile(dir+"/haa_"+grp+"-13TeV_ma5.root");
    TFile * file10 = new TFile(dir+"/haa_"+grp+"-13TeV_ma10.root");
    TFile * file15 = new TFile(dir+"/haa_"+grp+"-13TeV_ma15.root"); 
    for (auto sig : signals) {
      TH1D * h5  = (TH1D*)file5->Get(sig);
      TH1D * h10 = (TH1D*)file10->Get(sig);
      TH1D * h15 = (TH1D*)file15->Get(sig);
      hist5->Add(hist5,h5);
      hist10->Add(hist10,h10);
      hist15->Add(hist15,h15);
    }
  }

  // to get background model with proper uncertainties from
  // RooT file created by combine utility
  TFile * fileFit = new TFile(dir+"/fitDiagnosticsTest.root");

  gStyle->SetOptStat(0000);

  int nBins = 6;
  int nBinsNew = (nBins+1)*nBins/2;

  hist5->Scale(0.05);
  hist10->Scale(0.05);
  hist15->Scale(0.05);

  TString shapesDir("shapes_fit_b");
  if (prefit)
    shapesDir = "shapes_prefit";

  TH1D * histData = new TH1D("data","",21,0.,21.);
  TH1D * histBkgdX = new TH1D("bkgd","",21,0.,21.);
  isFirst = true;
  for (auto grp : group) {
    TFile * file  = new TFile(dir+"/haa_"+grp+"-13TeV_ma10.root");
    TH1D * hd = (TH1D*)file->Get("data_obs");
    TH1D * hb = (TH1D*)fileFit->Get(shapesDir+"/"+channel_map[grp]+"/total_background");
    histData->Add(histData,hd,1.,1.);
    histBkgdX->Add(histBkgdX,hb,1.,1.);
  }

  TH1D * histBkgd = new TH1D("histBkgd","",21,0.,21.0);

  TH1D * histBkgdErr = (TH1D*)histBkgdX->Clone("histBkgdErr");
  int new_idx = CreateTransparentColor(kGray,0.5);
  histBkgdErr->SetFillColor(new_idx);
  histBkgdErr->SetFillStyle(1001);
  histBkgdErr->SetMarkerStyle(0);
  histBkgdErr->SetMarkerSize(0);
  histBkgdErr->SetLineColor(4);
  histBkgdErr->SetLineWidth(2);

  histBkgd->GetXaxis()->SetLabelSize(0.06);
  histBkgd->GetXaxis()->SetLabelOffset(0.01);
  histBkgd->GetYaxis()->SetTitleOffset(0.77);
  histBkgd->GetYaxis()->SetTitle("Events / bin");
  histBkgd->GetYaxis()->SetTitleSize(0.06);
  histBkgd->GetYaxis()->SetLabelSize(0.06);
  histBkgd->GetYaxis()->SetTickLength(0.025);
  histBkgd->GetXaxis()->SetTickLength(0.025);
  histBkgd->GetYaxis()->SetTickSize(0.02);
  histBkgd->GetXaxis()->SetTickSize(0.02);

  histBkgd->GetYaxis()->SetRangeUser(0,1.1*histBkgd->GetMaximum());
  if (logY) histBkgd->GetYaxis()->SetRangeUser(0.2,50001);
  histBkgd->SetLineColor(4);
  histBkgd->SetLineWidth(2);
  histBkgd->SetLineStyle(1);
  
  hist5->SetLineColor(kRed);
  hist5->SetLineWidth(3);
  hist5->SetLineStyle(1);

  hist10->SetLineColor(kGreen+1);
  hist10->SetLineWidth(3);
  hist10->SetLineStyle(1);

  hist15->SetLineColor(kMagenta+2);
  hist15->SetLineWidth(3);
  hist15->SetLineStyle(1);

  TH1D * ratio = (TH1D*)histData->Clone("ratio");
  TH1D * ratioErr = (TH1D*)histData->Clone("ratioErr");
  ratioErr->SetFillColor(new_idx);
  ratioErr->SetFillStyle(1001);
  ratioErr->SetMarkerStyle(21);
  ratioErr->SetMarkerSize(0);
  ratioErr->SetLineColor(4);
  ratioErr->SetLineWidth(3);

  double x[21];
  double ex[21];
  double y[21];
  double eyhigh[21];
  double eylow[21];
  double yR[21];
  double erlow[21];
  double erhigh[21];

  int nBinsX = histBkgdErr->GetNbinsX();
  for (int iB=1; iB<=nBinsX; ++iB) {
    histBkgd->SetBinContent(iB,histBkgdX->GetBinContent(iB));
    histBkgd->SetBinError(iB,0);
    histData->SetBinError(iB,0);
    hist5->SetBinError(iB,0);
    hist10->SetBinError(iB,0);
    hist15->SetBinError(iB,0);
    double xD = histData->GetBinContent(iB);
    double eD = histData->GetBinError(iB);
    double xB = histBkgdErr->GetBinContent(iB);
    double eB = histBkgdErr->GetBinError(iB);
    double r = 1000;
    double er = 0;
    double erB = eB/xB;
    r = xD/xB;
    er = eD/xB;
    if (r<0.001) r = 0.001;
    if (erB>0.9) erB = 0.9;
    ratio->SetBinContent(iB,r);
    ratio->SetBinError(iB,0.);
    ratioErr->SetBinContent(iB,1);
    ratioErr->SetBinError(iB,erB);
    x[iB-1] = iB - 0.5;
    ex[iB-1] = 0.;
    y[iB-1] = histData->GetBinContent(iB);
    eylow[iB-1] = -0.5 + TMath::Sqrt(histData->GetBinContent(iB)+0.25);
    eyhigh[iB-1] = 0.5 + TMath::Sqrt(histData->GetBinContent(iB)+0.25);
    yR[iB-1] = xD/xB;
    erlow[iB-1] = eylow[iB-1]/xB;
    erhigh[iB-1] = eyhigh[iB-1]/xB;
  }
  TGraphAsymmErrors * data = new TGraphAsymmErrors(21,x,y,ex,ex,eylow,eyhigh);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.3);
  data->SetMarkerColor(1);
  data->SetLineWidth(2);

  TGraphAsymmErrors * ratioERR = new TGraphAsymmErrors(21,x,yR,ex,ex,erlow,erhigh);
  ratioERR->SetMarkerStyle(20);
  ratioERR->SetMarkerSize(1.3);
  ratioERR->SetMarkerColor(1);
  ratioERR->SetLineWidth(2);


  int binNumber = 1;
  for (int i=1; i<=nBins; ++i) {
    for (int j=i; j<=nBins;++j) {
      char charLabel[10];
      sprintf(charLabel,"(%1i,%1i)",i,j);
      TString label(charLabel);
      histBkgd->GetXaxis()->SetBinLabel(binNumber,"");
      ratioErr->GetXaxis()->SetBinLabel(binNumber,label);
      ratio->GetXaxis()->SetBinLabel(binNumber,label);
      binNumber++;
    }
  }

  ratio->GetXaxis()->LabelsOption("v");
  ratioErr->GetXaxis()->LabelsOption("v");

  TCanvas * canv = MakeCanvas("canv","",750,600);
  TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.09);
  upper->SetRightMargin(0.05);
  upper->SetBottomMargin(0.02);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);


  histBkgd->Draw("h");
  histBkgdErr->Draw("e2same");
  histBkgd->Draw("hsame");
  hist5->Draw("hsame");
  hist10->Draw("hsame");
  hist15->Draw("hsame");
  if (!blindData)  data->Draw("pe1");
  TLegend * leg = new TLegend(0.2,0.70,0.42,0.9);
  SetLegendStyle(leg);
  leg->SetTextSize(0.055);
  if (!blindData)  leg->AddEntry(data,"observed","lp");
  leg->AddEntry(histBkgdErr,"bkg(+unc)","lf");
  TLegend * leg1 = new TLegend(0.5,0.62,0.72,0.9);
  SetLegendStyle(leg1);
  leg1->SetTextSize(0.055);
  leg1->AddEntry(hist5,"m_{a_{1}} = 5 GeV","l");
  leg1->AddEntry(hist10,"m_{a_{1}} = 10 GeV","l");
  leg1->AddEntry(hist15,"m_{a_{1}} = 15 GeV","l");
  if (drawLeg) {
    leg->Draw(); leg1->Draw();
  }
  writeExtraText = true;
  CMS_lumi(upper,4,33); 
  float xLine = float(nBins);
  for (int i=1; i<=nBins; ++i) {
    TLine * line = new TLine(xLine,0.2,xLine,1000);
    line->SetLineWidth(1);
    line->SetLineStyle(3);
    line->Draw();
    xLine += nBins - i; 
  }

  if (logY) upper->SetLogy(true);

  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv->cd();


  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.09);
  lower->SetRightMargin(0.05);
  lower->SetTopMargin(0.026);
  lower->SetBottomMargin(0.35);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);

  ratioErr->GetYaxis()->SetRangeUser(0.0,3.2);
  ratioErr->GetXaxis()->SetNdivisions(210);
  ratioErr->GetYaxis()->SetNdivisions(505);
  ratioErr->GetXaxis()->SetLabelFont(42);
  ratioErr->GetXaxis()->SetLabelOffset(0.03);
  ratioErr->GetXaxis()->SetLabelSize(0.25);
  ratioErr->GetXaxis()->SetTitleSize(0.13);
  ratioErr->GetXaxis()->SetTitleOffset(1.2);
  ratioErr->GetYaxis()->SetTitle("obs/bkg");
  ratioErr->GetYaxis()->SetLabelFont(42);
  ratioErr->GetYaxis()->SetLabelOffset(0.015);
  ratioErr->GetYaxis()->SetLabelSize(0.14);
  ratioErr->GetYaxis()->SetTitleSize(0.14);
  ratioErr->GetYaxis()->SetTitleOffset(0.3);
  ratioErr->GetXaxis()->SetTickLength(0.025);
  ratioErr->GetYaxis()->SetTickLength(0.025);
  ratioErr->GetXaxis()->SetTickSize(0.02);
  ratioErr->GetYaxis()->SetTickSize(0.02);
  ratioErr->GetYaxis()->SetLabelOffset(0.01);
   
  ratioErr->Draw("e2");
  ratioERR->Draw("pe1");
  ratioErr->GetXaxis()->LabelsOption("v");

  TLine * line = new TLine(0,1,21,1);
  line->SetLineColor(4);
  line->SetLineWidth(2);
  line->Draw();
  xLine = float(nBins);
  for (int i=1; i<=nBins; ++i) {
    TLine * line = new TLine(xLine,0.,xLine,3.2);
    line->SetLineWidth(2);
    line->SetLineStyle(3);
    line->Draw();
    xLine += nBins - i;
  }
  //  ratioERR->Draw("pe1");
  //  ratio->Draw("pe1same");
  
  lower->Modified();
  lower->RedrawAxis();
  canv->cd();
  canv->SetSelected(canv);

  if (prefit) 
    canv->Print("mass2D_"+era+"_prefit.png");
  else
    canv->Print("mass2D_"+era+"_postfit.png");

}

//
  
