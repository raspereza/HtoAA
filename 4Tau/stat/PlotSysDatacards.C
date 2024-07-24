#include "HtoH.h"
#include "HttStylesNew.cc"

// -------------------------
// systematics (bkgd)
// -------------------------
// CMS_haa4t_uncCorr_FSR
// CMS_haa4t_uncCorr_ISR
// CMS_haa4t_uncCorr_nonQCD
// -------------------------

// ---------------------------------
// signal (ggh, vbf, vh, tth, mmtt)
// ---------------------------------
// CMS_haa4t_eff_trkiso_$ERA
// ---------------------------------

void PlotSysDatacards(
		      TString era = "2016",
		      //		      TString sysName = "CMS_haa4t_uncCorr_FSR",
		      TString sysName = "CMS_haa4t_eff_trkiso_2016",
		      TString histName = "ggh",
		      TString mass = "10"
		       ) {

  TString plotDir = ".";
  TString inputDir = ".";

  TString fileName = inputDir + "/haa_"+era+"-13TeV_ma"+mass;

  bool print = true;
  float YMin = 0.81;
  float YMax = 1.19;
  float upRange = -100;

  TString header = era+":"+histName+":"+sysName;
  TString SysLeg = sysName;
  TString ytitle("Events / bin");
  TString xtitle("bin"); 
  bool logX = false;
  bool logY = false;

  SetStyle();
  gStyle->SetErrorX(0);
  TFile * file = new TFile(fileName+".root");
  TH1D * histNominal = (TH1D*)file->Get(histName);
  TH1D * histUp = (TH1D*)file->Get(histName+"_"+sysName+"Up");
  TH1D * histDown = (TH1D*)file->Get(histName+"_"+sysName+"Down");
  //  std::cout << histNominal << " " <<  histUp << " " << histDown << std::endl;
  if (histNominal==NULL) {
    std::cout << "Histogram " << histName << "is absent in datacards RooT file" << std::endl;
    return;
  }
  if (histUp==NULL) {
    std::cout << "Histogram " << histName << "_"<< sysName << "Up is absent in datacards RooT file"	<< std::endl;
    return;
  }
  if (histDown==NULL) {
    std::cout << "Histogram " << histName << "_"<< sysName << "Down is absent in datacards RooT file"	<< std::endl;    
    return;
  }
  double xNominal = histNominal->GetSumOfWeights();
  double xUp = histUp->GetSumOfWeights();
  double xDown = histDown->GetSumOfWeights();

  std::cout << std::endl;
  if (histName=="bkgd")
    std::cout << "Era : " << era << " ; template : " << histName << " ; systematics : " << sysName << std::endl;
  else
    std::cout << "Era : " << era << " ; template : " << histName << "_ma" << mass << " ; systematics : " << sysName << std::endl;
  std::cout << sysName << "   lnN    " << xDown/xNominal << "/" << xUp/xNominal << std::endl;
  std::cout << std::endl;
  int nBins = histNominal->GetNbinsX(); 
  float xmax = histNominal->GetBinLowEdge(nBins+1)-0.01;

  if (print) {
    for (int iB=1; iB<=nBins; ++iB) {
      double x = histNominal->GetBinContent(iB);
      double ex = histNominal->GetBinError(iB);
      double xmin = histNominal->GetBinLowEdge(iB);
      double xmax = histNominal->GetBinLowEdge(iB+1);
      //      double x = histUp->GetBinContent(iB);
      //      double ex = histUp->GetBinError(iB);
      double ylower = histUp->GetBinContent(iB);
      double yupper = histDown->GetBinContent(iB);
      printf("%2i : %7.2f +/- %5.2f, %7.2f - %7.2f\n",iB,x,ex,ylower,yupper);

    }
  }
  std::cout << std::endl;
  std::cout << "Total = " << histNominal->GetSumOfWeights() << std::endl;
  std::cout << std::endl;

  InitData(histNominal);

  histNominal->GetXaxis()->SetTitleSize(0.0);
  histNominal->GetXaxis()->SetTitleOffset(1.2);

  histNominal->GetYaxis()->SetTitleSize(0.07);
  histNominal->GetYaxis()->SetTitleOffset(1.0);
  histNominal->GetYaxis()->SetLabelSize(0.045);

  double maximum = 0;
  for (int ib=1; ib<=nBins; ++ib) {
    double xM = histNominal->GetBinContent(ib)+histNominal->GetBinError(ib);
    if (xM>maximum) maximum = xM;
  }

  if (histUp->GetMaximum()>maximum) maximum = histUp->GetMaximum();
  if (histDown->GetMaximum()>maximum) maximum = histDown->GetMaximum();

  histNominal->GetYaxis()->SetRangeUser(0.01,1.2*maximum);
  histNominal->SetLineColor(1);
  histUp->SetLineColor(2);
  histDown->SetLineColor(4);
  histDown->SetLineStyle(3);
  histNominal->SetMarkerColor(1);
  histUp->SetMarkerColor(2);
  histDown->SetMarkerColor(4);
  histNominal->SetMarkerSize(1.4);
  histNominal->GetYaxis()->SetTitle(ytitle);
  histNominal->GetXaxis()->SetTitle(xtitle);
  histUp->GetYaxis()->SetTitle(ytitle);
  histUp->GetXaxis()->SetTitle(xtitle);
  histDown->GetYaxis()->SetTitle(ytitle);
  histDown->GetXaxis()->SetTitle(xtitle);
  histUp->SetLineWidth(2);
  histDown->SetLineWidth(2);
  TH1D * ratioUp = (TH1D*)histUp->Clone("ratioUp");
  TH1D * ratioDown = (TH1D*)histDown->Clone("ratioDown");
  TH1D * ratioCentral = (TH1D*)histNominal->Clone("ratioCentral");
  //  ratioCentral->SetFillStyle(3013);
  //  ratioCentral->SetFillColor(1);
  //  ratioCentral->SetMarkerStyle(21);
  //  ratioCentral->SetMarkerSize(0);


  for (int iB=1; iB<=nBins; ++iB) {
    histUp->SetBinError(iB,0); 
    histDown->SetBinError(iB,0); 
    float xUp = histUp->GetBinContent(iB);
    float xDown = histDown->GetBinContent(iB);
    float xCentral = histNominal->GetBinContent(iB);
    float xratioUp = 1;
    float xratioDown = 1;
    if (xCentral>0) {
      xratioUp   = xUp/xCentral;
      xratioDown = xDown/xCentral;
    }
    ratioUp->SetBinContent(iB,xratioUp);
    ratioDown->SetBinContent(iB,xratioDown);
    ratioUp->SetBinError(iB,0);
    ratioDown->SetBinError(iB,0);
    ratioCentral->SetBinContent(iB,1);
    ratioCentral->SetBinError(iB,0);
    if (histNominal->GetBinContent(iB)>0)
      ratioCentral->SetBinError(iB,histNominal->GetBinError(iB)/histNominal->GetBinContent(iB));
  }

  if (upRange>0) 
    histUp->GetYaxis()->SetRangeUser(0.1,upRange);

  histUp->GetYaxis()->SetTitleOffset(1.4);

  histNominal->SetTitle(header);

  TCanvas * canv1 = MakeCanvas("canv1", "", 600, 700);
  TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.17);
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

  histNominal->Draw("hpe");
  histUp->Draw("hsame");
  histDown->Draw("hsame");
  TLegend * leg = new TLegend(0.35,0.68,0.60,0.9);
  SetLegendStyle(leg);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(histNominal,"Central","ep");
  leg->AddEntry(histUp,"Up","l");
  leg->AddEntry(histDown,"Down","l");
  leg->Draw();
  upper->SetLogx(logX);
  upper->SetLogy(logY);
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  ratioUp->SetTitle("");
  ratioUp->GetYaxis()->SetRangeUser(YMin,YMax);
  ratioUp->GetYaxis()->SetNdivisions(505);
  ratioUp->GetXaxis()->SetLabelFont(42);
  ratioUp->GetXaxis()->SetLabelOffset(0.04);
  ratioUp->GetXaxis()->SetLabelSize(0.10);
  ratioUp->GetXaxis()->SetTitleSize(0.13);
  ratioUp->GetXaxis()->SetTitleOffset(1.2);
  ratioUp->GetYaxis()->SetTitle("ratio");
  ratioUp->GetYaxis()->SetLabelFont(42);
  ratioUp->GetYaxis()->SetLabelOffset(0.015);
  ratioUp->GetYaxis()->SetLabelSize(0.1);
  ratioUp->GetYaxis()->SetTitleSize(0.14);
  ratioUp->GetYaxis()->SetTitleOffset(0.5);
  ratioUp->GetXaxis()->SetTickLength(0.07);
  ratioUp->GetYaxis()->SetTickLength(0.04);
  ratioUp->GetYaxis()->SetLabelOffset(0.01);

  // ------------>Primitives in pad: lower
  TPad * lower = new TPad("lower", "pad",0,0,1,0.32);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.17);
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
  
  ratioUp->GetXaxis()->SetRangeUser(61,1499);
  
  ratioUp->Draw("h");
  ratioDown->Draw("hsame");
  ratioCentral->Draw("he1same");
  
  lower->SetLogx(logX);
  lower->Modified();
  lower->RedrawAxis();
  canv1->cd();
  canv1->Modified();
  canv1->cd();
  if (histName=="bkgd")
    canv1->Print(plotDir+"/"+era+"_"+histName+"_"+sysName+".png");
  else
    canv1->Print(plotDir+"/"+era+"_"+histName+"_ma"+mass+"_"+sysName+".png");
  
} 
