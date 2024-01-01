#include "HttStylesNew.cc"
#include "HtoH.h"
#include "CMS_lumi.C"

// mass_Selected_SR
// mass_Modelled_SR
// mass_Selected_Loose
// mass_Modelled_Loose

void PlotClosure(
		 TString era = "2018",
		 bool signalRegion = true,
		 bool normalize = true,
		 bool logY = true
		 ) {
  SetStyle();

  TString legend = "LooseIso";
  TString histName = "mass_Selected_LooseIso";
  TString histNameModel = "mass_Modelled_LooseIso";
  
  if (signalRegion) {
    legend = "SR";
    histName = "mass_Selected_SR";
    histNameModel = "mass_Modelled_SR";
  }

  int nBins = 6;
  double bins[7]     = {0,1,2,3,4,6,20};


  TFile * file = new TFile("/nfs/dust/cms/user/rasp/Run/HtoAA/QCDModel_"+era+".root");

  TH1D * histX = (TH1D*)file->Get(histName);
  TH1D * histModelX = (TH1D*)file->Get(histNameModel);

  std::cout << "Direct   = " << histX->GetSumOfWeights() <<  std::endl;
  std::cout << "Model   =  " << histModelX->GetSumOfWeights() <<  std::endl;

  TH1D * hist = (TH1D*)TH1DtoTH1D(histX,nBins,bins,true,"_rebinned");
  TH1D * histModel = (TH1D*)TH1DtoTH1D(histModelX,nBins,bins,true,"_rebinned");
    

  hist->SetLineColor(1);
  hist->SetMarkerColor(1);
  hist->SetMarkerSize(1.2);
  hist->SetMarkerStyle(20);
  hist->GetXaxis()->SetLabelSize(0.);
  hist->GetYaxis()->SetTitle("normalized");

  histModel->SetLineColor(2);
  histModel->SetMarkerColor(2);
  histModel->SetMarkerSize(0);
  histModel->SetMarkerStyle(0);

  hist->Scale(1.0/hist->GetSumOfWeights());
  histModel->Scale(1.0/histModel->GetSumOfWeights());
  if (logY) hist->GetYaxis()->SetRangeUser(0.0001,2.);

  TH1D * ratioH = (TH1D*)hist->Clone("ratioH");
  
  ratioH->GetYaxis()->SetRangeUser(0.001,1.999);

  for (int iB=1; iB<=nBins; ++iB) {
    float den = histModel->GetBinContent(iB);
    if (den<1e-6) {
      std::cout << "bin : " << iB << "  den = " << den << std::endl;
      den = 1;
    }
    float x = hist->GetBinContent(iB);
    float ratioX = x/den;
    float e = hist->GetBinError(iB);
    float ratioE = e/den;
    ratioH->SetBinContent(iB,ratioX);
    ratioH->SetBinError(iB,ratioE);
  }


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

  hist->GetYaxis()->SetRangeUser(1e-3,1.);
  hist->Draw();
  histModel->Draw("hsame");
  
  lumi_13TeV = "2017, 41.5 fb^{-1}";
  if (era=="2016")
    lumi_13TeV = "2016, 36.3 fb^{-1}";
  if (era=="2018")
    lumi_13TeV = "2018, 59.8 fb^{-1}";

  writeExtraText = true;
  extraText   = "Simulation";
  CMS_lumi(upper,4,33); 

  TLegend * leg = new TLegend(0.55,0.55,0.9,0.8);
  SetLegendStyle(leg);
  leg->SetHeader(legend);
  leg->SetTextSize(0.05);
  leg->AddEntry(hist,"actual selection","lp");
  leg->AddEntry(histModel,"model","l");
  leg->Draw();

  if (logY) upper->SetLogy(true);
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
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

  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.2);
  ratioH->SetLineColor(1);

  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleSize(0.13);
  ratioH->GetXaxis()->SetTitleOffset(1.2);
  ratioH->GetYaxis()->SetTitle("obs/exp");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.015);
  ratioH->GetYaxis()->SetLabelSize(0.13);
  ratioH->GetYaxis()->SetTitleSize(0.14);
  ratioH->GetYaxis()->SetTitleOffset(0.6);
  ratioH->GetXaxis()->SetTickLength(0.07);
  ratioH->GetYaxis()->SetTickLength(0.04);
  ratioH->GetYaxis()->SetLabelOffset(0.01);
  ratioH->GetXaxis()->SetTitle("m_{#mu,trk} [GeV]");

  ratioH->Draw("e1");
  TLine * line = new TLine(0.,1.,20.,1.);
  line->SetLineColor(2);
  line->Draw();
  ratioH->Draw("e1same");

  lower->Modified();
  lower->RedrawAxis();
  canv1->cd();
  canv1->SetSelected(canv1);
  canv1->Print(histName+"_"+era+".png");


}
