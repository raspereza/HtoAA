#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"

//TH1D * deltaRMuonPionH = new TH1D("deltaRMuonPionH","",200,0,2);
//TH1D * pionPtH = new TH1D("pionPtH","",100,0,100);

void PlotDistributions(TString histName = "deltaRMuonPionH",
		       TString dir = "./",
		       float xLower = 0,
		       float xUpper = 2,
		       TString xtitle = "#delta R(#mu,track)",
		       TString ytitle = "normalized to unity") {
  
  SetStyle();

  TFile * file4 = new TFile(dir+"/SUSYGluGluToHToAA_AToTauTau_M-125_M-4.root");
  TFile * file8 = new TFile(dir+"/SUSYGluGluToHToAA_AToTauTau_M-125_M-7.root");
  TFile * file15 = new TFile(dir+"/SUSYGluGluToHToAA_AToTauTau_M-125_M-15.root");
  
  TH1D * hist4 = (TH1D*)file4->Get(histName);
  TH1D * hist8 = (TH1D*)file8->Get(histName);
  TH1D * hist15 = (TH1D*)file15->Get(histName);

  int nBins = hist4->GetNbinsX();
  float xMin = hist4->GetBinLowEdge(1);
  float xMax = hist4->GetBinLowEdge(nBins+1);
  std::cout << histName << "   nBins = " << nBins << "   xmin = " << xMin << "   xmax = " << xMax << std::endl;

  int nBinsNew = nBins;
  double bins[200];

  hist4->Scale(1.0/hist4->GetSumOfWeights());
  hist8->Scale(1.0/hist8->GetSumOfWeights());
  hist15->Scale(1.0/hist15->GetSumOfWeights());

  std::cout << "new number of bins : "; std::cin >> nBinsNew;
  float binWidth = (xMax-xMin)/float(nBinsNew);
  for (int iB=0; iB<=nBinsNew; ++iB) 
    bins[iB] = xMin + float(iB)*binWidth;
  
  TH1D * histNew4  = TH1DtoTH1D(hist4,nBinsNew,bins,true,"_4");
  TH1D * histNew8  = TH1DtoTH1D(hist8,nBinsNew,bins,true,"_8");
  TH1D * histNew15 = TH1DtoTH1D(hist15,nBinsNew,bins,true,"_15");

  for (int iB=1; iB<=nBinsNew; ++iB) {
    histNew4->SetBinError(iB,0);
    histNew8->SetBinError(iB,0);
    histNew15->SetBinError(iB,0);
  }

  histNew4->SetLineColor(1);
  histNew8->SetLineColor(2);
  histNew15->SetLineColor(4);

  histNew4->SetMarkerColor(1);
  histNew8->SetMarkerColor(2);
  histNew15->SetMarkerColor(4);

  histNew4->SetMarkerSize(0);
  histNew8->SetMarkerSize(0);
  histNew15->SetMarkerSize(0);

  histNew4->SetLineWidth(2);
  histNew8->SetLineWidth(2);
  histNew15->SetLineWidth(2);

  float yMax = histNew4->GetMaximum();
  if (yMax<histNew8->GetMaximum()) yMax = histNew8->GetMaximum();
  if (yMax<histNew15->GetMaximum()) yMax = histNew15->GetMaximum();

  histNew4->GetXaxis()->SetRangeUser(xLower,xUpper);
  histNew4->GetYaxis()->SetRangeUser(0.,1.1*yMax);

  histNew4->GetXaxis()->SetTitle(xtitle);
  histNew4->GetYaxis()->SetTitle(ytitle);

  TCanvas * canv = MakeCanvas("canv","",600,600);
  histNew4->Draw();
  histNew8->Draw("same");
  histNew15->Draw("same");
  TLegend * leg = new TLegend(0.6,0.54,0.8,0.80);
  SetLegendStyle(leg);
  leg->SetTextSize(0.046);
  leg->AddEntry(histNew4,"m_{a} = 4 GeV","l");
  leg->AddEntry(histNew8,"m_{a} = 7 GeV","l");
  leg->AddEntry(histNew15,"m_{a} = 15 GeV","l");
  leg->Draw();
  writeExtraText = true;
  extraText   = "Simulation";
  CMS_lumi(canv,4,33); 

  canv->Update();
  canv->Print("figures/"+histName+".png");

  TFile * fileOutput = new TFile("deltaR_Haa.root","recreate");
  fileOutput->cd();
  histNew4->Write("dR_m4");
  histNew8->Write("dR_m7");
  histNew15->Write("dR_m10");
  fileOutput->Close();

}
