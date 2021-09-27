#include "CMS_lumi.C"
#include "HttStylesNew.cc"

void PlotSignal(TString mass = "4") {

  SetStyle();

  TString dir("/nfs/dust/cms/user/rasp/Run/Run2018/H2aa/");

  TFile * file = new TFile(dir+"/haa_2018-13TeV_ma"+mass+".root");
  TH1D * mmtt = (TH1D*)file->Get("mmtt");
  TH1D * tttt = (TH1D*)file->Get("ggh");

  TString signalNames[3] = {"vbf","vh","tth"}; 
  for (int iS=0; iS<3; ++iS) {
    TH1D * histS = (TH1D*)file->Get(signalNames[iS]);
    tttt->Add(tttt,histS);
  }
  tttt->Scale(0.05);
  mmtt->Scale(0.03);
  int nBins = tttt->GetNbinsX();
  for (int iB=1; iB<=nBins; ++iB) {
    tttt->SetBinError(iB,0.);
    mmtt->SetBinError(iB,0.);
  } 
  tttt->SetLineColor(2);
  mmtt->SetLineColor(4);
  tttt->SetLineWidth(2);
  mmtt->SetLineWidth(2);
  
  int nBins1D = 6;
  int binNumber = 1;
  for (int i=1; i<=nBins1D; ++i) {
    for (int j=i; j<=nBins1D;++j) {
      char charLabel[10];
      sprintf(charLabel," (%1i,%1i)",i,j);
      TString label(charLabel);
      tttt->GetXaxis()->SetBinLabel(binNumber,label);
      binNumber++;
    }
  }
  tttt->GetYaxis()->SetRangeUser(0.1,2000);
  tttt->GetYaxis()->SetTitleOffset(0.9);
  tttt->GetYaxis()->SetTitle("Events / bin");
  tttt->GetXaxis()->SetLabelSize(0.054);
  tttt->GetXaxis()->SetTickLength(0.025);
  tttt->GetYaxis()->SetTickLength(0.025);
  tttt->GetXaxis()->SetTickSize(0.02);
  tttt->GetYaxis()->SetTickSize(0.02);
  tttt->GetXaxis()->SetLabelOffset(0.01);
  tttt->GetXaxis()->SetLabelSize(0.09);
  tttt->GetXaxis()->LabelsOption("v");

  TCanvas * canv = MakeCanvas("canv","",900,600);
  canv->SetLeftMargin(0.1);
  tttt->Draw("h");
  mmtt->Draw("hsame");
  TLegend * leg = new TLegend(0.15,0.8,0.4,0.92);
  SetLegendStyle(leg);
  leg->SetTextSize(0.07);
  leg->SetHeader("m_{a_{1}} = "+mass+" GeV");
  leg->Draw();
  TLegend * leg1 = new TLegend(0.5,0.7,0.8,0.92);
  SetLegendStyle(leg1);
  leg1->SetTextSize(0.07);
  leg1->AddEntry(tttt,"4#tau","l");
  leg1->AddEntry(mmtt,"2#mu2#tau","l");
  leg1->Draw();
  extraText = "Simulation";
  writeExtraText = true;
  CMS_lumi(canv,4,33); 
  float xLine = float(nBins1D);
  for (int i=1; i<=nBins1D; ++i) {
    TLine * line = new TLine(xLine,0.1,xLine,100);
    line->SetLineWidth(2);
    line->SetLineStyle(3);
    line->Draw();
    xLine += nBins1D - i; 
  }


  canv->SetLogy(true);
  canv->Update();
  canv->Print("mass2D_m"+mass+".png");


}
