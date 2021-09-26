#include "HttStylesNew.cc"
#include "HtoH.h"

void CorrCoefficients() {

  TString dir("/nfs/dust/cms/user/rasp/Run/Run2018/H2aa");
  SetStyle();
  
  TFile * file = new TFile(dir+"/DoubleMuon_Run2018.root");
  TH1D * hist1Dold = (TH1D*)file->Get("InvMassTrackPlusMuon1D_ControlXH");
  TH2D * hist2Dold = (TH2D*)file->Get("InvMassTrackPlusMuon2D_ControlXH");

  std::cout << "Statistics : " << hist2Dold->GetSumOfWeights() << std::endl;

  int nBinsNew = 6;
  double bins[7]     = {0,1,2,3,4,6,20};
  double binsCorr[7] = {0,1,2,3,4,6,12};

  TH1D * hist1D = (TH1D*)TH1DtoTH1D(hist1Dold,nBinsNew,bins,true,"_new");
  TH2D * hist2D = (TH2D*)TH2DtoTH2D(hist2Dold,nBinsNew,bins,nBinsNew,bins,"_sigNew");
  TH2D * corrCoeff = new TH2D("corrCoeff","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  TH2D * corrCoeffX = new TH2D("corrCoeffX","",nBinsNew,bins,nBinsNew,bins);
  
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      printf("[%1i,%1i] = %5.0f\n",iB,jB,hist2D->GetBinContent(iB,jB));
    }
  }  
  std::cout << std::endl;
  hist1D->Scale(1/hist1D->GetSumOfWeights());
  hist2D->Scale(1/hist2D->GetSumOfWeights());

  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      float x = hist1D->GetBinContent(iB);
      float y = hist1D->GetBinContent(jB);
      float denominator = x*y;
      if (iB!=jB)
	denominator *= 2.0;
      float numerator  = hist2D->GetBinContent(iB,jB);
      float enumerator = hist2D->GetBinError(iB,jB);
      float corr  = numerator / denominator ;
      float ecorr = enumerator / denominator;
      corr = floor(100*corr+0.5)/100;
      ecorr = floor(100*ecorr+0.5)/100;
      printf("[%1i,%1i] : %5.2f +/- %5.2f \n",iB,jB,corr,ecorr);
      corrCoeff->SetBinContent(iB,jB,corr);
      corrCoeff->SetBinError(iB,jB,ecorr);
      corrCoeffX->SetBinContent(iB,jB,corr);
      corrCoeffX->SetBinError(iB,jB,ecorr);
    }
  }

  corrCoeff->GetXaxis()->SetNdivisions(207);
  corrCoeff->GetYaxis()->SetNdivisions(207);
  corrCoeff->GetYaxis()->SetTitleOffset(1.0);
  corrCoeff->SetMarkerSize(1.2);
  corrCoeff->GetXaxis()->SetTitle("m_{1} [GeV]");
  corrCoeff->GetYaxis()->SetTitle("m_{2} [GeV]");

  TCanvas * canv = MakeCanvas("canv","",700,700);
  corrCoeff->Draw("texte");
  
  for (int i=1; i<nBinsNew; ++i) {
    float yL = binsCorr[i];
    TLine * lineX = new TLine(0,yL,12,yL);
    TLine * lineY = new TLine(yL,0,yL,12);
    lineX->SetLineWidth(1);
    lineY->SetLineWidth(1);
    lineX->Draw();
    lineY->Draw();
  }
  
  canv->Update();
  canv->Print("CorrCoefficients_data.png");

  TFile * fileCorr = new TFile("CorrCoefficients_data.root","recreate");
  fileCorr->cd("");
  corrCoeffX->Write("corrCoeff");
  fileCorr->Close();

}
