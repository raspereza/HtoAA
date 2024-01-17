#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"

void CorrCoefficientsPoisson(TString era = "2017",
			     bool extendedSideBand = false) {

  //  gROOT->SetBatch(true);
  int nBinsNew = 6;
  double bins[7]     = {0, 1, 2, 3, 4, 5.2, 20};
  double binsCorr[7] = {0, 1, 2, 3, 4, 5.2, 12};

  std::map<TString, TString> LUMI_label = {
    {"2016"     ,"2016, 36.3 fb^{-1}"},
    {"2016_preVFP" ,"2016, preVFP, 19.5 fb^{-1}"},
    {"2016_postVFP","2016, postVFP, 16.8 fb^{-1}"},
    {"2017"     ,"2017, 41.5 fb^{-1}"},
    {"2018"     ,"2018, 59.8 fb^{-1}"}
  };
  
  lumi_13TeV = LUMI_label[era];

  TString dir = "/nfs/dust/cms/user/rasp/Run/HtoAA/" + era;
  SetStyle();
  
  TString histName1D("InvMassTrackPlusMuon1D_ControlXH");
  TString histName2D("InvMassTrackPlusMuon2D_ControlXH");
  if (extendedSideBand) {
    histName1D = "InvMassTrackPlusMuon1D_ControlYH";
    histName2D = "InvMassTrackPlusMuon2D_ControlYH";
  }

  TFile * file = new TFile(dir+"/Data.root");
  TH1D * hist1Dold = (TH1D*)file->Get(histName1D);
  TH2D * hist2Dold = (TH2D*)file->Get(histName2D);

  std::cout << "Statistics : " << hist2Dold->GetSumOfWeights() << std::endl;

  TH1D * hist1D = (TH1D*)TH1DtoTH1D(hist1Dold,nBinsNew,bins,true,"_new");
  TH2D * hist2D = (TH2D*)TH2DtoTH2D(hist2Dold,nBinsNew,bins,nBinsNew,bins,"_sigNew");
  TH2D * corrCoeff = new TH2D("corrCoeff","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  TH2D * corrCoeffH = new TH2D("corrCoeffH","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  TH2D * corrCoeffL = new TH2D("corrCoeffL","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  TH2D * corrCoeffX = new TH2D("corrCoeffX","",nBinsNew,bins,nBinsNew,bins);
  
  double norm2D = hist2D->GetSumOfWeights(); 

  TCanvas * canv1 = MakeCanvas("canv1","",1000,700);
  corrCoeff->GetXaxis()->SetNdivisions(207);
  corrCoeff->GetYaxis()->SetNdivisions(207);
  corrCoeff->GetYaxis()->SetTitleOffset(1.0);
  corrCoeff->SetMarkerSize(1.2);
  corrCoeff->GetXaxis()->SetTitle("m_{1} [GeV]");
  corrCoeff->GetYaxis()->SetTitle("m_{2} [GeV]");
  corrCoeff->Draw("texte");
  TLatex * latexBin1 = new TLatex();
  latexBin1->SetTextSize(0.039);
  latexBin1->SetTextFont(32);
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      char label[10];
      int x = int(hist2D->GetBinContent(iB,jB));
      sprintf(label,"%4i",x);
      TString Label(label);
      double binX = corrCoeffH->GetXaxis()->GetBinCenter(iB);
      double binWidthX = corrCoeffH->GetXaxis()->GetBinWidth(iB); 
      double binY = corrCoeffH->GetYaxis()->GetBinCenter(jB);
      double binWidthY = corrCoeffH->GetYaxis()->GetBinWidth(jB);       
      double textX = binX - 0.45;
      double textY = binY;
      latexBin1->SetTextAlign(12);
      latexBin1->DrawLatex(textX,textY,Label);
      //      textY = binY - 0.05*binWidthY;
      //      latexBin->DrawLatex(textX,textY,LabelError);
    }
  }
  
  for (int i=1; i<nBinsNew; ++i) {
    double yL = binsCorr[i];
    TLine * lineX = new TLine(0,yL,12,yL);
    TLine * lineY = new TLine(yL,0,yL,12);
    lineX->SetLineWidth(1);
    lineY->SetLineWidth(1);
    lineX->Draw();
    lineY->Draw();
  }
  
  writeExtraText = true;
  extraText = "Preliminary";
  CMS_lumi(canv1,4,33);

  canv1->Update();
  if (extendedSideBand)
    canv1->Print("stat_data_"+era+"_ext.png");
  else
    canv1->Print("stat_data_"+era+".png");

  
  std::cout << "Statistics in 2D distribution : " << std::endl;
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = hist2D->GetBinContent(iB,jB);
      printf("[%1i,%1i] = %4.0f\n",iB,jB,hist2D->GetBinContent(iB,jB));
      double elow = -0.5 + TMath::Sqrt(x+0.25);
      double ehigh = 0.5 + TMath::Sqrt(x+0.25);
      x = x/norm2D;
      elow = elow/norm2D;
      ehigh = ehigh/norm2D;
      corrCoeffL->SetBinContent(iB,jB,x);
      corrCoeffL->SetBinError(iB,jB,elow);
      corrCoeffH->SetBinContent(iB,jB,x);
      corrCoeffH->SetBinError(iB,jB,ehigh);
    }
  }  
  
  std::cout << std::endl;
  hist1D->Scale(1/hist1D->GetSumOfWeights());
  hist2D->Scale(1/hist2D->GetSumOfWeights());  

  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = hist1D->GetBinContent(iB);
      double y = hist1D->GetBinContent(jB);
      double denominator = x*y;
      if (iB!=jB)
	denominator *= 2.0;
      double numerator  = hist2D->GetBinContent(iB,jB);
      double enumerator = hist2D->GetBinError(iB,jB);
      double corr  = numerator / denominator ;
      double ecorr = enumerator / denominator;
      corr = floor(100*corr+0.5)/100;
      ecorr = floor(100*ecorr+0.5)/100;
      printf("[%1i,%1i] : %5.2f +/- %5.2f \n",iB,jB,corr,ecorr);
      double ehigh = corrCoeffH->GetBinError(iB,jB)/denominator;
      double elow = corrCoeffL->GetBinError(iB,jB)/denominator;
      corrCoeffH->SetBinContent(iB,jB,corr);
      corrCoeffH->SetBinError(iB,jB,ehigh);
      corrCoeffL->SetBinContent(iB,jB,corr);
      corrCoeffL->SetBinError(iB,jB,elow);      
      corrCoeffX->SetBinContent(iB,jB,corr);
      corrCoeffX->SetBinError(iB,jB,ecorr);
      if (iB<0) {
	corrCoeff->SetBinContent(iB,jB,corr);
	corrCoeff->SetBinError(iB,jB,ecorr);
      }
    }
  }

  corrCoeff->GetXaxis()->SetNdivisions(207);
  corrCoeff->GetYaxis()->SetNdivisions(207);
  corrCoeff->GetYaxis()->SetTitleOffset(1.0);
  corrCoeff->SetMarkerSize(1.2);
  corrCoeff->GetXaxis()->SetTitle("m_{1} [GeV]");
  corrCoeff->GetYaxis()->SetTitle("m_{2} [GeV]");

  TCanvas * canv = MakeCanvas("canv","",1000,700);
  corrCoeff->Draw("texte");
  TLatex * latexBin = new TLatex();
  latexBin->SetTextSize(0.029);
  latexBin->SetTextFont(32);
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = corrCoeffH->GetBinContent(iB,jB);
      double elow = corrCoeffL->GetBinError(iB,jB);
      double ehigh = corrCoeffH->GetBinError(iB,jB);
      double error = corrCoeffX->GetBinError(iB,jB);
      char label[10];
      sprintf(label,"%4.2f^{+%4.2f}_{ -%4.2f}",x,ehigh,elow);
      TString Label(label);
      double binX = corrCoeffH->GetXaxis()->GetBinCenter(iB);
      double binWidthX = corrCoeffH->GetXaxis()->GetBinWidth(iB); 
      double binY = corrCoeffH->GetYaxis()->GetBinCenter(jB);
      double binWidthY = corrCoeffH->GetYaxis()->GetBinWidth(jB);       
      double textX = binX - 0.45;
      double textY = binY;
      latexBin->SetTextAlign(12);
      latexBin->DrawLatex(textX,textY,Label);
      //      textY = binY - 0.05*binWidthY;
      //      latexBin->DrawLatex(textX,textY,LabelError);
    }
  }

  
  for (int i=1; i<nBinsNew; ++i) {
    double yL = binsCorr[i];
    TLine * lineX = new TLine(0,yL,12,yL);
    TLine * lineY = new TLine(yL,0,yL,12);
    lineX->SetLineWidth(1);
    lineY->SetLineWidth(1);
    lineX->Draw();
    lineY->Draw();
  }
  
  writeExtraText = true;
  extraText = "Preliminary";
  CMS_lumi(canv,4,33);

  canv->Update();
  if (extendedSideBand)
    canv->Print("CorrCoefficients_data_"+era+"_ext.png");
  else
    canv->Print("CorrCoefficients_data_"+era+".png");

  TFile * fileCorr = new TFile("CorrCoefficients_data_"+era+".root","recreate");
  fileCorr->cd("");
  corrCoeffX->Write("corrCoeff");
  fileCorr->Close();

}
