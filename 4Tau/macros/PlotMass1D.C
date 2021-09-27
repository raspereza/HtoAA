#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"

void PlotMass1D(float xLower = 0, // lower boundary in x axis
		float xUpper = 12, // upper boundary in x axis
		TString xtitle = "m_{#mu,trk} [GeV]",
		TString ytitle = "normalised to unity",
		bool drawLeg = true,
		bool logY = true,
		bool blindData = true) {
  
  SetStyle();

  TString dir("/nfs/dust/cms/user/rasp/Run/Run2018/H2aa/");

  TString samples[5] = {"SUSYGluGluToHToAA_AToTauTau_M-125",
			"SUSYVBFToHToAA_AToTauTau_M-125",
			"SUSYVH_HToAA_AToTauTau_M-125",
			"SUSYttH_HToAA_AToTauTau_M-125",
			"SUSYGluGluToHToAA_AToMuMu_AToTauTau"};

  double xsecGGH = 48.52;
  double xsecVBF = 3.779;
  double xsecVH  = 1.369 + 0.8824;
  double xsecTTH = 0.5065;
  double massTau = 1.777;
  double massMu  = 0.106;
  double massRatio = (massMu*massMu)/(massTau*massTau);  

  double xsec[5] = {xsecGGH,xsecVBF,xsecVH,xsecTTH,xsecGGH};
  TString massT[4] = {"4","7","10","15"};
  double massA[4] = {4,7,10,15};

  int nBinsNew = 6;
  double bins[7] = {0,1,2,3,4,6,12};

  TH1D * histMass[4];
  TH1D * histMass4Tau[4];
  
  TFile * file = new TFile(dir+"/DoubleMuon_Run2018.root");
  TH1D * histOld   = (TH1D*)file->Get("InvMassH");
  for (int iM=0; iM<4; ++iM) {
    double aF = 2*massTau/massA[iM];
    double SF = 2*massRatio/TMath::Sqrt(1-aF*aF);
    double xsecGGH_mmtt = xsecGGH * SF;
    histMass[iM] = new TH1D("hist"+massT[iM],"",nBinsNew,bins);
    histMass4Tau[iM] = new TH1D("hist4Tau"+massT[iM],"",nBinsNew,bins);
    xsec[4] = xsecGGH_mmtt;
    for (int iS=0; iS<5; ++iS) {
      TFile * fileSig = new TFile(dir+"/"+samples[iS]+"_M-"+massT[iM]+".root");
      TH1D * histWeightH = (TH1D*)fileSig->Get("histWeightsH");
      double ngen = histWeightH->GetSumOfWeights();
      double norm = xsec[iS]/ngen;
      TH1D * histSigOld = (TH1D*)fileSig->Get("InvMassH");
      TH1D * histSig = TH1DtoTH1D(histSigOld,nBinsNew,bins,true,"_"+samples[iS]+massT[iM]);
      histMass[iM]->Add(histMass[iM],histSig,1.0,norm); 
      if (iS<4) histMass4Tau[iM]->Add(histMass4Tau[iM],histSig,1.0,norm);
    }
  }


  TH1D * hist4  = histMass[0];
  TH1D * hist7  = histMass[1];
  TH1D * hist10 = histMass[2];
  TH1D * hist15 = histMass[3];

  TH1D * hist4Tau  = histMass4Tau[0];
  TH1D * hist7Tau  = histMass4Tau[1];
  TH1D * hist10Tau = histMass4Tau[2];
  TH1D * hist15Tau = histMass4Tau[3];

  double norm4  = hist4->GetSumOfWeights();
  double norm7  = hist7->GetSumOfWeights();
  double norm10 = hist10->GetSumOfWeights();
  double norm15 = hist15->GetSumOfWeights();

  std::cout << "norm4  = " << norm4 << std::endl;
  std::cout << "norm7  = " << norm7 << std::endl;
  std::cout << "norm10 = " << norm10 << std::endl;
  std::cout << "norm15 = " << norm15 << std::endl;

  TH1D * histN23Old = (TH1D*)file->Get("InvMassN23H");
  TH1D * histN45Old = (TH1D*)file->Get("InvMassN45H");

  TH1D * histN23 = (TH1D*)TH1DtoTH1D(histN23Old,nBinsNew,bins,true,"_new");
  TH1D * histN45 = (TH1D*)TH1DtoTH1D(histN45Old,nBinsNew,bins,true,"_new");
  TH1D * hist    = (TH1D*)TH1DtoTH1D(histOld,nBinsNew,bins,true,"_new");
  TH1D * histBkg     = (TH1D*)histN23->Clone("histBkg");
  TH1D * histBkgUp   = (TH1D*)histN23->Clone("histBkgUp");
  TH1D * histBkgDown = (TH1D*)histN23->Clone("histBkgDown");

  float normN23 = histN23->GetSumOfWeights();
  float normN45 = histN45->GetSumOfWeights();
  double norm      = hist->GetSumOfWeights();  

  for (int i=1; i<=nBinsNew; ++i) {
    float xN23 = histN23->GetBinContent(i)/normN23;
    float eN23 = histN23->GetBinError(i)/normN23;
    float xN45 = histN45->GetBinContent(i)/normN45;

    float xcorr = xN45/xN23;
    float xCentral = xN23;
    histBkg->SetBinContent(i,xCentral);
    histBkg->SetBinError(i,eN23);

    float xUp = xcorr*xCentral;
    float xDown = xCentral/xcorr;
    histBkgUp->SetBinContent(i,xUp);
    histBkgDown->SetBinContent(i,xDown);
  }


  double normBkg   = histBkg->GetSumOfWeights();
  double normBkgUp = histBkgUp->GetSumOfWeights();

  for (int i=1; i<=nBinsNew; ++i) {
    float binWidth = hist->GetBinLowEdge(i+1)-hist->GetBinLowEdge(i);
    float xC = histBkg->GetBinContent(i);
    float xU = histBkgUp->GetBinContent(i);
    float errStat = histBkg->GetBinError(i);
    float errSys = 2*(xU-xC);
    float eTot = TMath::Sqrt(errSys*errSys+errStat*errStat);
    hist->SetBinContent(i,hist->GetBinContent(i)/norm);
    hist->SetBinError(i,hist->GetBinError(i)/norm);
    histBkg->SetBinContent(i,xC);
    histBkg->SetBinError(i,eTot);
    hist4->SetBinContent(i,hist4->GetBinContent(i)/norm4);
    hist7->SetBinContent(i,hist7->GetBinContent(i)/norm7);
    hist10->SetBinContent(i,hist10->GetBinContent(i)/norm10);
    hist15->SetBinContent(i,hist15->GetBinContent(i)/norm15);
    hist4Tau->SetBinContent(i,hist4Tau->GetBinContent(i)/norm4);
    hist7Tau->SetBinContent(i,hist7Tau->GetBinContent(i)/norm7);
    hist10Tau->SetBinContent(i,hist10Tau->GetBinContent(i)/norm10);
    hist15Tau->SetBinContent(i,hist15Tau->GetBinContent(i)/norm15);

    hist4->SetBinError(i,0);
    hist7->SetBinError(i,0);
    hist10->SetBinError(i,0);
    hist15->SetBinError(i,0);
    hist4Tau->SetBinError(i,0);
    hist7Tau->SetBinError(i,0);
    hist10Tau->SetBinError(i,0);
    hist15Tau->SetBinError(i,0);
  }

  std::cout << "norm4  = " << hist4->GetSumOfWeights() << std::endl;
  std::cout << "norm7  = " << hist7->GetSumOfWeights() << std::endl;
  std::cout << "norm10 = " << hist10->GetSumOfWeights() << std::endl;
  std::cout << "norm15 = " << hist15->GetSumOfWeights() << std::endl;

  InitData(hist);
  hist->SetMarkerSize(1.3);

  TH1D * histErr = (TH1D*)histBkg->Clone("histErr");
  histErr->SetFillStyle(3013);
  histErr->SetFillColor(1);
  histErr->SetMarkerStyle(21);
  histErr->SetMarkerSize(0);
  histErr->SetLineColor(4);
  histErr->SetLineWidth(3);

  // ratio histograms
  TH1D * ratioH = (TH1D*)hist->Clone("ratioH");
  TH1D * ratioErrH = (TH1D*)histErr->Clone("ratioErr");
  for (int iB=1; iB<=nBinsNew; ++iB) {
    double data = hist->GetBinContent(iB);
    double dataE = hist->GetBinError(iB);
    double bkg   = histErr->GetBinContent(iB);
    double bkgE  = histErr->GetBinError(iB);    
    histBkg->SetBinError(iB,0);
    double ratio = data/bkg;
    double ratioE = dataE/bkg;
    double ratioBkgE = bkgE/bkg;
    ratioH->SetBinContent(iB,ratio);
    ratioH->SetBinError(iB,ratioE);
    ratioErrH->SetBinContent(iB,1);
    ratioErrH->SetBinError(iB,ratioBkgE);
  }
  ratioErrH->GetYaxis()->SetRangeUser(0.501,1.499);
  ratioErrH->GetXaxis()->SetNdivisions(210);
  ratioErrH->GetYaxis()->SetNdivisions(205);
  ratioErrH->GetXaxis()->SetLabelFont(42);
  ratioErrH->GetXaxis()->SetLabelOffset(0.04);
  ratioErrH->GetXaxis()->SetLabelSize(0.14);
  ratioErrH->GetXaxis()->SetTitleSize(0.13);
  ratioErrH->GetXaxis()->SetTitleOffset(1.2);
  ratioErrH->GetYaxis()->SetTitle("obs/bkg");
  ratioErrH->GetXaxis()->SetTitle(xtitle);
  ratioErrH->GetYaxis()->SetLabelFont(42);
  ratioErrH->GetYaxis()->SetLabelOffset(0.015);
  ratioErrH->GetYaxis()->SetLabelSize(0.13);
  ratioErrH->GetYaxis()->SetTitleSize(0.14);
  ratioErrH->GetYaxis()->SetTitleOffset(0.5);
  ratioErrH->GetXaxis()->SetTickLength(0.07);
  ratioErrH->GetYaxis()->SetTickLength(0.04);
  ratioErrH->GetYaxis()->SetLabelOffset(0.01);

  histBkg->GetYaxis()->SetRangeUser(0,0.5);
  if (logY) histBkg->GetYaxis()->SetRangeUser(0.5e-2,5.001);
  histBkg->SetLineColor(kBlue);
  histBkg->SetLineWidth(3);
  histBkg->SetLineStyle(1);
  histBkg->GetYaxis()->SetTitleOffset(1.2);
  histBkg->GetXaxis()->SetTitle(xtitle);
  histBkg->GetYaxis()->SetTitle(ytitle);
  histBkg->GetXaxis()->SetLabelSize(0.);

  histBkg->SetLineColor(kBlue);
  histBkg->SetLineWidth(3);
  histBkg->SetLineStyle(1);
  histBkg->GetXaxis()->SetNdivisions(210);

  hist4->SetLineColor(kMagenta);
  hist4->SetLineWidth(3);
  hist4->SetLineStyle(1);

  hist4Tau->SetLineColor(kMagenta);
  hist4Tau->SetLineWidth(3);
  hist4Tau->SetLineStyle(2);

  hist7->SetLineColor(kRed);
  hist7->SetLineWidth(3);
  hist7->SetLineStyle(1);

  hist7Tau->SetLineColor(kRed);
  hist7Tau->SetLineWidth(3);
  hist7Tau->SetLineStyle(2);

  hist10->SetLineColor(kGreen+1);
  hist10->SetLineWidth(3);
  hist10->SetLineStyle(1);

  hist10Tau->SetLineColor(kGreen+1);
  hist10Tau->SetLineWidth(3);
  hist10Tau->SetLineStyle(2);

  hist15->SetLineColor(kOrange);
  hist15->SetLineWidth(3);
  hist15->SetLineStyle(1);

  hist15Tau->SetLineColor(kOrange);
  hist15Tau->SetLineWidth(3);
  hist15Tau->SetLineStyle(2);

  if (blindData) {
    hist->SetBinContent(5,1e-13);
    hist->SetBinContent(6,1e-13);
    hist->SetBinError(5,1e-13);
    hist->SetBinError(6,1e-13);
    ratioH->SetBinContent(5,1e-13);
    ratioH->SetBinContent(6,1e-13);
    ratioH->SetBinError(5,1e-13);
    ratioH->SetBinError(6,1e-13);
  }    

  TCanvas * canv = MakeCanvas("canv","",600,700);
  TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(5);
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

  std::cout << "Background : " << histBkg->GetSumOfWeights() << std::endl;

  hist4->Scale(1.05);

  histBkg->Draw("hist");
  histErr->Draw("e2same");
  hist4->Draw("hsame");
  hist7->Draw("hsame");
  hist10->Draw("hsame");
  hist15->Draw("hsame");
  //  hist4Tau->Draw("hsame");
  //  hist7Tau->Draw("hsame");
  //  hist10Tau->Draw("hsame");
  //  hist15Tau->Draw("hsame");
  histBkg->Draw("hsame");
  hist->Draw("e1same");

  TLegend * leg = new TLegend(0.20,0.62,0.5,0.86);
  SetLegendStyle(leg);
  leg->SetTextSize(0.04);
  TLegend * leg1 = new TLegend(0.55,0.62,0.85,0.86);
  SetLegendStyle(leg1);
  leg1->SetTextSize(0.04);
  if (!blindData)
    leg->AddEntry(hist,"observed","lp");
  leg->AddEntry(histErr,"bkg(+unc.)","lf");
  leg->AddEntry(hist4,"m_{a_{1}} = 4 GeV","l");
  leg1->AddEntry(hist7,"m_{a_{1}} = 7 GeV","l");
  leg1->AddEntry(hist10,"m_{a_{1}} = 10 GeV","l");
  leg1->AddEntry(hist15,"m_{a_{1}} = 15 GeV","l");
  if (drawLeg) {leg->Draw(); leg1->Draw();}
  writeExtraText = false;
  extraText   = "Preliminary";
  CMS_lumi(upper,4,33); 

  if (logY) upper->SetLogy(true);
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv->cd();

  // ------------>Primitives in pad: lower                                                                                                                                                              
  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(5);
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

  ratioErrH->Draw("e2");
  ratioH->Draw("e1same");
  TLine * line = new TLine(0,1,12,1);
  line->SetLineColor(4);
  line->Draw();
  ratioH->Draw("e1same");

  lower->Modified();
  lower->RedrawAxis();
  canv->cd();
  canv->SetSelected(canv);

  canv->Update();
  canv->Print("figures/mass1D.png");
  //  canv->Print("mass1D.pdf","Portrait pdf");

}

  
