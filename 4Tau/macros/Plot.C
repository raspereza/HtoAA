#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"

// these histograms are filled after applying same-sign dimuon selection
// no further cuts are applied
//
// muon kinematics ->
//TH1D * ptLeadingMuH = new TH1D("ptLeadingMuH","",50,0,100);
//TH1D * ptTrailingMuH = new TH1D("ptTrailingMuH","",50,0,100);
//TH1D * etaLeadingMuH = new TH1D("etaLeadingMuH","",50,-2.5,2.5);
//TH1D * etaTrailingMuH = new TH1D("etaTrailingMuH","",50,-2.5,2.5);
//TH1D * dimuonMassH = new TH1D("dimuonMassH","",500,0,500);
//
// number of tracks within dR<0.4 around muons
// we select events where there is only one track accompanies each of muons   
//TH1D * nTracksLeadingMuH = new TH1D("nTracksLeadingMuH","",21,-0.5,20.5);
//TH1D * nTracksTrailingMuH = new TH1D("nTracksTrailingMuH","",21,-0.5,20.5); 


void Plot(TString histName = "ptTrailingMuH", // histogram name
	  TString xtitle = "p_{T}^{#mu2} [GeV]", // title of x axis
	  TString ytitle = "Events / 2 GeV", // title of y axis
	  TString Mass = "10", // ma [GeV]
	  float xLower = 10, // lower boundary of x axis 
	  float xUpper = 100, // upper boundary of x axis
	  float yLower = 1000, // lower boundary of y axis (in case when plotting y axis in log scale)
	  int nBinsNew = 50, // new number of bins
	  bool logY = true, // log or linear scale of Y axis
	  bool drawLeg = true) { // draw or not legend

  // directory where RooT files are located
  TString dir = "/nfs/dust/cms/user/rasp/Run/Run2018/H2aa";
  lumi_13TeV = "2018, 59.7 fb^{-1}";

  float qcdScale = 0.60; // add-hoc K-factor for QCD background

  // mu = sigma(pp->H+X)*BR(H->aa->4tau)/sigma(pp->H+X,SM) (can be > 1 to visualise signal)
  float mu = 10.0; 
  
  // Cross sections in pb are given at mH = 125.09 GeV
  // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
  float ggHXSec = 48.61 * mu;
  float VBFXSec = 3.766 * mu;
  float VHXSec = (1.358 + 0.880)*mu;
  float ttHXSec = 0.5065*mu;

  char muChar[5];
  sprintf(muChar,"%1i",int(mu));
  TString Mu(muChar);

  SetStyle();

  TFile * file = new TFile(dir+"/DoubleMuon_Run2018.root");  

  TString samples[26] = {"WW_13TeV-pythia8", // (0)
			 "WZ_13TeV-pythia8", // (1) 
			 "ZZ_13TeV-pythia8", // (2)
			 "WJetsToLNu",       // (3)
			 "ST_t-channel_antitop", // (4)
			 "ST_t-channel_top",     // (5)
			 "ST_tW_antitop",     // (6)
			 "ST_tW_top",         // (7)
			 "TTTo2L2Nu",         // (8)
			 "TTToSemiLeptonic",  // (9)
			 "TTToHadronic",      // (10)
			 "DYJetsToLL_M-50",   // (11)
			 "QCD_Pt-20to30_MuEnrichedPt5",    // (12)
			 "QCD_Pt-30to50_MuEnrichedPt5",    // (13)
			 "QCD_Pt-50to80_MuEnrichedPt5",    // (14)
			 "QCD_Pt-80to120_MuEnrichedPt5",   // (15)
			 "QCD_Pt-120to170_MuEnrichedPt5",  // (16)
			 "QCD_Pt-170to300_MuEnrichedPt5",   // (17)
			 "QCD_Pt-300to470_MuEnrichedPt5", // (18)
			 "QCD_Pt-470to600_MuEnrichedPt5", // (19)
			 "QCD_Pt-600to800_MuEnrichedPt5", // (20)
			 "QCD_Pt-800to1000_MuEnrichedPt5", // (21)
			 "SUSYGluGluToHToAA_AToTauTau_M-125_M-"+Mass, // (22)
			 "SUSYVBFToHToAA_AToTauTau_M-125_M-"+Mass,      // (23)
			 "SUSYVH_HToAA_AToTauTau_M-125_M-"+Mass, // (24)
			 "SUSYttH_HToAA_AToTauTau_M-125_M-"+Mass // (25)
  };

  float xsec[26] = {
    115.0,   // WW (0)
    47.13,   // WZ (1)
    16.52,   // ZZ (2)
    61526.7, // WJets (3)
    80.95,   // ST_t-channel_antitop (4)
    136.95,  // ST_t-channel_top     (5)
    35.6,    // ST_tW_antitop        (6)
    35.6,    // ST_tW_top            (7)
    88.29,   // TTbarToLN            (8)
    365.35,  // TTbarSemileptonic    (9)
    377.96,  // TTbarHadronic       (10)    
    5765,    // DYJetsToLL_M-50     (11)
    558528000*0.0053,  // (12)
    139803000*0.01182, // (13)
    19222500*0.02276,  // (14)
    2758420*0.03844,   // (15)
    469797*0.05362,    // (16)
    117989*0.07335,    // (17)
    7820.25*0.10196,   // (18)
    645.528*0.12242,   // (19)
    187.109*0.13412,   // (20)
    32.3486*0.14552,   // (21)
    ggHXSec,  // (22)
    VBFXSec,  // (23)
    VHXSec,   // (24)
    ttHXSec
  };

  for (int iQCD=12; iQCD<22; ++iQCD)
    xsec[iQCD] *= qcdScale;

  float lumi = 59740;

  TH1D * histDataOld = (TH1D*)file->Get(histName);

  int nBins = histDataOld->GetNbinsX();
  float xMin = histDataOld->GetBinLowEdge(1);
  float xMax = histDataOld->GetBinLowEdge(nBins+1);

  std::cout << std::endl;
  std::cout << "Histogram " << histName << " : " << "nbins = " << nBins
	    << " , min = " << xMin
	    << " , max = " << xMax << std::endl;
  std::cout << std::endl;

  // new binning
  if (nBins%nBinsNew!=0) { 
    std::cout << "new number of bins = " << nBinsNew 
	      << "  not multiple of " << nBins << std::endl;
    return;
  }
  double bins[300];
  float binWidth = (xMax-xMin)/float(nBinsNew);
  for (int iB=0; iB<=nBinsNew; ++iB)
    bins[iB] = xMin + float(iB)*binWidth;

  TH1D * histData = TH1DtoTH1D(histDataOld,nBinsNew,bins,true,"_Data_new");

  TH1D * ewkHist = new TH1D("ewkHist","",nBinsNew,bins);
  TH1D * ttHist  = new TH1D("ttHist", "",nBinsNew,bins);
  TH1D * qcdHist = new TH1D("qcdHist","",nBinsNew,bins);
  TH1D * zHist   = new TH1D("zHist",  "",nBinsNew,bins);
  TH1D * sigHist = new TH1D("sigHist","",nBinsNew,bins);

  int nSamples = 26;

  for (int iS=0; iS<nSamples; ++iS) {
    TFile * fileMC = new TFile(dir+"/"+samples[iS]+".root");
    TH1D * histOld = (TH1D*)fileMC->Get(histName);
    TH1D * hist = TH1DtoTH1D(histOld,nBinsNew,bins,true,"_new_"+samples[iS]);
    TH1D * eventCount = (TH1D*)fileMC->Get("histWeightsH");
    float nGen = eventCount->GetSumOfWeights();
    float norm = xsec[iS]*lumi/nGen;
    TH1D * tempHist = ewkHist;
    if (iS>7&&iS<11)
      tempHist = ttHist;
    if (iS==11)
      tempHist = zHist;
    if (iS>11&&iS<22)
      tempHist = qcdHist;
    if (iS>=22)
      tempHist = sigHist;
    tempHist->Add(tempHist,hist,1.,norm);

  }

  // Simple systematic uncertainties
  float lumiSys = 0.03;
  float lepSys  = 0.04;

  float ttSys  = 0.05;
  float ewkSys = 0.10;
  float qcdSys = 0.20;
  float zSys   = 0.05;

  float ttE2 = 0;
  float ewkE2 = 0;
  float qcdE2 = 0;
  float zE2 = 0;

  float ttTot = 0;
  float ewkTot = 0;
  float qcdTot = 0;
  float zTot = 0;

  for (int iB=1; iB<=nBinsNew; ++iB) {

    float qcdX = qcdHist->GetBinContent(iB);
    float qcdE = qcdHist->GetBinError(iB);
    
    float ewkX = ewkHist->GetBinContent(iB);
    float ewkE = ewkHist->GetBinError(iB);
    
    float ttX  = ttHist->GetBinContent(iB);
    float ttE  = ttHist->GetBinError(iB);
    
    float zX  = zHist->GetBinContent(iB);
    float zE  = zHist->GetBinError(iB);
    
    ttE2 += ttE*ttE;
    qcdE2 += qcdE*qcdE;
    ewkE2 += ewkE*ewkE;
    zE2 += zE*zE;

    ttTot += ttX;
    qcdTot += qcdX;
    ewkTot += ewkX;
    zTot += zX;

    if (zX<0) zX = 0;
    
    float ttErr   = ttX*ttSys;
    float ewkErr  = ewkX*ewkSys;
    float qcdErr  = qcdX*qcdSys;
    float zErr    = zX*zSys;
    
    ewkX += zX;
    ttX  += ewkX;
    qcdX += ttX;

    float lumiErr = qcdX*lumiSys;
    float lepErr  = qcdX*lepSys;

    float totErr = TMath::Sqrt(lumiErr*lumiErr+
			       lepErr*lepErr+
			       ttErr*ttErr+
			       ewkErr*ewkErr+
			       qcdErr*qcdErr+
			       ewkE*ewkE+
			       zE*zE+
			       ttE*ttE+
			       qcdE*qcdE);
    
    ewkHist->SetBinContent(iB,ewkX);
    ttHist->SetBinContent(iB,ttX);
    zHist->SetBinContent(iB,zX);
    qcdHist->SetBinContent(iB,qcdX);
    qcdHist->SetBinError(iB,totErr);
  }

  float ttE = TMath::Sqrt(ttE2);
  float ewkE = TMath::Sqrt(ewkE2);
  float qcdE = TMath::Sqrt(qcdE2);
  float zE = TMath::Sqrt(zE2);

  float dataYieldE = TMath::Sqrt(histData->GetSumOfWeights());
  
  std::cout << std::endl;
  std::cout << "QCD : " << qcdTot << " +/- " << qcdE << std::endl;
  std::cout << "EWK : " << ewkTot << " +/- " << ewkE << std::endl;
  std::cout << "TTJ : " << ttTot  << " +/- " << ttE << std::endl;
  std::cout << "DY  : " << zTot   << " +/- " << zE << std::endl;

  float totBkg = qcdTot + ewkTot + ttTot + zTot;
  float eTotBkg = TMath::Sqrt(qcdE*qcdE+ewkE*ewkE+ttE*ttE+zE*zE);
  std::cout << std::endl;
  std::cout << "Bkgd   : " << totBkg << " +/- " << eTotBkg << std::endl;
  std::cout << "Data   : " << histData->GetSumOfWeights() << std::endl;
  std::cout << "Signal : " << sigHist->GetSumOfWeights() << std::endl;

  TH1D * bkgdErr = (TH1D*)qcdHist->Clone("bkgdErr");
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);
  bkgdErr->SetFillStyle(3444);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);  
  
  for (int iB=1; iB<=nBinsNew; ++iB) {
    qcdHist->SetBinError(iB,0);
    ewkHist->SetBinError(iB,0);
    ttHist->SetBinError(iB,0);
    zHist->SetBinError(iB,0);
    sigHist->SetBinError(iB,0);
  }

  InitData(histData);
  InitSignal(sigHist);
  sigHist->SetLineColor(kBlue);
  sigHist->SetLineStyle(2);
  sigHist->SetLineWidth(2);
  InitHist(qcdHist,"","",TColor::GetColor("#FFCCFF"),1001);
  InitHist(ewkHist,"","",TColor::GetColor("#DE5A6A"),1001);
  InitHist(ttHist,"","",TColor::GetColor("#9999CC"),1001);
  //  InitHist(wHist,"","",TColor::GetColor("#6F2D35"),1001);
  InitHist(zHist,"","",TColor::GetColor("#FFCC66"),1001);
  histData->GetXaxis()->SetTitle(xtitle);
  histData->GetYaxis()->SetTitle(ytitle);
  histData->GetYaxis()->SetTitleOffset(1.3);
  histData->GetYaxis()->SetTitleSize(0.06);
  histData->GetXaxis()->SetRangeUser(xLower,xUpper);
  float yUpper = histData->GetMaximum();
  if (logY)
    histData->GetYaxis()->SetRangeUser(yLower,20*yUpper);
  else
    histData->GetYaxis()->SetRangeUser(0,1.2*yUpper);

  histData->SetMarkerSize(1.2);
  histData->GetXaxis()->SetLabelSize(0);
  histData->GetYaxis()->SetLabelSize(0.06);

  //  nData = histData->GetSum();
  //  float nMC   = ttHist->GetSum();
  //  float eData = TMath::Sqrt(nData);

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

  histData->Draw("e1");
  qcdHist->Draw("sameh");
  ttHist->Draw("sameh");
  ewkHist->Draw("sameh");
  histData->Draw("e1same");
  bkgdErr->Draw("e2same");
  sigHist->Draw("same");
  histData->Draw("e1same");
  float chi2 = 0;
  for (int iB=1; iB<=nBinsNew; ++iB) {
    float xData = histData->GetBinContent(iB);
    float xMC = qcdHist->GetBinContent(iB);
    if (xMC>1e-1) {
      float diff2 = (xData-xMC)*(xData-xMC);
      chi2 += diff2/xMC;
    }
  }
  std::cout << std::endl;
  std::cout << "Chi2 = " << chi2 << std::endl;
  std::cout << std::endl;

  TLegend * leg = new TLegend(0.52,0.5,0.84,0.78);
  SetLegendStyle(leg);
  leg->SetTextSize(0.04);
  leg->AddEntry(histData,"Observed","lp");
  leg->AddEntry(qcdHist,"QCD multijets","f");
  //  leg->AddEntry(zHist,"Drell-Yan","f");
  leg->AddEntry(ttHist,"t#bar{t}+single top","f");
  leg->AddEntry(ewkHist,"electroweak","f");
  leg->AddEntry(sigHist,"H#rightarrow aa#rightarrow 4#tau(#mu="+Mu+",m_{a}="+Mass+")","l");
  if (drawLeg) leg->Draw();
  writeExtraText = true;
  extraText   = "Internal";
  CMS_lumi(upper,4,33); 
  //  plotchannel("same-sign #mu#mu");

  if (logY) 
    upper->SetLogy(true);
    
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

  TH1D * ratioH = (TH1D*)histData->Clone("ratioH");
  TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.2);
  ratioH->SetLineColor(1);
  ratioH->GetYaxis()->SetRangeUser(0.0001,1.999);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleSize(0.13);
  ratioH->GetXaxis()->SetTitleOffset(1.2);
  ratioH->GetYaxis()->SetTitle("obs/bkg");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.015);
  ratioH->GetYaxis()->SetLabelSize(0.13);
  ratioH->GetYaxis()->SetTitleSize(0.14);
  ratioH->GetYaxis()->SetTitleOffset(0.55);
  ratioH->GetXaxis()->SetTickLength(0.07);
  ratioH->GetYaxis()->SetTickLength(0.04);
  ratioH->GetYaxis()->SetLabelOffset(0.01);
  ratioH->GetYaxis()->SetNdivisions(505);

  for (int iB=1; iB<=nBinsNew; ++iB) {
    float x1 = histData->GetBinContent(iB);
    float x2 = qcdHist->GetBinContent(iB);
    ratioErrH->SetBinContent(iB,1.0);
    ratioErrH->SetBinError(iB,0.0);
    float xBkg = bkgdErr->GetBinContent(iB);
    float errBkg = bkgdErr->GetBinError(iB);
    if (xBkg>0) {
      float relErr = errBkg/xBkg;
      float err = TMath::Sqrt(relErr*relErr);
      ratioErrH->SetBinError(iB,err);

    }
    if (x1>0&&x2>0) {
      float e1 = histData->GetBinError(iB);
      float ratio = x1/x2;
      float eratio = e1/x2;
      ratioH->SetBinContent(iB,ratio);
      ratioH->SetBinError(iB,eratio);
    }
    else {
      ratioH->SetBinContent(iB,1000);
    }
  }


  // ------------>Primitives in pad: lower
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

  ratioH->Draw("e1");
  ratioErrH->Draw("e2same");
  TLine * line = new TLine(xLower,1,xUpper,1);
  line->SetLineColor(2);
  line->SetLineStyle(3);
  line->Draw();
  ratioH->Draw("e1same");

  lower->Modified();
  lower->RedrawAxis();
  canv1->cd();
  canv1->SetSelected(canv1);

  canv1->Print(histName+".png");

}
