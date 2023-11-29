#include "HttStylesNew.cc"
#include "HtoH.h"
#include "CMS_lumi.C"
//  TH1D * hist1Dold = (TH1D*)file->Get("InvMassTrackPlusMuon1D_"+Suffix);
//  TH2D * hist2Dold = (TH2D*)file->Get("InvMassTrackPlusMuon2D_"+Suffix);
//  TH1D * hist1Dold = (TH1D*)file->Get("ModelInvMassH");
//  TH2D * hist2Dold = (TH2D*)file->Get("ModelInvMass2DH");

void CorrCoefficientsMC_Era(
			     TString era = "2018",
			     bool signalRegion = false,
			     bool btagVeto = true
			     ) {


  std::map<TString,double> eraLumi = {
    {"2016_preVFP", 19520},
    {"2016_postVFP",16810},
    {"2017",        41480},
    {"2018",        59830}
  };



  std::map<TString,TString> eraLabel = {
    {"2016_preVFP", "2016 (preVFP), 19.5 fb^{-1}"},
    {"2016_postVFP","2016 (postVFP), 16.8 fb^{-1}"},
    {"2017",        "2017 41.5 fb^{-1}"},
    {"2018",        "2018 59.8 fb^{-1}"}
  };

  lumi_13TeV = eraLabel[era];

  TString dir("/nfs/dust/cms/user/rasp/Run/Run"+era+"_UL/mutrk");
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  SetStyle();

  TString baseName1D = "InvMass";
  TString baseName2D = "InvMass2D";
  TString Suffix = "_ControlYH";
  if (signalRegion) 
    Suffix = "H";    

  //  TFile * file = new TFile(dir+"/QCD_Pt-20toInf_MuEnrichedPt15.root");
  //  TFile * file = new TFile(dir+"/QCD_Pt-15To20_MuEnrichedPt5.root");
  TFile * file = new TFile(dir+"/QCD_Pt-20To30_MuEnrichedPt5.root");
  TH1D * hist1Dold = (TH1D*)file->Get(baseName1D+Suffix);
  TH2D * hist2Dold = (TH2D*)file->Get(baseName2D+Suffix);
  TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
  std::cout << hist1Dold << " " << hist2Dold << " " << " " << histWeightsH << std::endl;

  //  double norm = 720648000 * 0.00042/histWeightsH->GetSumOfWeights(); // QCD_Pt-20toInf_MuEnrichedPt15
  //  double norm = 1273190000 * 0.003/histWeightsH->GetSumOfWeights(); // QCD_Pt-15To20_MuEnrichedPt5
  double norm = 558528000 * 0.0053/histWeightsH->GetSumOfWeights(); // QCD_Pt-20To30_MuEnrichedPt5
  int nBins = hist1Dold->GetNbinsX();
  for (int iB=1; iB<=nBins; ++iB) {
    hist1Dold->SetBinContent(iB,norm*hist1Dold->GetBinContent(iB));
    hist1Dold->SetBinError(iB,norm*hist1Dold->GetBinError(iB));
    for (int jB=1; jB<=nBins; ++jB) {
      hist2Dold->SetBinContent(iB,jB,norm*hist2Dold->GetBinContent(iB,jB));    
      hist2Dold->SetBinError(iB,jB,norm*hist2Dold->GetBinError(iB,jB));
    }
  }
  
  TString SamplesMuEnrichedQCD[11] = {
    "QCD_Pt-15To20_MuEnrichedPt5",
    "QCD_Pt-20To30_MuEnrichedPt5",
    "QCD_Pt-30To50_MuEnrichedPt5",
    "QCD_Pt-50To80_MuEnrichedPt5",
    "QCD_Pt-80To120_MuEnrichedPt5",
    "QCD_Pt-120To170_MuEnrichedPt5",
    "QCD_Pt-170To300_MuEnrichedPt5",
    "QCD_Pt-300To470_MuEnrichedPt5",
    "QCD_Pt-470To600_MuEnrichedPt5",
    "QCD_Pt-600To800_MuEnrichedPt5",
    "QCD_Pt-800To1000_MuEnrichedPt5"
  };
  
  double xsecMuEnrichedQCD[11] = {
    1273190000*0.003,
    558528000*0.0053,
    139803000*0.01182,
    19222500*0.02276,
    2758420*0.03844,
    469797*0.05362,
    117989*0.07335,
    7820.25*0.10196,
    645.528*0.12242,
    187.109*0.13412,
    32.3486*0.14552
  };

  int nBinsNew;
  double *bins;
  double *binsCorr;

  if (btagVeto) {
    nBinsNew = 5;
    double tempBins[6] = {0, 1, 2, 3, 4, 20};
    double tempBinsCorr[6] = {0, 1, 2, 3, 4, 12};
    bins = tempBins;
    binsCorr = tempBinsCorr;
  } else {
    nBinsNew = 6;
    double tempBins[7] = {0, 1, 2, 3, 4, 6, 20};
    double tempBinsCorr[7] = {0, 1, 2, 3, 4, 6, 12};
    bins = tempBins;
    binsCorr = tempBinsCorr;
  }

  TH1D * hist1D = (TH1D*)TH1DtoTH1D(hist1Dold,nBinsNew,bins,true,"_new");
  TH2D * hist2D = (TH2D*)TH2DtoTH2D(hist2Dold,nBinsNew,bins,nBinsNew,bins,"_new");

  for (int iS=2; iS<11; ++iS) {
    TFile * fileSample = new TFile(dir+"/"+SamplesMuEnrichedQCD[iS]+".root");
    TH1D * hist1DoldSample = (TH1D*)fileSample->Get(baseName1D+Suffix);
    TH2D * hist2DoldSample = (TH2D*)fileSample->Get(baseName2D+Suffix);
    TH1D * histWeightsSampleH = (TH1D*)fileSample->Get("histWeightsH");
    double normSample = xsecMuEnrichedQCD[iS]/histWeightsSampleH->GetSumOfWeights();
    for (int iB=1; iB<=nBins; ++iB) {
      hist1DoldSample->SetBinContent(iB,normSample*hist1DoldSample->GetBinContent(iB));
      hist1DoldSample->SetBinError(iB,normSample*hist1DoldSample->GetBinError(iB));
      for (int jB=1; jB<=nBins; ++jB) {
	hist2DoldSample->SetBinContent(iB,jB,normSample*hist2DoldSample->GetBinContent(iB,jB));
	hist2DoldSample->SetBinError(iB,jB,normSample*hist2DoldSample->GetBinError(iB,jB));
      }
    }
    TH1D * hist1DSample = (TH1D*)TH1DtoTH1D(hist1DoldSample,nBinsNew,bins,true,SamplesMuEnrichedQCD[iS]+"_new");
    TH2D * hist2DSample = (TH2D*)TH2DtoTH2D(hist2DoldSample,nBinsNew,bins,nBinsNew,bins,SamplesMuEnrichedQCD[iS]+"_new");
    hist1D->Add(hist1D,hist1DSample);
    hist2D->Add(hist2D,hist2DSample);

  }

  TH2D * corrCoeff = new TH2D("corrCoeff","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  TH2D * corrCoeffX = new TH2D("corrCoeffX","",nBinsNew,bins,nBinsNew,bins);
  
  hist1D->Scale(1/hist1D->GetSum());
  hist2D->Scale(1/hist2D->GetSum());

  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = hist1D->GetBinContent(iB);
      double y = hist1D->GetBinContent(jB);
      double denominator = x*y;
      if (iB!=jB)
	denominator *= 2.0;
      double numerator  = hist2D->GetBinContent(jB,iB);
      double enumerator = hist2D->GetBinError(jB,iB);
      if (iB!=jB) { 
	numerator += hist2D->GetBinContent(iB,jB);
	double enumerator2 = hist2D->GetBinError(iB,jB);
	enumerator = TMath::Sqrt(enumerator*enumerator+enumerator2*enumerator2);
      }
      double corr  = numerator / denominator ;
      double ecorr = enumerator / denominator;
      if (iB==5&&jB==5&&!signalRegion) {
	corr *= 1.15;
	ecorr *= 1.07;
      }
      corr = floor(100*corr+0.5)/100;
      ecorr = floor(100*ecorr+0.5)/100;
      printf("[%1i,%1i] = %5.2f +/- %5.2f \n",iB,jB,corr,ecorr);
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
    double yL = binsCorr[i];
    TLine * lineX = new TLine(0,yL,12,yL);
    TLine * lineY = new TLine(yL,0,yL,12);
    lineX->SetLineWidth(1);
    lineY->SetLineWidth(1);
    lineX->Draw();
    lineY->Draw();
  }

  writeExtraText = "Simulation";
  CMS_lumi(canv,4,33);
  
  canv->Update();

  TString suffix("control");
  if (signalRegion)
    suffix = "signal";
  suffix += "_mc";

  canv->Print("CorrCoefficients_"+suffix+"_"+era+".png");

  TFile * fileCorr = new TFile("CorrCoefficients_"+suffix+"_"+era+".root","recreate");
  fileCorr->cd("");
  corrCoeffX->Write("corrCoeff");
  fileCorr->Close();


}
