#include "HttStylesNew.cc"
#include "HtoH.h"
#include "CMS_lumi.C"
//  TH1D * hist1Dold = (TH1D*)file->Get("InvMassTrackPlusMuon1D_"+Suffix);
//  TH2D * hist2Dold = (TH2D*)file->Get("InvMassTrackPlusMuon2D_"+Suffix);
//  TH1D * hist1Dold = (TH1D*)file->Get("ModelInvMassH");
//  TH2D * hist2Dold = (TH2D*)file->Get("ModelInvMass2DH");

// Systematic uncertainty (parton shower scale) 
// parameterized as a function of (mu,trk) mass
double sysUnc(double mass,
	      TString era) {
  double alpha = 0.017;
  double beta = 10.2;
  if (era=="2017") {
    alpha = 0.018;
    beta = 11.5;
  }
  else if (era=="2018") {
    alpha = 0.013;
    beta = 13.0;
  }
  return alpha*TMath::Exp(mass/beta);    

} 

void CorrCoefficientsMC_Era(
			    TString era = "2018",
			    TString subfolder = "bin5p2",
			    bool signalRegion = false,
			    bool extendedSideband = false,
			    bool applySystematics = false
			    ) {
  
  gROOT->SetBatch(true);
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  SetStyle();

  TString baseName1D = "InvMass";
  TString baseName2D = "InvMass2D";
  TString Suffix = "_ControlXH";
  if (extendedSideband) {
    Suffix = "_ControlYH";
  }
  if (signalRegion) {
    Suffix = "H";        
  }

  int nBinsNew = 6;
  double bins[7]     = {0, 1, 2, 3, 4, 5.2, 20};
  double binsCorr[7] = {0, 1, 2, 3, 4, 5.2, 12.};

  // lumi in invfb
  std::map<TString,double> eraLumi = {
    {"2016_preVFP", 19.520},
    {"2016_postVFP",16.810},
    {"2017",        14.480},
    {"2018",        59.830}
  };

  std::map<TString,TString> eraLabel = {
    {"2016_preVFP", "2016 (preVFP), 19.5 fb^{-1}"},
    {"2016_postVFP","2016 (postVFP), 16.8 fb^{-1}"},
    {"2016",        "2016 36.3 fb^{-1}"},
    {"2017",        "2017 41.5 fb^{-1}"},
    {"2018",        "2018 59.8 fb^{-1}"}
  };

  std::map<TString,vector<TString>> eraGroup = {
    {"2016_preVFP",{"2016_preVFP"}},
    {"2016_postVFP",{"2016_postVFP"}},
    {"2016",{"2016_preVFP","2016_postVFP"}},
    {"2017",{"2017"}},
    {"2018",{"2018"}},
  };

  vector<TString> samples = eraGroup[era];
  TH1D * hist1D;
  TH2D * hist2D;
  bool isFirst = true;
  
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

  for (auto sample : samples) {
    TString dir = "/nfs/dust/cms/user/rasp/Run/QCDModel/"+sample+"/"+subfolder;
    std::cout << std::endl;
    double lumi = eraLumi[sample];
    for (int iS=1; iS<11; ++iS) {
      TString name = SamplesMuEnrichedQCD[iS];
      TFile * fileSample = new TFile(dir+"/"+name+".root");
      if (fileSample->IsZombie() || fileSample==NULL) {
	std::cout << "file " << name << ".root is not found in folder " << std::endl;
	std::cout << dir << std::endl;
	return;
      }
      TH1D * hist1DoldSample = (TH1D*)fileSample->Get(baseName1D+Suffix);
      if (hist1DoldSample==NULL) {
	std::cout << "histogram named " << baseName1D << Suffix 
		  << " not found in file " << name << ".root" << std::endl;
	return;
      } 
      TH2D * hist2DoldSample = (TH2D*)fileSample->Get(baseName2D+Suffix);
      if (hist2DoldSample==NULL) {
	std::cout << "histogram named " << baseName2D << Suffix 
		  << " not found in file " << name << ".root" << std::endl;
	return;
      } 
      TH1D * histWeightsSampleH = (TH1D*)fileSample->Get("histWeightsH");
      double nevents = histWeightsSampleH->GetSumOfWeights();
      std::cout << sample << " : " 
		<< SamplesMuEnrichedQCD[iS] << " : " 
		<< int(nevents) << std::endl;
      double normSample = xsecMuEnrichedQCD[iS]*lumi/nevents;
      hist1DoldSample->Scale(normSample);
      hist2DoldSample->Scale(normSample);
      TH1D * hist1DSample = (TH1D*)TH1DtoTH1D(hist1DoldSample,nBinsNew,bins,true,sample+name);
      TH2D * hist2DSample = (TH2D*)TH2DtoTH2D(hist2DoldSample,nBinsNew,bins,nBinsNew,bins,sample+name);
      if (isFirst) {
	hist1D = (TH1D*)hist1DSample->Clone("hist1D");
	hist2D = (TH2D*)hist2DSample->Clone("hist2D");
	isFirst = false;
      }
      else {
	hist1D->Add(hist1D,hist1DSample);
	hist2D->Add(hist2D,hist2DSample);
      }
    }
    std::cout << std::endl;
  }

  TH2D * corrCoeff = new TH2D("corrCoeff","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  TH2D * corrCoeffX = new TH2D("corrCoeffX","",nBinsNew,bins,nBinsNew,bins);
  TH2D * corrCoeffH = new TH2D("corrCoeffH","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  TH2D * corrCoeffL = new TH2D("corrCoeffL","",nBinsNew,binsCorr,nBinsNew,binsCorr);

  std::cout << std::endl;
  hist1D->Scale(1/hist1D->GetSumOfWeights());
  hist2D->Scale(1/hist2D->GetSumOfWeights());

  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = hist1D->GetBinContent(iB);
      double y = hist1D->GetBinContent(jB);
      double denominator = x*y;
      double ex = hist1D->GetBinError(iB);
      double ey = hist1D->GetBinError(jB);
      double rx = ex/x;
      double ry = ey/y;
      double rdenominator = 2*rx;
      if (iB!=jB) {
	denominator *= 2.0;
	rdenominator = 2.0*TMath::Sqrt(rx*rx+ry*ry);
      }
      double numerator  = hist2D->GetBinContent(jB,iB);
      double enumerator = hist2D->GetBinError(jB,iB);
      if (iB!=jB) {
	numerator += hist2D->GetBinContent(iB,jB);
	enumerator *= 2.0;
      }
      double renumenator = enumerator/numerator;
      double rcorr = TMath::Sqrt(renumenator*renumenator+rdenominator*rdenominator);
      double corr  = numerator / denominator ;
      double stat = rcorr*corr;
      double center1 = hist1D->GetXaxis()->GetBinWidth(iB);
      double center2 = hist1D->GetXaxis()->GetBinWidth(jB);
      double sys1 = corr*sysUnc(center1,era);
      double sys2 = corr*sysUnc(center2,era);
      //      double ecorr = TMath::Sqrt(stat*stat+sys1*sys1+sys2*sys2);
      double ecorr = stat;
      corr = floor(100*corr+0.5)/100;
      ecorr = floor(100*ecorr+0.5)/100;
      printf("[%1i,%1i] = %5.2f +/- %5.2f \n",iB,jB,corr,ecorr);
      //      corrCoeff->SetBinContent(iB,jB,corr);
      //      corrCoeff->SetBinError(iB,jB,ecorr);
      corrCoeffX->SetBinContent(iB,jB,corr);
      corrCoeffX->SetBinError(iB,jB,ecorr);
      double ehigh = corrCoeffH->GetBinError(iB,jB)/denominator;
      double elow = corrCoeffL->GetBinError(iB,jB)/denominator;
      corrCoeffH->SetBinContent(iB,jB,corr);
      corrCoeffH->SetBinError(iB,jB,ecorr);
      corrCoeffL->SetBinContent(iB,jB,corr);
      corrCoeffL->SetBinError(iB,jB,ecorr);      
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

  lumi_13TeV = eraLabel[era];
  writeExtraText = true;
  extraText = "Simulation";  
  CMS_lumi(canv,4,33);  
  canv->Update();

  TString suffix("control");
  if (extendedSideband)
    suffix += "_ext";
  if (signalRegion)
    suffix = "signal";
  suffix += "_mc";

  canv->Print("CorrCoefficients_"+suffix+"_"+"_"+era+".png");

  TFile * fileCorr = new TFile("CorrCoefficients_"+suffix+"_"+"_"+era+".root","recreate");
  fileCorr->cd("");
  corrCoeffX->Write("corrCoeff");
  fileCorr->Close();


}
