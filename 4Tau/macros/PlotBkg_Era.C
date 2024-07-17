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

void GetCorrCoeff(TH1D * hist1D, TH2D * hist2D, int iB, int jB, double * output) {

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
  //  systematics (will be added later)
  //  double center1 = hist1D->GetXaxis()->GetBinWidth(iB);
  //  double center2 = hist1D->GetXaxis()->GetBinWidth(jB);
  //  double sys1 = corr*sysUnc(center1,era);
  //  double sys2 = corr*sysUnc(center2,era);
  //  double ecorr = TMath::Sqrt(stat*stat+sys1*sys1+sys2*sys2);
  double ecorr = stat;
  output[0] = floor(100*corr+0.5)/100;
  output[1] = floor(100*ecorr+0.5)/100;

}

void PlotBkg_Era(
		 TString era = "2018",
		 bool signalRegion = false,
		 bool extendedSideband = false,
		 bool applySystematics = false
		 ) {
  
  TString folderQCD = "/nfs/dust/cms/user/rasp/Run/QCDModel";
  TString folderNonQCD = "/nfs/dust/cms/user/sreelatl/Analyses/H2aa_4tau/Run2/Jul24";
  // new definition of LooseIso : Nsig==1 AND (Niso==3 OR Niso==4)
  TString subfolder = "bin5p2_v3"; // for binning [0,1,2,3,4,5.2,12] 

  std::cout << std::endl;
  if (signalRegion) {
    std::cout << "+++++ Signal region +++++" << std::endl;
  }
  else {
    if (extendedSideband) {
      std::cout << "+++++ Extended Loose-Iso region +++++" << std::endl;
    }
    else {
      std::cout << "+++++ Loose-Iso region +++++" << std::endl;
    }
  }
  std::cout << std::endl;
  
  gROOT->SetBatch(true);
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  SetStyle();

  TString baseName1D = "InvMass";
  TString baseName2D = "InvMass2D";
  TString Suffix = "_ControlXH";
  //  TString histNameNonQCD_1D = "InvMassTrackPlusMuon1D_ControlXH";
  //  TString histNameNonQCD_2D = "InvMassTrackPlusMuon2D_ControlXH";
  //  new definition of Loose-Iso sideband
  TString histNameNonQCD_1D = "InvMassTrackPlusMuon1D_ControlCH";
  TString histNameNonQCD_2D = "InvMassTrackPlusMuon2D_ControlCH";
  if (extendedSideband) {
    Suffix = "_ControlYH";
    //    histNameNonQCD_1D = "InvMassTrackPlusMuon1D_ControlYH";
    //    histNameNonQCD_2D = "InvMassTrackPlusMuon2D_ControlYH";
    //    new definition of Loose-Iso sideband
    histNameNonQCD_1D = "InvMassTrackPlusMuon1D_ControlDH";
    histNameNonQCD_2D = "InvMassTrackPlusMuon2D_ControlDH";
  }
  if (signalRegion) {
    Suffix = "H";
    histNameNonQCD_1D = "InvMassH";
    histNameNonQCD_2D = "InvMass2DH";
  }

  int nBinsNew = 6;
  double bins[7]     = {0, 1, 2, 3, 4, 5.2, 20};
  double binsCorr[7] = {0, 1, 2, 3, 4, 5.2, 12.};

  // lumi in invfb
  std::map<TString,double> eraLumi = {
    {"2016_preVFP", 19.520},
    {"2016_postVFP",16.810},
    {"2017",        41.480},
    {"2018",        59.830}
  };

  // QCD k-factor
  std::map<TString,double> eraKFactor = {
    {"2016_preVFP",  0.55},
    {"2016_postVFP", 0.55},
    {"2017",         0.60},
    {"2018",         0.70}
  };

  std::map<TString,TString> eraLabel = {
    {"2016_preVFP", "2016 (preVFP), 19.5 fb^{-1}"},
    {"2016_postVFP","2016 (postVFP), 16.8 fb^{-1}"},
    {"2016",        "2016 36.3 fb^{-1}"},
    {"2017",        "2017 41.5 fb^{-1}"},
    {"2018",        "2018 59.8 fb^{-1}"},
    {"Run2",        "138 fb^{-1}"}
  };

  std::map<TString,vector<TString>> eraGroup = {
    {"2016_preVFP",{"2016_preVFP"}},
    {"2016_postVFP",{"2016_postVFP"}},
    {"2016",{"2016_preVFP","2016_postVFP"}},
    {"2017",{"2017"}},
    {"2018",{"2018"}},
    {"Run2",{"2016_preVFP","2016_postVFP","2017","2018"}}
  };

  vector<TString> samples = eraGroup[era];
  TH1D * hist1D; // total 1D histogram
  TH2D * hist2D; // total 2D histogram
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

  std::cout << std::endl;
  std::cout << "folder with QCD samples : " << folderQCD << std::endl;
  std::cout << std::endl;
  
  for (auto sample : samples) {
    std::cout << std::endl;
    double lumi = eraLumi[sample];
    double kfactor = eraKFactor[sample];
    TString dir = folderQCD+"/"+sample+"/"+subfolder;
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
      double normSample = kfactor*xsecMuEnrichedQCD[iS]*lumi/nevents;
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

  TH1D * hist1D_qcdonly = (TH1D*)hist1D->Clone("hist1D_qcdonly");
  TH2D * hist2D_qcdonly = (TH2D*)hist2D->Clone("hist2D_qcdonly");

  TString SamplesNonQCD[11] = {
    "DYJetsToLL_M-10to50",
    "DYJetsToLL_M-50",
    "TTTo2L2Nu",
    "TTToSemiLeptonic",
    "ST_t-channel_top",
    "ST_t-channel_antitop",
    "ST_tW_top",
    "ST_tW_antitop",
    "WW_13TeV-pythia8",
    "WZ_13TeV-pythia8",
    "ZZ_13TeV-pythia8",
  };
  
  double xsecNonQCD[11] = {
    21610.0,
    6077.22,
    88.29,
    365.35,
    136.02,
    80.95,
    35.85,
    35.85,
    118.7,
    27.68,
    12.19
  };

  std::cout << std::endl;
  std::cout << "folder with non-QCD samples : " << folderNonQCD << std::endl;
  std::cout << std::endl;
  
  TH1D * hist1D_nonQCD;
  TH2D * hist2D_nonQCD;
  isFirst = true;
  for (auto sample : samples) {
    TString dir = folderNonQCD+"/"+sample;
    std::cout << std::endl;
    double lumi = eraLumi[sample];
    for (int iS=0; iS<11; ++iS) {
      TString name = SamplesNonQCD[iS];
      TFile * fileSample = new TFile(dir+"/"+name+".root");
      if (fileSample->IsZombie() || fileSample==NULL) {
	std::cout << "file " << name << ".root is not found in folder " << std::endl;
	std::cout << dir << std::endl;
	return;
      }
      TH1D * hist1DoldSample = (TH1D*)fileSample->Get(histNameNonQCD_1D);
      if (hist1DoldSample==NULL) {
	std::cout << "histogram named " << histNameNonQCD_1D
		  << " not found in file " << name << ".root" << std::endl;
	return;
      } 
      TH2D * hist2DoldSample = (TH2D*)fileSample->Get(histNameNonQCD_2D);
      if (hist2DoldSample==NULL) {
	std::cout << "histogram named " << histNameNonQCD_2D
		  << " not found in file " << name << ".root" << std::endl;
	return;
      } 
      TH1D * histWeightsSampleH = (TH1D*)fileSample->Get("histWeightsH");
      double nevents = histWeightsSampleH->GetSumOfWeights();
      std::cout << sample << " : " 
		<< SamplesNonQCD[iS] << " : " 
		<< " Ngen = " << nevents
		<< "   Nsel " << hist2DoldSample->GetSumOfWeights() << std::endl;
      double normSample = xsecNonQCD[iS]*lumi/nevents;
      hist1DoldSample->Scale(normSample);
      hist2DoldSample->Scale(normSample);
      TH1D * hist1DSample = (TH1D*)TH1DtoTH1D(hist1DoldSample,nBinsNew,bins,true,sample+name);
      TH2D * hist2DSample = (TH2D*)TH2DtoTH2D(hist2DoldSample,nBinsNew,bins,nBinsNew,bins,sample+name);
      if (isFirst) {
	hist1D_nonQCD = (TH1D*)hist1DSample->Clone("hist1D_nonQCD");
	hist2D_nonQCD = (TH2D*)hist2DSample->Clone("hist2D_nonQCD");
	isFirst = false;
      }
      else {
	hist1D_nonQCD->Add(hist1D_nonQCD,hist1DSample);
	hist2D_nonQCD->Add(hist2D_nonQCD,hist2DSample);
      }
      hist1D->Add(hist1D,hist1DSample);
      hist2D->Add(hist2D,hist2DSample);
    }
  }

  std::cout << std::endl;
  std::cout << "----------------------------" << std::endl;
  std::cout << std::endl;
  std::cout << "Checksum 1D distribution    " << std::endl;
  std::cout << std::endl;
  for (int iB=1; iB<=nBinsNew; ++iB) {
    double x_QCD = hist1D_qcdonly->GetBinContent(iB);
    double x_nonQCD = hist1D_nonQCD->GetBinContent(iB);
    double x_total = hist1D->GetBinContent(iB);
    double f_QCD = x_QCD/x_total;
    double f_nonQCD = x_nonQCD/x_total;
    double f_total = f_QCD + f_nonQCD;
    printf("%1i : %5.3f + %5.3f = %5.3f\n",iB,f_QCD,f_nonQCD,f_total);
  }

  
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "----------------------------" << std::endl;
  std::cout << "Checksum 2D distribution    " << std::endl;
  std::cout << std::endl;
  TH2D * purity2D  = new TH2D("purity2D","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x_QCD = hist2D_qcdonly->GetBinContent(jB,iB);
      double x_nonQCD = hist2D_nonQCD->GetBinContent(jB,iB);
      double x_total = hist2D->GetBinContent(jB,iB);
      double e_nonQCD = hist2D_nonQCD->GetBinError(jB,iB);
      if (iB!=jB) {
	x_QCD += hist2D_qcdonly->GetBinContent(iB,jB); 
	x_nonQCD += hist2D_nonQCD->GetBinContent(iB,jB);
	x_total += hist2D->GetBinContent(iB,jB);
	double error = hist2D_nonQCD->GetBinError(iB,jB);
	double e2 = e_nonQCD*e_nonQCD + error*error;
	e_nonQCD = TMath::Sqrt(e2);
      }
      double f_QCD = x_QCD/x_total;
      double f_nonQCD = x_nonQCD/x_total;
      e_nonQCD = e_nonQCD/x_total;
      double f_total = f_QCD + f_nonQCD;
      purity2D->SetBinContent(iB,jB,f_nonQCD);
      purity2D->SetBinError(iB,jB,e_nonQCD);
      printf("[%1i,%1i] : %5.3f + %5.3f = %5.3f\n",iB,jB,f_QCD,f_nonQCD,f_total);
    }
  }
  
  TH2D * corrCoeff = new TH2D("corrCoeffH","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  TH2D * corrCoeffX = new TH2D("corrCoeffX","",nBinsNew,bins,nBinsNew,bins);
  TH2D * corrCoeffX_qcdonly = new TH2D("corrCoeffX_qcdonly","",nBinsNew,bins,nBinsNew,bins);

  hist1D->Scale(1/hist1D->GetSumOfWeights());
  hist2D->Scale(1/hist2D->GetSumOfWeights());

  hist1D_qcdonly->Scale(1/hist1D_qcdonly->GetSumOfWeights());
  hist2D_qcdonly->Scale(1/hist2D_qcdonly->GetSumOfWeights());
  
  std::cout << std::endl;
  std::cout << "C(i,j)      total bkg    :      QCD only" << std::endl;
  std::cout << std::endl;
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double output[2];
      GetCorrCoeff(hist1D,hist2D,iB,jB,output);
      double corr  = output[0];
      double ecorr = output[1];
      GetCorrCoeff(hist1D_qcdonly,hist2D_qcdonly,iB,jB,output);
      double corr_qcdonly  = output[0];
      double ecorr_qcdonly = output[1];
      printf("[%1i,%1i] = %5.2f +/- %5.2f  :  %5.2f +/- %5.2f\n",
	     iB,jB,corr,ecorr,corr_qcdonly,ecorr_qcdonly);
      corrCoeffX->SetBinContent(iB,jB,corr);
      corrCoeffX->SetBinError(iB,jB,ecorr);
      corrCoeffX->SetBinContent(iB,jB,corr);
      corrCoeffX->SetBinError(iB,jB,ecorr);
      corrCoeffX_qcdonly->SetBinContent(iB,jB,corr_qcdonly);
      corrCoeffX_qcdonly->SetBinError(iB,jB,ecorr_qcdonly);
    }
  }

  std::cout << std::endl;
  std::cout << "-------------------------------------" << std::endl;
  std::cout << std::endl;
  std::cout << "Fraction of non-QCD bkg " << std::endl;
  std::cout << std::endl;
  
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double purity = purity2D->GetBinContent(iB,jB);
      double epurity = purity2D->GetBinError(iB,jB);
      printf("[%1i,%1i] = %5.3f +/- %5.3f\n",iB,jB,purity,epurity);
    }
  }  
  std::cout << std::endl;

  //
  // Correlation coefficients
  //
  corrCoeff->GetXaxis()->SetNdivisions(207);
  corrCoeff->GetYaxis()->SetNdivisions(207);
  corrCoeff->GetYaxis()->SetTitleOffset(1.0);
  corrCoeff->SetMarkerSize(1.2);
  corrCoeff->GetXaxis()->SetTitle("m_{1} [GeV]");
  corrCoeff->GetYaxis()->SetTitle("m_{2} [GeV]");

  TLatex * latexBin = new TLatex();
  latexBin->SetTextSize(0.03);
  latexBin->SetTextFont(32);

  //
  // total background
  // 
  TCanvas * canv0 = MakeCanvas("canv0","",1000,700);
  corrCoeff->Draw("texte");

  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = corrCoeffX->GetBinContent(iB,jB);
      double error = corrCoeffX->GetBinError(iB,jB);
      char label[30];
      sprintf(label,"%4.2f^{+%4.2f}_{-%4.2f}",x,error,error);
      TString Label(label);
      double binX = corrCoeff->GetXaxis()->GetBinCenter(iB);
      double binWidthX = corrCoeff->GetXaxis()->GetBinWidth(iB); 
      double binY = corrCoeff->GetYaxis()->GetBinCenter(jB);
      double binWidthY = corrCoeff->GetYaxis()->GetBinWidth(jB);       
      double textX = binX - 0.46;
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
  CMS_lumi(canv0,4,33);  
  canv0->Update();

  TString suffix("control");
  if (extendedSideband)
    suffix += "_ext";
  if (signalRegion)
    suffix = "signal";
  suffix += "_mc";
  
  canv0->Print("CorrCoefficients_"+suffix+"_"+era+".png");

  
  //
  // total background
  // 
  TCanvas * canv1 = MakeCanvas("canv1","",1000,700);
  corrCoeff->Draw("texte");

  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = corrCoeffX_qcdonly->GetBinContent(iB,jB);
      double error = corrCoeffX_qcdonly->GetBinError(iB,jB);
      char label[30];
      sprintf(label,"%4.2f^{+%4.2f}_{-%4.2f}",x,error,error);
      TString Label(label);
      double binX = corrCoeff->GetXaxis()->GetBinCenter(iB);
      double binWidthX = corrCoeff->GetXaxis()->GetBinWidth(iB); 
      double binY = corrCoeff->GetYaxis()->GetBinCenter(jB);
      double binWidthY = corrCoeff->GetYaxis()->GetBinWidth(jB);       
      double textX = binX - 0.46;
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
  CMS_lumi(canv1,4,33);  
  canv1->Update();
  canv1->Print("CorrCoefficients_"+suffix+"_qcdonly_"+era+".png");

  
  //
  // purity
  //
  purity2D->GetXaxis()->SetNdivisions(207);
  purity2D->GetYaxis()->SetNdivisions(207);
  purity2D->GetYaxis()->SetTitleOffset(1.0);
  purity2D->SetMarkerSize(1.2);
  purity2D->GetXaxis()->SetTitle("m_{1} [GeV]");
  purity2D->GetYaxis()->SetTitle("m_{2} [GeV]");
    
  TCanvas * canv2 = MakeCanvas("canv2","",1000,700);
  corrCoeff->Draw("texte");
  TLatex * latexBin2 = new TLatex();
  latexBin2->SetTextSize(0.037);
  latexBin2->SetTextFont(32);
  
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = purity2D->GetBinContent(iB,jB);
      double ex = purity2D->GetBinError(iB,jB);
      char label[10];
      if (x<0.005)
	sprintf(label,"<0.01");
      else
	sprintf(label,"%5.2f",x);
      TString Label(label);
      double binX = purity2D->GetXaxis()->GetBinCenter(iB);
      double binWidthX = purity2D->GetXaxis()->GetBinWidth(iB); 
      double binY = purity2D->GetYaxis()->GetBinCenter(jB);
      double binWidthY = purity2D->GetYaxis()->GetBinWidth(jB);       
      double textX = binX - 0.45;
      double textY = binY;
      latexBin2->SetTextAlign(12);
      latexBin2->DrawLatex(textX,textY,Label);
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
  CMS_lumi(canv2,4,33);  
  canv2->Update();
  canv2->Print("purityQCD_"+suffix+"_"+era+".png");

  TString outputName = "CorrCoefficients_"+suffix+"_"+era+".root";
  TFile * fileCorr = new TFile(outputName,"recreate");
  fileCorr->cd("");
  corrCoeffX->Write("corrCoeff");
  corrCoeffX_qcdonly->Write("corrCoeff_qcdonly");
  fileCorr->Close();


}
