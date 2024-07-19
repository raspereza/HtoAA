#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"

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
  double ecorr = stat;
  output[0] = floor(100*corr+0.5)/100;
  output[1] = floor(100*ecorr+0.5)/100;
  output[2] = denominator;
}

void CorrCoefficientsContamination(
				   TString era = "2017", // 2016, 2016_preVFP, 2016_postVFP, 2017, 2018
				   int mass = 10, // ma = 4-15 (GeV)
				   bool extendedSideBand = false // extended sideband
				   ) {

  SetStyle();
  gROOT->SetBatch(true);

  TString folder = "/nfs/dust/cms/user/sreelatl/Analyses/H2aa_4tau/Run2/Jul24";
  TString workdir = std::getenv("PWD");
  
  int nBinsNew = 6;
  double bins[7]     = {0, 1, 2, 3, 4, 5.2, 20};
  double binsCorr[7] = {0, 1, 2, 3, 4, 5.2, 12};
  
  TString histName1D("InvMassTrackPlusMuon1D_ControlCH");
  TString histName2D("InvMassTrackPlusMuon2D_ControlCH");
  if (extendedSideBand) {
    histName1D = "InvMassTrackPlusMuon1D_ControlDH";
    histName2D = "InvMassTrackPlusMuon2D_ControlDH";
  }

  char massChar[3];
  if (mass<4||mass>15) {
    std::cout << std::endl;
    std::cout << "invalid specified mass " << mass << std::endl;
    std::cout << "mass should be in the range [4,15]" << std::endl;
    std::cout << std::endl;
    return;
  }
  else {
    std::cout << std::endl;
    std::cout << "Estimating contamination for signal with ma = " << mass << std::endl;
  }
  
  if (mass>9)
    sprintf(massChar,"%2i",mass);
  else {
    sprintf(massChar,"%1i",mass);
  }

  TString Mass(massChar);
  
  std::map<TString, TString> LUMI_label = {
    {"2016"        ,"2016, 36.3 fb^{-1}"},
    {"2016_preVFP" ,"2016, preVFP, 19.5 fb^{-1}"},
    {"2016_postVFP","2016, postVFP, 16.8 fb^{-1}"},
    {"2017"        ,"2017, 41.5 fb^{-1}"},
    {"2018"        ,"2018, 59.8 fb^{-1}"},
    {"Run2"        ,"138 fb^{-1}"}
  };

  // lumi in invpb
  std::map<TString,double> eraLumi = {
    {"2016_preVFP", 19520},
    {"2016_postVFP",16810},
    {"2017",        41480},
    {"2018",        59830}
  };

  std::map<TString,vector<TString>> eraGroup = {
    {"2016_preVFP",{"2016_preVFP"}},
    {"2016_postVFP",{"2016_postVFP"}},
    {"2016",{"2016_preVFP","2016_postVFP"}},
    {"2017",{"2017"}},
    {"2018",{"2018"}},
    {"Run2",{"2016_preVFP","2016_postVFP","2017","2018"}}
  };

  std::vector<TString> eras = eraGroup[era];

  std::cout << std::endl;
  std::cout << "Eras included : ";
  for (auto currentEra : eras) {
    std::cout << currentEra << " ";
  }
  std::cout << std::endl;
  
  double xsecGGH = 48.61;
  double xsecVBF = 3.766;
  double xsecVH  = 1.358 + 0.880;
  double massD = double(mass);
  
  double massTau = 1.777;
  double massMu  = 0.106;
  double massRatio = (massMu * massMu) / (massTau * massTau);
  double aF = 2 * massTau / massD;
  double SF = 2 * massRatio / TMath::Sqrt(1 - aF * aF);
  double xsecMMTT = (xsecGGH + xsecVBF + xsecVH) * SF;
  double BR = 0.16;
    
  map<TString, double> signalXSec;
  map<TString, TString> signalProcess = {
    {"GGH", "SUSYGluGluToHToAA_AToTauTau_M-125_M-"},
    {"VBF", "SUSYVBFToHToAA_AToTauTau_M-125_M-"},
    {"VH", "SUSYVH_HToAA_AToTauTau_M-125_M-"},
    {"MMTT", "SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-"}
  };

  map<TString, TString> signalLabel = {
    {"GGH", "gg->H 4tau   "},
    {"VBF", "qqH   4tau   "},
    {"VH",  "VH    4tau   "},
    {"MMTT","incl  2mu2tau"}
  };

  vector<TString> signals = {"GGH","VBF","VH","MMTT"};
  
  map<TString, TString> eraLabel = {
    {"2016_preVFP" ,"2016APV"},
    {"2016_postVFP","2016   "},
    {"2017"        ,"2017   "},
    {"2018"        ,"2018   "}
  };
  
  signalXSec["GGH"] = xsecGGH*BR;
  signalXSec["VBF"] = xsecVBF*BR;
  signalXSec["VH"] = xsecVH*BR;
  signalXSec["MMTT"] = xsecMMTT*BR;
  
  lumi_13TeV = LUMI_label[era];

  //
  // frame histogram
  // 
  TH2D * corrCoeff = new TH2D("corrCoeff","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  corrCoeff->GetXaxis()->SetNdivisions(207);
  corrCoeff->GetYaxis()->SetNdivisions(207);
  corrCoeff->GetYaxis()->SetTitleOffset(1.0);
  corrCoeff->SetMarkerSize(1.2);
  corrCoeff->GetXaxis()->SetTitle("m_{1} [GeV]");
  corrCoeff->GetYaxis()->SetTitle("m_{2} [GeV]");

  TH2D * corrCoeffX = new TH2D("corrCoeffX","",nBinsNew,bins,nBinsNew,bins);
  TH2D * corrCoeffH = new TH2D("corrCoeffH","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  TH2D * corrCoeffL = new TH2D("corrCoeffL","",nBinsNew,binsCorr,nBinsNew,binsCorr);

  TH2D * corrCoeffDataX = new TH2D("corrCoeffDataX","",nBinsNew,bins,nBinsNew,bins);
  TH2D * corrCoeffDataH = new TH2D("corrCoeffDataH","",nBinsNew,binsCorr,nBinsNew,binsCorr);
  TH2D * corrCoeffDataL = new TH2D("corrCoeffDataL","",nBinsNew,binsCorr,nBinsNew,binsCorr);

  TH1D * hist1D = new TH1D("hist1D","",nBinsNew,bins);
  TH2D * hist2D = new TH2D("hist2D","",nBinsNew,bins,nBinsNew,bins);

  std::cout << std::endl;
  for (auto currentEra : eras) {
    TString dir = folder + "/" + currentEra;
    TString filename = dir+"/DoubleMuon_Run"+currentEra+".root";
    TFile * file = new TFile(filename);
    if (file->IsZombie() || file==NULL) {
      std::cout << "File " << filename << " cannot be opened" << std::endl;
      return;
    }
    TH1D * hist1Dold = (TH1D*)file->Get(histName1D);
    TH2D * hist2Dold = (TH2D*)file->Get(histName2D);
    if (hist1Dold==NULL) {
      std::cout << "Histogram " << histName1D << " is not found in file " << filename << std::endl;
      return;
    }
    if (hist2Dold==NULL) {
      std::cout << "Histogram " << histName2D << " is not found in file " << filename << std::endl;
      return;
    }
    char yield[20];
    sprintf(yield,"%6i",int(hist2Dold->GetSumOfWeights()));
    std::cout << eraLabel[currentEra] << " : data          = " << yield << std::endl;
    TH1D * hist1Dsample = (TH1D*)TH1DtoTH1D(hist1Dold,nBinsNew,bins,true,"_new");
    TH2D * hist2Dsample = (TH2D*)TH2DtoTH2D(hist2Dold,nBinsNew,bins,nBinsNew,bins,"_2DNew");
    hist1D->Add(hist1D,hist1Dsample,1.,1.);
    hist2D->Add(hist2D,hist2Dsample,1.,1.);
  }

  TH1D * hist1D_data = (TH1D*)hist1D->Clone("hist1D_data");
  TH2D * hist2D_data = (TH2D*)hist2D->Clone("hist2D_data");

  TH1D * hist1D_signal = new TH1D("hist1D_signal","",nBinsNew,bins);
  TH2D * hist2D_signal = new TH2D("hist2D_signal","",nBinsNew,bins,nBinsNew,bins);

  for (auto currentEra : eras) {
    TString dir = folder + "/" + currentEra; 
    for (auto signal : signals) {
      TString signalName = signal;
      TString signalSample = signalProcess[signal];
      TString filename = dir+"/"+signalSample+Mass+".root";
      TFile * file = new TFile(filename);
      if (file->IsZombie() || file==NULL) {
	std::cout << "File " << filename << " cannot be opened" << std::endl;
	return;
      }
      TH1D * hist1Dold = (TH1D*)file->Get(histName1D);
      TH2D * hist2Dold = (TH2D*)file->Get(histName2D);
      if (hist1Dold==NULL) {
	std::cout << "Histogram " << histName1D << " is not found in file " << filename << std::endl;
	return;
      }
      if (hist2Dold==NULL) {
	std::cout << "Histogram " << histName2D << " is not found in file " << filename << std::endl;
	return;
      }
      char yield[20];
      sprintf(yield,"%6.2f",hist2Dold->GetSumOfWeights());
      std::cout << eraLabel[currentEra] << " : " << signalLabel[signalName]
		<< " = " << yield << std::endl;
      TH1D * histWeightsSampleH = (TH1D*)file->Get("histWeightsH");
      double nevents = histWeightsSampleH->GetSumOfWeights();
      double xsec = signalXSec[signalName];
      double lumi = eraLumi[currentEra];
      double norm = xsec*lumi/nevents;
      TH1D * hist1Dsample = (TH1D*)TH1DtoTH1D(hist1Dold,nBinsNew,bins,true,"_new");
      TH2D * hist2Dsample = (TH2D*)TH2DtoTH2D(hist2Dold,nBinsNew,bins,nBinsNew,bins,"_sigNew");
      hist1Dsample->Scale(norm);
      hist2Dsample->Scale(norm);
      hist1D->Add(hist1D,hist1Dsample,1.,-1.);
      hist2D->Add(hist2D,hist2Dsample,1.,-1.);
      hist1D_signal->Add(hist1D_signal,hist1Dsample,1.,1.);
      hist2D_signal->Add(hist2D_signal,hist2Dsample,1.,1.);
    }
  }

  std::cout << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << "Checksum 1D distribution" << std::endl;
  std::cout << std::endl;
  for (int iB=1; iB<=nBinsNew; ++iB) {
    double x_data = hist1D_data->GetBinContent(iB);
    double x_signal = hist1D_signal->GetBinContent(iB);
    double x_total = hist1D->GetBinContent(iB);
    double f_data   = x_data/x_total;
    double f_signal = x_signal/x_total;
    double f_total = f_data - f_signal;
    double f = 100.*x_signal/x_data;
    printf("%1i : %5.3f - %5.3f = %5.3f -> Nsig/Ndata = %4.1f%%\n",iB,f_data,f_signal,f_total,f);
  }
  
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << "Checksum 2D distribution" << std::endl;
  std::cout << std::endl;
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x_data = hist2D_data->GetBinContent(jB,iB);
      double x_signal = hist2D_signal->GetBinContent(jB,iB);
      double x_total = hist2D->GetBinContent(jB,iB);
      if (iB!=jB) {
	x_data += hist2D_data->GetBinContent(iB,jB); 
	x_signal += hist2D_signal->GetBinContent(iB,jB);
	x_total += hist2D->GetBinContent(iB,jB);
      }
      double f_data = x_data/x_total;
      double f_signal = x_signal/x_total;
      double f_total = f_data - f_signal;
      double f = 100.*x_signal/x_data;
      printf("[%1i,%1i] : %5.3f - %5.3f = %5.3f -> Nsig/Ndata = %4.1f%%\n",iB,jB,f_data,f_signal,f_total,f);
    }
  }

  std::cout << std::endl;
  
  double norm2D = hist2D->GetSumOfWeights(); 
  double norm2D_data = hist2D_data->GetSumOfWeights();

  TLatex * latexBin = new TLatex();
  latexBin->SetTextSize(0.038);
  latexBin->SetTextFont(32);
  
  //
  // Plot statistics in data
  //
  TCanvas * canvStatData = MakeCanvas("canvStatData","",1000,700);
  corrCoeff->Draw("text");
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      char label[10];
      int x = int(hist2D_data->GetBinContent(iB,jB));
      sprintf(label,"%5i",x);
      TString Label(label);
      double binX = corrCoeffH->GetXaxis()->GetBinCenter(iB);
      double binWidthX = corrCoeffH->GetXaxis()->GetBinWidth(iB); 
      double binY = corrCoeffH->GetYaxis()->GetBinCenter(jB);
      double binWidthY = corrCoeffH->GetYaxis()->GetBinWidth(jB);       
      double textX = binX - 0.43;
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
  CMS_lumi(canvStatData,4,33);

  canvStatData->Update();
  canvStatData->RedrawAxis();
  if (extendedSideBand)
    canvStatData->Print(workdir+"/figures/stat_data_"+era+"_ext.png");
  else
    canvStatData->Print(workdir+"/figures/stat_data_"+era+".png");

  std::cout << std::endl;
  //
  // Plot statistics in data with subtracted signal
  // 
  TCanvas * canvSubtrSig = MakeCanvas("canvSubtrSig","",1000,700);
  corrCoeff->Draw("text");
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      char label[10];
      double x = hist2D->GetBinContent(iB,jB);
      if (x>=1000.)
	sprintf(label,"%5.0f",x);
      else 
	sprintf(label,"%5.1f",x);
      TString Label(label);
      double binX = corrCoeffH->GetXaxis()->GetBinCenter(iB);
      double binWidthX = corrCoeffH->GetXaxis()->GetBinWidth(iB); 
      double binY = corrCoeffH->GetYaxis()->GetBinCenter(jB);
      double binWidthY = corrCoeffH->GetYaxis()->GetBinWidth(jB);       
      double textX = binX - 0.43;
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
  CMS_lumi(canvSubtrSig,4,33);

  canvSubtrSig->Update();
  canvSubtrSig->RedrawAxis();
  if (extendedSideBand)
    canvSubtrSig->Print(workdir+"/figures/stat_data_ma"+Mass+"_"+era+"_ext.png");
  else
    canvSubtrSig->Print(workdir+"/figures/stat_data_ma"+Mass+"_"+era+".png");

  //
  // Computing statistics in 2D distribution
  //  
  std::cout << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "Statistics in 2D distribution" << std::endl;
  std::cout << std::endl;
  std::cout << " bin  |  data | data-signal  " << std::endl;
  std::cout << "------|-------|--------------" << std::endl;
  //           "[1,1] |  5995 |    5995.0 
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = hist2D->GetBinContent(iB,jB);
      double x_data = hist2D_data->GetBinContent(iB,jB);
      printf("[%1i,%1i] | %5i |   %7.1f\n",iB,jB,int(x_data),x);
      
      double elow = TMath::Max(double(0.0),double(-0.5 + TMath::Sqrt(x+0.25)));
      double ehigh = TMath::Max(double(0.0),double(0.5 + TMath::Sqrt(x+0.25)));
      x = x/norm2D;
      elow = elow/norm2D;
      ehigh = ehigh/norm2D;
      
      corrCoeffL->SetBinContent(iB,jB,x);
      corrCoeffL->SetBinError(iB,jB,elow);
      corrCoeffH->SetBinContent(iB,jB,x);
      corrCoeffH->SetBinError(iB,jB,ehigh);

      double elow_data = TMath::Max(double(0.0),double(-0.5 + TMath::Sqrt(x_data+0.25)));
      double ehigh_data = TMath::Max(double(0.0),double(0.5 + TMath::Sqrt(x_data+0.25)));
      x_data = x_data/norm2D_data;
      elow_data = elow_data/norm2D_data;
      ehigh_data = ehigh_data/norm2D_data;

      corrCoeffDataL->SetBinContent(iB,jB,x_data);
      corrCoeffDataL->SetBinError(iB,jB,elow_data);
      corrCoeffDataH->SetBinContent(iB,jB,x_data);
      corrCoeffDataH->SetBinError(iB,jB,ehigh_data);
    }
  }  
  
  std::cout << std::endl;
  hist1D->Scale(1.0/hist1D->GetSumOfWeights());
  hist2D->Scale(1.0/hist2D->GetSumOfWeights());  
  hist1D_data->Scale(1.0/hist1D_data->GetSumOfWeights());
  hist2D_data->Scale(1.0/hist2D_data->GetSumOfWeights());
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "     Correlation coefficients C(i,j)" << std::endl;
  std::cout << std::endl;
  std::cout << " bin  |      data       |    data-signal   " << std::endl;
  std::cout << "------|-----------------|------------------" << std::endl;
  //           "[1,1] |  0.99 +/-  0.02 |  0.99 +/-  0.02
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {

      double output[3];
      
      GetCorrCoeff(hist1D,hist2D,iB,jB,output);
      double corr = output[0];
      double ecorr = output[1];
      double denominator = output[2];
      double ehigh = corrCoeffH->GetBinError(iB,jB)/denominator;
      double elow = corrCoeffL->GetBinError(iB,jB)/denominator;

      GetCorrCoeff(hist1D_data,hist2D_data,iB,jB,output);
      double corr_data = output[0];
      double ecorr_data = output[1];
      double denominator_data = output[2];
      double ehigh_data = corrCoeffDataH->GetBinError(iB,jB)/denominator_data;
      double elow_data = corrCoeffDataL->GetBinError(iB,jB)/denominator_data;

      printf("[%1i,%1i] | %5.2f +/- %5.2f | %5.2f +/- %5.2f\n",iB,jB,corr_data,ecorr_data,corr,ecorr);

      corrCoeffH->SetBinContent(iB,jB,corr);
      corrCoeffH->SetBinError(iB,jB,ehigh);
      corrCoeffL->SetBinContent(iB,jB,corr);
      corrCoeffL->SetBinError(iB,jB,elow);      
      corrCoeffX->SetBinContent(iB,jB,corr);
      corrCoeffX->SetBinError(iB,jB,ecorr);

      corrCoeffDataH->SetBinContent(iB,jB,corr_data);
      corrCoeffDataH->SetBinError(iB,jB,ehigh_data);
      corrCoeffDataL->SetBinContent(iB,jB,corr_data);
      corrCoeffDataL->SetBinError(iB,jB,elow_data);      
      corrCoeffDataX->SetBinContent(iB,jB,corr_data);
      corrCoeffDataX->SetBinError(iB,jB,ecorr_data);

    }
  }

  //
  // Plotting correlation coefficients in data with signal subtracted
  //
  std::cout << std::endl;
  latexBin->SetTextSize(0.03);
  latexBin->SetTextFont(32);
  TCanvas * canv1 = MakeCanvas("canv1","",1000,700);
  corrCoeff->Draw("texte");
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = corrCoeffH->GetBinContent(iB,jB);
      double elow = corrCoeffL->GetBinError(iB,jB);
      double ehigh = corrCoeffH->GetBinError(iB,jB);
      double error = corrCoeffX->GetBinError(iB,jB);
      char label[30];
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
  CMS_lumi(canv1,4,33);

  canv1->RedrawAxis();
  canv1->Update();
  if (extendedSideBand)
    canv1->Print(workdir+"/figures/CorrCoefficients_data_ma"+Mass+"_"+era+"_ext.png");
  else
    canv1->Print(workdir+"/figures/CorrCoefficients_data_ma"+Mass+"_"+era+".png");

  //
  // Plotting correlation coefficients in data without signal subtraction
  //
  std::cout << std::endl;
  TCanvas * canv2 = MakeCanvas("canv2","",1000,700);
  corrCoeff->Draw("texte");
  for (int iB=1; iB<=nBinsNew; ++iB) {
    for (int jB=iB; jB<=nBinsNew; ++jB) {
      double x = corrCoeffDataH->GetBinContent(iB,jB);
      double elow = corrCoeffDataL->GetBinError(iB,jB);
      double ehigh = corrCoeffDataH->GetBinError(iB,jB);
      double error = corrCoeffDataX->GetBinError(iB,jB);
      char label[30];
      sprintf(label,"%4.2f^{+%4.2f}_{ -%4.2f}",x,ehigh,elow);
      TString Label(label);
      double binX = corrCoeffDataH->GetXaxis()->GetBinCenter(iB);
      double binWidthX = corrCoeffDataH->GetXaxis()->GetBinWidth(iB); 
      double binY = corrCoeffDataH->GetYaxis()->GetBinCenter(jB);
      double binWidthY = corrCoeffDataH->GetYaxis()->GetBinWidth(jB);       
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
  CMS_lumi(canv2,4,33);

  canv2->RedrawAxis();
  canv2->Update();
  if (extendedSideBand)
    canv2->Print(workdir+"/figures/CorrCoefficients_data_"+era+"_ext.png");
  else
    canv2->Print(workdir+"/figures/CorrCoefficients_data_"+era+".png");

  std::cout << std::endl;
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  TFile * fileCorr = new TFile(workdir+"/CorrCoefficients_data_ma"+Mass+"_"+era+".root","recreate");
  fileCorr->cd("");
  corrCoeffX->Write("corrCoeff");
  corrCoeffDataX->Write("corrCoeffDataOnly");
  fileCorr->Close();

  TFile * fileCorrData = new TFile(workdir+"/CorrCoefficients_data_"+era+".root","recreate");
  fileCorrData->cd("");
  corrCoeffDataX->Write("corrCoeff");
  fileCorrData->Close();

}
