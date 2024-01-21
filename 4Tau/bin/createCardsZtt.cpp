#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <string> 
#include <map>

#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "HtoAA/Utilities/interface/Config.h"
#include "HtoAA/Utilities/src/Config.cc"
#include <iomanip>

using namespace std;

std::map<TString,std::vector<TString>> groups = {
  {"EWK",{"WW_13TeV-pythia8",
	  "WZ_13TeV-pythia8",
	  "ZZ_13TeV-pythia8",
	  "WJetsToLNu"}},
  
  {"TT",{"TTTo2L2Nu",
	 "TTToSemiLeptonic",
	 "TTToHadronic",
	 "ST_t-channel_top",
	 "ST_t-channel_antitop",
	 "ST_tW_top",
	 "ST_tW_antitop" }},
};

std::map<TString,double> eraLumi = {
  {"2016_preVFP", 19520.},
  {"2016_postVFP",16810.},
  {"2017",        41480.},
  {"2018",        59830.}
};

std::map<TString, TString> LUMI_label = {
  {"2016"     ,"2016, 36.3 fb^{-1}"},
  {"2016_preVFP" ,"2016, preVFP, 19.5 fb^{-1}"},
  {"2016_postVFP","2016, postVFP, 16.8 fb^{-1}"},
  {"2017"     ,"2017, 41.5 fb^{-1}"},
  {"2018"     ,"2018, 59.8 fb^{-1}"}
};


std::map<TString, double> sample_xsec_2016pre = {
  {"WW_13TeV-pythia8",118.7},
  {"WZ_13TeV-pythia8",27.68},
  {"ZZ_13TeV-pythia8",12.19},
  {"WJetsToLNu",61526.7},
  {"ST_t-channel_top",136.02},
  {"ST_t-channel_antitop",80.95},
  {"ST_tW_top",35.85}, 
  {"ST_tW_antitop",35.85},
  {"TTTo2L2Nu",88.29},
  {"TTToSemiLeptonic",365.35},
  {"TTToHadronic",377.96},
  {"DYJetsToLL_M-10to50",21610.0},
  {"DYJetsToLL_M-50",6077.22},
  {"DYJetsToTT_M-50",6077.22},
  {"DYJetsToLL_M-50_amcatnlo",6077.22},
  {"DYJetsToTT_M-50_amcatnlo",6077.22}
  
};

std::map<TString, double> sample_xsec_2016post = {
  {"WW_13TeV-pythia8",118.7},
  {"WZ_13TeV-pythia8",27.68},
  {"ZZ_13TeV-pythia8",12.19},
  {"WJetsToLNu",61526.7},
  {"ST_t-channel_top",136.02},
  {"ST_t-channel_antitop",80.95},
  {"ST_tW_top",35.85}, 
  {"ST_tW_antitop",35.85},
  {"TTTo2L2Nu",88.29},
  {"TTToSemiLeptonic",365.35},
  {"TTToHadronic",377.96},
  {"DYJetsToLL_M-10to50",21610.0},
  {"DYJetsToLL_M-50",6077.22},
  {"DYJetsToTT_M-50",6077.22},
  {"DYJetsToLL_M-50_amcatnlo",6077.22},
  {"DYJetsToTT_M-50_amcatnlo",6077.22}
  
};

std::map<TString, double> sample_xsec_2017 = {
  {"WW_13TeV-pythia8",118.7},
  {"WZ_13TeV-pythia8",27.68},
  {"ZZ_13TeV-pythia8",12.19},
  {"WJetsToLNu",61526.7},
  {"ST_t-channel_top",136.02},
  {"ST_t-channel_antitop",80.95},
  {"ST_tW_top",35.85}, 
  {"ST_tW_antitop",35.85},
  {"TTTo2L2Nu",88.29},
  {"TTToSemiLeptonic",365.35},
  {"TTToHadronic",377.96},
  {"DYJetsToLL_M-10to50",21610.0},
  {"DYJetsToLL_M-50",6077.22}, 
  {"DYJetsToTT_M-50",6077.22}, 
  {"DYJetsToLL_M-50_amcatnlo",6077.22},
  {"DYJetsToTT_M-50_amcatnlo",6077.22}
  
};

std::map<TString, double> sample_xsec_2018 = {
  {"WW_13TeV-pythia8",118.7},
  {"WZ_13TeV-pythia8",27.68},
  {"ZZ_13TeV-pythia8",12.19},
  {"WJetsToLNu",61526.7},
  {"ST_t-channel_top",136.02},
  {"ST_t-channel_antitop",80.95},
  {"ST_tW_top",35.85}, 
  {"ST_tW_antitop",35.85},
  {"TTTo2L2Nu",88.29},
  {"TTToSemiLeptonic",365.35},
  {"TTToHadronic",377.96},
  {"DYJetsToLL_M-10to50",21610.0},
  {"DYJetsToLL_M-50",6077.22}, 
  {"DYJetsToTT_M-50",6077.22}, 
  {"DYJetsToLL_M-50_amcatnlo",6077.22},
  {"DYJetsToTT_M-50_amcatnlo",6077.22}
  
};

std::map<TString, std::map<TString,double> > sample_xsec_ERA = {
  {"2016_preVFP",sample_xsec_2016pre},
  {"2016_postVFP",sample_xsec_2016post},
  {"2017",sample_xsec_2017},
  {"2018",sample_xsec_2018}
};

// systematic uncertainties
// applied for plotting
// somewhat arbitrary
std::map<TString,double> sysSamplesZTT = {
  {"ZLL",1.5},
  {"ZTT",0.15},
  {"EWK",0.2},
  {"TT" ,0.2},
  {"QCD",0.2}
};

std::map<TString,double> sysSamplesZMM = {
  {"ZLL",0.07},
  {"ZTT",0.07},
  {"EWK",0.2},
  {"TT" ,0.2},
};

std::map<TString, TString > xtitleZmm = {
  {"massMuMuH","dimuon mass [GeV]"},
  {"ptMuMuH","dimuon p_{T} [GeV]"},
  {"ptLeadingMuH","leading muon p_{T} [GeV]"},
  {"ptTrailingMuH","trailing muon p_{T} [GeV]"},
  {"etaLeadingMuH","leading muon #eta"},
  {"etaTrailingMuH","trailing muon #eta"}
};

std::map<TString,std::vector<float> > histRangesZmm = {
  {"massMuMuH",{30,60,120}},
  {"ptMuMuH",{50,0,100}},
  {"ptLeadingMuH",{50,0,100}},
  {"ptTrailingMuH",{50,0,100}},
  {"etaLeadingMuH",{24,-2.4,2.4}},
  {"etaTrailingMuH",{24,-2.4,2.4}},
};
  
// color map
std::map<TString,TString> colorSamples = {
  {"ZLL","#4496C8"},
  {"ZTT","#FFCC66"},
  {"EWK","#DE5A6A"},
  {"TT" ,"#9999CC"},
  {"QCD","#FFCCFF"}
};

std::vector<TString> eras = {"2016","2017","2018"};
  
//  ***** TREE from which datacards are produced ***********  
//  TTree * mutrkTree = new TTree("mutrkTree","Muon Track Tree");
//  mutrkTree->Branch("os",&os_tree,"os/O");
//  mutrkTree->Branch("mt",&mt_tree,"mt/F");
//  mutrkTree->Branch("met",&met_tree,"met/F");
//  mutrkTree->Branch("dzeta",&dzeta_tree,"dzeta/F");
//  mutrkTree->Branch("njets",&njets_tree,"njets/i");
//  mutrkTree->Branch("massmutrk",&massmutrk_tree,"massmutrk/F");
//  mutrkTree->Branch("dphimutrk",&dphimutrk_tree,"dphimutrk/F");
//  mutrkTree->Branch("trktyp",&trktyp_tree,"trktyp/I");
//  mutrkTree->Branch("trkpt",&trkpt_tree,"trkpt/F");
//  mutrkTree->Branch("trketa",&trketa_tree,"trketa/F");
//  mutrkTree->Branch("muiso",&muiso_tree,"muiso/F");
//  mutrkTree->Branch("mupt",&mupt_tree,"mupt/F");
//  mutrkTree->Branch("mueta",&mueta_tree,"mueta/F");
//  mutrkTree->Branch("weight",&weight_tree,"weight/F");

std::map<TString,TString> legend_pt = {
  {"5to10","p_{T}=[5,10] GeV"},
  {"5to15","p_{T}=[5,15] GeV"},
  {"10to15","p_{T}=[10,15] GeV"},
  {"15to20","p_{T}=[15,20] GeV"},
  {"20toInf","p_{T}>20 GeV"},
};

std::map<TString,TString> legend_mt = {
  {"lowMT","low m_{T}"},
  {"highMT","high m_{T}"},
};

std::map<TString,TString> legend_cone = {
  {"cone0","#DeltaR=0"},
  {"cone0p15","#DeltaR=0.15"},
  {"cone0p30","#DeltaR=0.3"},
  {"cone0p45","#DeltaR=0.45"}
};

// accessing files 
std::map<TString, TFile*> accessFiles(TString dir,
				      TString era) {

  std::map<TString, TFile*> files;
  TString filename = dir+"/"+era+"/Data.root";
  files["Data"] = new TFile(filename);
  if (files["Data"]==NULL || files["Data"]->IsZombie()) {
    std::cout << "file " << filename << "does not exist " << std::endl;
    exit(-1);
  }
  for (auto group : groups) {
    TString name = group.first;
    std::vector<TString> samples = group.second;
    for (auto sample : samples) {
      filename = dir+"/"+era+"/"+sample+".root";
      files[sample] = new TFile(filename);
      if (files[sample]==NULL || files[sample]->IsZombie()) {
	std::cout << "file " << filename << "does not exist " << std::endl;
	exit(-1);
      }
    }
  }

  return files;

}


// adding histograms 
void addHistograms(std::map<TString, TH1D*> histograms,
		   std::map<TString, TH1D*> added_histograms) {

  for (auto histogram : histograms) {
    TString name = histogram.first;
    TH1D * hist = histogram.second;
    TH1D * added_hist = added_histograms[name];
    hist->Add(hist,added_hist,1.,1.);
  }

}

// removing negative bins
void fixNegativeBins(TH1D * hist) {

  int nbins = hist->GetNbinsX();
  for (int ib=1; ib<=nbins; ++ib) {
    double x = hist->GetBinContent(ib);
    if (x<0) {
  hist->SetBinContent(ib,0.0);
      hist->SetBinError(ib,0.0);
    }      
  }    

}
// =====================================================
// ========== Getting histograms (ZLL) =================
// =====================================================

std::map<TString, TH1D*> GetHistos(TString histName,
				   TString era,
				   TString region,
				   std::map<TString, TFile*> files,
				   int nbins,
				   double * bins) {
  
  std::map<TString, TH1D*> histograms;
  std::map<TString, double> sample_xsec = sample_xsec_ERA[era];
  std::cout << "-- extracting histogram " << histName << std::endl;
  std::cout << "-- era : " << era << "    region : " << region << std::endl;
  std::cout << std::endl;
   
  double lumi = eraLumi[era];

  //  for (auto file : files)
  //    std::cout << file.first << " " << file.second << std::endl;


  TH1D * histDataOld = (TH1D*)files["Data"]->Get(histName);

  //  std::cout << histDataOld << " " << histDataSSOld << std::endl;

  int nBins = histDataOld->GetNbinsX();
  double xMin = histDataOld->GetBinLowEdge(1);
  double xMax = histDataOld->GetBinLowEdge(nBins+1);

  std::cout << "histogram " << histName << " : " << "nbins = " << nBins
	    << "   xmin = " << xMin
	    << "   xmax = " << xMax << std::endl;
  std::cout << std::endl;

  //  double bins[300];
  //  int nbins = nBins;
  //  std::cout << "New number of bins : ";
  //  std::cin >> nbins;
  
  histograms["Data"] = TH1DtoTH1D(histDataOld,nbins,bins,true,"_Data_"+region+"_"+era);

  for (auto group : groups) {
    TString name = group.first; 
    histograms[name] = new TH1D(name+"_"+region+"_"+era,"",nbins,bins);
    vector<TString> samples = group.second;
    for (auto sample : samples) {
      TH1D * histOld = (TH1D*)files[sample]->Get(histName);
      TH1D * hist = TH1DtoTH1D(histOld,nbins,bins,true,"_"+sample+"_"+region+"_"+era);
      TH1D * eventCount = (TH1D*)files[sample]->Get("histWeightsH");
      double nGen = eventCount->GetSumOfWeights();
      double norm = sample_xsec[sample]*lumi/nGen;
      hist->Scale(norm);
      double yield = hist->GetSumOfWeights();
      histograms[name]->Add(histograms[name],hist,1,1);
    }
  }
  std::cout << std::endl;
  double yield_data = histograms["Data"]->GetSumOfWeights();
  double yield_total = 0;
  for (auto group : groups) {
    TString name = group.first;
    TH1D * hist = histograms[name];
    double yield_sample = hist->GetSumOfWeights();
    printf("%7s : %8.0f\n",name.Data(),yield_sample);
    yield_total += yield_sample;
  }
  std::cout << std::endl;
  TString name="Total";
  printf("%7s : %8.0f\n",name.Data(),yield_total);
  name="Data";
  printf("%7s : %8.0f\n",name.Data(),yield_data);
  std::cout << std::endl;

  return histograms;

}


// =========================================================
// ========= Getting histograms from Trees =================
// =========================================================
std::map<TString, TH1D*> GetHistosFromTree(TString treeName,
					   TString cuts,
					   TString var,
					   TString era,
					   TString region,
					   map<TString, TFile*> files,
					   int nbins,
					   double * bins) {
  
  TCanvas * canv = new TCanvas("canv","",500,500);

  std::map<TString, TH1D*> histograms;
  std::map<TString, double> sample_xsec = sample_xsec_ERA[era];

  double lumi = eraLumi[era];

  std::cout << "-- filling histograms from " << treeName << std::endl;
  std::cout << "-- era : " << era << "   region : " << region << std::endl;
  std::cout << std::endl;

  TString cutsOS  = "weight*("+cuts+"&&os>0.5)";
  TString cutsSS  = "weight*("+cuts+"&&os<0.5)";
  TString cutsJet = "weight*("+cuts+"&&os>0.5&&trktyp==0)"; // jet->tau fake

  TString nameOS = "Data_OS_"+region+"_"+era;
  TString nameSS = "Data_SS_"+region+"_"+era;
  histograms["Data"] = new TH1D(nameOS,"",nbins,bins);
  histograms["Data_SS"] = new TH1D(nameSS,"",nbins,bins);
  TTree * tree = (TTree*)files["Data"]->Get(treeName);
  TString draw=var+">>"+nameOS;
  std::cout << draw << "  " << cutsOS << std::endl;
  tree->Draw(draw,cutsOS);
  draw=var+">>"+nameSS;
  std::cout << draw << "  " << cutsSS <<std::endl;
  tree->Draw(draw,cutsSS);
  std::cout << std::endl;

  nameOS = "DataSubtr_OS_"+region+"_"+era;
  nameSS = "DataSubtr_SS_"+region+"_"+era;
  histograms["DataSubtr"] = (TH1D*)histograms["Data"]->Clone(nameOS);
  histograms["DataSubtr_SS"] = (TH1D*)histograms["Data_SS"]->Clone(nameSS);

  TString nameJet;
  for (auto group : groups) {

    TString name = group.first; 
    histograms[name] = new TH1D(name+"_"+region+"_"+era,"",nbins,bins);
    nameSS = name+"_SS";
    histograms[nameSS] = new TH1D(name+"_"+region+"_"+era+"_SS","",nbins,bins);
    nameJet = name+"_Jet";
    histograms[nameJet] = new TH1D(name+"_"+region+"_"+era+"_Jet","",nbins,bins);
    
    vector<TString> samples = group.second;

    for (auto sample : samples) {

      TString nameSample = sample+"_"+region+"_"+era;
      TH1D * hist = new TH1D(nameSample,"",nbins,bins);
      TString nameSampleSS = nameSample+"_SS";
      TH1D * histSS = new TH1D(nameSampleSS,"",nbins,bins);
      TString nameSampleJet = nameSample+"_Jet";
      TH1D * histJet = new TH1D(nameSampleJet,"",nbins,bins);

      TTree * tree = (TTree*)files[sample]->Get(treeName);
      draw=var+">>"+nameSample;
      tree->Draw(draw,cutsOS);
      draw=var+">>"+nameSampleSS;
      tree->Draw(draw,cutsSS);	
      draw=var+">>"+nameSampleJet;
      tree->Draw(draw,cutsJet);	
      TH1D * eventCount = (TH1D*)files[sample]->Get("histWeightsH");
      double nGen = eventCount->GetSumOfWeights();
      double norm = sample_xsec[sample]*lumi/nGen;
      hist->Scale(norm);
      histSS->Scale(norm);
      histJet->Scale(norm);
      double yield = hist->GetSumOfWeights();
      double yieldSS = histSS->GetSumOfWeights();
      //      printf("%20s: %6.0f  %6.0f\n",sample.Data(),yield,yieldSS);
      histograms[name]->Add(histograms[name],hist,1.,1.);
      histograms[nameSS]->Add(histograms[nameSS],histSS,1.,1.);
      histograms[nameJet]->Add(histograms[nameJet],histJet,1.,1.);
      
    }

    histograms["DataSubtr"]->Add(histograms["DataSubtr"],histograms[name],1.,-1.);
    histograms["DataSubtr_SS"]->Add(histograms["DataSubtr_SS"],histograms[nameSS],1.,-1.);

  }
  printf("            total :    jets :      SS |\n");
  printf("--------------------------------------+\n");
  double yield_data = histograms["Data"]->GetSumOfWeights();
  double yield_data_ss = histograms["Data_SS"]->GetSumOfWeights();
  double yield_total = 0.0;
  double yield_total_ss = 0.0;
  double yield_total_jet = 0.0;
  for (auto group : groups) {
    nameOS = group.first;
    nameSS = nameOS + "_SS";
    nameJet = nameOS + "_Jet";
    double yield_sample = histograms[nameOS]->GetSumOfWeights();
    double yield_sample_ss = histograms[nameSS]->GetSumOfWeights();
    double yield_sample_jet = histograms[nameJet]->GetSumOfWeights();
    printf("%7s = %7.0f : %7.0f : %7.0f | \n",
	   nameOS.Data(),yield_sample,yield_sample_jet,yield_sample_ss);
    yield_total += yield_sample;
    yield_total_ss += yield_sample_ss;
    yield_total_jet += yield_sample_jet;
  }

  printf("--------------------------------------+\n");
  TString name="Total";
  printf("%7s = %7.0f : %7.0f : %7.0f |\n",name.Data(),yield_total,yield_total_jet,yield_total_ss);
  name="Data";
  printf("%7s = %7.0f :         : %7.0f |\n",name.Data(),yield_data,yield_data_ss);
  std::cout << std::endl;

  delete canv;

  return histograms;

}

// ===================================================
// === Getting QCD OS/SS scale factor ================
// ===================================================

std::map<TString,double> GetNormQCD(
				    std::map<TString,TH1D*> histograms,
				    TString region,
				    TString era,
				    TString plotdir
				    ) {

  double norm[2];
  double xOS = histograms["DataSubtr"]->GetSumOfWeights();
  double xSS = histograms["DataSubtr_SS"]->GetSumOfWeights();
  // uncertainty is 50% of subtracted non-QCD MC 
  double eOS = 0.2*(histograms["Data"]->GetSumOfWeights() - histograms["DataSubtr"]->GetSumOfWeights());
  double eSS = 0.2*(histograms["Data_SS"]->GetSumOfWeights() - histograms["DataSubtr_SS"]->GetSumOfWeights());

  TString filename = plotdir+"/QCD_"+region+"_"+era+".root";
  TFile * file = new TFile(filename,"recreate");
  file->cd("");
  histograms["DataSubtr"]->Write("QCD_OS");
  histograms["DataSubtr_SS"]->Write("QCD_SS");
  file->Close();

  double scale = xOS/xSS;
  double rOS = eOS/xOS;
  double rSS = eSS/xSS;
  // limit uncertainty from above by 100%
  double rscale = TMath::Min(double(1.0),double(TMath::Sqrt(rOS*rOS+rSS*rSS)));
  double escale = rscale;
  
  std::map<TString,double> results;
  results["scale"]=scale;
  results["error"]=escale;

  std::cout << "-- era = " << era << "  region = " << region << std::endl; 
  printf("--> QCD OS/SS ratio = %4.2f +/- %4.2f\n",scale,escale);
  std::cout << std::endl;

  return results;

}

//**************************************************

TH1D * GetTemplateQCD(
		      TString era,
		      TString region,
		      std::map<TString,TH1D*> histograms,
		      std::map<TString,double> norm
		      ) {

  TString name = "QCD_"+region+"_"+era;
  TH1D * QCD = (TH1D*)histograms["DataSubtr_SS"]->Clone(name);
  int nbins = QCD->GetNbinsX();
  double scale = norm["scale"];
  double error = norm["error"];
  QCD->Scale(scale);
  for (int ib=1; ib<=nbins; ++ib) {
    double x = QCD->GetBinContent(ib);
    double estat = QCD->GetBinError(ib);
    double esys = x*error;
    double etot = TMath::Sqrt(estat*estat+esys*esys);
    QCD->SetBinContent(ib,x*scale);
    QCD->SetBinError(ib,etot);
  }

  return QCD;

}

// =====================================================
// ==================== Plotting =======================
// =====================================================

void Plot(std::map<TString, TH1D*> histograms,
	  TString xtitle,
	  TString region,
	  TString legend,
	  TString era,
	  bool plotLegend,
	  TString plotdir) {

  std::cout << "-- plotting distribution " << std::endl;
  std::cout << "-- era : " << era << "  region : " << region << std::endl;
  std::cout << std::endl;

  SetStyle();

  TString postfix = "_" + region + "_" + era + "_plot";
  TH1D * Data = (TH1D*)histograms["Data"]->Clone("Data"+postfix);
  int nBins = Data->GetNbinsX();

  std::map<TString,TH1D*> plot_sample;
  std::map<TString,double> sysSamples = sysSamplesZTT;
  if (region.Contains("zmm")) sysSamples = sysSamplesZMM;
  for (auto sysSample : sysSamples) {
    TString name = sysSample.first;
    double sys = sysSample.second;
    TH1D * hist = (TH1D*)histograms[name]->Clone(name+postfix);
    plot_sample[name] = hist;
    for (int iB=0; iB<nBins; ++iB) {
      double x = hist->GetBinContent(iB);
      double e = hist->GetBinError(iB);
      double error = TMath::Sqrt(e*e+x*x*sys*sys);
      hist->SetBinError(iB,error);
    }
  }

  TH1D * ZLL = plot_sample["ZLL"];
  TH1D * ZTT = plot_sample["ZTT"];
  TH1D * TT  = plot_sample["TT"];
  TH1D * EWK = plot_sample["EWK"];
  TH1D * QCD = NULL;
  if (region.Contains("zmm")) {
    QCD = (TH1D*)EWK->Clone("QCD_zmm"); 
    int nbins = QCD->GetNbinsX();
    for (int bin=1; bin<=nbins; ++bin) {
      QCD->SetBinContent(bin,0.);
      QCD->SetBinError(bin,0.);
    }
    plot_sample["QCD"] = QCD;
  }
  else {
    QCD = plot_sample["QCD"];
  }

  printf("ZLL   : %8.0f\n",ZLL->GetSumOfWeights());
  printf("EWK   : %8.0f\n",EWK->GetSumOfWeights());
  printf("TT    : %8.0f\n",TT->GetSumOfWeights());
  printf("QCD   : %8.0f\n",QCD->GetSumOfWeights());
  printf("ZTT   : %8.0f\n",ZTT->GetSumOfWeights());

  TT->Add(TT,QCD,1.,1.);
  EWK->Add(EWK,TT,1.,1.);
  if (region.Contains("zmm")) {
    ZTT->Add(ZTT,EWK,1.,1.);
    ZLL->Add(ZLL,ZTT,1.,1.);
  } 
  else {
    ZLL->Add(ZLL,EWK,1.,1.);
    ZTT->Add(ZTT,ZLL,1.,1.);
  }
  TH1D * bkgdErr = NULL;
  if (region.Contains("zmm")) {
    bkgdErr = (TH1D*)ZLL->Clone("bkgdErr");
  }
  else {
    bkgdErr = (TH1D*)ZTT->Clone("bkgdErr");
  }
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);
  printf("----------------\n");
  printf("Total : %8.0f\n",bkgdErr->GetSumOfWeights());
  printf("Data  : %8.0f\n",Data->GetSumOfWeights());
  std::cout << std::endl;
  
  for (int iB=1; iB<=nBins; ++iB) {
    TT->SetBinError(iB,0);
    EWK->SetBinError(iB,0);
    ZLL->SetBinError(iB,0);
    ZTT->SetBinError(iB,0);
    QCD->SetBinError(iB,0);
  }

  InitData(Data);  

  for (auto colorSample : colorSamples) {
    TString name = colorSample.first;
    TString color = colorSample.second;
    TH1D * hist = plot_sample[name];
    InitHist(hist,xtitle,"Events",TColor::GetColor(color),1001);
  }

  double ymax = Data->GetMaximum();
  if (bkgdErr->GetMaximum()>ymax) ymax = bkgdErr->GetMaximum();
  ZTT->GetYaxis()->SetRangeUser(0.,1.2*ymax);
  ZLL->GetYaxis()->SetRangeUser(0.,1.2*ymax);
  //  QCD->GetYaxis()->SetRangeUser(0.,1.2*ymax);

  TCanvas * canv = MakeCanvas("canv","",600,600);

  if (region.Contains("zmm")) {
    ZLL->Draw("h");
  }
  else {
    ZTT->Draw("h");
    ZLL->Draw("hsame");
    EWK->Draw("hsame");
    TT->Draw("hsame");
    QCD->Draw("hsame");
  }
  bkgdErr->Draw("e2same");
  Data->Draw("e1same");
    
  TLegend * leg;
  if (region.Contains("zmm")) {
    leg = new TLegend(0.66,0.57,0.90,0.75);
    leg->SetTextSize(0.04);    
  }
  else {
    leg = new TLegend(0.66,0.4,0.90,0.75);
    leg->SetTextSize(0.03);    
  }
  SetLegendStyle(leg);
  leg->SetHeader(legend);
  leg->AddEntry(Data,"Data","lp");
  if (region.Contains("zmm")) {
    leg->AddEntry(ZLL,"Simulation","f");
    //    leg->AddEntry(EWK,"other","f");
  }
  else {
    leg->AddEntry(ZTT,"Z#rightarrow#tau#tau","f");
    leg->AddEntry(ZLL,"Z#rightarrow#mu#mu","f");
    leg->AddEntry(EWK,"electroweak","f");
    leg->AddEntry(TT,"t#bar{t}","f");
    leg->AddEntry(QCD,"QCD","f");
  }
  if (plotLegend)
    leg->Draw();

  lumi_13TeV = LUMI_label[era];
  writeExtraText = true;
  extraText = "Preliminary";
  CMS_lumi(canv,4,33); 

  canv->RedrawAxis();
  canv->Update();
  canv->Print(plotdir+"/"+region+"_"+era+".png");
  std::cout << std::endl;
  delete canv;

}

// =====================================================
// ==== Writing datacards (Z->mumu) ====================
// =====================================================
void WriteDatacardsZMM(std::map<TString, TH1D*> histograms, 
		       TString region,
		       TString era,
		       TString folder) {

  double data = histograms["Data"]->GetSumOfWeights();
  double zll = histograms["ZLL"]->GetSumOfWeights();
  double ztt = histograms["ZTT"]->GetSumOfWeights();
  double tt = histograms["TT"]->GetSumOfWeights();
  double ewk = histograms["EWK"]->GetSumOfWeights();

  TH1D * Data = new TH1D("Data_zmm","",1,0.,1.); 
  Data->SetBinContent(1,data);
  TH1D * ZLL = new TH1D("ZLL_zmm","",1,0.,1.);
  ZLL->SetBinContent(1,zll+ztt);
  TH1D * TT = new TH1D("TT_zmm","",1,0.,1.);
  TT->SetBinContent(1,tt);
  TH1D * EWK = new TH1D("EWK_zmm","",1,0.,1.);
  EWK->SetBinContent(1,ewk);
  
  TString baseName = folder+"/"+region + "_" + era;
  // Creating datacards ->
  TString fileName = region + "_" + era + ".root";

  std::cout << "Writing shapes to " << fileName << std::endl;

  TFile * fileOutput = new TFile(folder+"/"+fileName,"recreate");
  fileOutput->cd("");
  Data->Write("data_obs");
  EWK->Write("EWK");
  TT->Write("TT");
  ZLL->Write("ZLL");
  fileOutput->Close();

  TString cardName = baseName + ".txt";
  std::cout << "Writing cards to " << cardName << std::endl;

  ostringstream str;
  str << cardName ;
  string nn = str.str();
  const char * p = nn.c_str();

  std::ofstream textFile(p);
  textFile << "imax 1   number of channels" << std::endl;
  textFile << "jmax *   number of backgrounds" << std::endl;
  textFile << "kmax *   number of nuisance parameters" << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "observation " << Data->GetSumOfWeights() << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "shapes * * "  << fileName << "     $PROCESS       $PROCESS_$SYSTEMATIC " << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "bin  ";
  for (int i=0; i<3; ++i)
    textFile << "  " << region;
  textFile << std::endl;
  textFile << "process              ZLL     EWK      TT" << std::endl;
  textFile << "process                1       2       3" << std::endl;
  textFile << "rate  ";
  textFile << ZLL->GetSumOfWeights() << "  "
	   << EWK->GetSumOfWeights() << "  "
	   << TT->GetSumOfWeights() << std::endl;
  textFile << "eff_m        lnN    1.03    1.03    1.03" << std::endl;
  textFile << "ewkXsec      lnN       -    1.05       -" << std::endl;
  textFile << "tjXsec       lnN       -       -    1.07" << std::endl;
  textFile << "ZNorm rateParam " << region << " ZLL  1.0 [0.5,1.5]" << std::endl; 
  textFile << std::endl;

}

// =====================================================
// ===== Writing datacards (simple model) ==============
// =====================================================
void WriteDatacardsSimple(std::map<TString, TH1D*> histograms, 
			  std::map<TString, double> norm,
			  TString region,
			  TString era,
			  TString folder) {

  TH1D * Data = histograms["Data"];

  TH1D * ZLL = histograms["ZLL"];
  fixNegativeBins(ZLL);

  TH1D * ZTT = histograms["ZTT"];
  fixNegativeBins(ZTT);

  TH1D * TT = histograms["TT"];
  fixNegativeBins(TT);

  TH1D * EWK = histograms["EWK"];
  fixNegativeBins(EWK);

  TH1D * QCD = histograms["QCD"];
  fixNegativeBins(QCD);

  TString baseName = folder+"/"+region + "_" + era + "_simpleModel";
  // Creating datacards ->
  TString fileName = region + "_" + era + "_simpleModel.root";

  std::cout << std::endl;
  std::cout << "Writing shapes to " << fileName << std::endl;

  TFile * fileOutput = new TFile(folder+"/"+fileName,"recreate");
  fileOutput->cd("");

  histograms["Data"]->Write("data_obs");

  EWK->Write("EWK");
  TT->Write("TT");
  ZTT->Write("ZTT");
  ZLL->Write("ZLL");
  QCD->Write("QCD");

  fileOutput->Close();

  double qcdErr = 1.0 + norm["error"];
  if (qcdErr>2.0) qcdErr = 2.0; // no greater than 100%
  char qcdErr_str[20];
  sprintf(qcdErr_str,"%4.2f",qcdErr);
  TString QCDErr(qcdErr_str);

  TString cardName = baseName + ".txt";
  std::cout << "Writing cards to " << cardName << std::endl;

  ostringstream str;
  str << cardName ;
  string nn = str.str();
  const char * p = nn.c_str();

  std::ofstream textFile(p);
  textFile << "imax 1   number of channels" << std::endl;
  textFile << "jmax *   number of backgrounds" << std::endl;
  textFile << "kmax *   number of nuisance parameters" << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "observation " << histograms["Data"]->GetSumOfWeights() << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "shapes * * "  << fileName << "     $PROCESS       $PROCESS_$SYSTEMATIC " << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "bin  ";
  for (int i=0; i<5; ++i)
    textFile << "  " << region;
  textFile << std::endl;
  textFile << "process              ZTT     EWK      TT     ZLL     QCD" << std::endl;
  textFile << "process                0       1       2       3       4" << std::endl;
  textFile << "rate  ";

  textFile << ZTT->GetSumOfWeights() << "  "
	   << EWK->GetSumOfWeights() << "  "
	   << TT->GetSumOfWeights() << "  "
	   << ZLL->GetSumOfWeights() << "  "
	   << QCD->GetSumOfWeights() << std::endl;  

  textFile << "eff_m        lnN    1.03    1.03    1.03    1.06      -" << std::endl;
  textFile << "jetFakeEWK   lnN       -    1.30       -       -      -" << std::endl;
  textFile << "jetFakeTTJ   lnN       -       -    1.30       -      -" << std::endl;
  textFile << "ewkXsec      lnN       -    1.05       -       -      -" << std::endl;
  textFile << "tjXsec       lnN       -       -    1.07       -      -" << std::endl;
  textFile << "muonVeto     lnN       -       -       -     2.0      -" << std::endl;
  textFile << "qcdNorm      lnN       -       -       -       -   " << QCDErr << std::endl;
  textFile << "lumi         lnN    1.02    1.02    1.02    1.02      -" << std::endl;
  textFile << "zjXsec       lnN    1.05       -       -       -      -" << std::endl;
  textFile << "* autoMCStats 0" << std::endl;
  
  std::cout << "Datacards " << cardName << " are created " << std::endl;
  std::cout << std::endl;

}

// =====================================================
// == Writing datacards (advanced model) ===============
// =====================================================

void WriteDatacards(std::map<TString, TH1D*> histograms, 
		    std::map<TString, double> norm,
		    bool floatFakes,
		    bool isZnorm,
		    TString region,
		    TString era,
		    TString folder) {

  TH1D * Data = histograms["Data"];

  TH1D * ZLL_Jet = histograms["ZLL_Jet"];
  TH1D * ZLL_Mu  = (TH1D*)histograms["ZLL"]->Clone("ZLL_Mu");
  ZLL_Mu->Add(ZLL_Mu,ZLL_Jet,1.,-1.);
  fixNegativeBins(ZLL_Mu);
  fixNegativeBins(ZLL_Jet);

  TH1D * ZTT_Jet = histograms["ZTT_Jet"];
  TH1D * ZTT_Tau = (TH1D*)histograms["ZTT"]->Clone("ZTT_Tau");
  ZTT_Tau->Add(ZTT_Tau,ZTT_Jet,1.,-1.);
  fixNegativeBins(ZTT_Tau);
  fixNegativeBins(ZTT_Jet);

  TH1D * TT_Jet = histograms["TT_Jet"];
  TH1D * TT_Tau = (TH1D*)histograms["TT"]->Clone("TT_Tau");
  TT_Tau->Add(TT_Tau,TT_Jet,1.,-1.);
  fixNegativeBins(TT_Tau);
  fixNegativeBins(TT_Jet);

  TH1D * EWK_Jet = histograms["EWK_Jet"];
  TH1D * EWK_Tau = (TH1D*)histograms["EWK"]->Clone("EWK_Tau");
  EWK_Tau->Add(EWK_Tau,EWK_Jet,1.,-1.);
  fixNegativeBins(EWK_Tau);
  fixNegativeBins(EWK_Jet);

  TH1D * QCD = histograms["QCD"];
  fixNegativeBins(QCD);

  TString baseName = folder+"/"+region + "_" + era;
  // Creating datacards ->
  TString fileName = region + "_" + era + ".root";

  std::cout << std::endl;
  std::cout << "Writing shapes to " << fileName << std::endl;

  TFile * fileOutput = new TFile(folder+"/"+fileName,"recreate");
  fileOutput->cd("");

  histograms["Data"]->Write("data_obs");

  EWK_Jet->Write("EWK_Jet");
  EWK_Tau->Write("EWK_Tau");

  TT_Jet->Write("TT_Jet");
  TT_Tau->Write("TT_Tau");

  ZTT_Jet->Write("ZTT_Jet");
  ZTT_Tau->Write("ZTT_Tau");

  ZLL_Jet->Write("ZLL_Jet");
  ZLL_Mu->Write("ZLL_Mu");

  QCD->Write("QCD");
  
  fileOutput->Close();

  double qcdErr = 1 + norm["error"];
  if (qcdErr>2.0) qcdErr = 2.0; // no greater than 100%
  char qcdErr_str[20];
  sprintf(qcdErr_str,"%4.2f",qcdErr);
  TString QCDErr(qcdErr_str);

  TString cardName = baseName;
  //  if (floatFakes) 
  //    cardName += "_floatFakes";
  //  if (isZnorm)
  //    cardName += "_znorm";

  cardName += ".txt";

  std::cout << "Writing cards to " << cardName << std::endl;

  ostringstream str;
  str << cardName ;
  string nn = str.str();
  const char * p = nn.c_str();

  bool lowpt = 
    region.Contains("_5to15") || 
    region.Contains("_5to10") || 
    region.Contains("_10to15");

  std::ofstream textFile(p);
  textFile << "imax 1   number of channels" << std::endl;
  textFile << "jmax *   number of backgrounds" << std::endl;
  textFile << "kmax *   number of nuisance parameters" << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "observation " << histograms["Data"]->GetSumOfWeights() << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "shapes * * "  << fileName << "     $PROCESS       $PROCESS_$SYSTEMATIC " << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "bin  ";
  for (int i=0; i<9; ++i)
    textFile << "  " << region;
  textFile << std::endl;
  textFile << "process  ZTT_Tau EWK_Tau TT_Tau ZTT_Jet EWK_Jet TT_Jet ZLL_Jet ZLL_Mu   QCD" << std::endl;
  textFile << "process       -2      -1      0       1       2      3       4      5     6" << std::endl;
  textFile << "rate "
	   << ZTT_Tau->GetSumOfWeights() << " "
	   << EWK_Tau->GetSumOfWeights() << " "
	   << TT_Tau->GetSumOfWeights() << " "
	   << ZTT_Jet->GetSumOfWeights() << " "
	   << EWK_Jet->GetSumOfWeights() << " "
	   << TT_Jet->GetSumOfWeights() << " "
	   << ZLL_Jet->GetSumOfWeights() << " "
	   << ZLL_Mu->GetSumOfWeights() << " "
	   << QCD->GetSumOfWeights() << std::endl;
  if (isZnorm) {
    textFile << "eff_m     lnN     -   1.03   1.03      -   1.03   1.03   1.03      -     -" << std::endl;
    textFile << "lumi      lnN     -   1.02   1.02      -   1.02   1.02   1.02      -     -" << std::endl;
    textFile << "zjXSec    lnN     -      -      -      -      -      -      -      -     -" << std::endl;
  }
  else {
    textFile << "eff_m     lnN  1.03   1.03   1.03   1.03   1.03   1.03   1.03   1.06     -" << std::endl;
    textFile << "lumi      lnN     -   1.02   1.02      -   1.02   1.02   1.02      -     -" << std::endl;
    textFile << "zjXSec    lnN  1.05      -      -   1.05      -      -      -   1.05     -" << std::endl;
  }
  textFile << "ttjXsec   lnN     -      -   1.07      -      -   1.07      -      -     -" << std::endl;
  textFile << "vvXsec    lnN     -   1.10      -      -      -      -      -      -     -" << std::endl;
  textFile << "wjXsec    lnN     -      -      -      -   1.05      -      -      -     -" << std::endl;
  textFile << "qcdNorm   lnN     -      -      -      -      -      -      -      -  " << QCDErr << std::endl;
  // znorm
  if (isZnorm) {
    textFile << "ZNorm rateParam " << region << " ZTT_Tau 1.0 [0.5,1.5]" << std::endl;
    textFile << "ZNorm rateParam " << region << " ZTT_Jet 1.0 [0.5,1.5]" << std::endl;
    textFile << "ZNorm rateParam " << region << " ZLL_Jet 1.0 [0.5,1.5]" << std::endl;
    textFile << "ZNorm rateParam " << region << " ZLL_Mu  1.0 [0.5,1.5]" << std::endl;
  }
  // float fakes
  if (floatFakes) {
    textFile << "JFake rateParam " << region << " ZTT_Jet 1.0 [0.0,5.0]" << std::endl;
    textFile << "JFake rateParam " << region << " ZLL_Jet 1.0 [0.0,5.0]" << std::endl;
    textFile << "JFake rateParam " << region << " TT_Jet  1.0 [0.0,5.0]" << std::endl;
    textFile << "JFake rateParam " << region << " EWK_Jet 1.0 [0.0,5.0]" << std::endl;
    if (lowpt)
      textFile << "MuonVeto  lnN     -      -      -      -      -      -      -   1.30     -" << std::endl;
    else 
      textFile << "MuonVeto rateParam " << region << " ZLL_Mu 1.0 [0.0,5.0]" << std::endl;
  }
  else {
    textFile << "JFakeEWK  lnN     -      -      -   1.50   1.50     -    1.50      -     -" << std::endl;
    textFile << "JFakeTTJ  lnN     -      -      -      -      -   1.50      -      -     -" << std::endl;
    if (lowpt)
      textFile << "MuonVeto  lnN     -      -      -      -      -      -      -   1.30     -" << std::endl; 
    else 
      textFile << "MuonVeto  lnN     -      -      -      -      -      -      -   2.00     -" << std::endl;
  }
  textFile << "* autoMCStats 0" << std::endl;
  
  std::cout << "Datacards " << cardName << " are created " << std::endl;
  std::cout << std::endl;

}

int main (int argc, char * argv[])
{

  if (argc!=3) {
    std::cout << "Usage : createCardsZtt [config] [era]" << std::endl;
    exit(-1);
  }

  if (fopen(argv[1],"r")==NULL) {
    std::cout << "Configuration file " << argv[1] << " does not exist" << std::endl;
    exit(-1);
  }

  // ************************
  // **** CONFIGURATION *****
  // ************************
  Config cfg(argv[1]);
  const string inputdir = cfg.get<string>("InputFolder");
  const string plotdir = cfg.get<string>("PlotFolder"); 
  const string carddir = cfg.get<string>("DatacardsFolder");
  const float xmin = cfg.get<float>("XMin");
  const float xmax = cfg.get<float>("XMax");
  const int nbins = cfg.get<int>("NBins");
  const string histName = cfg.get<string>("HistNameZmm");
  const float xminZmm = cfg.get<float>("XMinZmm");
  const float xmaxZmm = cfg.get<float>("XMaxZmm");
  const int nbinsZmm = cfg.get<int>("NBinsZmm");
  const string lowMT = cfg.get<string>("CutLowMT");
  const string highMT = cfg.get<string>("CutHighMT");
  const string dPhi = cfg.get<string>("DeltaPhiCut");
  const string trkPtCut = cfg.get<string>("TrkPtCutForCones");
  const bool useZmumu = cfg.get<bool>("UseZMuMu");
  const bool floatFakes = cfg.get<bool>("FloatFakes");
  const bool isNLO = cfg.get<bool>("IsNLO");
  const bool controlPlots = cfg.get<bool>("ControlPlotsZmm"); 
  // **** end configuration ****

  // define regions and measurement bins
  vector<TString> pTbins = {"5to10","5to15","10to15","15to20","20toInf"};
  vector<TString> dRcones = {"cone0","cone0p15","cone0p30","cone0p45"};
  vector<TString> mT_regions = {"lowMT", "highMT"};

  // define cuts ->
  TString TrkPtCut = TString("trkpt>")+TString(trkPtCut);
  TString DeltaPhiCut = TString("dphimutrk>")+TString(dPhi);
  TString CutLowMT = TString("mt<")+TString(lowMT);
  TString CutHighMT = TString("mt>")+TString(highMT);
  TString DatacardsDir(carddir);

  // plotting and datacards folders ->
  TString PlotDir(plotdir);
  TString dir(inputdir);
  
  // LO or NLO DT samples ->
  if (isNLO) {
    groups["ZLL"] = {"DYJetsToLL_M-10to50","DYJetsToLL_M-50_amcatnlo"};
    groups["ZTT"] = {"DYJetsToTT_M-50_amcatnlo"};
  }
  else {
    groups["ZLL"] = {"DYJetsToLL_M-10to50","DYJetsToLL_M-50"};
    groups["ZTT"] = {"DYJetsToTT_M-50"};
  }

  std::map<TString,TString> ptcuts = {
    {"5to15","trkpt>5.0&&trkpt<15.0"},
    {"5to10","trkpt>5.0&&trkpt<10.0"},
    {"10to15","trkpt>10.0&&trkpt<15.0"},
    {"15to20","trkpt>15.0&&trkpt<20.0"},
    {"20toInf","trkpt>20.0"},
  };

  std::map<TString,TString> mtcuts = {
    {"lowMT",CutLowMT},
    {"highMT",CutHighMT},
  };

  TString era(argv[2]);

  if (era=="2016"||era=="2017"||era=="2018"||era=="2016_preVFP"||era=="2016_postVFP") {
    std::cout << std::endl;
    std::cout << "Era : " << era << std::endl;
    std::cout << std::endl;
  }
  else {
    std::cout << "Uknown era : " << era << std::endl;
    exit(-1);
  }

  // tree name
  TString treeName = "mutrkTree";

  SetStyle();
  gROOT->SetBatch(true);
  
  // accessing files 
  std::map<TString,TFile*> files;
  std::map<TString,TFile*> filesX;
  TString era_current = era;
  if (era=="2016") era_current = "2016_preVFP";
  files = accessFiles(dir,era_current);
  if (era=="2016") {
  era_current = "2016_postVFP";
    filesX = accessFiles(dir,era_current);
  }

  // binning of m(mu,trk) histograms
  double x0 = xmin;
  double width = (xmax-xmin)/double(nbins);
  double bins[100];
  for (int iB=0; iB<=nbins; ++iB)
    bins[iB] = x0 + width*double(iB);


  // ZLL histogram name -->
  TString HistNameDY(histName);

  std::map<TString, std::vector<float> > histogramsDY;
  if (controlPlots) 
    histogramsDY = histRangesZmm;
  else 
    histogramsDY[HistNameDY] = {float(nbinsZmm),xminZmm,xmaxZmm};

  for (auto histoDY : histogramsDY) {
  
    TString HistName = histoDY.first;
    std::vector<float> bn = histoDY.second;

    int nbinsDY = int(bn[0]+0.1);
    double xminDY = bn[1];
    double xmaxDY = bn[2];
    double binsDY[100];

    double widthDY = (xmaxDY-xminDY)/double(nbinsDY);
    for (int bin=0; bin<=nbinsDY; ++bin) 
      binsDY[bin] = xminDY + double(bin)*widthDY;

    TString regionZmm("zmm");
    if (controlPlots) { 
      regionZmm += "_"+HistName;
    }

    era_current = era;
    if (era=="2016") era_current = "2016_preVFP";
    std::map<TString, TH1D*> histsZmm = GetHistos(HistName,
						  era_current,
						  regionZmm,
						  files,
						  nbinsDY,
						  binsDY);
    if (era=="2016") {
      era_current = "2016_postVFP";
      std::map<TString, TH1D*> histsZmmX = GetHistos(HistName,
						     era_current,
						     regionZmm,
						     filesX,
						     nbinsDY,
						     binsDY);
      addHistograms(histsZmm,histsZmmX);
      
    }

    TString legend("Z#rightarrow#mu#mu");
    bool plotLegend = true;
    if (HistName=="etaLeadingMuH"||HistName=="etaTrailingMuH")
      plotLegend = false;

    if (controlPlots) {
      std::cout << std::endl;
      std::cout << "-- Control plot Z->mumu (rebinned) -> " << HistName << " : nbins = " << nbinsDY
		<< "   xmin = " << xminDY
		<< "   xmax = " << xmaxDY << std::endl;
    }
  
    Plot(histsZmm,xtitleZmm[HistName],regionZmm,legend,era,plotLegend,PlotDir);
    if (!controlPlots) 
      WriteDatacardsZMM(histsZmm,regionZmm,era,DatacardsDir);

    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

  }
  if (controlPlots)
    return 0;


  // *********************
  // *** track pT bins ***
  // *********************
  for (auto pt : pTbins) { // track pT bins
    for (auto mt : mT_regions) { // low and high MT regions

      // antiisolated muon region
      // computing OS/SS extrapolation factors for QCD
      TString region = "ztt_" + mt + "_" + pt + "_antiiso";
      TString baseCuts = "muiso>0.2&&"+mtcuts[mt]+"&&"+ptcuts[pt]+"&&"+DeltaPhiCut;
      TString var("massmutrk");
      era_current = era;
      if (era=="2016") era_current = "2016_preVFP";
      std::map<TString,TH1D*> histosInvIso = 
	GetHistosFromTree(treeName,
			  baseCuts,
			  var,
			  era_current,
			  region,
			  files,
			  nbins,
			  bins);
      if (era=="2016") {
	era_current = "2016_postVFP";
	std::map<TString,TH1D*> histosInvIsoX = 
	  GetHistosFromTree(treeName,
			    baseCuts,
			    var,
			    era_current,
			    region,
			    filesX,
			    nbins,
			    bins);
	addHistograms(histosInvIso,histosInvIsoX);
      }

      // computing OS/SS extrapolation
      // factors for QCD background
      std::map<TString, double> norm = 
	GetNormQCD(histosInvIso,
		   region,
		   era,
		   PlotDir);

      // isolated muon region 
      // building templates for datacards
      region = "ztt_" + mt + "_" + pt;
      baseCuts = "muiso<0.15&&"+mtcuts[mt]+"&&"+ptcuts[pt]+"&&"+DeltaPhiCut;
      era_current = era;
      if (era=="2016") era_current = "2016_preVFP";
      std::map<TString,TH1D*> histos = 
	GetHistosFromTree(
			  treeName,
			  baseCuts,
			  var,
			  era_current,
			  region,
			  files,
			  nbins,
			  bins);
      if (era=="2016") {
	era_current = "2016_postVFP";
	std::map<TString,TH1D*> histosX = 
	  GetHistosFromTree(
			    treeName,
			    baseCuts,
			    var,
			    era_current,
			    region,
			    filesX,
			    nbins,
			    bins);
      
	addHistograms(histos,histosX);
      }

      TH1D * histQCD = GetTemplateQCD(era,
				      region,
				      histos,
				      norm);
      histos["QCD"] = histQCD;      
      TString legend = legend_mt[mt] + " " + legend_pt[pt];
      TString xtitle("m_{#mu,trk} [GeV]");

      bool isFakes = floatFakes;
      bool isZnorm = useZmumu;
      WriteDatacards(histos,norm,isFakes,isZnorm,region,era,DatacardsDir);
      bool plotLegend = true;
      if (mt=="highMT"&&(pt=="20toInf"))
	plotLegend = false;
      Plot(histos,xtitle,region,legend,era,plotLegend,PlotDir);

    }
  }

  // ********************
  // *** random cones ***
  // ********************
  for (auto cone : dRcones) { // cone
    for (auto mt : mT_regions) { // low and high MT regions

      bool isNonZeroCone = cone=="cone0p15" || cone=="cone0p30" || cone=="cone0p45";
      // antiisolated muon region
      // computing OS/SS extrapolation factors for QCD
      TString region = "ztt_" + mt + "_" + cone + "_antiiso";
      TString baseCuts = "muiso>0.2&&"+mtcuts[mt]+"&&"+DeltaPhiCut+"&&"+TrkPtCut;
      if (isNonZeroCone) 
	baseCuts += "&&"+cone+">0.5";
      TString var("massmutrk");
      era_current = era;
      if (era=="2016") era_current = "2016_preVFP";
      std::map<TString,TH1D*> histosInvIso = 
	GetHistosFromTree(treeName,
			  baseCuts,
			  var,
			  era_current,
			  region,
			  files,
			  nbins,
			  bins);
      if (era=="2016") {
	era_current = "2016_postVFP";
	std::map<TString,TH1D*> histosInvIsoX = 
	  GetHistosFromTree(treeName,
			    baseCuts,
			    var,
			    era_current,
			    region,
			    filesX,
			    nbins,
			    bins);
	addHistograms(histosInvIso,histosInvIsoX);
      }

      // computing OS/SS extrapolation
      // factors for QCD background
      std::map<TString, double> norm = 
	GetNormQCD(histosInvIso,
		   region,
		   era,
		   PlotDir);

      // isolated muon region 
      // building templates for datacards
      region = "ztt_" + mt + "_" + cone;
      baseCuts = "muiso<0.15&&"+mtcuts[mt]+"&&"+DeltaPhiCut+"&&"+TrkPtCut;
      if (isNonZeroCone)
	baseCuts += "&&"+cone+">0.5";
      era_current = era;
      if (era=="2016") era_current = "2016_preVFP";
      std::map<TString,TH1D*> histos = 
	GetHistosFromTree(
			  treeName,
			  baseCuts,
			  var,
			  era_current,
			  region,
			  files,
			  nbins,
			  bins);
      if (era=="2016") {
	era_current = "2016_postVFP";
	std::map<TString,TH1D*> histosX = 
	  GetHistosFromTree(
			    treeName,
			    baseCuts,
			    var,
			    era_current,
			    region,
			    filesX,
			    nbins,
			    bins);
      
	addHistograms(histos,histosX);
      }

      TH1D * histQCD = GetTemplateQCD(era,
				      region,
				      histos,
				      norm);
      histos["QCD"] = histQCD;      
      TString legend = legend_mt[mt] + " " + legend_cone[cone];
      TString xtitle("m_{#mu,trk} [GeV]");

      bool isFakes = floatFakes;
      bool isZnorm = useZmumu;
      WriteDatacards(histos,norm,isFakes,isZnorm,region,era,DatacardsDir);
      bool plotLegend = true;
      Plot(histos,xtitle,region,legend,era,plotLegend,PlotDir);

    }
  }


}
