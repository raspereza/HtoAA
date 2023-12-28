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

  {"ZLL",{"DYJetsToLL_M-10to50",
	  "DYJetsToLL_M-50"}},

  {"ZTT",{"DYJetsToTT_M-50"}}
  
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
  //  {"DYJetsToLL_M-50",1.16*6077.22}, // scale factor from Z->mumu
  //  {"DYJetsToTT_M-50",1.16*6077.22}, // scale factor from Z->mumu
  {"DYJetsToLL_M-50",6077.22},
  {"DYJetsToTT_M-50",6077.22}
  
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
  //  {"DYJetsToLL_M-50",1.16*6077.22}, // scale factor from Z->mumu
  //  {"DYJetsToTT_M-50",1.16*6077.22}, // scale factor from Z->mumu
  {"DYJetsToLL_M-50",6077.22},
  {"DYJetsToTT_M-50",6077.22},
  
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
  //  {"DYJetsToLL_M-50",0.95*6077.22}, // scale factor from Z->mumu
  //  {"DYJetsToTT_M-50",0.95*6077.22}, // scale factor from Z->mumu
  {"DYJetsToLL_M-50",6077.22}, // scale factor from Z->mumu
  {"DYJetsToTT_M-50",6077.22}, // scale factor from Z->mumu
  
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
  //  {"DYJetsToLL_M-50",1.07*6077.22}, // scale factor from Z->mumu
  //  {"DYJetsToTT_M-50",1.07*6077.22}, // scale factor from Z->mumu
  {"DYJetsToLL_M-50",6077.22}, // scale factor from Z->mumu
  {"DYJetsToTT_M-50",6077.22}, // scale factor from Z->mumu
  
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
std::map<TString,double> sysSamples = {
  {"ZLL",1.5},
  {"ZTT",0.15},
  {"EWK",0.3},
  {"TT" ,0.3},
  {"QCD",0.2}
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
  std::cout << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl;
  std::cout << " Extracting histogram " << histName << std::endl;
  std::cout << " Region " << region << std::endl;
  std::cout << std::endl;
   
  double lumi = eraLumi[era];

  //  for (auto file : files)
  //    std::cout << file.first << " " << file.second << std::endl;


  TH1D * histDataOld = (TH1D*)files["Data"]->Get(histName);

  //  std::cout << histDataOld << " " << histDataSSOld << std::endl;

  int nBins = histDataOld->GetNbinsX();
  double xMin = histDataOld->GetBinLowEdge(1);
  double xMax = histDataOld->GetBinLowEdge(nBins+1);

    std::cout << "Histogram " << histName << " : " << "nbins = " << nBins
	    << " , min = " << xMin
	    << " , max = " << xMax << std::endl;
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

  std::cout << "+++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl;
  std::cout << "Filling histograms from " << treeName << std::endl;
  std::cout << "Era : " << era << "   Region " << region << std::endl;
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
  TFile * file = new TFile(filename+".root","recreate");
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

  std::cout << std::endl;
  std::cout << "Era = " << era << "  Region = " << region << std::endl; 
  printf("QCD OS/SS ratio = %4.2f +/- %4.2f\n",scale,escale);
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
	  TString plotdir) {

  std::cout << std::endl;
  std::cout << "Plotting distributions " << region << std::endl;
  std::cout << "Era = " << era << "  Region " << region << std::endl;
  std::cout << std::endl;

  SetStyle();

  TString postfix = "_" + region + "_" + era + "_plot";
  TH1D * Data = (TH1D*)histograms["Data"]->Clone("Data"+postfix);
  int nBins = Data->GetNbinsX();

  std::map<TString,TH1D*> plot_sample;

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
  TH1D * QCD = plot_sample["QCD"];

  printf("ZLL   : %7.0f\n",ZLL->GetSumOfWeights());
  printf("EWK   : %7.0f\n",EWK->GetSumOfWeights());
  printf("TT    : %7.0f\n",TT->GetSumOfWeights());
  printf("QCD   : %7.0f\n",QCD->GetSumOfWeights());
  printf("ZTT   : %7.0f\n",ZTT->GetSumOfWeights());

  TT->Add(TT,QCD,1.,1.);
  EWK->Add(EWK,TT,1.,1.);
  ZLL->Add(ZLL,EWK,1.,1.);
  ZTT->Add(ZTT,ZLL,1.,1.);

  TH1D * bkgdErr = (TH1D*)ZTT->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);
  printf("----------------\n");
  printf("Total : %7.0f\n",bkgdErr->GetSumOfWeights());
  printf("Data  : %7.0f\n",histograms["Data"]->GetSumOfWeights());
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
  if (ZTT->GetMaximum()>ymax) ymax = ZTT->GetMaximum();
  ZTT->GetYaxis()->SetRangeUser(0.,1.2*ymax);

  TCanvas * canv = MakeCanvas("canv","",600,600);

  if (region=="zll") {
    ZLL->Draw("h");
    EWK->Draw("hsame");
    TT->Draw("hsame");
  }
  else {
    ZTT->Draw("h");
    ZLL->Draw("hsame");
  }
  if (region=="zll") {
    TT->Draw("hsame");
  }
  else {
    TT->Draw("hsame");
    QCD->Draw("hsame");
  }
  bkgdErr->Draw("e2same");
  Data->Draw("e1same");
    
  TLegend * leg = new TLegend(0.63,0.4,0.90,0.75);
  SetLegendStyle(leg);
  leg->SetHeader(legend);
  leg->SetTextSize(0.03);
  leg->AddEntry(Data,"Data","lp");
  if (region=="zll") {
    leg->AddEntry(ZLL,"Z#rightarrow#mu#mu","f");
    leg->AddEntry(EWK,"electroweak","f");
    leg->AddEntry(TT,"t#bar{t}","f");    
  }
  else {
    leg->AddEntry(ZTT,"Z#rightarrow#tau#tau","f");
    leg->AddEntry(ZLL,"Z#rightarrow#mu#mu","f");
    leg->AddEntry(EWK,"electroweak","f");
    leg->AddEntry(TT,"t#bar{t}","f");
    leg->AddEntry(QCD,"QCD","f");
  }
  leg->Draw();

  lumi_13TeV = LUMI_label[era];
  writeExtraText = true;
  extraText = "Preliminary";
  CMS_lumi(canv,4,33); 

  canv->RedrawAxis();
  canv->Update();
  std::cout << std::endl;
  canv->Print(plotdir+"/"+region+"_"+era+".png");
  std::cout << std::endl;
  delete canv;

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

  TString baseName = folder+"/"+region + "_" + era + "_simple";
  // Creating datacards ->
  TString fileName = baseName + ".root";

  std::cout << std::endl;
  std::cout << "Writing shapes to " << fileName << std::endl;

  TFile * fileOutput = new TFile(fileName,"recreate");
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
  textFile << "process           ZTT     EWK      TT     ZLL     QCD" << std::endl;
  textFile << "process             0       1       2       3       4" << std::endl;
  textFile << "rate  ";

  textFile << ZTT->GetSumOfWeights() << "  "
	   << EWK->GetSumOfWeights() << "  "
	   << TT->GetSumOfWeights() << "  "
	   << ZLL->GetSumOfWeights() << "  "
	   << QCD->GetSumOfWeights() << std::endl;  

  textFile << "eff_m     lnN    1.03    1.03    1.03    1.03       -" << std::endl;
  textFile << "jetFake   lnN       -     1.1     1.1       -       -" << std::endl;
  textFile << "ewkXsec   lnN       -    1.05       -       -       -" << std::endl;
  textFile << "tjXsec    lnN       -       -    1.07       -       -" << std::endl;
  textFile << "zllNorm   lnN       -       -       -     2.0       -" << std::endl;
  textFile << "qcdNorm   lnN       -       -       -       -   " << QCDErr << std::endl;
  textFile << "lumi      lnN    1.02    1.02    1.02    1.02       -" << std::endl;
  textFile << "zjXsec    lnN    1.05       -       -       -       -" << std::endl;
  textFile << "* autoMCStats 0" << std::endl;
  
  std::cout << "Datacards " << cardName << " are created " << std::endl;
  std::cout << std::endl;

}

// =====================================================
// == Writing datacards (advanced model) ===============
// =====================================================

void WriteDatacards(std::map<TString, TH1D*> histograms, 
		    std::map<TString, double> norm,
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

  TString cardName = baseName + ".txt";
  std::cout << "Writing cards to " << cardName << std::endl;

  ostringstream str;
  str << cardName ;
  string nn = str.str();
  const char * p = nn.c_str();

  bool lowpt = (era=="2018" || era=="2016") &&
    (region.Contains("_5to15") || region.Contains("_5to10") || region.Contains("_10to15"));


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
  textFile << "eff_m     lnN  1.03   1.03   1.03   1.03   1.03   1.03   1.03   1.06     -" << std::endl;
  textFile << "ttjXsec   lnN     -      -   1.07      -      -   1.07      -      -     -" << std::endl;
  textFile << "ewkXsec   lnN     -   1.10      -      -   1.10      -      -      -     -" << std::endl;
  textFile << "zjXSec    lnN  1.05      -      -   1.05      -      -      -      -     -" << std::endl;
  // to make fit more stable in low pt bin for eras 2016 and 2018
  if (lowpt) 
    textFile << "qcdNorm   lnN     -      -      -      -      -      -      -      -  1.01" << std::endl;
  else 
    textFile << "qcdNorm   lnN     -      -      -      -      -      -      -      - " << QCDErr << std::endl;
  textFile << "lumi      lnN     -   1.02   1.02      -   1.02   1.02   1.02      -     -" << std::endl;
  // to make fit more stable in low pt bin for eras 2016 and 2018
  if (lowpt) {
    textFile << "JFake rateParam " << region << " ZTT_Jet 1.0 [0.0,5.0]" << std::endl;
    textFile << "JFake rateParam " << region << " TT_Jet  1.0 [0.0,5.0]" << std::endl;
    textFile << "JFake rateParam " << region << " EWK_Jet 1.0 [0.0,5.0]" << std::endl;
  }
  else {
    textFile << "JFake rateParam " << region << " ZTT_Jet 1.0 [0.0,5.0]" << std::endl;
    textFile << "JFake rateParam " << region << " TT_Jet  1.0 [0.0,5.0]" << std::endl;
    textFile << "JFake rateParam " << region << " EWK_Jet 1.0 [0.0,5.0]" << std::endl;
    textFile << "JFake rateParam " << region << " ZLL_Jet 1.0 [0.0,5.0]" << std::endl;
    textFile << "MuonVeto rateParam " << region << " ZLL_Mu 1.0 [0.0,5.0]" << std::endl;
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


  Config cfg(argv[1]);
  const string inputdir = cfg.get<string>("InputFolder");
  const string plotdir = cfg.get<string>("PlotFolder"); 
  const string carddir = cfg.get<string>("DatacardsFolder");
  const float xmin = cfg.get<float>("XMin");
  const float xmax = cfg.get<float>("XMax");
  const int nbins = cfg.get<int>("NBins");
  const string lowMT = cfg.get<string>("CutLowMT");
  const string highMT = cfg.get<string>("CutHighMT");
  const bool simpleCards = cfg.get<bool>("SimpleCards");

  TString CutLowMT = TString("mt<")+TString(lowMT);
  TString CutHighMT = TString("mt>")+TString(highMT);
  TString DatacardsDir(carddir);
  TString PlotDir(plotdir);
  TString dir(inputdir);

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

  if (era=="2016"||era=="2017"||era=="2018") {
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

  std::vector<TString> pT_bins = {"5to15"};
  //  std::vector<TString> pT_bins = {"5to15","15to20","20toInf"};
  std::vector<TString> mT_regions = {"lowMT", "highMT"};
  //  std::vector<TString> mT_regions = {"lowMT"};

  // Checking ZLL sample
  
  int nbins_Zmm = 1;
  double bins_Zmm[2] = {80.,102.};
  TString region_Zmm = "zmm";
  TString histName("massMuMuH");
  era_current = era;
  if (era=="2016") era_current = "2016_preVFP";
  std::map<TString, TH1D*> histsZmm = GetHistos(histName,
						era_current,
						region_Zmm,
						files,
						nbins_Zmm,
						bins_Zmm);
  if (era=="2016") {
    era_current = "2016_postVFP";
    std::map<TString, TH1D*> histsZmmX = GetHistos(histName,
						   era_current,
						   region_Zmm,
						   filesX,
						   nbins_Zmm,
						   bins_Zmm);
  }
  exit(1);

  //
  // muon+track selection
  //
  for (auto pt : pT_bins) { // track pT bins
    for (auto mt : mT_regions) { // low and high MT regions

      // antiisolated muon region
      // computing OS/SS extrapolation factors for QCD
      TString region = "ztt_" + mt + "_" + pt + "_antiiso";
      TString baseCuts = "muiso>0.2&&"+mtcuts[mt]+"&&"+ptcuts[pt];
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
      baseCuts = "muiso<0.15&&"+mtcuts[mt]+"&&"+ptcuts[pt];
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
      //      WriteDatacardsSimple(histos,norm,region,era,DatacardsDir);
      WriteDatacards(histos,norm,region,era,DatacardsDir);
      Plot(histos,xtitle,region,legend,era,PlotDir);

    }
  }


}
