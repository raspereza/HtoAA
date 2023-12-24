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
  {"2016_pre" ,"2016, preVFP, 19.5 fb^{-1}"},
  {"2016_post","2016, postVFP, 16.8 fb^{-1}"},
  {"2017"     ,"2017, 41.5 fb^{-1}"},
  {"2018"     ,"2018, 59.8 fb^{-1}"}
};

std::map<TString,std::vector<TString>> groups = {
  {"EWK",{"WW_13TeV-pythia8",
	  "WZ_13TeV-pythia8",
	  "ZZ_13TeV-pythia8",
	  "ST_t-channel_top",
	  "ST_t-channel_antitop",
	  "ST_tW_top",
	  "ST_tW_antitop"}},

  {"WJ",{"WJetsToLNu"}},

  {"TT",{"TTTo2L2Nu",
	 "TTToSemiLeptonic",
	 "TTToHadronic"}},

  {"ZLL",{"DYJetsToLL_M-10to50",
	  "DYJetsToLL_M-50"}},

  {"ZTT",{"DYJetsToTT_M-50"}}
  
};

std::map<TString, float> sample_xsec = {
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
  
};

// systematic uncertainties
std::map<TString,double> sysSamples = {
  {"ZLL",1.0},
  {"ZTT",0.1},
  {"EWK",0.2},
  {"WJ",1.0},
  {"TT",0.2},
  {"QCD",0.1}
};

std::map<TString,TString> colorSamples = {
  {"ZLL","#4496C8"},
  {"ZTT","#FFCC66"},
  {"EWK","#DE5A6A"},
  {"WJ","#c6f74a"},
  {"TT","#9999CC"},
  {"QCD","#FFCCFF"}
};

//  signal region (to measure track id/iso scale factors)
//  TH1D * invMassMuTrkH         = new TH1D("invMassMuTrk","",500,0,500);
//  TH1D * invMassMuTrkH_2p5to5  = new TH1D("invMassMuTrk_2p5to5","",500,0,500);
//  TH1D * invMassMuTrkH_5to10   = new TH1D("invMassMuTrk_5to10","",500,0,500);
//  TH1D * invMassMuTrkH_10to15  = new TH1D("invMassMuTrk_10to15","",500,0,500);
//  TH1D * invMassMuTrkH_15to20  = new TH1D("invMassMuTrk_15to20","",500,0,500);
//  TH1D * invMassMuTrkH_20toInf = new TH1D("invMassMuTrk_20toInf","",500,0,500);
  
//  high mT control region (to constrain backgrounds)
//  TH1D * invMassMuTrkH_highMT         = new TH1D("invMassMuTrk_highMT","",500,0,500);
//  TH1D * invMassMuTrkH_2p5to5_highMT  = new TH1D("invMassMuTrk_highMT_2p5to5","",500,0,500);
//  TH1D * invMassMuTrkH_5to10_highMT   = new TH1D("invMassMuTrk_highMT_5to10","",500,0,500);
//  TH1D * invMassMuTrkH_10to15_highMT  = new TH1D("invMassMuTrk_highMT_10to15","",500,0,500);
//  TH1D * invMassMuTrkH_15to20_highMT  = new TH1D("invMassMuTrk_highMT_15to20","",500,0,500);
//  TH1D * invMassMuTrkH_20toInf_highMT = new TH1D("invMassMuTrk_highMT_20toInf","",500,0,500);

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

//  Z->mumu region (to constrain DY normalisation) 
//  TH1D * TH1D * massMuMuH    = new TH1D("massMuMuH","",22,80.,102.);

std::map<TString,TString> legend_pt = {
  {"5to10","5<p_{T}<10 GeV"},
  {"10to15","10<p_{T}<15 GeV"},
  {"15to20","15<p_{T}<20 GeV"},
  {"20toInf","p_{T}>20 GeV"},
};

std::map<TString,TString> legend_mt = {
  {"lowMT","low m_{T}"},
  {"highMT","high m{T}"},
};

// definition of cuts

std::map<TString,TString> ptcuts = {
  {"5to10","trkpt>5.0&&trkpt<10.0"},
  {"10to15","trkpt>10.0&&trkpt<15.0"},
  {"15to20","trkpt>15.0&&trkpt<20.0"},
  {"20toInf","trkpt>20.0"},
};

std::map<TString,TString> mtcuts = {
  {"lowMT","mt<30.0"},
  {"highMT","mt>40.0"},
};

// folder
TString dir = "/nfs/dust/cms/user/rasp/Run/MuTrk";
// tree name
TString treeName = "mutrkTree";

// =====================================================
// ========== Getting histograms =======================
// =====================================================

std::map<TString, TH1D*> GetHistos(TString histName,
				   TString era,
				   TString region,
				   std::map<TString, TFile*> files,
				   int nbins,
				   double * bins) {
  
  std::map<TString, TH1D*> histograms;

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
  TH1D * histDataSSOld = (TH1D*)files["Data_SS"]->Get(histName);

  //  std::cout << histDataOld << " " << histDataSSOld << std::endl;

  int nBins = histDataOld->GetNbinsX();
  double xMin = histDataOld->GetBinLowEdge(1);
  double xMax = histDataOld->GetBinLowEdge(nBins+1);

  std::cout << std::endl;
  std::cout << "Histogram " << histName << " : " << "nbins = " << nBins
	    << " , min = " << xMin
	    << " , max = " << xMax << std::endl;
  std::cout << std::endl;

  //  double bins[300];
  //  int nbins = nBins;
  //  std::cout << "New number of bins : ";
  //  std::cin >> nbins;
  
  histograms["Data"] = TH1DtoTH1D(histDataOld,nbins,bins,true,"_Data_"+region+"_"+era);
  histograms["Data_SS"] = TH1DtoTH1D(histDataSSOld,nbins,bins,true,"_Data_SS_"+region+"_"+era);
  histograms["QCD"] = (TH1D*)histograms["Data_SS"]->Clone("QCD_"+region+"_"+era);

  std::cout << std::endl;
  for (auto group : groups) {
    TString name = group.first; 
    histograms[name] = new TH1D(name+"_"+region+"_"+era,"",nbins,bins);
    histograms[name+"_SS"] = new TH1D(name+"_SS_"+region+"_"+era,"",nbins,bins);
    vector<TString> samples = group.second;
    for (auto sample : samples) {
      TH1D * histOld = (TH1D*)files[sample]->Get(histName);
      TH1D * hist = TH1DtoTH1D(histOld,nbins,bins,true,"_"+sample+"_"+region+"_"+era);
      TH1D * histSSOld = (TH1D*)files[sample+"_SS"]->Get(histName);
      TH1D * histSS = TH1DtoTH1D(histSSOld,nbins,bins,true,"_"+sample+"_SS_"+region+"_"+era);
      TH1D * eventCount = (TH1D*)files[sample]->Get("histWeightsH");
      double nGen = eventCount->GetSumOfWeights();
      double norm = sample_xsec[sample]*lumi/nGen;
      hist->Scale(norm);
      histSS->Scale(norm);
      double yield = hist->GetSumOfWeights();
      double yieldSS = histSS->GetSumOfWeights(); 
      printf("%20s: %7.0f %7.0f\n",sample.Data(),yield,yieldSS);
      histograms[name]->Add(histograms[name],hist,1,1);
      histograms[name+"_SS"]->Add(histograms[name+"_SS"],histSS,1,1);
    }
    histograms["QCD"]->Add(histograms["QCD"],histograms[name+"_SS"],1,-1);
  }
  std::cout << std::endl;
  std::cout << std::endl;
  double yield_data = histograms["Data"]->GetSumOfWeights();
  double yield_data_ss = histograms["Data_SS"]->GetSumOfWeights();
  double yieldQCD = histograms["QCD"]->GetSumOfWeights();
  double yield_total = yieldQCD;
  double yield_total_ss = yieldQCD;
  TString name("QCD");
  printf("%7s : OS = %7.0f  SS = %7.0f\n",name.Data(),yieldQCD,yieldQCD);  
  for (auto group : groups) {
    TString name = group.first;
    TH1D * hist = histograms[name];
    TH1D * hist_ss = histograms[name+"_SS"];
    double yield_sample = hist->GetSumOfWeights();
    double yield_sample_ss = hist_ss->GetSumOfWeights();
    printf("%7s : OS = %7.0f  SS = %7.0f\n",name.Data(),yield_sample,yield_sample_ss);
    yield_total += yield_sample;
    yield_total_ss += yield_sample_ss;
  }
  std::cout << std::endl;
  name="Total";
  printf("%7s : OS = %7.0f  SS = %7.0f\n",name.Data(),yield_total,yield_total_ss);
  name="Data";
  printf("%7s : OS = %7.0f  SS = %7.0f\n",name.Data(),yield_data,yield_data_ss);
  std::cout << std::endl;

  return histograms;

}

// =====================================================
// ===== Getting histograms from Trees =================
// =====================================================
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

  double lumi = eraLumi[era];

  std::cout << "+++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl;
  std::cout << "Plotting histograms from " << treeName << std::endl;
  std::cout << "Region " << region << std::endl;
  std::cout << std::endl;

  TString cutsOS = "weight*("+cuts+"&&os>0.5)";
  TString cutsSS = "weight*("+cuts+"&&os<0.5)";
  //  std::cout << "cutsOS = " << cutsOS << std::endl;
  //  std::cout << "cutsSS = " << cutsSS << std::endl;


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

  for (auto group : groups) {

    TString name = group.first; 
    histograms[name] = new TH1D(name+"_"+region+"_"+era,"",nbins,bins);
    nameSS = name+"_SS";
    histograms[nameSS] = new TH1D(name+"_"+region+"_"+era+"_SS","",nbins,bins);

    vector<TString> samples = group.second;

    for (auto sample : samples) {

      TString nameSample = sample+"_"+region+"_"+era;
      TH1D * hist = new TH1D(nameSample,"",nbins,bins);
      TString nameSampleSS = nameSample+"_SS";
      TH1D * histSS = new TH1D(nameSampleSS,"",nbins,bins);

      TTree * tree = (TTree*)files[sample]->Get(treeName);
      draw=var+">>"+nameSample;
      tree->Draw(draw,cutsOS);
      draw=var+">>"+nameSampleSS;
      tree->Draw(draw,cutsSS);	
      TH1D * eventCount = (TH1D*)files[sample]->Get("histWeightsH");
      double nGen = eventCount->GetSumOfWeights();
      double norm = sample_xsec[sample]*lumi/nGen;
      hist->Scale(norm);
      histSS->Scale(norm);
      double yield = hist->GetSumOfWeights();
      double yieldSS = histSS->GetSumOfWeights();
      printf("%20s: %6.0f  %6.0f\n",sample.Data(),yield,yieldSS);
      histograms[name]->Add(histograms[name],hist,1.,1.);
      histograms[nameSS]->Add(histograms[nameSS],histSS,1.,1.);
      
    }

    histograms["DataSubtr"]->Add(histograms["DataSubtr"],histograms[name],1.,-1.);
    histograms["DataSubtr_SS"]->Add(histograms["DataSubtr_SS"],histograms[nameSS],1.,-1.);

  }

  std::cout << std::endl;
  double yield_data = histograms["Data"]->GetSumOfWeights();
  double yield_data_ss = histograms["Data_SS"]->GetSumOfWeights();
  double yield_total = 0.0;
  double yield_total_ss = 0.0;
  for (auto group : groups) {
    nameOS = group.first;
    nameSS = nameOS + "_SS";
    double yield_sample = histograms[nameOS]->GetSumOfWeights();
    double yield_sample_ss = histograms[nameSS]->GetSumOfWeights();
    printf("%7s : %6.0f  %6.0f\n",nameOS.Data(),yield_sample,yield_sample_ss);
    yield_total += yield_sample;
    yield_total_ss += yield_sample_ss;
  }
  std::cout << std::endl;
  //  name="Total";
  //  printf("%7s : %6.0f  %6.0f\n",name.Data(),yield_total,yield_total_ss);
  TString name="Data";
  printf("%7s : %6.0f  %6.0f\n",name.Data(),yield_data,yield_data_ss);
  std::cout << std::endl;

  delete canv;

  return histograms;

}

// ===================================================
// === Getting QCD OS/SS scale factor ================
// ===================================================

std::map<TString,double> GetNormQCD(std::map<TString,TH1D*> histograms) {

  double norm[2];
  double xOS = histograms["DataSubtr"]->GetSumOfWeights();
  double xSS = histograms["DataSubtr_SS"]->GetSumOfWeights();
  double eOS = histograms["Data"]->GetSumOfWeights() - histograms["DataSubtr"]->GetSumOfWeights();
  double eSS = histograms["Data_SS"]->GetSumOfWeights() - histograms["DataSubtr_SS"]->GetSumOfWeights();

  double scale = xOS/xSS;
  double rOS = eOS/xOS;
  double rSS = eSS/xSS;
  double rscale = TMath::Sqrt(rOS*rOS+rSS*rSS);
  double escale = scale*rscale;
  
  std::map<TString,double> results;
  results["scale"]=scale;
  results["error"]=escale;

  std::cout << std::endl;
  std::cout << "QCD OS/SS ratio = " << scale << " +/- " << escale << std::endl;
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
  double error = norm["error"]/norm["scale"];
  QCD->Scale(scale);
  for (int ib=1; ib<=nbins; ++ib) {
    double x = QCD->GetBinContent(ib);
    double estat = QCD->GetBinError(ib);
    double esys = x*error;
    double etot = TMath::Sqrt(estat*estat+esys*esys);
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

  SetStyle();


  TH1D * Data = histograms["Data"];
  int nBins = Data->GetNbinsX();

  for (auto sysSample : sysSamples) {
    TString name = sysSample.first;
    double sys = sysSample.second;
    TH1D * hist = histograms[name];
    for (int iB=0; iB<nBins; ++iB) {
      double x = hist->GetBinContent(iB);
      double e = hist->GetBinError(iB);
      double error = TMath::Sqrt(e*e+x*x*sys*sys);
      hist->SetBinError(iB,error);
    }
  }

  TH1D * ZLL = histograms["ZLL"];
  TH1D * ZTT = histograms["ZTT"];
  TH1D * EWK = histograms["EWK"];
  TH1D * TT = histograms["TT"];
  TH1D * W  = histograms["WJ"];
  TH1D * QCD = histograms["QCD"];

  std::cout << std::endl;
  std::cout << "Plotting distributions in region " << region << std::endl;
  std::cout << std::endl;
  printf("ZLL   : %6.0f\n",ZLL->GetSumOfWeights());
  printf("EWK   : %6.0f\n",EWK->GetSumOfWeights());
  printf("TT    : %6.0f\n",TT->GetSumOfWeights());
  printf("W     : %6.0f\n",W->GetSumOfWeights());
  printf("QCD   : %6.0f\n",QCD->GetSumOfWeights());
  printf("ZTT   : %6.0f\n",ZTT->GetSumOfWeights());

  TT->Add(TT,EWK,1.,1.);
  QCD->Add(QCD,TT,1.,1.);
  W->Add(W,QCD,1.,1.);
  ZLL->Add(ZLL,W,1.,1.);
  ZTT->Add(ZTT,ZLL,1.,1.);

  TH1D * bkgdErr = (TH1D*)ZTT->Clone("bkgdErr");
  bkgdErr->SetFillStyle(3013);
  bkgdErr->SetFillColor(1);
  bkgdErr->SetMarkerStyle(21);
  bkgdErr->SetMarkerSize(0);

  std::cout << std::endl;
  printf("Total : %6.0f\n",bkgdErr->GetSumOfWeights());
  printf("Data  : %6.0f\n",histograms["Data"]->GetSumOfWeights());
  std::cout << std::endl;

  
  for (int iB=1; iB<=nBins; ++iB) {
    TT->SetBinError(iB,0);
    EWK->SetBinError(iB,0);
    W->SetBinError(iB,0);
    ZLL->SetBinError(iB,0);
    ZTT->SetBinError(iB,0);
    QCD->SetBinError(iB,0);
  }
  
  InitData(Data);

  InitHist(ZTT,xtitle,"Events",TColor::GetColor("#FFCC66"),1001);
  InitHist(ZLL,"","",TColor::GetColor("#4496C8"),1001);
  InitHist(QCD,"","",TColor::GetColor("#FFCCFF"),1001);
  InitHist(TT,"","",TColor::GetColor("#9999CC"),1001);
  InitHist(W,"","",TColor::GetColor("#c6f74a"),1001);
  InitHist(EWK,"","",TColor::GetColor("#DE5A6A"),1001);

  double ymax = Data->GetMaximum();
  if (ZTT->GetMaximum()>ymax) ymax = ZTT->GetMaximum();
  ZTT->GetYaxis()->SetRangeUser(0.,1.2*ymax);

  TCanvas * canv = MakeCanvas("canv","",600,600);

  ZTT->Draw("h");
  ZLL->Draw("hsame");
  W->Draw("hsame");
  QCD->Draw("hsame");
  TT->Draw("hsame");
  EWK->Draw("hsame");
  bkgdErr->Draw("e2same");
  Data->Draw("e1same");
    
  TLegend * leg = new TLegend(0.61,0.42,0.90,0.77);
  SetLegendStyle(leg);
  leg->SetHeader(legend);
  leg->SetTextSize(0.04);
  leg->AddEntry(Data,"Data","lp");
  leg->AddEntry(ZTT,"Z#rightarrow#tau#tau","f");
  leg->AddEntry(ZLL,"Z#rightarrow#mu#mu","f");
  leg->AddEntry(W,"W+jets","f");
  leg->AddEntry(QCD,"QCD","f");
  leg->AddEntry(TT,"t#bar{t}","f");
  leg->AddEntry(EWK,"electroweak","f");
  leg->Draw();

  lumi_13TeV = LUMI_label[era];
  writeExtraText = true;
  extraText = "Preliminary";
  CMS_lumi(canv,4,33); 

  canv->RedrawAxis();
  canv->Update();
  canv->Print(plotdir+"/"+region+"_"+era+".png");

  delete canv;

}

// =====================================================
// ============== Writing datacards ====================
// =====================================================

void WriteDatacards(std::map<TString, TH1D*> histograms, 
		    TString region,
		    TString era,
		    TString folder) {

  // Creating datacards ->
  TString fileName = folder+"/"+region + "_" + era + ".root";
  TFile * fileOutput = new TFile(fileName,"recreate");
  fileOutput->cd("");
  histograms["Data"]->Write("data_obs");
  histograms["EWK"]->Write("EWK");
  histograms["TT"]->Write("TT");
  histograms["WJ"]->Write("WJ");
  histograms["QCD"]->Write("QCD");
  histograms["ZLL"]->Write("ZLL");
  if (region.Contains("ztt_"))
    histograms["ZTT"]->Write("ZTT");
  else
    histograms["ZTT"]->Write("ZTT_mm");

  ostringstream str;
  TString Datacards = folder+"/"+region + "_" + era + ".txt";
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
  for (int i=0; i<6; ++i)
    textFile << "  " << region;
  textFile << std::endl;
  if (region.Contains("ztt_")) {
    textFile << "process                 WJ      TT    QCD     EWK    ZLL    ZTT" << std::endl;
    textFile << "process                  1       2      3       4      5      0" << std::endl;
  }
  else {
    textFile << "process                 WJ      TT    QCD     EWK    ZLL   ZTT_mm" << std::endl;
    textFile << "process                  1       2      3       4      5      6  " << std::endl;
  }    
  textFile << "rate                     -1      -1     -1     -1     -1    -1   " << std::endl;
  if (region.Contains("ztt_")) {      
    textFile << "eff_m          lnN    1.03    1.03      -    1.03   1.03   1.03" << std::endl;
    textFile << "ttjNorm        lnN       -    1.10      -       -      -      -" << std::endl;
    textFile << "ewkNorm        lnN       -       -      -    1.10      -      -" << std::endl;
    textFile << "qcdNorm        lnN       -       -    2.0       -      -      -" << std::endl;
    textFile << "wjNorm         lnN     2.0       -      -       -      -      -" << std::endl;
    textFile << "zllNorm        lnN       -       -      -       -    2.0      -" << std::endl;
    textFile << "lumi           lnN       -    1.02      -    1.02      -   1.02" << std::endl;
    textFile << "ZNorm rateParam "   << region << " ZTT 1.0 [0.5,1.5]" << std::endl;
  }
  else {
    textFile << "eff_m          lnN    1.06    1.06      -    1.06   1.06   1.06" << std::endl;
    textFile << "ttjNorm        lnN       -    1.10      -       -      -      -" << std::endl;
    textFile << "ewkNorm        lnN       -       -      -    1.10      -      -" << std::endl;
    textFile << "qcdNorm_mm     lnN       -       -   1.50       -      -      -" << std::endl;
    textFile << "wjNorm_mm      lnN    1.50       -      -       -      -      -" << std::endl;
    textFile << "lumi           lnN       -    1.02      -    1.02      -   1.02" << std::endl;
    textFile << "ZNorm rateParam  " << region << " ZLL 1.0 [0.0,5.0]" << std::endl;
    textFile << "ZNorm rateParam  " << region << " ZTT_mm 1.0 [0.0,5.0]" << std::endl;
  }
  textFile << "* autoMCStats 0" << std::endl;
  

}

int main (int argc, char * argv[])
{

  Config cfg(argv[1]);
  const string inputdir = cfg.get<string>("InputFolder");
  const string plotdir = cfg.get<string>("PlotFolder"); 
  const string carddir = cfg.get<string>("DatacardsFolder");

  TString era(argv[2]);
  TString dir(inputdir);

  SetStyle();
  gROOT->SetBatch(true);
  

  // accessing files 
  std::map<TString,TFile*> files;
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

  // binning of histograms
  int nbins = 18;
  double bins[100];
  for (int iB=0; iB<=nbins; ++iB)
    bins[iB] = 20 + 10.0*double(iB);


  //  std::vector<TString> pT_bins = {"5to10", "10to15", "15to20", "20toInf"};
  std::vector<TString> pT_bins = {"20toInf"};
  std::vector<TString> mT_regions = {"lowMT", "highMT"};
  //  std::vector<TString> mT_regions = {"lowMT"};
  
  for (auto pt : pT_bins) {
    for (auto mt : mT_regions) {
      TString region = "ztt_" + mt + "_" + pt + "_antiiso";
      TString baseCuts = "muiso>0.2&&"+mtcuts[mt]+"&&"+ptcuts[pt]+"&&TMath::Abs(weight)<100";
      TString var("massmutrk");
      std::map<TString,TH1D*> histosInvIso = GetHistosFromTree(
							       treeName,
							       baseCuts,
							       var,
							       era,
							       region,
							       files,
							       nbins,
							       bins);
      
      std::map<TString, double> norm = GetNormQCD(histosInvIso);

      region = "ztt_" + mt + "_" + pt;
      baseCuts = "muiso<0.15&&"+mtcuts[mt]+"&&"+ptcuts[pt];
      std::map<TString,TH1D*> histos = GetHistosFromTree(
							 treeName,
							 baseCuts,
							 var,
							 era,
							 region,
							 files,
							 nbins,
							 bins);
      TH1D * histQCD = GetTemplateQCD(era,
				      region,
				      histos,
				      norm);
      histos["QCD"] = histQCD;      
      TString legend = legend_mt[mt] + " " + legend_pt[pt];
      Plot(histos,"m_{#mu,trk} [GeV]",region,legend,era,plotdir);
    }
  }

  
//  TH1D * invMassMuTrkH_20toInf_highMT = new TH1D("invMassMuTrk_highMT_20toInf","",500,0,500);

//  Z->mumu region (to constrain DY normalisation) 
//  TH1D * TH1D * massMuMuH    = new TH1D("massMuMuH","",22,80.,102.);



//	TH1D * massMuMuH    = new TH1D("massMuMuH","",22,80.,102.);
//	TH1D * massMuMuH    = new TH1D("massMuMuH","",22,80.,102.);

//	for (int iB = 0; iB < npTbins; ++iB) CreateCards(pT_bins[iB]);
}
