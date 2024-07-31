#include "HttStylesNew.cc"
#include "HtoH.h"
#include "CMS_lumi.C"

// [Model]InvMassHighMuLooseIso
// [Model]InvMassHighMuIso

void CreateClosure(TString era = "2018") {

  SetStyle();

  TString dir("/nfs/dust/cms/user/rasp/Run/Run2016/H2aa/QCD_Model");
  TString outdir("/nfs/dust/cms/user/rasp/Run/HtoAA");

  TString code_sel_sig   = "sel_"+era+"_sig";
  TString code_model_sig = "model_"+era+"_sig";
  
  TString code_sel_loose   = "sel_"+era+"_loose";
  TString code_model_loose = "model_"+era+"_loose";

  std::map<TString, TString> inputHistNames;
  inputHistNames[code_sel_sig] = "InvMassHighMuIsoAllH";
  inputHistNames[code_model_sig] = "ModelInvMassHighMuIsoH";
  inputHistNames[code_sel_loose] = "InvMassHighMuLooseIsoAllH";
  inputHistNames[code_model_loose] = "ModelInvMassHighMuLooseIsoH";

  std::map<TString, TString> outputHistNames;
  outputHistNames[code_sel_sig]     = "mass_Selected_SR";
  outputHistNames[code_model_sig]   = "mass_Modelled_SR";
  outputHistNames[code_sel_loose]   = "mass_Selected_LooseIso";
  outputHistNames[code_model_loose] = "mass_Modelled_LooseIso";

  TString SamplesQCD[10] = {"QCD_Pt_15to30_pythia8",
			    "QCD_Pt_30to50_pythia8",
			    "QCD_Pt_50to80_pythia8",
			    "QCD_Pt_80to120_pythia8",
			    "QCD_Pt_120to170_pythia8",
			    "QCD_Pt_170to300_pythia8",
			    "QCD_Pt_300to470_pythia8",
			    "QCD_Pt_470to600_pythia8",
			    "QCD_Pt_600to800_pythia8",
			    "QCD_Pt_800to1000_pythia8"};

  double xsecQCD[10] = {1837410000,
		       140932000,
		       19204300,
		       2762530,
		       471100,
		       117276,
		       7823,
		       648.2,
		       186.9,
		       32.29};


  int nBins = 8;
  double binsX[9] = {0,1,2,3,4,5,6,8,20};
  double bins[9] =  {0,1,2,3,4,5,6,8,12}; 

  int nBinsOut = 6;
  double binsOut[7] = {0,1,2,3,4,5.2,12};
  
  map<TString,double> scaleFactor = {
    {"model_2016_sig"  ,1.01},
    {"sel_2016_sig"    ,1.02},
    {"model_2016_loose",1.08},
    {"sel_2016_loose"  ,0.98},

    {"model_2017_sig"  ,1.01},
    {"sel_2017_sig"    ,1.02},
    {"model_2017_loose",1.07},
    {"sel_2017_loose"  ,0.99},

    {"model_2018_sig"  ,1.03},
    {"sel_2018_sig"    ,1.03},
    {"model_2018_loose",1.08},
    {"sel_2018_loose"  ,1.00},

  };

  map<TString, vector<double>> shift = {
    //                      1    2    3    4    5    6    8   12
    {"model_2016_sig",  {1.05,0.93,0.80,0.95,1.23,4.04,4.85,2.6}},
    {"sel_2016_sig",    {1.05,0.98,0.87,1.02,1.22,6.00,3.22,5.1}},
    {"model_2016_loose",{1.09,0.95,0.75,0.82,0.99,1.41,1.53,4.5}},
    {"sel_2016_loose",  {0.92,1.07,1.01,1.95,1.29,1.56,1.00,9.0}},
    //                      1    2    3    4    5    6    8   12
    {"model_2017_sig",  {1.02,1.01,0.77,1.00,1.39,4.20,3.89,3.4}},
    {"sel_2017_sig",    {1.00,1.00,0.90,1.41,0.70,4.22,4.00,7.0}},
    {"model_2017_loose",{1.10,0.90,0.80,0.80,0.93,1.19,2.52,5.0}},
    {"sel_2017_loose",  {1.00,1.00,1.21,1.71,1.31,2.70,1.00,5.0}},
    //                      1    2    3    4    5    6    8   12
    {"model_2018_sig",  {0.99,0.93,0.83,0.97,1.42,4.02,4.00,3.5}},
    {"sel_2018_sig",    {1.01,1.05,0.97,1.02,0.72,2.00,6.70,9.0}},
    {"model_2018_loose",{1.19,0.89,0.77,0.82,0.95,1.30,2.10,4.5}},
    {"sel_2018_loose",  {1.02,0.97,1.03,1.84,2.29,1.05,0.60,9.0}},

  };

  std::map<TString, TH1D*> histograms;

  for (int iS=1; iS<10; ++iS) {

    TFile * fileSample = new TFile(dir+"/"+SamplesQCD[iS]+".root");
    TH1D * histWeightsH = (TH1D*)fileSample->Get("histWeightsH");
    double norm = xsecQCD[iS]/histWeightsH->GetSumOfWeights();

    std::map<TString, TH1D*> histos;
    for (auto inputHist : inputHistNames) {
      TString codename = inputHist.first;
      TString histname = inputHist.second;
      TH1D * hist = (TH1D*)fileSample->Get(histname);
      hist->Scale(norm);
      if (iS==1) {
	histograms[codename]=hist;
      }
      else {
	TH1D * histogram = histograms[codename];
	histogram->Add(histogram,hist);
      }
    }
  }

  TString outname = outdir + "/QCDModel_"+era+".root";
  TFile * outfile = new TFile(outname,"recreate");
  outfile->cd("");

  for (auto outputHist : outputHistNames) {

    TString codename = outputHist.first;
    TString histname = outputHist.second;

    TH1D * histX = histograms[codename];
    TH1D * histY = (TH1D*)TH1DtoTH1D(histX,nBins,binsX,true,"_rebinned");
    TH1D * hist = new TH1D(histname,"",nBins,bins);

    std::vector<double> shifts = shift[codename];
    double sf = scaleFactor[codename];

    for (int ib=1; ib<=nBins; ++ib) {

      double x = histY->GetBinContent(ib);
      double e = histY->GetBinError(ib);
      hist->SetBinContent(ib,sf*shifts[ib-1]*x);
      hist->SetBinError(ib,sf*shifts[ib-1]*e);      
      if (ib==nBins)
	hist->SetBinError(ib,1.4*hist->GetBinError(ib));

    }
    TH1D * histOut = (TH1D*)TH1DtoTH1D(hist,nBinsOut,binsOut,true,"_histOut");
    outfile->cd("");
    histOut->Write(histname);

  }


  outfile->Close();



}
