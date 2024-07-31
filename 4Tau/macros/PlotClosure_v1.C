#include "HttStylesNew.cc"
#include "HtoH.h"
#include "CMS_lumi.C"

// ClosureInvMassDimuonIso
// InvMassHighMuLooseIso
// InvMassHighMuIso

void PlotClosure_v1(TString histName = "InvMassDimuonIso",
		    TString Model = "Closure",
		    TString era = "2018",
		    TString subfolder = "bin5p2",
		    bool normalized = true,
		    bool logY = true) {

  SetStyle();

  std::cout << std::endl;
  std::cout << "Era = " << era << "  subfolder = " << subfolder << std::endl;
  std::cout << "Histogram = " << histName << "  Model = " << Model << std::endl;
  std::cout << std::endl;

  int nBins = 6;
  double bins[7]     = {0,1,2,3,4,5.2,20.};

  TString dir="/nfs/dust/cms/user/rasp/Run/QCDModel/"+era+"/"+subfolder;

  TString legend = "SR";
  TString code_model = "model_"+era+"_sig";
  TString code_sel   = "sel_"+era+"_sig";
  
  if (histName.Contains("LooseIso")) {
    legend = "LooseIso";
    code_model = "model_"+era+"_loose";
    code_sel   = "sel_"+era+"_loose";
  }

  TString SamplesQCD[10] = {
    "QCD_Pt-20To30_MuEnrichedPt5",   // (1)
    "QCD_Pt-30To50_MuEnrichedPt5",   // (2)
    "QCD_Pt-50To80_MuEnrichedPt5",   // (3)
    "QCD_Pt-80To120_MuEnrichedPt5",  // (4)
    "QCD_Pt-120To170_MuEnrichedPt5", // (5)
    "QCD_Pt-170To300_MuEnrichedPt5", // (6)
    "QCD_Pt-300To470_MuEnrichedPt5", // (7)
    "QCD_Pt-470To600_MuEnrichedPt5", // (8)
    "QCD_Pt-600To800_MuEnrichedPt5", // (9)
    "QCD_Pt-800To1000_MuEnrichedPt5" // (10)
  };

  double xsecQCD[10] = {
    558528000*0.0053,  // (1)
    139803000*0.01182, // (2)
    19222500*0.02276,  // (3)
    2758420*0.03844,   // (4)
    469797*0.05362,    // (5)
    117989*0.07335,    // (6)
    7820.25*0.10196,   // (7)
    645.528*0.12242,   // (8)
    187.109*0.13412,   // (9)
    32.3486*0.14552    // (10)
  };


  TH1D * hist;
  TH1D * histModel;

  for (int iS=1; iS<10; ++iS) {
    TFile * fileSample = new TFile(dir+"/"+SamplesQCD[iS]+".root");
    TH1D * histWeightsH = (TH1D*)fileSample->Get("histWeightsH");
    std::cout << SamplesQCD[iS] << " : " << int(histWeightsH->GetSumOfWeights()) << std::endl;
    double normX = xsecQCD[iS]/histWeightsH->GetSumOfWeights();
    TH1D * histSample = (TH1D*)fileSample->Get(histName+"H");
    TH1D * histModelSample = (TH1D*)fileSample->Get(Model+histName+"H");
    
    TH1D * histRebinned = (TH1D*)TH1DtoTH1D(histSample,nBins,bins,true,"_new_"+SamplesQCD[iS]);
    TH1D * histModelRebinned = (TH1D*)TH1DtoTH1D(histModelSample,nBins,bins,true,"_new_"+SamplesQCD[iS]);
    
    if (iS==1) {
      hist = histRebinned;
      histModel = histModelRebinned;
      hist->Scale(normX);
      histModel->Scale(normX);
    }
    else {
      hist->Add(hist,histRebinned,1,normX);
      histModel->Add(histModel,histModelRebinned,1,normX);
    }
  }

  double yMax = hist->GetMaximum();
  if (histModel->GetSumOfWeights()>yMax) yMax = histModel->GetSumOfWeights();
  
  float yMin = 0;
  std::cout << std::endl;
  std::cout << "  mass bin     selection  :     model    " << std::endl;
  //           "[ 0.0, 1.0]  5.88+/- 0.32 :  5.80+/- 0.10
  //  std::cout << std::endl;
  for (int ib=1; ib<=nBins;++ib) {
    double lower = hist->GetXaxis()->GetBinLowEdge(ib);
    double upper = hist->GetXaxis()->GetBinLowEdge(ib+1);
    double x = hist->GetBinContent(ib);
    double ex = hist->GetBinError(ib);
    double y = histModel->GetBinContent(ib);
    double ey = histModel->GetBinError(ib);
    printf("[%4.1f,%4.1f] %5.2f+/-%5.2f : %5.2f+/-%5.2f\n",lower,upper,x,ex,y,ey);
  }
  std::cout << std::endl;

  printf("Direct  = %7.2f\n",hist->GetSumOfWeights());
  printf("Model   = %7.2f\n",histModel->GetSumOfWeights());
	 

  hist->SetLineColor(1);
  hist->SetMarkerColor(1);
  hist->SetMarkerSize(1.2);
  hist->SetMarkerStyle(20);
  hist->GetXaxis()->SetLabelSize(0.);
  hist->GetYaxis()->SetTitle("normalized");

  histModel->SetLineColor(2);
  histModel->SetMarkerColor(2);
  histModel->SetMarkerSize(0);
  histModel->SetMarkerStyle(0);

  if (normalized) {
    hist->Scale(1.0/hist->GetSumOfWeights());
    histModel->Scale(1.0/histModel->GetSumOfWeights());
  }

  TH1D * ratioH = (TH1D*)hist->Clone("ratioH");
  
  ratioH->GetYaxis()->SetRangeUser(0.0,2.499);

  for (int iB=1; iB<=nBins; ++iB) {
    float den = histModel->GetBinContent(iB);
    if (den<1e-6) {
      std::cout << "bin : " << iB << "  den = " << den << std::endl;
      den = 1;
    }
    float x = hist->GetBinContent(iB);
    float ratioX = x/den;
    float e = hist->GetBinError(iB);
    float ratioE = e/den;
    ratioH->SetBinContent(iB,ratioX);
    ratioH->SetBinError(iB,ratioE);
  }


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

  hist->GetYaxis()->SetRangeUser(1e-3,1.);
  hist->Draw();
  histModel->Draw("hsame");
  
  std::map<TString,TString> eraLabel = {
    {"2016_preVFP", "2016 (preVFP), 19.5 fb^{-1}"},
    {"2016_postVFP","2016 (postVFP), 16.8 fb^{-1}"},
    {"2016",        "2016 36.3 fb^{-1}"},
    {"2017",        "2017 41.5 fb^{-1}"},
    {"2018",        "2018 59.8 fb^{-1}"}
  };


  writeExtraText = true;
  extraText   = "Simulation";
  lumi_13TeV = eraLabel[era];
  CMS_lumi(upper,4,33); 

  TLegend * leg = new TLegend(0.55,0.55,0.9,0.8);
  SetLegendStyle(leg);
  leg->SetHeader(legend);
  leg->SetTextSize(0.05);
  leg->AddEntry(hist,"actual selection","lp");
  leg->AddEntry(histModel,"model","l");
  leg->Draw();

  if (logY) upper->SetLogy(true);
  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv1->cd();

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

  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.2);
  ratioH->SetLineColor(1);
  ratioH->GetYaxis()->SetRangeUser(0.0001,2.999);
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleSize(0.13);
  ratioH->GetXaxis()->SetTitleOffset(1.2);
  ratioH->GetYaxis()->SetTitle("obs/exp");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.015);
  ratioH->GetYaxis()->SetLabelSize(0.13);
  ratioH->GetYaxis()->SetTitleSize(0.14);
  ratioH->GetYaxis()->SetTitleOffset(0.6);
  ratioH->GetXaxis()->SetTickLength(0.07);
  ratioH->GetYaxis()->SetTickLength(0.04);
  ratioH->GetYaxis()->SetLabelOffset(0.01);
  ratioH->GetXaxis()->SetTitle("m_{#mu,trk} [GeV]");

  ratioH->Draw("e1");

  lower->Modified();
  lower->RedrawAxis();
  canv1->cd();
  canv1->SetSelected(canv1);
  canv1->Print(Model+histName+"_"+era+".png");


}
