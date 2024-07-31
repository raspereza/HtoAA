#include "HttStylesNew.cc"
#include "HtoH.h"
#include "CMS_lumi.C"

std::map<TString, double> eraLumi = {
  {"2016_preVFP", 19520},
  {"2016_postVFP", 16810},
  {"2017", 41480},
  {"2018", 59830}
};

std::map<TString, double> nonQCDsamples= {
  {"DYJetsToLL_M-10to50",21610.0},
  {"DYJetsToLL_M-50",6077.22},
  {"TTTo2L2Nu",88.29},
  {"TTToSemiLeptonic",365.35},
  {"ST_t-channel_top",136.02},
  {"ST_t-channel_antitop",80.95},
  {"ST_tW_top",35.85},
  {"ST_tW_antitop",35.85},
  {"WW_13TeV-pythia8",118.7},
  {"WZ_13TeV-pythia8",27.68},
  {"ZZ_13TeV-pythia8",12.19}
};

std::map<TString, vector<TString> > groups = {
  {"2016",{"2016_postVFP","2016_preVFP"}},
  {"2016_preVFP",{"2016_preVFP"}},
  {"2016_postVFP",{"2016_postVFP"}},
  {"2017",{"2017"}},
  {"2018",{"2018"}},
  {"Run2",{"2016_postVFP","2016_preVFP","2017","2018"}}
};

// mass_Selected_SR
// mass_Modelled_SR
// mass_Selected_Loose
// mass_Modelled_Loose

void PlotClosure2D(
		   TString era = "Run2",
		   bool signalRegion = false,
		   bool logY = true,
		   double yminRatio = 0.201,
		   double ymaxRatio = 1.799
		   ) {
  SetStyle();

  TString legend("#bf{LooseIso}");
  TString histName("mass_Selected_LooseIso");
  TString histNameModel("mass_Modelled_LooseIso");
  TString histNameData("InvMassTrackPlusMuon2D_ControlCH");
  TString label("control");
  
  if (signalRegion) {
    legend = "#bf{SR}";
    histName = "mass_Selected_SR";
    histNameModel = "mass_Modelled_SR";
    histNameData = "InvMass2DH";
    label = "signal";
  }

  std::cout << std::endl;
  std::cout << "Era : " << era << "   Region : " << label << std::endl;
  std::cout << std::endl;

  int nBins = 6;
  double bins[7]     = {0,1,2,3,4,5.2,12};
  double binsX[7]    = {0,1,2,3,4,5.2,20};
  
  TString folder("/nfs/dust/cms/user/sreelatl/Analyses/H2aa_4tau/Run2/Jul24");
  TFile * file = new TFile("/nfs/dust/cms/user/rasp/Run/HtoAA/QCDModel_"+era+".root");
  TString pwd = TString(std::getenv("CMSSW_BASE")) + "/src/HtoAA/4Tau/macros";
  TString fileNameCorr = pwd + "/CorrCoefficients_" + label + "_mc_" + era + ".root";
  TFile * fileCorr = new TFile(fileNameCorr);
  
  TH1D * histModel = (TH1D*)file->Get(histNameModel);

  TH2D * histData_2D = new TH2D("histData_2D","",nBins,binsX,nBins,binsX);
  TH2D * histNonQCD_2D = new TH2D("histNonQCD_2D","",nBins,binsX,nBins,binsX);

  vector<TString> eras = groups[era];
  
  for (auto grp : eras) {
    TString dir = folder + "/" + grp;
    double luminosity = eraLumi[grp];
    for (auto sample : nonQCDsamples) {
      TString sampleName = sample.first;
      double xsec = sample.second;
      TString fileName = dir + "/" + sampleName + ".root";
      TFile * file = new TFile(fileName);
      TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
      TH2D * histOld = (TH2D*)file->Get(histNameData);
      TH2D * histNew = (TH2D*)TH2DtoTH2D(histOld,nBins,binsX,nBins,binsX,"_rebinned");
      double ngen = histWeightsH->GetSumOfWeights();
      double norm = xsec*luminosity/ngen;
      histNew->Scale(norm);
      histNonQCD_2D->Add(histNonQCD_2D,histNew);
    }
    TString dataSampleName = "DoubleMuon_Run" + grp;
    TString fileName = dir + "/" + dataSampleName + ".root";
    TFile * file =  new TFile(fileName);
    TH2D * histOld = (TH2D*)file->Get(histNameData);
    TH2D * histNew = (TH2D*)TH2DtoTH2D(histOld,nBins,binsX,nBins,binsX,"_rebinned");
    histData_2D->Add(histData_2D,histNew);
  }

  int nBins1D = (nBins+1)*nBins/2;
  TH1D * histData = new TH1D("histData","",nBins1D,0.,double(nBins1D));
  TH1D * histNonQCD = new TH1D("histNonQCD","",nBins1D,0.,double(nBins1D));
  TH1D * histBkgd = new TH1D("histBkgd","",nBins1D,0.,double(nBins1D));
  TH2D * corrCoeff = (TH2D*)fileCorr->Get("corrCoeff");
  
  histModel->Scale(1.0/histModel->GetSumOfWeights());

  int iBin = 0;
  for (int iB=1; iB<=nBins; ++iB) {
    for (int jB=iB; jB<=nBins; ++jB) {
      iBin++;
      double x = histModel->GetBinContent(iB);
      double y = histModel->GetBinContent(jB);
      double prod = x*y;
      double xData = histData_2D->GetBinContent(iB,jB);
      double eData = histData_2D->GetBinError(iB,jB);
      double xNonQCD = histNonQCD_2D->GetBinContent(iB,jB);
      if (iB!=jB) {
	xNonQCD += histNonQCD_2D->GetBinContent(jB,iB);
	xData += histData_2D->GetBinContent(jB,iB);
	double err = histData_2D->GetBinError(jB,iB);
	eData = TMath::Sqrt(eData*eData+err*err);
	prod *= 2;
      }
      double xcorr = corrCoeff->GetBinContent(iB,jB);
      double ecorr = corrCoeff->GetBinError(iB,jB);
      double xModel = prod*xcorr;
      double eModel = prod*ecorr;
      histData->SetBinContent(iBin,xData);
      histData->SetBinError(iBin,eData);
      histNonQCD->SetBinContent(iBin,xNonQCD);
      histNonQCD->SetBinError(iBin,0.);
      histBkgd->SetBinContent(iBin,xModel);
      histBkgd->SetBinError(iBin,eModel);
    }
  }
  
  double normData = histData->GetSumOfWeights();
  double normBkgd = histBkgd->GetSumOfWeights();
  histBkgd->Scale(normData/normBkgd);
  
  TH1D * histBkgdErr = (TH1D*)histBkgd->Clone("histBkgdErr");
  histBkgdErr->SetMarkerStyle(21);
  histBkgdErr->SetMarkerSize(0);
  histBkgdErr->SetFillStyle(3444);
  histBkgdErr->SetFillColor(1);
  histBkgdErr->SetLineWidth(2);
  histBkgdErr->SetLineColor(4);  

  histBkgd->GetXaxis()->SetLabelSize(0.06);
  histBkgd->GetXaxis()->SetLabelOffset(0.01);
  histBkgd->GetYaxis()->SetTitleOffset(0.77);
  histBkgd->GetYaxis()->SetTitle("Events / bin");
  histBkgd->GetYaxis()->SetTitleSize(0.06);
  histBkgd->GetYaxis()->SetLabelSize(0.06);
  histBkgd->GetYaxis()->SetTickLength(0.025);
  histBkgd->GetXaxis()->SetTickLength(0.025);
  histBkgd->GetYaxis()->SetTickSize(0.02);
  histBkgd->GetXaxis()->SetTickSize(0.02);

  histBkgd->GetYaxis()->SetRangeUser(0,1.2*histBkgd->GetMaximum());
  if (logY) histBkgd->GetYaxis()->SetRangeUser(0.5*histNonQCD->GetMinimum(),10.*histBkgd->GetMaximum());
  histBkgd->SetLineColor(4);
  histBkgd->SetLineWidth(2);
  histBkgd->SetLineStyle(1);
  
  TH1D * ratio = (TH1D*)histData->Clone("ratio");
  TH1D * ratioErr = (TH1D*)histData->Clone("ratioErr");
  ratioErr->SetMarkerStyle(21);
  ratioErr->SetMarkerSize(0);
  ratioErr->SetFillStyle(3444);
  ratioErr->SetFillColor(1);
  ratioErr->SetMarkerStyle(21);
  ratioErr->SetMarkerSize(0);  

  double x[21];
  double ex[21];
  double y[21];
  double eyhigh[21];
  double eylow[21];
  double yR[21];
  double erlow[21];
  double erhigh[21];

  int nBinsX = histBkgdErr->GetNbinsX();
  for (int iB=1; iB<=nBinsX; ++iB) {
    histBkgd->SetBinContent(iB,histBkgd->GetBinContent(iB));
    histBkgd->SetBinError(iB,0);
    histData->SetBinError(iB,0);
    double xD = histData->GetBinContent(iB);
    double eD = histData->GetBinError(iB);
    double xB = histBkgdErr->GetBinContent(iB);
    double eB = histBkgdErr->GetBinError(iB);
    double r = 1000;
    double er = 0;
    double erB = eB/xB;
    double xNonQCD = histNonQCD->GetBinContent(iB);
    if (xB>100.)
      printf("Bin %2i -> D = %5.0f  :  B = %6.0f +/- %4.0f  :  non-QCD = %5.1f\n",iB,xD,xB,eB,xNonQCD);
    else
      printf("Bin %2i -> D = %5.0f  :  B = %6.1f +/- %4.1f  :  non-QCD = %5.1f\n",iB,xD,xB,eB,xNonQCD);
    r = xD/xB;
    er = eD/xB;
    if (r<0.001) r = 0.001;
    if (erB>0.9) erB = 0.9;
    ratio->SetBinContent(iB,r);
    ratio->SetBinError(iB,0.);
    ratioErr->SetBinContent(iB,1);
    ratioErr->SetBinError(iB,erB);
    x[iB-1] = iB - 0.5;
    ex[iB-1] = 0.;
    y[iB-1] = histData->GetBinContent(iB);
    eylow[iB-1] = -0.5 + TMath::Sqrt(histData->GetBinContent(iB)+0.25);
    eyhigh[iB-1] = 0.5 + TMath::Sqrt(histData->GetBinContent(iB)+0.25);
    yR[iB-1] = xD/xB;
    erlow[iB-1] = eylow[iB-1]/xB;
    erhigh[iB-1] = eyhigh[iB-1]/xB;
  }
  
  double normNonQCD = histNonQCD->GetSumOfWeights();
  double normTot = histBkgd->GetSumOfWeights();
  double fractionNonQCD = 100.*normNonQCD/normTot;
  std::cout << std::endl;
  printf("Total background = %4.0f   non-QCD = %3.0f   fraction = %4.1f%%\n",normTot,normNonQCD,fractionNonQCD);
  std::cout << std::endl;

  TGraphAsymmErrors * data = new TGraphAsymmErrors(21,x,y,ex,ex,eylow,eyhigh);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.3);
  data->SetMarkerColor(1);
  data->SetLineWidth(2);

  TGraphAsymmErrors * ratioERR = new TGraphAsymmErrors(21,x,yR,ex,ex,erlow,erhigh);
  ratioERR->SetMarkerStyle(20);
  ratioERR->SetMarkerSize(1.3);
  ratioERR->SetMarkerColor(1);
  ratioERR->SetLineWidth(2);

  InitHist(histNonQCD,"","",kCyan,1001);
  
  int binNumber = 1;
  for (int i=1; i<=nBins; ++i) {
    for (int j=i; j<=nBins;++j) {
      char charLabel[10];
      sprintf(charLabel,"(%1i,%1i)",i,j);
      TString label(charLabel);
      if (binNumber<=nBins1D) {
	histBkgd->GetXaxis()->SetBinLabel(binNumber,"");
	ratioErr->GetXaxis()->SetBinLabel(binNumber,label);
	ratio->GetXaxis()->SetBinLabel(binNumber,label);
      }
      binNumber++;
    }
  }

  ratio->GetXaxis()->LabelsOption("v");
  ratioErr->GetXaxis()->LabelsOption("v");

  TCanvas * canv = MakeCanvas("canv","",750,600);
  TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.09);
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

  histBkgd->Draw("h");
  histNonQCD->Draw("sameh");
  histBkgdErr->Draw("e2same");
  histBkgd->Draw("hsame");
  data->Draw("pe1");

  TLegend * leg = new TLegend(0.65,0.66,0.85,0.9);
  SetLegendStyle(leg);
  leg->SetTextSize(0.055);
  leg->AddEntry(data,"observed","lp");
  leg->AddEntry(histBkgdErr,"bkg(+unc) (MC)","lf");
  leg->AddEntry(histNonQCD,"non-QCD (MC)","f");
  leg->Draw();

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextAngle(0);
  latexLabel.SetTextColor(kBlack);
  latexLabel.SetTextSize(0.075);
  latexLabel.DrawLatex(0.45,0.77,legend);
  
  writeExtraText = true;
  CMS_lumi(upper,4,0,-0.03); 

  float xLine = float(nBins);
  for (int i=1; i<nBins; ++i) {
    TLine * line = new TLine(xLine,0.2,xLine,1000);
    line->SetLineWidth(1);
    line->SetLineStyle(3);
    line->Draw();
    xLine += nBins - i; 
  }

  if (logY) upper->SetLogy(true);

  upper->Draw("SAME");
  upper->RedrawAxis();
  upper->Modified();
  upper->Update();
  canv->cd();

  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.09);
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

  ratioErr->GetYaxis()->SetRangeUser(yminRatio,ymaxRatio);
  ratioErr->GetXaxis()->SetNdivisions(410);
  ratioErr->GetYaxis()->SetNdivisions(505);
  ratioErr->GetXaxis()->SetLabelFont(42);
  ratioErr->GetXaxis()->SetLabelOffset(0.03);
  ratioErr->GetXaxis()->SetLabelSize(0.25);
  ratioErr->GetXaxis()->SetTitleSize(0.13);
  ratioErr->GetXaxis()->SetTitleOffset(1.2);
  ratioErr->GetYaxis()->SetTitle("obs/bkg");
  ratioErr->GetYaxis()->SetLabelFont(42);
  ratioErr->GetYaxis()->SetLabelOffset(0.015);
  ratioErr->GetYaxis()->SetLabelSize(0.14);
  ratioErr->GetYaxis()->SetTitleSize(0.14);
  ratioErr->GetYaxis()->SetTitleOffset(0.3);
  ratioErr->GetXaxis()->SetTickLength(0.025);
  ratioErr->GetYaxis()->SetTickLength(0.025);
  ratioErr->GetXaxis()->SetTickSize(0.02);
  ratioErr->GetYaxis()->SetTickSize(0.02);
  ratioErr->GetYaxis()->SetLabelOffset(0.01);
   
  ratioErr->Draw("e2");
  ratioERR->Draw("pe1");
  ratioErr->GetXaxis()->LabelsOption("v");

  TLine * line = new TLine(0,1,nBins1D,1);
  line->SetLineColor(4);
  line->SetLineWidth(2);
  line->Draw();
  xLine = float(nBins1D);
  
  lower->Modified();
  lower->RedrawAxis();
  xLine = float(nBins);
  for (int i=1; i<nBins; ++i) {
    TLine * line = new TLine(xLine,yminRatio,xLine,ymaxRatio);
    line->SetLineWidth(1);
    line->SetLineStyle(3);
    line->Draw();
    xLine += nBins - i; 
  }

  ratioERR->Draw("pe1");
  //  ratio->Draw("pe1same");
  canv->cd();
  canv->SetSelected(canv);

  std::cout << std::endl;
  canv->Print("figures/mass2D_"+era+"_"+label+".png");
  std::cout << std::endl;


}
