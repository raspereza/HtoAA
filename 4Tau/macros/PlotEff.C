#include "HttStylesNew.cc"
#include "CMS_lumi.C"
// EtaLt0p9
// Eta0p9to1p2
// Eta1p2to2p1
// EtaGt2p1
void PlotEff(TString era = "2016",
	     TString leg = "highPt",
	     TString EtaBin = "EtaGt2p1") {

  int nMax = 15;
 
  SetStyle();

  map<TString, TString> legnameMap = {
    {"2016_highPt","Mu17"},
    {"2016_lowPt","Mu8"},
    {"2017_highPt","Mu17"},
    {"2017_lowPt","Mu8"},
    {"2018_highPt","Mu18"},
    {"2018_lowPt","Mu9"}
  };

  map<TString, TString> filenameMap = {
    //    {"2016_highPt","HLT_Mu17Mu8_Muon17"},
    //    {"2016_lowPt","HLT_Mu17Mu8_Muon8"},
    {"2016_highPt","DoubleMuon_Mu17_UL_Iso"},
    {"2016_lowPt","DoubleMuon_Mu8_UL_Iso"},
    {"2017_highPt","HLT_Mu17Mu8_2017_Mu17_coarse"},
    {"2017_lowPt","HLT_Mu17Mu8_2017_Mu8_coarse"},
    {"2018_highPt","HLT_DoubleMu_Mu18leg"},
    {"2018_lowPt","HLT_DoubleMu_Mu9leg"}
  };

  map<TString, TString> etaBinMap = {
    {"EtaLt0p9","|#eta|<0.9"},
    {"Eta0p9to1p2","0.9<|#eta|<1.2"},
    {"Eta1p2to2p1","1.2<|#eta|<2.1"},
    {"EtaGt2p1","2.1<|#eta|<2.4"}
  };
  
  TString codeName = era + "_" + leg;
  TString Header = legnameMap[codeName]+", "+etaBinMap[EtaBin];

  TString FileName = filenameMap[codeName];

  double ptcut = 10;
  if (leg=="highPt")
    ptcut = 19;

  TString cmsswBase = TString(getenv("CMSSW_BASE"));


  TString folder = cmsswBase + "/src/HtoAA/data";
  TFile * file = new TFile(folder+"/"+FileName+".root");

  TString histNameData = "ZMass"+EtaBin+"_Data";
  TString histNameMC = "ZMass"+EtaBin+"_MC";

  TGraphAsymmErrors * effData = (TGraphAsymmErrors*)file->Get(histNameData);
  TGraphAsymmErrors * effMC = (TGraphAsymmErrors*)file->Get(histNameMC);

  effData->SetLineColor(2);
  effData->SetMarkerColor(2);
  effData->SetMarkerSize(1.5);
  effData->SetMarkerStyle(20);

  effMC->SetLineColor(4);
  effMC->SetMarkerColor(4);
  effMC->SetMarkerSize(1.5);
  effMC->SetMarkerStyle(21);

  int nBins = effData->GetN();
  float bins[100];
  std::cout << "nBins = " << nBins << std::endl;
  std::cout << "efficiencies -> " << std::endl;
  for (int iB=0; iB<nBins; ++iB) {
    bins[iB] = effData->GetX()[iB] - effData->GetErrorX(iB);
    float lower = bins[iB];
    float upper = effData->GetX()[iB] + effData->GetErrorX(iB);
    printf("[%2i,%3i] : %5.3f - %5.3f\n",int(lower),int(upper),effData->GetY()[iB],effMC->GetY()[iB]);

  }
  bins[nBins] = effData->GetX()[nBins-1] + effData->GetErrorX(nBins-1);
  std::cout << bins[nBins] << std::endl;

  TH1D * SF = new TH1D("SF","",nBins,bins);
  SF->GetXaxis()->SetTitle("muon p_{T} [GeV]");
  SF->SetMarkerColor(1);
  SF->SetMarkerStyle(20);
  SF->SetMarkerSize(1.2);
  SF->SetLineColor(1);
  SF->GetYaxis()->SetRangeUser(0.801,1.199);
  SF->GetYaxis()->SetNdivisions(505);
  SF->GetXaxis()->SetLabelFont(42);
  SF->GetXaxis()->SetLabelOffset(0.04);
  SF->GetXaxis()->SetLabelSize(0.14);
  SF->GetXaxis()->SetTitleSize(0.13);
  SF->GetXaxis()->SetTitleOffset(1.2);
  SF->GetYaxis()->SetTitle("obs/exp");
  SF->GetYaxis()->SetLabelFont(42);
  SF->GetYaxis()->SetLabelOffset(0.015);
  SF->GetYaxis()->SetLabelSize(0.13);
  SF->GetYaxis()->SetTitleSize(0.14);
  SF->GetYaxis()->SetTitleOffset(0.5);
  SF->GetXaxis()->SetTickLength(0.07);
  SF->GetYaxis()->SetTickLength(0.04);
  SF->GetYaxis()->SetLabelOffset(0.01);

  SF->GetXaxis()->SetRangeUser(6.001,99.9999);
    
  std::cout << std::endl;
  std::cout << "Scale factors -> " << std::endl;
  for (int iB=0; iB<nBins; ++iB) {
    float num = effData->GetY()[iB];
    float eNum = effData->GetErrorY(iB);
    float den = effMC->GetY()[iB];    
    float sf = 0;
    float sfE = 0; 
    if (den>0) { 
      sf = num/den; 
      sfE = eNum/den;
    }
    SF->SetBinContent(iB+1,sf);
    SF->SetBinError(iB+1,sfE);
    float lower = SF->GetBinLowEdge(iB+1);
    float upper = SF->GetBinLowEdge(iB+2);
    printf("[%2i,%3i] : %5.3f +/- %5.3f\n",int(lower),int(upper),sf,sfE);
  }
  std::cout << std::endl;

  TH2F * frame = new TH2F("frame","",2,6,100,2,0.,1.2);
  frame->GetXaxis()->SetTitle();
  frame->GetYaxis()->SetTitle("efficiency");
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetXaxis()->SetLabelSize(0);

  TCanvas * canv = MakeCanvas("canv","",600,700);
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

  //  upper->SetGridx();
  //  upper->SetGridy();

  frame->Draw();
  effMC->Draw("PEsame");
  effData->Draw("PEsame");
  TLegend * legend = new TLegend(0.4,0.17,0.94,0.40);
  SetLegendStyle(legend);
  legend->SetHeader(Header);
  legend->SetTextSize(0.06);
  legend->AddEntry(effData,"Data","lp");
  legend->AddEntry(effMC,"MC (Z#rightarrow#mu#mu)","lp");
  legend->Draw();
  TLine * line = new TLine(ptcut,0,ptcut,1.1);
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->Draw();
  lumi_13TeV = "2016, 36.3 fb^{-1}";
  if (era=="2017") lumi_13TeV = "2017, 41.5 fb^{-1}";
  if (era=="2018") lumi_13TeV = "2018, 59.8 fb^{-1}";
  writeExtraText = true;
  extraText   = "Preliminary";
  CMS_lumi(upper,4,33); 


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
  SF->Draw("e1");
  TLine * lineSF = new TLine(ptcut,0.801,ptcut,1.199);
  lineSF->SetLineColor(2);
  lineSF->SetLineStyle(2);
  lineSF->Draw();
  lower->Modified();
  lower->RedrawAxis();
  canv->cd();
  canv->SetSelected(canv);
  canv->Update();
  canv->Print("TriggerEff_"+leg+"_"+EtaBin+"_"+era+".png");
  //  canv->Print(FileName+EtaBin+"_eff.pdf");
  //  canv->Print(FileName+EtaBin+"_eff.pdf","Portrait pdf");
  //  canv->Print(FileName+EtaBin+"_eff.C");
  //  TFile * canvasFile = new TFile(FileName+"_"+EtaBin+"_eff.root","recreate");
  //  canvasFile->cd();
  //  canv->Write("canvas");
  //  canvasFile->Close();
  //  canv->Update();


}
