#include "HttStylesNew.cc"
void CreateBTagEff(TString era = "2018") {

  SetStyle();
  gStyle->SetOptStat(0000);

  TString baseFolder = "/nfs/dust/cms/user/rasp/BTagReshaping";
  TString fileName = baseFolder + "/" + era+"/TTToHadronic.root";

  TFile * file = new TFile(fileName);
  std::vector<TString> flavors = {"b","c","oth"};

  std::vector<TH2D*> effv; effv.clear();
  for (auto flavor : flavors) {
    TH2D * num = (TH2D*)file->Get(flavor+"_pass");
    TH2D * den = (TH2D*)file->Get(flavor+"_all");
    TH2D * eff = (TH2D*)num->Clone("btag_eff_"+flavor);
    int nBinsX = eff->GetNbinsX();
    int nBinsY = eff->GetNbinsY();
    for (int xB=1; xB<=nBinsX; ++xB) {
      for (int yB=1; yB<=nBinsY; ++yB) {
	double xnum = num->GetBinContent(xB,yB);
	double xden = den->GetBinContent(xB,yB);
	double xeff = 0.0;
	if (xden>0) xeff = xnum/xden;
	eff->SetBinContent(xB,yB,xeff);
      }
    }
    effv.push_back(eff);
  }

  TString output = "beff_Tight_" + era + ".root"; 
  TFile * outfile = new TFile(output,"recreate");
  outfile->cd("");
  for (auto eff : effv) {
    TString name = eff->GetName();
    eff->SetTitle(name+" "+"  CMS simulation "+era);
    eff->GetXaxis()->SetTitle("jet p_{T} (GeV)");
    eff->GetYaxis()->SetTitle("jet |#eta|");
    eff->Write(name);
    TCanvas * canv = MakeCanvas(name+"canv","",700,600);
    eff->Draw("colz");
  }
  outfile->Close();
  


}
