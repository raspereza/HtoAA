#include "HttStylesNew.cc"

void Compatibility(
		   TString folder = "GoF_Run2_bonly", // folder with RooT files
		   TString Algo = "saturated", // algorithm
		   TString legend = "H(125)#rightarrowaa#rightarrow4#tau (Run2)",
		   int bins = 50 // number of bins in the histogram of toys
		   ) {

  SetStyle();
  gStyle->SetOptStat(0);

  TFile * fileObs = new TFile(folder+"/gof_obs.root");
  if (fileObs==NULL||fileObs->IsZombie()) {
    TString fullpath = folder+"/gof_obs.root";
    std::cout << "file does not exist : " << fullpath << std::endl;
    std::cout << "Make sure that you have run GoF test and saved output in the folder " << std::endl; 
    std::cout << folder << std::endl;
    return;
  }
  double obs;
  float xMax = -1;
  float xMin = 1e+10;
  TTree * treeObs = (TTree*)fileObs->Get("limit");
  treeObs->SetBranchAddress("limit",&obs);
  treeObs->GetEntry(0);
  if (obs<xMin) xMin = obs;
  if (obs>xMax) xMax = obs;

  TFile * file = new TFile(folder+"/gof_exp.root");  
  if (file==NULL||file->IsZombie()) {
    TString fullpath = folder+"/gof_exp.root";
    std::cout << "file does not exist : " << fullpath << std::endl;
    std::cout << "Make sure that you have run GoF test and saved output in the folder " << std::endl; 
    std::cout << folder << std::endl;
    return;
  }

  double limit;

  TTree * tree = (TTree*)file->Get("limit");
  tree->SetBranchAddress("limit",&limit);


  int entries = tree->GetEntries();
  float count = 0;
  float Entries = 0;
  for (int i=0; i<entries; ++i) {
    tree->GetEntry(i);
    if (limit>xMax) xMax = limit;
    if (limit<xMin) xMin = limit;
    Entries += 1.0;
    if (limit>obs) {
      count += 1.0;
    }
  }

  float xlower = 0.8*xMin;
  float xupper = 1.2*xMax;

  TH1F * chi2 = new TH1F("chi2",legend,bins,xlower,xupper);
  InitSignal(chi2);
  chi2->SetLineStyle(1);
  chi2->SetLineWidth(2);

  for (int i=0; i<entries; ++i) {
    tree->GetEntry(i);
    chi2->Fill(float(limit));
  }

  std::cout << "Minimal value of q        = " << xMin << std::endl;
  std::cout << "Maximal value of q        = " << xMax << std::endl;
  std::cout << "Observed value of q       = " << obs << std::endl;
  std::cout << "Number of toys with q>obs = " << count << std::endl;
  std::cout << "Total number of toys      = " << entries << std::endl;
  float prob = count / float(Entries);
  std::cout << "p-value                   = " << prob << std::endl;
  char probChar[5];
  sprintf(probChar,"%4.2f",prob);
  TString probStr = TString("p-value = ") + TString(probChar);

  float yMax = 0.6*chi2->GetMaximum();
  float yArrow = 0.05*chi2->GetMaximum();

  chi2->SetNdivisions(505,"X");

  TCanvas * canv = new TCanvas("canv","",700,600);
  if (Algo.Contains("saturated"))
    chi2->GetXaxis()->SetTitle("#chi^{2}");
  if (Algo.Contains("KS"))
    chi2->GetXaxis()->SetTitle("q_{KS}");
  if (Algo.Contains("AD"))
    chi2->GetXaxis()->SetTitle("q_{AD}");
  chi2->GetYaxis()->SetTitle("toys");
  chi2->GetYaxis()->SetTitleOffset(1.2);
  chi2->GetXaxis()->SetTitleOffset(1.2);
  chi2->Draw();
  float arrW = 0.02*(chi2->GetBinLowEdge(bins+1)-chi2->GetBinLowEdge(1));
  TLine * line = new TLine(obs,0,obs,yMax);
  line->SetLineWidth(3);
  line->SetLineColor(4);
  line->Draw();
  TLine * line1 = new TLine(obs,0,obs-arrW,yArrow);
  line1->SetLineWidth(3);
  line1->SetLineColor(4);
  line1->Draw();
  TLine * line2 = new TLine(obs,0,obs+arrW,yArrow);
  line2->SetLineWidth(3);
  line2->SetLineColor(4);
  line2->Draw();  
  TLatex * tex = new TLatex(0.62,0.82,probStr);
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetLineWidth(2);
  tex->Draw();

  canv->Update();
  canv->Print(folder+"_"+Algo+".png");

}
