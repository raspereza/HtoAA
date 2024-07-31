#include "HttStylesNew.cc"

void HtoAA_Comparison() {

  SetStyle();

  // expected limits
  std::map<int, double> sigma_4t_ma = {
    {4,0.068},
    {5,0.046},
    {6,0.038},
    {7,0.026},
    {8,0.018},
    {9,0.014},
    {10,0.011},
    {11,0.011},
    {12,0.012},
    {13,0.013},
    {14,0.022},
    {15,0.043}
  };

  double m_tau = 1.777;
  double m_mu = 0.106;

  TH1D * hist = new TH1D("hist","",11,4.5,15.5);
  //ma  | sigma
  // 4  |  1.05

  std::cout << " ma | sigma " << std::endl;
  std::cout << "----+-------" << std::endl;
  for (auto p: sigma_4t_ma) {        
    double ma = p.first;
    double sigma_4t = p.second;
    double denominator = m_tau*m_tau*TMath::Sqrt(1-4*m_tau*m_tau/(ma*ma));
    double numerator = m_mu*m_mu*TMath::Sqrt(1-4*m_mu*m_mu/(ma*ma));
    double ratio = 2.*numerator/denominator;
    double sigma_mmtt = sigma_4t * ratio;
    printf(" %2i |  %4.2f\n",p.first,1000*sigma_mmtt);
    int bin = hist->FindBin(ma);
    hist->SetBinContent(bin,sigma_mmtt);
  }

  InitSignal(hist);
  hist->SetLineColor(4);
  hist->SetMarkerColor(4);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->SetLineStyle(1);
  hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitle("m_{a_{1}} [GeV]");
  hist->GetYaxis()->SetTitle("#sigma(pp#rightarrowH)#timesB(H#rightarrowa_{1}a_{1}#rightarrow2#mu2#tau)/#sigma_{SM}");
  hist->GetYaxis()->SetNdivisions(006);
  hist->GetXaxis()->SetNdivisions(020);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetRangeUser(0.,0.5e-3);


  TCanvas * canv = MakeCanvas("canv","",700,600);
  hist->Draw("lp");
  canv->SetGridx(true);
  canv->SetGridy(true);
  canv->RedrawAxis();
  canv->Update();
  canv->Print("SUS-23-005_comp.png");

}
