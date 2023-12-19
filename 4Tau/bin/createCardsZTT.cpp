#include "HtoH.h"
void CreateCards(TString mass="5", // mass of pseudoscalar
		 bool Azimov = true, // replace data by background expectations 
		 bool correlation = true // apply correlations
                 ) {

  TString dir = "./";
  
  double massD = 4.0;
  if (mass=="4")
    massD = 4.0;
  if (mass=="5")
    massD = 5.0;
  if (mass=="6")
    massD = 6.0;
  if (mass=="7")
    massD = 7.0;
  if (mass=="8")
    massD = 8.0;
  if (mass=="9")
    massD = 9.0;
  if (mass=="10")
    massD = 10.0;
  if (mass=="11")
    massD = 11.0;
  if (mass=="12")
    massD = 12.0;
  if (mass=="13")
    massD = 13.0;
  if (mass=="14")
    massD = 14.0;
  if (mass=="15")
    massD = 15.0;

  bool crude = false;

  std::cout << std::endl;
  std::cout << "Mass = " << mass << std::endl;
  std::cout << std::endl;

  double lumi = 59740;
  // Cross sections in pb are given at mH = 125.09 GeV
  // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
  double xsecGGH = 48.61;
  double xsecVBF = 3.766;
  double xsecVH  = 1.358+0.880;
  double xsecTTH = 0.5065;
  
  double massTau = 1.777;
  double massMu  = 0.106;
  double massRatio = (massMu*massMu)/(massTau*massTau);
  double aF = 2*massTau/massD;
  double SF = 2*massRatio/TMath::Sqrt(1-aF*aF);

  // contribution of H->aa->(2mu)(2tau)
  double xsecMMTT = (xsecGGH+xsecVBF+xsecVH+xsecTTH) * SF;

  TFile * file     = new TFile(dir+"/DoubleMuon_Run2018.root");
  TFile * fileGGH  = new TFile(dir+"/SUSYGluGluToHToAA_AToTauTau_M-125_M-"+mass+".root");
  TFile * fileVBF  = new TFile(dir+"/SUSYVBFToHToAA_AToTauTau_M-125_M-"+mass+".root");
  TFile * fileVH   = new TFile(dir+"/SUSYVH_HToAA_AToTauTau_M-125_M-"+mass+".root");
  TFile * fileTTH  = new TFile(dir+"/SUSYttH_HToAA_AToTauTau_M-125_M-"+mass+".root");
  TFile * fileMMTT = new TFile(dir+"/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-"+mass+".root");

  TH2D * histOld    = (TH2D*)file->Get("InvMass2DH");

  TH2D * histGGHOld  = (TH2D*)fileGGH->Get("InvMass2DH");
  TH2D * histVBFOld  = (TH2D*)fileVBF->Get("InvMass2DH");
  TH2D * histVHOld   = (TH2D*)fileVH ->Get("InvMass2DH");
  TH2D * histTTHOld  = (TH2D*)fileTTH->Get("InvMass2DH");
  TH2D * histMMTTOld = (TH2D*)fileMMTT->Get("InvMass2DH");

  TH2D * histGGH_btagUp_Old = (TH2D*)fileGGH->Get("InvMass2DH_btagUp");
  TH2D * histGGH_btagDown_Old = (TH2D*)fileGGH->Get("InvMass2DH_btagDown");
  TH2D * histVBF_btagUp_Old = (TH2D*)fileVBF->Get("InvMass2DH_btagUp");
  TH2D * histVBF_btagDown_Old = (TH2D*)fileVBF->Get("InvMass2DH_btagDown");
  TH2D * histVH_btagUp_Old = (TH2D*)fileVH->Get("InvMass2DH_btagUp");
  TH2D * histVH_btagDown_Old = (TH2D*)fileVH->Get("InvMass2DH_btagDown");
  TH2D * histTTH_btagUp_Old = (TH2D*)fileTTH->Get("InvMass2DH_btagUp");
  TH2D * histTTH_btagDown_Old = (TH2D*)fileTTH->Get("InvMass2DH_btagDown");
  TH2D * histMMTT_btagUp_Old = (TH2D*)fileMMTT->Get("InvMass2DH_btagUp");
  TH2D * histMMTT_btagDown_Old = (TH2D*)fileMMTT->Get("InvMass2DH_btagDown");

  TH1D * histWeightsGGH = (TH1D*)fileGGH->Get("histWeightsH");
  double nGenGGH = histWeightsGGH->GetSumOfWeights();
  
  TH1D * histWeightsVBF = (TH1D*)fileVBF->Get("histWeightsH");
  double nGenVBF = histWeightsVBF->GetSumOfWeights();
  
  TH1D * histWeightsVH = (TH1D*)fileVH->Get("histWeightsH");
  double nGenVH = histWeightsVH->GetSumOfWeights();
  
  TH1D * histWeightsTTH = (TH1D*)fileTTH->Get("histWeightsH");
  double nGenTTH = histWeightsTTH->GetSumOfWeights();

  TH1D * histWeightsMMTT = (TH1D*)fileMMTT->Get("histWeightsH");
  double nGenMMTT = histWeightsMMTT->GetSumOfWeights();
  
  double gghNorm = xsecGGH*lumi/nGenGGH;
  double vbfNorm = xsecVBF*lumi/nGenVBF;
  double vhNorm  = xsecVH*lumi/nGenVH;
  double tthNorm = xsecTTH*lumi/nGenTTH;
  double mmttNorm = xsecMMTT*lumi/nGenMMTT;

  TFile * fileCorr = new TFile(dir+"/CorrCoefficients_data.root");
  TH2D * corrCoeff = (TH2D*)fileCorr->Get("corrCoeff");

  TFile * fileCorrCR = new TFile(dir+"/CorrCoefficients_control_mc.root");
  TH2D * corrCoeffCR = (TH2D*)fileCorrCR->Get("corrCoeff");

  TFile * fileCorrSR = new TFile(dir+"/CorrCoefficients_signal_mc.root");
  TH2D * corrCoeffSR = (TH2D*)fileCorrSR->Get("corrCoeff");

  TH1D * hist1dN23Old  = (TH1D*)file->Get("InvMassN23H");
  TH1D * hist1dN45Old  = (TH1D*)file->Get("InvMassN45H");

  int nBinsNew = 6;
  double bins[7] = {0,1,2,3,4,6,20};

  TH1D * hist1d    = (TH1D*)TH1DtoTH1D(hist1dN23Old,nBinsNew,bins,true,"_new");
  TH1D * hist1dN45 = (TH1D*)TH1DtoTH1D(hist1dN45Old,nBinsNew,bins,true,"_new");
  hist1d->Scale(1.0/hist1d->GetSumOfWeights());
  hist1dN45->Scale(1.0/hist1dN45->GetSumOfWeights());

  TH1D * uncertMass1DH = new TH1D("uncertMass1DH","",nBinsNew,bins);
  for (int iB=1; iB<=nBinsNew; ++iB) {
    double num = hist1d->GetBinContent(iB);
    double den = hist1dN45->GetBinContent(iB);
    double ratio = num/den;
    uncertMass1DH->SetBinContent(iB,ratio);
  }

  TH2D * histGGH = (TH2D*)TH2DtoTH2D(histGGHOld,nBinsNew,bins,nBinsNew,bins,"_gghNew");
  TH2D * histVBF = (TH2D*)TH2DtoTH2D(histVBFOld,nBinsNew,bins,nBinsNew,bins,"_vbfNew");
  TH2D * histVH  = (TH2D*)TH2DtoTH2D(histVHOld,nBinsNew,bins,nBinsNew,bins,"_vhNew");
  TH2D * histTTH = (TH2D*)TH2DtoTH2D(histTTHOld,nBinsNew,bins,nBinsNew,bins,"_tthNew");
  TH2D * histMMTT = (TH2D*)TH2DtoTH2D(histMMTTOld,nBinsNew,bins,nBinsNew,bins,"_mmttNew");

  TH2D * histGGH_btagUp = (TH2D*)TH2DtoTH2D(histGGH_btagUp_Old,nBinsNew,bins,nBinsNew,bins,"_ggh_btagUp_New");
  TH2D * histGGH_btagDown = (TH2D*)TH2DtoTH2D(histGGH_btagDown_Old,nBinsNew,bins,nBinsNew,bins,"_ggh_btagDown_New");
  TH2D * histVBF_btagUp = (TH2D*)TH2DtoTH2D(histVBF_btagUp_Old,nBinsNew,bins,nBinsNew,bins,"_vbf_btagUp_New");
  TH2D * histVBF_btagDown = (TH2D*)TH2DtoTH2D(histVBF_btagDown_Old,nBinsNew,bins,nBinsNew,bins,"_vbf_btagDown_New");
  TH2D * histVH_btagUp = (TH2D*)TH2DtoTH2D(histVH_btagUp_Old,nBinsNew,bins,nBinsNew,bins,"_vh_btagUp_New");
  TH2D * histVH_btagDown = (TH2D*)TH2DtoTH2D(histVH_btagDown_Old,nBinsNew,bins,nBinsNew,bins,"_vh_btagDown_New");
  TH2D * histTTH_btagUp = (TH2D*)TH2DtoTH2D(histTTH_btagUp_Old,nBinsNew,bins,nBinsNew,bins,"_tth_btagUp_New");
  TH2D * histTTH_btagDown = (TH2D*)TH2DtoTH2D(histTTH_btagDown_Old,nBinsNew,bins,nBinsNew,bins,"_tth_btagDown_New");
  TH2D * histMMTT_btagUp = (TH2D*)TH2DtoTH2D(histMMTT_btagUp_Old,nBinsNew,bins,nBinsNew,bins,"_mmtt_btagUp_New");
  TH2D * histMMTT_btagDown = (TH2D*)TH2DtoTH2D(histMMTT_btagDown_Old,nBinsNew,bins,nBinsNew,bins,"_mmtt_btagDown_New");

  TH2D * histData = (TH2D*)TH2DtoTH2D(histOld,nBinsNew,bins,nBinsNew,bins,"_dataNew");

  double bkgNorm = histData->GetSumOfWeights();
  double sideBandNorm = 1/hist1d->GetSumOfWeights();

  std::cout << "Bkg  Norm = " << bkgNorm << std::endl;
  std::cout << "Signal expectations for BR(H->aa->4tau)=1.0" << std::endl;
  std::cout << "ggH  Norm = " << gghNorm*histGGH->GetSumOfWeights() << std::endl;
  std::cout << "VBF  Norm = " << vbfNorm*histVBF->GetSumOfWeights() << std::endl;
  std::cout << "VH   Norm = " << vhNorm*histVH->GetSumOfWeights() << std::endl;
  std::cout << "TTH  Norm = " << tthNorm*histTTH->GetSumOfWeights() << std::endl;
  std::cout << "MMTT Norm = " << mmttNorm*histMMTT->GetSumOfWeights() << std::endl;
  std::cout << std::endl;

  int nBins1D = nBinsNew * (nBinsNew + 1) / 2;

  // unrolled 1D signal templates
  TH1D * ggh = new TH1D("ggh","",nBins1D,0.,float(nBins1D));
  TH1D * vbf = new TH1D("vbf","",nBins1D,0.,float(nBins1D));
  TH1D * vh  = new TH1D("vh","",nBins1D,0.,float(nBins1D));
  TH1D * tth  = new TH1D("tth","",nBins1D,0.,float(nBins1D));
  TH1D * mmtt = new TH1D("mmtt","",nBins1D,0.,float(nBins1D));

  //unrolled 1D signal systematic templates
  TH1D * ggh_btagUp = new TH1D("ggh_btagUp","",nBins1D,0.,float(nBins1D));
  TH1D * ggh_btagDown = new TH1D("ggh_btagDown","",nBins1D,0.,float(nBins1D));
  TH1D * vbf_btagUp = new TH1D("vbf_btagUp","",nBins1D,0.,float(nBins1D));
  TH1D * vbf_btagDown = new TH1D("vbf_btagDown","",nBins1D,0.,float(nBins1D));
  TH1D * vh_btagUp = new TH1D("vh_btagUp","",nBins1D,0.,float(nBins1D));
  TH1D * vh_btagDown = new TH1D("vh_btagDown","",nBins1D,0.,float(nBins1D));  
  TH1D * tth_btagUp = new TH1D("tth_btagUp","",nBins1D,0.,float(nBins1D));
  TH1D * tth_btagDown = new TH1D("tth_btagDown","",nBins1D,0.,float(nBins1D));
  TH1D * mmtt_btagUp = new TH1D("mmtt_btagUp","",nBins1D,0.,float(nBins1D));
  TH1D * mmtt_btagDown = new TH1D("mmtt_btagDown","",nBins1D,0.,float(nBins1D));

  // central template and systematic background templates
  // 1d_up - up variation of 1D background pdf (mu-trk mass)
  // 2_down - down variation of 1D background pdf
  TH1D * bkgd[3];
  TString bkg_name[3] = {"","_1d_up","_1d_down"};
  for (int i = 0; i<3; ++i) {
    bkgd[i] = new TH1D("bkgd"+bkg_name[i],"",nBins1D,0.,float(nBins1D));
  }
  TH1D * bkgdCorrUp = new TH1D("bkgdCorrUp","",nBins1D,0.,float(nBins1D));
  TH1D * bkgdCorrDown = new TH1D("bkgdCorrDown","",nBins1D,0.,float(nBins1D));

  TH1D * data   = new TH1D("data","",nBins1D,0.,float(nBins1D));

  int iBin = 0;
  int dataAll = 0;
  double productAll = 0;
  // *******************************************
  // ********* Unrolling 2D templates **********
  // ***** constructing central bkg template ***
  // ***** and systematic bkgd templates ******* 
  // *******************************************
  for (int i=1; i<=nBinsNew; ++i) {
    for (int j=i; j<=nBinsNew; ++j) {

      iBin++;

      double xcorr = uncertMass1DH->GetBinContent(i);
      double ycorr = uncertMass1DH->GetBinContent(j);
      double xBkgd[3];
      double yBkgd[3];
      double exBkgd[3];
      double eyBkgd[3];
      xBkgd[0] = sideBandNorm*hist1d->GetBinContent(i);
      yBkgd[0] = sideBandNorm*hist1d->GetBinContent(j);
      exBkgd[0] = sideBandNorm*hist1d->GetBinError(i);
      eyBkgd[0] = sideBandNorm*hist1d->GetBinError(j);

      xBkgd[1] = xBkgd[0]*xcorr;
      yBkgd[1] = yBkgd[0]*ycorr;
      exBkgd[1] = exBkgd[0]*xcorr;
      eyBkgd[1] = eyBkgd[0]*ycorr;

      xBkgd[2] = xBkgd[0]/xcorr;
      yBkgd[2] = yBkgd[0]/ycorr;
      exBkgd[2] = exBkgd[0]/xcorr;
      eyBkgd[2] = eyBkgd[0]/ycorr;

      for (int itempl=0; itempl<3; ++itempl) {
	double product = xBkgd[itempl] * yBkgd[itempl];
	double err = exBkgd[itempl] * yBkgd[itempl] + eyBkgd[itempl] * xBkgd[itempl];

        double m1 = hist1d->GetBinCenter(i);
        double m2 = hist1d->GetBinCenter(j);

        int corrBin = corrCoeff->FindBin(m1, m2);
      
        double corrData   = corrCoeff->GetBinContent(corrBin);
        double corrDataE  = corrCoeff->GetBinError(corrBin);
        double corrDataR  = corrDataE/corrData; 
        double corrCR   = corrCoeffCR->GetBinContent(corrBin);
        double corrCRE  = corrCoeffCR->GetBinError(corrBin);
        double corrCRR  = corrCRE/corrCR;
        double corrSR   = corrCoeffSR->GetBinContent(corrBin);
        double corrSRE  = corrCoeffSR->GetBinError(corrBin);
        double corrSRR  = corrSRE/corrSR;
	
        double corrX = corrSR * corrData / corrCR;
	double corrUp = corrData;
	double corrDown = corrData * (corrSR * corrSR) / (corrCR * corrCR);
	double corrR = TMath::Sqrt(corrDataR*corrDataR+corrCRR*corrCRR+corrSRR*corrSRR);
	double corrE = corrX * corrR;

	double productUp = product;
	double productDown = product;
	if (i!=j) {
	  product *= 2.0;
	  err *= 2.0;
	  productUp *= 2.0;
	  productDown *= 2.0;
	}
	if (correlation) {
	  product *= corrX;
	  productUp *= corrUp;
	  productDown *= corrDown;
	  err *= corrX;
	}
	double errCor = product * corrE / corrX;

	err = TMath::Sqrt(err*err+errCor*errCor);
	bkgd[itempl]->SetBinContent(iBin,product);
	bkgd[itempl]->SetBinError(iBin,err);
	if (itempl==0) {
	  bkgdCorrUp->SetBinContent(iBin,productUp);
	  bkgdCorrDown->SetBinContent(iBin,productDown);
	  bkgdCorrUp->SetBinError(iBin,err);
	  bkgdCorrDown->SetBinError(iBin,err);
	  productAll += product;
	} 
      }
      double xGGH = histGGH->GetBinContent(i,j);
      double eGGH = histGGH->GetBinError(i,j);
      double xVBF = histVBF->GetBinContent(i,j);
      double eVBF = histVBF->GetBinError(i,j);
      double xVH  = histVH->GetBinContent(i,j);
      double eVH  = histVH->GetBinError(i,j);
      double xTTH = histTTH->GetBinContent(i,j);
      double eTTH = histTTH->GetBinError(i,j);
      double xMMTT = histMMTT->GetBinContent(i,j);
      double eMMTT = histMMTT->GetBinError(i,j);

      double xGGH_btagUp = histGGH_btagUp->GetBinContent(i,j);
      double eGGH_btagUp = histGGH_btagUp->GetBinError(i,j);
      double xGGH_btagDown = histGGH_btagDown->GetBinContent(i,j);
      double eGGH_btagDown = histGGH_btagDown->GetBinError(i,j);
      double xVBF_btagUp = histVBF_btagUp->GetBinContent(i,j);
      double eVBF_btagUp = histVBF_btagUp->GetBinError(i,j);
      double xVBF_btagDown = histVBF_btagDown->GetBinContent(i,j);
      double eVBF_btagDown = histVBF_btagDown->GetBinError(i,j);
      double xVH_btagUp = histVH_btagUp->GetBinContent(i,j);
      double eVH_btagUp = histVH_btagUp->GetBinError(i,j);
      double xVH_btagDown = histVH_btagDown->GetBinContent(i,j);
      double eVH_btagDown = histVH_btagDown->GetBinError(i,j);
      double xTTH_btagUp = histTTH_btagUp->GetBinContent(i,j);
      double eTTH_btagUp = histTTH_btagUp->GetBinError(i,j);
      double xTTH_btagDown = histTTH_btagDown->GetBinContent(i,j);
      double eTTH_btagDown = histTTH_btagDown->GetBinError(i,j);
      double xMMTT_btagUp = histMMTT_btagUp->GetBinContent(i,j);
      double eMMTT_btagUp = histMMTT_btagUp->GetBinError(i,j);
      double xMMTT_btagDown = histMMTT_btagDown->GetBinContent(i,j);
      double eMMTT_btagDown = histMMTT_btagDown->GetBinError(i,j);

      if (i!=j) {

	xGGH += histGGH->GetBinContent(j,i);
	double err1 = eGGH;
	double err2 = histGGH->GetBinError(j,i);
	eGGH = TMath::Sqrt(err1*err1+err2*err2);

        xGGH_btagUp += histGGH_btagUp->GetBinContent(j,i);
        err1 = eGGH_btagUp;
        err2 = histGGH_btagUp->GetBinError(j,i);
        eGGH_btagUp = TMath::Sqrt(err1*err1+err2*err2); 

        xGGH_btagDown += histGGH_btagDown->GetBinContent(j,i);
        err1 = eGGH_btagDown;
        err2 = histGGH_btagDown->GetBinError(j,i);
        eGGH_btagDown = TMath::Sqrt(err1*err1+err2*err2); 
	  
	xVBF += histVBF->GetBinContent(j,i);
	err1  = eVBF;
	err2  = histVBF->GetBinError(j,i);
	eVBF  = TMath::Sqrt(err1*err1+err2*err2);
	
        xVBF_btagUp += histVBF_btagUp->GetBinContent(j,i);
        err1 = eVBF_btagUp;
        err2 = histVBF_btagUp->GetBinError(j,i);
        eVBF_btagUp = TMath::Sqrt(err1*err1+err2*err2);

        xVBF_btagDown += histVBF_btagDown->GetBinContent(j,i);
        err1 = eVBF_btagDown;
        err2 = histVBF_btagDown->GetBinError(j,i);
        eVBF_btagDown = TMath::Sqrt(err1*err1+err2*err2);        

	xVH += histVH->GetBinContent(j,i);
	err1 = eVH;
	err2 = histVH->GetBinError(j,i);
	eVH = TMath::Sqrt(err1*err1+err2*err2);
        
        xVH_btagUp += histVH_btagUp->GetBinContent(j,i);
        err1 = eVH_btagUp;
        err2 = histVH_btagUp->GetBinError(j,i);
        eVH_btagUp = TMath::Sqrt(err1*err1+err2*err2);

        xVH_btagDown += histVH_btagDown->GetBinContent(j,i);
        err1 = eVH_btagDown;
        err2 = histVH_btagDown->GetBinError(j,i);
        eVH_btagDown = TMath::Sqrt(err1*err1+err2*err2);

	xTTH += histTTH->GetBinContent(j,i);
	err1 = eTTH;
	err2 = histTTH->GetBinError(j,i);
	eTTH = TMath::Sqrt(err1*err1+err2*err2);

        xTTH_btagUp += histTTH_btagUp->GetBinContent(j,i);
        err1 = eTTH_btagUp;
        err2 = histTTH_btagUp->GetBinError(j,i);
        eTTH_btagUp = TMath::Sqrt(err1*err1+err2*err2);

        xTTH_btagDown += histTTH_btagDown->GetBinContent(j,i);
        err1 = eTTH_btagDown;
        err2 = histTTH_btagDown->GetBinError(j,i);
        eTTH_btagDown = TMath::Sqrt(err1*err1+err2*err2);
	  
	xMMTT += histMMTT->GetBinContent(j,i);
	err1 = eMMTT;
	err2 = histMMTT->GetBinError(j,i);
	eMMTT = TMath::Sqrt(err1*err1+err2*err2);
        
        xMMTT_btagUp += histMMTT_btagUp->GetBinContent(j,i);
        err1 = eMMTT_btagUp;
        err2 = histMMTT_btagUp->GetBinError(j,i);
        eMMTT_btagUp = TMath::Sqrt(err1*err1+err2*err2);

        xMMTT_btagDown += histMMTT_btagDown->GetBinContent(j,i);
        err1 = eMMTT_btagDown;
        err2 = histMMTT_btagDown->GetBinError(j,i);
        eGGH_btagDown = TMath::Sqrt(err1*err1+err2*err2);	
  
      }

      ggh->SetBinContent(iBin,gghNorm*xGGH);
      ggh->SetBinError(iBin,gghNorm*eGGH);

      ggh_btagUp->SetBinContent(iBin,gghNorm*xGGH_btagUp);
      ggh_btagUp->SetBinError(iBin,gghNorm*eGGH_btagUp); 

      ggh_btagDown->SetBinContent(iBin,gghNorm*xGGH_btagDown);
      ggh_btagDown->SetBinError(iBin,gghNorm*eGGH_btagDown);

      vbf->SetBinContent(iBin,vbfNorm*xVBF);
      vbf->SetBinError(iBin,vbfNorm*eVBF);
     
      vbf_btagUp->SetBinContent(iBin,vbfNorm*xVBF_btagUp);
      vbf_btagUp->SetBinError(iBin,vbfNorm*eVBF_btagUp);

      vbf_btagDown->SetBinContent(iBin,vbfNorm*xVBF_btagDown);
      vbf_btagDown->SetBinError(iBin,vbfNorm*eVBF_btagDown);

      vh->SetBinContent(iBin,vhNorm*xVH);
      vh->SetBinError(iBin,vhNorm*eVH);

      vh_btagUp->SetBinContent(iBin,vhNorm*xVH_btagUp);
      vh_btagUp->SetBinError(iBin,vhNorm*eVH_btagUp);

      vh_btagDown->SetBinContent(iBin,vhNorm*xVH_btagDown);
      vh_btagDown->SetBinError(iBin,vhNorm*eVH_btagDown);

      tth->SetBinContent(iBin,tthNorm*xTTH);
      tth->SetBinError(iBin,tthNorm*eTTH);

      tth_btagUp->SetBinContent(iBin,tthNorm*xTTH_btagUp);
      tth_btagUp->SetBinError(iBin,tthNorm*eTTH_btagUp);

      tth_btagDown->SetBinContent(iBin,tthNorm*xTTH_btagDown);
      tth_btagDown->SetBinError(iBin,tthNorm*eTTH_btagDown);

      mmtt->SetBinContent(iBin,mmttNorm*xMMTT);
      mmtt->SetBinError(iBin,mmttNorm*eMMTT);

      mmtt_btagUp->SetBinContent(iBin,mmttNorm*xMMTT_btagUp);
      mmtt_btagUp->SetBinError(iBin,mmttNorm*eMMTT_btagUp);

      mmtt_btagDown->SetBinContent(iBin,mmttNorm*xMMTT_btagDown);
      mmtt_btagDown->SetBinError(iBin,mmttNorm*eMMTT_btagDown);

      double xData = histData->GetBinContent(i,j);
      double eData = histData->GetBinError(i,j);
      if (i!=j) {
	xData += histData->GetBinContent(j,i);
	double err1 = eData;
	double err2 = histData->GetBinError(j,i);
	eData = TMath::Sqrt(err1*err1);
      }
      data->SetBinContent(iBin,xData);
      data->SetBinError(iBin,eData);
      dataAll = dataAll + xData;
    }
  }

  float scaleUp = bkgNorm/bkgdCorrUp->GetSumOfWeights();
  float scaleDown = bkgNorm/bkgdCorrDown->GetSumOfWeights();
  for (int iB=1; iB<=nBins1D; ++iB) {
    double xBkg = bkgdCorrUp->GetBinContent(iB);
    double eBkg = bkgdCorrUp->GetBinError(iB);
    bkgdCorrUp->SetBinContent(iB,scaleUp*xBkg);
    bkgdCorrUp->SetBinError(iB,scaleUp*eBkg);
    xBkg = bkgdCorrDown->GetBinContent(iB);
    eBkg = bkgdCorrDown->GetBinError(iB);
    bkgdCorrDown->SetBinContent(iB,scaleDown*xBkg);
    bkgdCorrDown->SetBinError(iB,scaleDown*eBkg);
  }

  for (int i=0; i<3; ++i) {
    double scale = bkgNorm/bkgd[i]->GetSumOfWeights();
    for (int iB=1; iB<=nBins1D; ++iB) {
      double xBkg =  bkgd[i]->GetBinContent(iB);
      double eBkg =  bkgd[i]->GetBinError(iB);
      bkgd[i]->SetBinContent(iB,scale*xBkg);
      bkgd[i]->SetBinError(iB,scale*eBkg);
    }
  }

  TString BaseName = "haa_2018-13TeV_ma" + mass;
  TString rootFileName = BaseName+".root";
  TFile * fileInputs = new TFile(rootFileName,"recreate"); 
  if (Azimov)
    bkgd[0]->Write("data_obs");
  else 
    data->Write("data_obs");
  
  TString sysBkgd[3] = {"","_CMS_unc1d_2018Up","_CMS_unc1d_2018Down"}; // 
  for (int i=0; i<3; ++i) {
    bkgd[i]->Write("bkgd"+sysBkgd[i]);
  }
  
  ggh->Write("ggh");
  ggh_btagUp->Write("ggh_CMS_btag_2018Up");
  ggh_btagDown->Write("ggh_CMS_btag_2018Down");
  vbf->Write("vbf");
  vbf_btagUp->Write("vbf_CMS_btag_2018Up");
  vbf_btagDown->Write("vbf_CMS_btag_2018Down");
  vh->Write("vh");
  vh_btagUp->Write("vh_CMS_btag_2018Up");
  vh_btagDown->Write("vh_CMS_btag_2018Down");
  tth->Write("tth");
  tth_btagUp->Write("tth_CMS_btag_2018Up");
  tth_btagDown->Write("tth_CMS_btag_2018Down");
  mmtt->Write("mmtt");
  mmtt_btagUp->Write("mmtt_CMS_btag_2018Up");
  mmtt_btagDown->Write("mmtt_CMS_btag_2018Down");
  bkgdCorrUp->Write("bkgd_CMS_uncCorr_Up");
  bkgdCorrDown->Write("bkgd_CMS_uncCorr_Down");

  fileInputs->Close();

  ostringstream str;
  str << BaseName << ".txt";
  string nn = str.str();
  const char * p = nn.c_str();

  std::ofstream textFile(p);
  textFile << "imax 1   number of channels" << std::endl;
  textFile << "jmax *   number of backgrounds" << std::endl;
  textFile << "kmax *   number of nuisance parameters" << std::endl;
  textFile << "-----------------" << std::endl;
  if (Azimov)
    textFile << "observation " << bkgd[0]->GetSumOfWeights() << std::endl;
  else
    textFile << "observation " << data->GetSumOfWeights() << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "shapes * * " << rootFileName << "  $PROCESS    $PROCESS_$SYSTEMATIC " << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "bin    ";
  textFile << "    haa_2018"
	   << "    haa_2018"
	   << "    haa_2018"
	   << "    haa_2018"
	   << "    haa_2018"
	   << "    haa_2018" << std::endl;
  textFile << "process                mmtt      tth       vh      vbf      ggh    bkgd" << std::endl;
  textFile << "process                  -4       -3       -2       -1        0       1" << std::endl;
  textFile << "rate     " 
	   << mmtt->GetSumOfWeights() << "  " 
	   << tth->GetSumOfWeights() << "  "
	   << vh->GetSumOfWeights() << "  "
	   << vbf->GetSumOfWeights() << "  " 
	   << ggh->GetSumOfWeights() << "  " 
	   << bkgd[0]->GetSumOfWeights() << std::endl;
  textFile << "-----------------------------" << std::endl;
  textFile << "CMS_lumi                    lnN   1.016  1.016   1.016   1.016   1.016      -" << std::endl;
  textFile << "CMS_lumi_2018               lnN   1.025  1.025   1.025   1.025   1.025      -" << std::endl;
  textFile << "CMS_eff_m_2018              lnN   1.02    1.02    1.02    1.02    1.02      -" << std::endl;
  textFile << "CMS_trkiso_2018             lnN   1.12    1.12    1.12    1.12    1.12      -" << std::endl;
  textFile << "CMS_unc1d_2018            shape      -      -       -       -       -    1.00" << std::endl;
  textFile << "CMS_uncCorr_              shape      -      -       -       -       -    1.00" << std::endl;
  textFile << "CMS_btag_2018             shape   1.00    1.00    1.00     1.00   1.00       -" << std::endl;

  textFile << "QCDScale_ggH                lnN   1.046/0.933   -       -       -  1.046/0.933  -" << std::endl;
  textFile << "QCDScale_vbf                lnN      -          -       -  1.004/0.997  -       -" << std::endl;
  textFile << "QCDScale_vh                 lnN      -          - 1.018/0.983   -       -       -" << std::endl;
  textFile << "QCDScale_ttH                lnN      -    1.058/0.908   -    -       -       -" << std::endl;

  textFile << "PDF_ggh                     lnN   1.032      -       -       -    1.032      -" << std::endl;
  textFile << "PDF_vbf                     lnN      -       -       -    1.021      -       -" << std::endl;
  textFile << "PDF_vh                      lnN      -       -    1.018      -       -       -" << std::endl;
  textFile << "PDF_tth                     lnN      -    1.036      -       -       -       -" << std::endl;

  textFile << "acc_ggh                     lnN      -       -       -       -    1.025       -" << std::endl;
  textFile << "acc_ggh_mmtt                lnN   1.030      -       -       -    1.025       -" << std::endl;
  textFile << "acc_vbf                     lnN      -       -       -    1.02       -        -" << std::endl;
  textFile << "acc_vh                      lnN      -       -    1.021      -       -        -" << std::endl;
  textFile << "acc_tth                     lnN      -     1.02       -      -       -        -" << std::endl;

  textFile << "bkgNorm_2018   rateParam  haa_2018  bkgd  1  [0.5,1.5]" << std::endl;
  textFile << "* autoMCStats 0" << std::endl;
  textFile << std::endl;
  std::cout << "Datacards production completed for mass ma=" << mass << std::endl; 


}
