#include "HtoH.h"

std::map<TString, double> eraLumi = {
  {"2016_preVFP", 19520},
  {"2016_postVFP", 16810},
  {"2017", 41480},
  {"2018", 59830}
};

std::map<TString, double> CalculateAcceptances(const TString& dir, const std::map<TString, TString>& signalProcess, const TString& mass, const std::vector<TString>& group) {
    std::map<TString, double> acceptances;

    for (const auto& process : signalProcess) {
        double sumOfWeights = 0;
        double ngen = 0;
        
        for (const auto& era : group) {
            TString filePath = dir + "/" + era + "/" + process.second + mass + ".root";
            TFile* file = TFile::Open(filePath, "READ");
            if (!file || file->IsZombie()) {
                std::cerr << "Failed to open file: " << filePath << std::endl;
                continue;
            }

            TH1D* histWeights = dynamic_cast<TH1D*>(file->Get("histWeightsH"));
            TH2D* histSignal = dynamic_cast<TH2D*>(file->Get("InvMass2DH"));
            if (!histWeights || !histSignal) {
                std::cerr << "Failed to extract histogram from file: " << filePath << std::endl;
                file->Close();
                delete file;
                continue;
            }

            ngen += histWeights->GetSumOfWeights();
            sumOfWeights += histSignal->GetSumOfWeights();

            file->Close();
            delete file;
        }

        if (ngen > 0) {
            acceptances[process.first] = sumOfWeights / ngen;
        } else {
            acceptances[process.first] = 0;
        }
    }

    return acceptances;
}

std::map<TString, double> CalculateNormFactors(const TString& dir, const std::map<TString, TString>& signalProcess, const TString& mass, const TString& era, const std::map<TString, double>& acceptance4tau) {

        double xsecGGH = 48.61;
	double xsecVBF = 3.766;
	double xsecVH  = 1.358 + 0.880;
	double xsecTTH = 0.5065;
	double massD = mass.Atof();

	double massTau = 1.777;
	double massMu  = 0.106;
	double massRatio = (massMu * massMu) / (massTau * massTau);
	double aF = 2 * massTau / massD;
	double SF = 2 * massRatio / TMath::Sqrt(1 - aF * aF);

	double combinedMMTTXsec = xsecGGH +
	                          xsecVBF * (acceptance4tau.at("VBF") / acceptance4tau.at("GGH")) +
        	                  xsecVH * (acceptance4tau.at("VH") / acceptance4tau.at("GGH")) +
                	          xsecTTH * (acceptance4tau.at("TTH") / acceptance4tau.at("GGH"));

	double xsecMMTT = combinedMMTTXsec * SF;

	std::map<TString, double> normFactors;

	double lumi = eraLumi[era];

	for (const auto& process : signalProcess) {
        	TString filePath = dir + process.second + mass + ".root";
        	TFile* file = TFile::Open(filePath, "READ");
        	if (!file || file->IsZombie()) {
            		std::cerr << "Failed to open file: " << filePath << std::endl;
            		continue;
       		 }

		TH1D* histWeights = dynamic_cast<TH1D*>(file->Get("histWeightsH"));
		double nGen = histWeights ? histWeights->GetSumOfWeights() : 0;
		if (!histWeights) {
    		 std::cerr << "Failed to extract histogram from file: " << filePath << std::endl;
		}
        	file->Close();
        	delete file;

        	double xsec = 0;
        	if (process.first == "GGH") xsec = xsecGGH;
        	else if (process.first == "VBF") xsec = xsecVBF;
        	else if (process.first == "VH") xsec = xsecVH;
        	else if (process.first == "TTH") xsec = xsecTTH;
        	else if (process.first == "MMTT") xsec = xsecMMTT;

		if (nGen > 0) {
            		normFactors[process.first] = (xsec * lumi) / nGen;
       		 } else {
           		 normFactors[process.first] = 0;
       		 }
    	}
    	return normFactors;
}

double signalXSec(const TString& process, const TString& mass, const std::map<TString, double>& acceptance4tau) {
   
  double xsecGGH = 48.61;
  double xsecVBF = 3.766;
  double xsecVH  = 1.358 + 0.880;
  double xsecTTH = 0.5065;
  double massD = mass.Atof();

  double massTau = 1.777;
  double massMu  = 0.106;
  double massRatio = (massMu * massMu) / (massTau * massTau);
  double aF = 2 * massTau / massD;
  double SF = 2 * massRatio / TMath::Sqrt(1 - aF * aF);

  double combinedMMTTXsec = xsecGGH +
                            xsecVBF * (acceptance4tau.at("VBF") / acceptance4tau.at("GGH")) +
                            xsecVH * (acceptance4tau.at("VH") / acceptance4tau.at("GGH")) +
                            xsecTTH * (acceptance4tau.at("TTH") / acceptance4tau.at("GGH"));

  double xsecMMTT = combinedMMTTXsec * SF;

  double xsec = 0.;
  if (process=="GGH")
    xsec = xsecGGH;
  else if (process=="VBF")
    xsec = xsecVBF;
  else if (process=="VH")
    xsec = xsecVH;
  else if (process=="TTH")
    xsec = xsecTTH;
  else if (process=="MMTT")
    xsec = xsecMMTT;

  return xsec;

}

std::map<TString, TH2D*> GetHistograms(const TString& dir, const std::map<TString, TString>& signalProcess, std::vector<TString> variations, const TString& mass, const std::vector<TString>& group, const std::map<TString, double>& acceptance4tau) {

  std::map<TString, TH2D*> histogramsOld;
  
  for (const auto & process : signalProcess) {
    bool isFirst = true;
    double xsec = signalXSec(process.first,mass,acceptance4tau);
    for (const auto & era : group) {	    
      TString filePath = dir + "/" + era + "/" + process.second + mass + ".root";
      TFile* file = TFile::Open(filePath, "READ");
      double lumi = eraLumi[era];
      TH1D * histWeights = (TH1D*)file->Get("histWeightsH");
      double ngen = histWeights->GetSumOfWeights();
      double norm = 0; 
      if (ngen>0) norm = lumi*xsec/ngen;
      if (!file || file->IsZombie()) {
	std::cerr << "Failed to open file: " << filePath << std::endl;
	continue;
      }
      
      for (const auto& var : variations) {
	TString histName = "InvMass2DH" + var;
	TH2D* hist = (TH2D*)file->Get(histName);
	if (!hist) {
	  std::cerr << "Histogram " << histName << " not found in " << filePath << std::endl;
	}
	hist->Scale(norm);
	TString uniqueName = "hist" + process.first + var +"_Old" ;
	if (isFirst) {
	  hist->SetDirectory(0);	  
	  histogramsOld[uniqueName] = hist;
	}
	else {	  
	  histogramsOld[uniqueName]->Add(histogramsOld[uniqueName],hist);
	}
      }
      isFirst = false;
      file->Close();
      delete file;	
    }
  }
  return histogramsOld;
}

std::map<TString, TH2D*> RebinHistograms(const std::map<TString, TH2D*>& histSignalOld, int nBinsNew, double* bins) {
  std::map<TString, TH2D*> rebinnedHistograms;

  for (const auto& histPair : histSignalOld) {
    const TString& histNameOld = histPair.first;
    TH2D* originalHist = histPair.second;
    
    TString histName = histNameOld;
    histName.Resize(histName.Length() - 4);
    
    TH2D* rebinnedHist = (TH2D*)TH2DtoTH2D(originalHist, nBinsNew, bins, nBinsNew, bins, histName);
    rebinnedHist->SetDirectory(0);
    rebinnedHistograms[histName] = rebinnedHist;
    
  }
  
  return rebinnedHistograms;
}

void UnrollHistograms(const std::map<TString, TString>& signalProcess, std::vector<TString> variations, std::map<TString, TH2D*>& rebinnedHistograms, std::map<TString, TH1D*>& unrolledSigHistograms, int nBinsNew) {
    	

  int iBin = 0;
  
  for (int i=1; i<=nBinsNew; ++i) {
    for (int j=i; j<=nBinsNew; ++j) {
      
      iBin++;
      
      for (const auto& processPair : signalProcess) {
	TString processName = processPair.first;
	TString processLowerCase = processName;
	processLowerCase.ToLower();
	
	//	double normFactor = normFactors[processName];
	
	for (const auto& var : variations) {
	  TString rebinnedHistName = "hist" + processName + var;
	  TH2D* hist2D = rebinnedHistograms[rebinnedHistName];
	  
	  double content = hist2D->GetBinContent(i, j);
	  double error = hist2D->GetBinError(i, j);
	  
	  if (i != j) {
	    content += hist2D->GetBinContent(j, i);
	    double err1 = error;
	    double err2 = hist2D->GetBinError(j, i);
	    error = TMath::Sqrt(err1 * err1 + err2 * err2);
	  }
	  
	  // normalization is done in GetHistograms
	  //	  content *= normFactor;
	  //	  error *= normFactor;
	  
	  TString unrolledHistName = processLowerCase + var;
	  TH1D* unrolledHist = unrolledSigHistograms[unrolledHistName];
	  
	  unrolledHist->SetBinContent(iBin, unrolledHist->GetBinContent(iBin) + content);
	  double currentError = unrolledHist->GetBinError(iBin);
	  unrolledHist->SetBinError(iBin, sqrt(currentError * currentError + error * error));
	}
	
      }      
    }
  }   
}

std::map<TString, TString> CalculateJESUncertainties(const TString& dir, const std::map<TString, TString>& signalProcess, const TString& mass, const std::vector<TString>& group) {

    std::map<TString, TString> jesUncertaintyStrings;

    for (const auto& process : signalProcess) {

      double centralCounter = 0;
      double upCounter = 0;
      double downCounter = 0;
	
      for (const auto& era : group) {
	TString filePath = dir + "/" + era + "/" + process.second + mass + ".root";
	TFile* file = TFile::Open(filePath, "READ");
	double lumi = eraLumi[era];
	TH1D * histCentral = (TH1D*)file->Get("counter_btagH");
	TH1D * histUp = (TH1D*)file->Get("counter_btag_jesUpH");
	TH1D * histDown = (TH1D*)file->Get("counter_btag_jesDownH");
	TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
	double nGen = histWeightsH->GetSumOfWeights();
	double scale = lumi/nGen;
	centralCounter += scale*histCentral->GetSumOfWeights();
	upCounter += scale*histUp->GetSumOfWeights();
	downCounter += scale*histDown->GetSumOfWeights();
	file->Close();
	delete file;
      }      
      double jesUpRatio = 1;
      double jesDownRatio = 1;
      if (centralCounter>0) {
	if (upCounter>0) jesUpRatio = upCounter/centralCounter;
	if (downCounter) jesDownRatio = downCounter/centralCounter;
      }
      TString jesString = TString::Format("%.3f", jesUpRatio);
      jesUncertaintyStrings[process.first] = jesString;

    }
    return jesUncertaintyStrings;
}

std::map<TString, std::map<TString, double>> CalculateAdditionalUncertainties(const TString& dir, const std::map<TString, TString>& signalProcess, const TString& mass, const std::vector<TString>& group) {

  std::map<TString, std::map<TString, double>> signalUncertainties;

  std::vector<TString> uncTypes = {"btagCorrUp", "btagUncorrUp", "mistagCorrUp", "mistagUncorrUp", "prefireUp"};

  for (const auto& process : signalProcess) {

    double nFinalEvents = 0;
    std::map<TString, double> nUncEvents;
    for (const auto& uncType : uncTypes)
      nUncEvents[uncType] = 0;

    for (const auto& era : group) {

      TString filePath = dir + "/" + era + "/" + process.second + mass + ".root";
      TFile* file = TFile::Open(filePath, "READ");
      TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
      double lumi = eraLumi[era];
      double nGen = histWeightsH->GetSumOfWeights();
      double scale = lumi/nGen;
    
      TH1D* counter_FinalEventsH = dynamic_cast<TH1D*>(file->Get("counter_FinalEventsH"));
      nFinalEvents += scale*counter_FinalEventsH->GetSumOfWeights();
      
      if (nFinalEvents == 0) {
	file->Close();
	delete file;
	continue;
      }
    
      for (const auto& uncType : uncTypes) {
	TString histName = "counter_" + uncType;
	TH1D* counterHist = dynamic_cast<TH1D*>(file->Get(histName));
	nUncEvents[uncType] += scale*counterHist->GetSumOfWeights();
	
      }
      file->Close();
      delete file;
    }
    // Calculate uncertainty
    for (const auto& uncType : uncTypes) {
      double unc = 1.0;
      if (nFinalEvents>0&&nUncEvents[uncType]>0) unc = nUncEvents[uncType] / nFinalEvents;
      unc = std::round(unc * 1000.0) / 1000.0;
      signalUncertainties[process.first][uncType] = unc;
    }
  } 
  return signalUncertainties;
}

void CreateCards(TString mass="5", // mass of pseudoscalar
		 TString era="2018",// "2018", "2017", "2016_preVFP", "2016_postVFP" 
		 bool Azimov = true, // replace data by background expectations 
		 bool correlation = true, // apply mass correlation coefficients
		 bool MassUncertPerEras = true, // decorrelate mass uncert between eras
		 bool MassUncertPerBins = true, // decorrelate 1D pdf uncert across bins
		 TString bkgUnc = "0", // background uncertainty, if set 0, bkg norm is floated freely
		 bool multiProc = false // if set to true background normalizatin is treated as an additional POI by Combine
		 ) {

    TH1::SetDefaultSumw2(true);
    TH2::SetDefaultSumw2(true);

    std::map<TString, vector<TString> > groups = {
      {"2016",{"2016_postVFP","2016_preVFP"}},
      {"2016_preVFP",{"2016_preVFP"}},
      {"2016_postVFP",{"2016_postVFP"}},
      {"2017",{"2017"}},
      {"2018",{"2018"}},
    };
   
    vector<TString> group = groups[era];
 
    TString dir = "/nfs/dust/cms/user/sreelatl/Analyses/H2aa_4tau/Run2/Jul24";

    std::map<TString, TString> signalProcess = {
        {"GGH", "SUSYGluGluToHToAA_AToTauTau_M-125_M-"},
        {"VBF", "SUSYVBFToHToAA_AToTauTau_M-125_M-"},
        {"VH", "SUSYVH_HToAA_AToTauTau_M-125_M-"},
        {"TTH", "SUSYttH_HToAA_AToTauTau_M-125_M-"},
        {"MMTT", "SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-"}
    };

    std::map<TString, TString> signalProcess4tau = {
        {"GGH", signalProcess["GGH"]},
        {"VBF", signalProcess["VBF"]},
        {"VH", signalProcess["VH"]},
        {"TTH", signalProcess["TTH"]}
    };

    auto acceptance4tau = CalculateAcceptances(dir, signalProcess4tau, mass, group);

    std::vector<TString> variations = {"", "_btagUp", "_btagDown", "_mistagUp", "_mistagDown", "_trkIsoUp", "_trkIsoDown"};

    std::map<TString, TFile*> files;
    for (const auto& grp : group )
      files[grp] = new TFile(dir+"/"+grp+"/DoubleMuon_Run"+grp+".root");

    TH2D * histOld = NULL;
    
    std::vector<TString> histnames = {
      "InvMassN23H",
      "InvMassN45H",
      "InvMassHardestNtrk23H",
      "InvMassSoftestNtrk23H",
      "InvMassHardestNtrk1H",
      "InvMassSoftestNtrk1H"
    };
    std::map<TString,TH1D*> hist1D; 
    bool isFirst = true;
    for (const auto& grp : group ) {
      TFile * file = files[grp];
      TH2D * histOld_ = (TH2D*)file->Get("InvMass2DH");
      if (isFirst) {	  
	  histOld = histOld_;
	  histOld->SetDirectory(0);
      }
      else {
	histOld->Add(histOld,histOld_,1.,1.);
      }
      for (const auto& histname : histnames ) {
	TH1D * hist1D_ = (TH1D*)file->Get(histname);
	if (isFirst) {
	  hist1D[histname] = hist1D_;
	  hist1D[histname]->SetDirectory(0);
	}
	else {
	  hist1D[histname]->Add(hist1D[histname],hist1D_,1.,1.);
	}
      }
      file->Close();
      delete file;
      isFirst = false;
    }

    // Fetch signal root files and open histograms
    auto histSignalOld = GetHistograms(dir, signalProcess, variations, mass, group, acceptance4tau);
    //    auto normFactors = CalculateNormFactors(dir, signalProcess, mass, era);

    TFile * fileCorr;
    if (era == "2016_preVFP" || era == "2016_postVFP") {
      fileCorr = new TFile("CorrCoefficients_data_2016.root");
    } else {
	fileCorr = new TFile("CorrCoefficients_data_"+era+".root");
    } 
    TH2D * corrCoeff = (TH2D*)fileCorr->Get("corrCoeff");

    TFile * fileCorrCR;
    if (era == "2016_preVFP" || era == "2016_postVFP"){
        fileCorrCR = new TFile("CorrCoefficients_control_mc_2016.root");
    } else {
	fileCorrCR = new TFile("CorrCoefficients_control_mc_"+era+".root");
    }
    TH2D * corrCoeffCR = (TH2D*)fileCorrCR->Get("corrCoeff");
    TH2D * corrCoeffCR_nonQCDUp = (TH2D*)fileCorrCR->Get("corrCoeff_nonQCDUp");
    TH2D * corrCoeffCR_nonQCDDown = (TH2D*)fileCorrCR->Get("corrCoeff_nonQCDDown");
    TH2D * corrCoeffCR_ISRUp = (TH2D*)fileCorrCR->Get("corrCoeff_ISRUp");
    TH2D * corrCoeffCR_ISRDown = (TH2D*)fileCorrCR->Get("corrCoeff_ISRDown");
    TH2D * corrCoeffCR_FSRUp = (TH2D*)fileCorrCR->Get("corrCoeff_FSRUp");
    TH2D * corrCoeffCR_FSRDown = (TH2D*)fileCorrCR->Get("corrCoeff_FSRDown");

    TFile * fileCorrSR;
    if (era == "2016_preVFP" || era == "2016_postVFP"){
        fileCorrSR = new TFile("CorrCoefficients_signal_mc_2016.root");
    } else {
	fileCorrSR = new TFile("CorrCoefficients_signal_mc_"+era+".root");
    }
    TH2D * corrCoeffSR = (TH2D*)fileCorrSR->Get("corrCoeff");
    TH2D * corrCoeffSR_nonQCDUp = (TH2D*)fileCorrSR->Get("corrCoeff_nonQCDUp");
    TH2D * corrCoeffSR_nonQCDDown = (TH2D*)fileCorrSR->Get("corrCoeff_nonQCDDown");
    TH2D * corrCoeffSR_ISRUp = (TH2D*)fileCorrSR->Get("corrCoeff_ISRUp");
    TH2D * corrCoeffSR_ISRDown = (TH2D*)fileCorrSR->Get("corrCoeff_ISRDown");
    TH2D * corrCoeffSR_FSRUp = (TH2D*)fileCorrSR->Get("corrCoeff_FSRUp");
    TH2D * corrCoeffSR_FSRDown = (TH2D*)fileCorrSR->Get("corrCoeff_FSRDown");

    int nBinsNew = 6;
    double bins[7] = {0,1,2,3,4,5.2,20};

    std::map<TString, TH1D* > hist1D_rebinned;
    for (const auto& histname : histnames) {
      hist1D_rebinned[histname] = (TH1D*)TH1DtoTH1D(hist1D[histname],nBinsNew,bins,true,"_new");
    }

    TH1D * hist1d     = hist1D_rebinned["InvMassN23H"];
    TH1D * hist1dHardestUp   = hist1D_rebinned["InvMassHardestNtrk1H"];
    TH1D * hist1dHardestDown = hist1D_rebinned["InvMassHardestNtrk23H"];
    TH1D * hist1dSoftestUp   = hist1D_rebinned["InvMassSoftestNtrk1H"];
    TH1D * hist1dSoftestDown = hist1D_rebinned["InvMassSoftestNtrk23H"];
    
    hist1d->Scale(1.0/hist1d->GetSumOfWeights());
    hist1dHardestUp->Scale(1.0/hist1dHardestUp->GetSumOfWeights());
    hist1dHardestDown->Scale(1.0/hist1dHardestDown->GetSumOfWeights());
    hist1dSoftestUp->Scale(1.0/hist1dSoftestUp->GetSumOfWeights());
    hist1dSoftestDown->Scale(1.0/hist1dSoftestDown->GetSumOfWeights());

    TH1D * uncertMass1DH = new TH1D("uncertMass1DH","",nBinsNew,bins);
    for (int iB=1; iB<=nBinsNew; ++iB) {

      double numHardest = hist1dHardestUp->GetBinContent(iB);
      double denHardest = hist1dHardestDown->GetBinContent(iB);
      double numHardestE = hist1dHardestUp->GetBinError(iB);
      double denHardestE = hist1dHardestDown->GetBinError(iB);
      double ratioHardest = numHardest/denHardest;
      double numHardestR = numHardestE/numHardest;
      double denHardestR = denHardestE/denHardest;
      double ratioHardestR = TMath::Sqrt(numHardestR*numHardestR+denHardestR*denHardestR);
      double ratioHardestE = ratioHardest*ratioHardestR;

      double numSoftest = hist1dSoftestUp->GetBinContent(iB);
      double denSoftest = hist1dSoftestDown->GetBinContent(iB);
      double numSoftestE = hist1dSoftestUp->GetBinError(iB);
      double denSoftestE = hist1dHardestDown->GetBinError(iB);
      double ratioSoftest = numSoftest/denSoftest;
      double numSoftestR = numSoftestE/numSoftest;
      double denSoftestR = denSoftestE/denSoftest;
      double ratioSoftestR = TMath::Sqrt(numSoftestR*numSoftestR+denSoftestR*denSoftestR);
      double ratioSoftestE = ratioSoftest*ratioSoftestR;

      double ratio = ratioSoftest;
      double ratioE = ratioSoftestE;
      if (ratioHardest>ratio) {
	ratio = ratioHardest;
	ratioE = ratioHardestE;
      }

      uncertMass1DH->SetBinContent(iB,ratio);
      uncertMass1DH->SetBinError(iB,ratioE);
    }

    // Rebin signal histograms
    auto rebinnedSigHistograms = RebinHistograms(histSignalOld, nBinsNew, bins);

    TH2D * histData = (TH2D*)TH2DtoTH2D(histOld,nBinsNew,bins,nBinsNew,bins,"_dataNew");

    double bkgNorm = histData->GetSumOfWeights();
    double sideBandNorm = 1/hist1d->GetSumOfWeights();

    std::cout << "Bkg  Norm = " << bkgNorm << std::endl;
    std::cout << "Signal expectations for BR(H->aa->4tau)=1.0" << std::endl;
    std::cout << "ggH  Norm = " << rebinnedSigHistograms["histGGH"]->GetSumOfWeights() << std::endl;
    std::cout << "VBF  Norm = " << rebinnedSigHistograms["histVBF"]->GetSumOfWeights() << std::endl; 
    std::cout << "VH  Norm = " << rebinnedSigHistograms["histVH"]->GetSumOfWeights() << std::endl;
    std::cout << "TTH  Norm = " << rebinnedSigHistograms["histTTH"]->GetSumOfWeights() << std::endl;
    std::cout << "MMTT  Norm = " << rebinnedSigHistograms["histMMTT"]->GetSumOfWeights() << std::endl;
    std::cout << std::endl;  

    int nBins1D = nBinsNew * (nBinsNew + 1) / 2;

    // Create new histograms to store unrolled signal
    std::map<TString, TH1D*> unrolledSigHistograms;

    for (const auto& processPair : signalProcess) {
        const TString& process = processPair.first;
        TString processLowerCase = process;
        processLowerCase.ToLower();

        for (const auto& variation : variations) {
            TString histName = processLowerCase + variation;
            if (variation == "") {
            histName = processLowerCase;
            }

            unrolledSigHistograms[histName] = new TH1D(histName, "", nBins1D, 0., double(nBins1D));

        }   
    } 

    // central template and systematic background templates
    // 1d_up - up variation of 1D background pdf (mu-trk mass)
    // 1d_down - down variation of 1D background pdf (mu-trk mass)
    // per-bin variations are also considered when 
    // effect is decorrelated across the bins of 1D pdf f1D(m)
    TH1D * bkgd[15];
    TString bkg_name[15] = {"",
			    "_1d_up","_1d_down",
			    "_1d_bin1_up","_1d_bin1_down",
			    "_1d_bin2_up","_1d_bin2_down",
			    "_1d_bin3_up","_1d_bin3_down",
			    "_1d_bin4_up","_1d_bin4_down",
			    "_1d_bin5_up","_1d_bin5_down",
			    "_1d_bin6_up","_1d_bin6_down"};

    for (int i = 0; i<15; ++i) {
        bkgd[i] = new TH1D("bkgd"+bkg_name[i],"",nBins1D,0.,float(nBins1D));
    }

    TH1D *bkgd_FSRUp = new TH1D("bkgd_FSRUp", "", nBins1D, 0., float(nBins1D));
    TH1D *bkgd_FSRDown = new TH1D("bkgd_FSRDown", "", nBins1D, 0., float(nBins1D));
    TH1D *bkgd_ISRUp = new TH1D("bkgd_ISRUp", "", nBins1D, 0., float(nBins1D));
    TH1D *bkgd_ISRDown = new TH1D("bkgd_ISRDown", "", nBins1D, 0., float(nBins1D));
    TH1D *bkgd_nonQCDUp = new TH1D("bkgd_nonQCDUp", "", nBins1D, 0., float(nBins1D));
    TH1D *bkgd_nonQCDDown = new TH1D("bkgd_nonQCDDown", "", nBins1D, 0., float(nBins1D));

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
      double excorr = uncertMass1DH->GetBinError(i);
      double eycorr = uncertMass1DH->GetBinError(j);

      double xBkgd[20];
      double yBkgd[20];
      double exBkgd[20];
      double eyBkgd[20];

      // apply 100% of difference as uncertainty !
      // b-only fit largely constrain 
      // corresponding nuisance parameters  
      // 
      float xshift = 1.0 + 1.0*(xcorr-1.0);
      float yshift = 1.0 + 1.0*(ycorr-1.0);

      xBkgd[0] = hist1d->GetBinContent(i)*xshift;
      yBkgd[0] = hist1d->GetBinContent(j)*yshift;
      exBkgd[0] = hist1d->GetBinError(i)*xshift;
      eyBkgd[0] = hist1d->GetBinError(j)*yshift;

      for (unsigned int itempl=1; itempl<15; ++itempl) {
	xBkgd[itempl] = xBkgd[0];
	exBkgd[itempl] = exBkgd[0];
	yBkgd[itempl] = yBkgd[0];
	eyBkgd[itempl] = eyBkgd[0];
      }


      xBkgd[1] = xBkgd[0]*xshift;
      yBkgd[1] = yBkgd[0]*yshift;

      xBkgd[2] = xBkgd[0]/xshift;
      yBkgd[2] = yBkgd[0]/yshift;

      xBkgd[1+2*i] = xBkgd[0]*xshift;
      yBkgd[1+2*j] = yBkgd[0]*yshift;

      xBkgd[2+2*i] = xBkgd[0]/xshift;
      yBkgd[2+2*j] = yBkgd[0]/yshift;

      for (int itempl=0; itempl<15; ++itempl) {
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
	double corrR = TMath::Sqrt(corrDataR*corrDataR+corrCRR*corrCRR+corrSRR*corrSRR);
	double corrE = corrX * corrR;

	double corrCR_FSRUp = corrCoeffCR_FSRUp->GetBinContent(corrBin);
        double corrCR_FSRDown = corrCoeffCR_FSRDown->GetBinContent(corrBin);
        double corrCR_ISRUp = corrCoeffCR_ISRUp->GetBinContent(corrBin);
        double corrCR_ISRDown = corrCoeffCR_ISRDown->GetBinContent(corrBin);
        double corrCR_nonQCDUp = corrCoeffCR_nonQCDUp->GetBinContent(corrBin);
        double corrCR_nonQCDDown = corrCoeffCR_nonQCDDown->GetBinContent(corrBin);

	double corrSR_FSRUp = corrCoeffSR_FSRUp->GetBinContent(corrBin);
        double corrSR_FSRDown = corrCoeffSR_FSRDown->GetBinContent(corrBin);
        double corrSR_ISRUp = corrCoeffSR_ISRUp->GetBinContent(corrBin);
        double corrSR_ISRDown = corrCoeffSR_ISRDown->GetBinContent(corrBin);
        double corrSR_nonQCDUp = corrCoeffSR_nonQCDUp->GetBinContent(corrBin);
        double corrSR_nonQCDDown = corrCoeffSR_nonQCDDown->GetBinContent(corrBin);

	double corr_FSRUp = corrSR_FSRUp * corrData / corrCR_FSRUp;
	double corr_FSRDown = corrSR_FSRDown * corrData / corrCR_FSRDown;;
	double corr_ISRUp = corrSR_ISRUp * corrData / corrCR_ISRUp;
	double corr_ISRDown = corrSR_ISRDown * corrData / corrCR_ISRDown;
	double corr_nonQCDUp = corrSR_nonQCDUp * corrData / corrCR_nonQCDUp;
	double corr_nonQCDDown = corrSR_nonQCDDown * corrData / corrCR_nonQCDDown;

	double product_FSRUp = product;
        double product_FSRDown = product;
        double product_ISRUp = product;
        double product_ISRDown = product;
        double product_nonQCDUp = product;
        double product_nonQCDDown = product;

        if (i!=j) {
	  product *= 2.0;
	  err *= 2.0;
	  product_FSRUp *= 2.0;
          product_FSRDown *= 2.0;
          product_ISRUp *= 2.0;
          product_ISRDown *= 2.0;
          product_nonQCDUp *= 2.0;
          product_nonQCDDown *= 2.0;
	}
	if (correlation) {
	  product *= corrX;
	  err *= corrX;
	  product_FSRUp *= corr_FSRUp;
	  product_FSRDown *= corr_FSRDown;
	  product_ISRUp *= corr_ISRUp;
	  product_ISRDown *= corr_ISRDown;
	  product_nonQCDUp *= corr_nonQCDUp;
	  product_nonQCDDown *= corr_nonQCDDown;
	}
	double errCor = product * corrE / corrX;

        err = TMath::Sqrt(err*err+errCor*errCor);
	bkgd[itempl]->SetBinContent(iBin,product);
	bkgd[itempl]->SetBinError(iBin,err);
	if (itempl==0) {
	  bkgd_FSRUp->SetBinContent(iBin,product_FSRUp);
	  bkgd_FSRDown->SetBinContent(iBin,product*product/product_FSRUp);
	  bkgd_FSRUp->SetBinError(iBin,err);
	  bkgd_FSRDown->SetBinError(iBin,err);

	  bkgd_ISRUp->SetBinContent(iBin,product_ISRUp);
          bkgd_ISRDown->SetBinContent(iBin,product*product/product_ISRUp);
          bkgd_ISRUp->SetBinError(iBin,err);
          bkgd_ISRDown->SetBinError(iBin,err);

	  bkgd_nonQCDUp->SetBinContent(iBin,product_nonQCDUp);
          bkgd_nonQCDDown->SetBinContent(iBin,product*product/product_nonQCDUp);
          bkgd_nonQCDUp->SetBinError(iBin,err);
          bkgd_nonQCDDown->SetBinError(iBin,err);

	  productAll += product;
	} 
      }

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

    // Unrolled signal templates
    UnrollHistograms(signalProcess, variations, rebinnedSigHistograms, unrolledSigHistograms, nBinsNew);

    TH1D* bkgdSystematicHistogramsUp[] = {bkgd_FSRUp, bkgd_ISRUp, bkgd_nonQCDUp};
    TH1D* bkgdSystematicHistogramsDown[] = {bkgd_FSRDown, bkgd_ISRDown, bkgd_nonQCDDown};

    for (auto& histUp : bkgdSystematicHistogramsUp) {
        float scaleUp = bkgNorm / histUp->GetSumOfWeights();
        for (int iB = 1; iB <= nBins1D; ++iB) {
            double xBkg = histUp->GetBinContent(iB);
            double eBkg = histUp->GetBinError(iB);
            histUp->SetBinContent(iB, scaleUp * xBkg);
            histUp->SetBinError(iB, scaleUp * eBkg);
   	}
    }

    for (auto& histDown : bkgdSystematicHistogramsDown) {
        float scaleDown = bkgNorm / histDown->GetSumOfWeights();
        for (int iB = 1; iB <= nBins1D; ++iB) {
            double xBkg = histDown->GetBinContent(iB);
            double eBkg = histDown->GetBinError(iB);
            histDown->SetBinContent(iB, scaleDown * xBkg);
            histDown->SetBinError(iB, scaleDown * eBkg);
        }
    }

    for (int i=0; i<15; ++i) {
        double scale = bkgNorm / bkgd[i]->GetSumOfWeights();
        for (int iB=1; iB<=nBins1D; ++iB) {
            double xBkg =  bkgd[i]->GetBinContent(iB);
            double eBkg =  bkgd[i]->GetBinError(iB);
            bkgd[i]->SetBinContent(iB,scale * xBkg);
            bkgd[i]->SetBinError(iB,scale * eBkg);
        }
    }

    TString BaseName = "haa_"+era+"-13TeV_ma" + mass;
    TString rootFileName = BaseName+".root";
    TFile * fileInputs = new TFile(rootFileName,"recreate"); 
    if (Azimov)
        bkgd[0]->Write("data_obs");
    else 
        data->Write("data_obs");
  
    TString sysBkgd[7];
    sysBkgd[0] = TString("CMS_haa4t_unc1d");
    for (unsigned int ibin=1; ibin<7; ++ibin) {
      char uncerName[30];
      sprintf(uncerName,"CMS_haa4t_unc1d_bin%1i",ibin);
      sysBkgd[ibin] = TString(uncerName);
    }

    TString massCorrUncName("CMS_haa4t_uncCorr");

    // decorrelation between eras
    if (MassUncertPerEras) {
      //massCorrUncName += "_"+era;
      for (unsigned int ibin=0; ibin<7; ++ibin)
	sysBkgd[ibin] += "_"+era;
    }

    bkgd[0]->Write("bkgd");
    for (unsigned int i=0; i<7; ++i) {
      bkgd[1+2*i]->Write("bkgd_"+sysBkgd[i]+"Up");
      bkgd[2+2*i]->Write("bkgd_"+sysBkgd[i]+"Down");
    }

    bkgd_FSRUp->Write("bkgd_" + massCorrUncName + "_FSRUp");
    bkgd_FSRDown->Write("bkgd_" + massCorrUncName + "_FSRDown");
    bkgd_ISRUp->Write("bkgd_" + massCorrUncName + "_ISRUp");
    bkgd_ISRDown->Write("bkgd_" + massCorrUncName + "_ISRDown");
    bkgd_nonQCDUp->Write("bkgd_" + massCorrUncName + "_nonQCDUp");
    bkgd_nonQCDDown->Write("bkgd_" + massCorrUncName + "_nonQCDDown");

    // Write unrolled signal templates to the root file
    for (const auto& pair : unrolledSigHistograms) {
        const TString& originalName = pair.first; 
        TH1D* hist = pair.second;                 
 
        TString processName = originalName;
        for (const auto& var : variations) {
            processName.ReplaceAll(var, "");
        } 
 
        TString sysVariationSuffix;
        for (const auto& var : variations) {
            if (originalName.Contains(var) && var != "") {
                TString newVar = var;
                newVar.ReplaceAll("_btag", "_CMS_btag_");
                newVar.ReplaceAll("_mistag", "_CMS_mistag_");
                newVar.ReplaceAll("_trkIso", "_CMS_haa4t_eff_trkiso_");
                newVar.ReplaceAll("Up", era + "Up");
                newVar.ReplaceAll("Down", era + "Down");
                sysVariationSuffix = newVar;
                break;
            }
       }

    	TString writeName = processName + sysVariationSuffix;
    
    	if (hist) {
        	hist->SetName(writeName);
        	hist->Write();
    	}

    }

    fileInputs->Close();

    ostringstream str;
    str << BaseName << ".txt";
    string nn = str.str();
    const char * p = nn.c_str();

    TString lumiUnc_era;
    TString lumiUnc_corr;
    TString lumiUnc_1718_corr;

    if (era == "2016_preVFP" || era == "2016_postVFP" || era == "2016") {
        lumiUnc_era = "1.01";
	lumiUnc_corr = "1.006";
	lumiUnc_1718_corr = "-";
    } else if (era == "2017") {
        lumiUnc_era = "1.02";
	lumiUnc_corr = "1.009";
	lumiUnc_1718_corr = "1.006";
    } else if (era == "2018") {
        lumiUnc_era = "1.025";
	lumiUnc_corr = "1.02";
	lumiUnc_1718_corr = "1.002";
    } 

    auto jesUnc = CalculateJESUncertainties(dir, signalProcess, mass, group);
    auto Unc = CalculateAdditionalUncertainties(dir, signalProcess, mass, group);

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
    textFile << "    haa_"+era
	    << "    haa_"+era
	    << "    haa_"+era
	    << "    haa_"+era
	    << "    haa_"+era
	    << "    haa_"+era << std::endl;
    textFile << "process                mmtt      tth       vh      vbf      ggh    bkgd" << std::endl;
    if (multiProc)
      textFile << "process                  -5       -4       -3       -2       -1       0" << std::endl;
    else
      textFile << "process                  -4       -3       -2       -1        0       1" << std::endl;
    textFile << "rate     " 
	    << unrolledSigHistograms["mmtt"]->GetSumOfWeights() << "  " 
	    << unrolledSigHistograms["tth"]->GetSumOfWeights() << "  "
	    << unrolledSigHistograms["vh"]->GetSumOfWeights() << "  "
	    << unrolledSigHistograms["vbf"]->GetSumOfWeights() << "  " 
	    << unrolledSigHistograms["ggh"]->GetSumOfWeights() << "  " 
	    << bkgd[0]->GetSumOfWeights() << std::endl;
    textFile << "-----------------------------" << std::endl;
    if (bkgUnc!="0"&&!multiProc)
      textFile << "bkgNorm_" << era << "         lnN     -     -     -     -     -     " << bkgUnc << std::endl;
    textFile << "lumi_"+era+"                    lnN   " << lumiUnc_era << "  " << lumiUnc_era << "   " << lumiUnc_era << "   " << lumiUnc_era << "   " << lumiUnc_era << "      -" << std::endl;
    textFile << "lumi_13TeV_correlated        lnN   " << lumiUnc_corr << "  " << lumiUnc_corr << "  " << lumiUnc_corr << "  " << lumiUnc_corr << "  " << lumiUnc_corr << "     -" << std::endl;
    textFile << "lumi_13TeV_1718              lnN   " << lumiUnc_1718_corr << "  " << lumiUnc_1718_corr << "  " << lumiUnc_1718_corr << "  " << lumiUnc_1718_corr << "  " << lumiUnc_1718_corr << "     -" << std::endl;
    textFile << "CMS_eff_m                    lnN   1.03   1.03    1.03    1.03    1.03      -" << std::endl;
    textFile << "CMS_haa4t_eff_trkiso_"+era+"  shape   1.00   1.00    1.00    1.00    1.00      -" << std::endl;
    
    if (MassUncertPerBins) {
      for (unsigned int ibin=1; ibin<=6; ++ibin)
	textFile << sysBkgd[ibin] << "  shape    -      -       -       -       -    1.0" << std::endl;
    }
    else {
      textFile << sysBkgd[0] << "               shape      -      -       -       -       -    1.0" << std::endl;
    }

    textFile << massCorrUncName + "_FSR" << "    shape    -      -       -       -       -    1.0" << std::endl;
    textFile << massCorrUncName + "_ISR" << "    shape    -      -       -       -       -    1.0" << std::endl;
    textFile << massCorrUncName + "_nonQCD" << " shape    -      -       -       -       -    1.0" << std::endl;

    textFile << "CMS_scale_j_"+era+"             lnN  "  << jesUnc["MMTT"] << " " << jesUnc["TTH"] << " " << jesUnc["VH"] << " " << jesUnc["VBF"] << " " << jesUnc["GGH"] << "     -" << std::endl;      
    textFile << "CMS_btag_heavy_corr          lnN    "  << Unc["MMTT"]["btagCorrUp"] << " " << Unc["TTH"]["btagCorrUp"] << " " <<  Unc["VH"]["btagCorrUp"] << " " <<  Unc["VBF"]["btagCorrUp"] << " " <<  Unc["GGH"]["btagCorrUp"] <<  "   -" << std::endl;
    textFile << "CMS_btag_heavy_"+era+"          lnN    "  << Unc["MMTT"]["btagUncorrUp"] << " " << Unc["TTH"]["btagUncorrUp"] << " " << Unc["VH"]["btagUncorrUp"] << " " << Unc["VBF"]["btagUncorrUp"] << " " << Unc["GGH"]["btagUncorrUp"] <<  "   -" << std::endl;
    textFile << "CMS_btag_light_corr          lnN    "  << Unc["MMTT"]["mistagCorrUp"] << " " << Unc["TTH"]["mistagCorrUp"] << " " << Unc["VH"]["mistagCorrUp"] << " " <<  Unc["VBF"]["mistagCorrUp"] << " " << Unc["GGH"]["mistagCorrUp"] <<  "   -" << std::endl;
    textFile << "CMS_btag_light_"+era+"          lnN    "  << Unc["MMTT"]["mistagUncorrUp"] << " " << Unc["TTH"]["mistagUncorrUp"] << " " << Unc["VH"]["mistagUncorrUp"] << " " << Unc["VBF"]["mistagUncorrUp"] << " " << Unc["GGH"]["mistagUncorrUp"] << "   -" << std::endl; 
    textFile << "CMS_l1_ecal_prefiring_"+era+"   lnN     "  << Unc["MMTT"]["prefireUp"] << " " << Unc["TTH"]["prefireUp"] << " " << Unc["VH"]["prefireUp"] << " " << Unc["VBF"]["prefireUp"] << " " << Unc["GGH"]["prefireUp"] << "   -" << std::endl;

    textFile << "QCDscale_ggH                 lnN   1.046/0.933   -       -       -  1.046/0.933  -" << std::endl;
    textFile << "QCDscale_qqH                 lnN      -          -       -  1.004/0.997  -       -" << std::endl;
    textFile << "QCDscale_VH                  lnN      -          - 1.018/0.983   -       -       -" << std::endl;
    textFile << "QCDscale_ttH                 lnN      -    1.058/0.908   -    -       -       -" << std::endl;

    textFile << "pdf_Higgs_gg                 lnN   1.032      -       -       -    1.032      -" << std::endl;
    textFile << "pdf_Higgs_qqbar              lnN      -       -       -    1.021      -       -" << std::endl;
    textFile << "pdf_Higgs_gq                 lnN      -       -    1.018      -       -       -" << std::endl;
    textFile << "pdf_Higgs_ttH                lnN      -    1.036      -       -       -       -" << std::endl;

    textFile << "CMS_haa4t_acc_ggH            lnN      -       -       -       -    1.025       -" << std::endl;
    textFile << "CMS_haa4t_acc_ggH_mmtt       lnN   1.030      -       -       -    1.025       -" << std::endl;
    textFile << "CMS_haa4t_acc_VBF            lnN      -       -       -    1.02       -        -" << std::endl;
    textFile << "CMS_haa4t_acc_VH             lnN      -       -    1.021      -       -        -" << std::endl;
    textFile << "CMS_haa4t_acc_ttH            lnN      -     1.02       -      -       -        -" << std::endl;

    if (bkgUnc=="0"&&!multiProc)
      textFile << "CMS_haa4t_bkgNorm_"+era+" rateParam  haa_"+era+"  bkgd  1  [0.5,1.5]" << std::endl;
    textFile << "* autoMCStats 0 1" << std::endl;
    textFile << std::endl;
    std::cout << "Datacards production completed for mass ma=" << mass << std::endl; 


}
