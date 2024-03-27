#include "HtoH.h"

std::map<TString, TH2D*> GetHistograms(const TString& dir, const std::map<TString, TString>& signalProcess, std::vector<TString> variations, const TString& mass) {
	std::map<TString, TH2D*> histogramsOld;

	for (const auto& process : signalProcess) {
		TString filePath = dir + process.second + mass + ".root";
 		TFile* file = TFile::Open(filePath, "READ");
		if (!file || file->IsZombie()) {
         	   std::cerr << "Failed to open file: " << filePath << std::endl;
            	   continue;
       		}
		
		for (const auto& var : variations) {
		   TString histName = "InvMass2DH" + var;
	           TH2D* hist = dynamic_cast<TH2D*>(file->Get(histName));
		   if (!hist) {
                	std::cerr << "Histogram " << histName << " not found in " << filePath << std::endl;
           	   }
		hist->SetDirectory(0);

		TString uniqueName = "hist" + process.first + var +"_Old" ;
		histogramsOld[uniqueName] = hist;
		}
		
		file->Close();
        	delete file;	
	}
	
	return histogramsOld;
}

std::map<TString, double> CalculateNormFactors(const TString& dir, const std::map<TString, TString>& signalProcess, const TString& mass, const TString& era) {

	std::map<TString, double> eraLumi = {
       	 {"2016_preVFP", 19520},
       	 {"2016_postVFP", 16810},
       	 {"2017", 41480},
       	 {"2018", 59830}
        };
	
	double lumi = eraLumi[era];

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
	double xsecMMTT = (xsecGGH + xsecVBF + xsecVH + xsecTTH) * SF;

	std::map<TString, double> normFactors;

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

void UnrollHistograms(const std::map<TString, TString>& signalProcess, std::vector<TString> variations, std::map<TString, TH2D*>& rebinnedHistograms, std::map<TString, double>& normFactors, std::map<TString, TH1D*>& unrolledSigHistograms, int nBinsNew) {
    	

	int iBin = 0;
  
    	for (int i=1; i<=nBinsNew; ++i) {
	  for (int j=i; j<=nBinsNew; ++j) {

               iBin++;
    
	       for (const auto& processPair : signalProcess) {
                   TString processName = processPair.first;
                   TString processLowerCase = processName;
                   processLowerCase.ToLower();

        	   double normFactor = normFactors[processName];
	
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

                    	content *= normFactor;
                    	error *= normFactor;

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

std::map<TString, TString> CalculateJESUncertainties(const TString& dir, const std::map<TString, TString>& signalProcess, const TString& mass) {
    std::map<TString, TString> jesUncertaintyStrings;

    for (const auto& process : signalProcess) {
        TString filePath = dir + process.second + mass + ".root";
        TFile* file = TFile::Open(filePath, "READ");

        TH1D* counter_btagH = dynamic_cast<TH1D*>(file->Get("counter_btagH"));
        TH1D* counter_btag_jesUpH = dynamic_cast<TH1D*>(file->Get("counter_btag_jesUpH"));
        TH1D* counter_btag_jesDownH = dynamic_cast<TH1D*>(file->Get("counter_btag_jesDownH"));

        double jesUpRatio = counter_btag_jesUpH->GetSumOfWeights() / counter_btagH->GetSumOfWeights();
        double jesDownRatio = counter_btag_jesDownH->GetSumOfWeights() / counter_btagH->GetSumOfWeights();
        TString jesString = TString::Format("%.3f/%.3f", jesUpRatio, jesDownRatio);
        jesUncertaintyStrings[process.first] = jesString;

        file->Close();
        delete file;
    }
    return jesUncertaintyStrings;
}

std::map<TString, std::map<TString, double>> CalculateAdditionalUncertainties(const TString& dir, const std::map<TString, TString>& signalProcess, const TString& mass) {
    std::map<TString, std::map<TString, double>> signalUncertainties;

    for (const auto& process : signalProcess) {
        TString filePath = dir + process.second + mass + ".root";
        TFile* file = TFile::Open(filePath, "READ");

        TH1D* counter_FinalEventsH = dynamic_cast<TH1D*>(file->Get("counter_FinalEventsH"));
        double nFinalEvents = counter_FinalEventsH->GetSumOfWeights();
	if (nFinalEvents == 0) {
            file->Close();
            delete file;
            continue;
        }

        // Names of histograms for different uncertainties
        std::vector<TString> uncTypes = {"btagCorrUp", "btagUncorrUp", "mistagCorrUp", "mistagUncorrUp", "prefireUp", "prefireDown"};
        for (const auto& uncType : uncTypes) {
            TString histName = "counter_" + uncType;
            TH1D* counterHist = dynamic_cast<TH1D*>(file->Get(histName));
            double nUncEvents = counterHist->GetSumOfWeights();

            // Calculate uncertainty
            double unc = nUncEvents / nFinalEvents;
	    unc = std::round(unc * 1000.0) / 1000.0;
            signalUncertainties[process.first][uncType] = unc;
        }

        file->Close();
        delete file;
    }

    return signalUncertainties;
}
	
void CreateCards(TString mass="5", // mass of pseudoscalar
                    TString era="2018",// "2018", "2017", "2016_preVFP", "2016_postVFP" 
        	    bool Azimov = true, // replace data by background expectations 
	            bool correlation = true // apply correlations
                    ) {

    TString dir = "./";
    if  (era == "2016_preVFP" || era == "2016_postVFP") {
	dir = "./" + era + "/";
    }

    std::map<TString, TString> signalProcess = {
        {"GGH", "SUSYGluGluToHToAA_AToTauTau_M-125_M-"},
        {"VBF", "SUSYVBFToHToAA_AToTauTau_M-125_M-"},
        {"VH", "SUSYVH_HToAA_AToTauTau_M-125_M-"},
        {"TTH", "SUSYttH_HToAA_AToTauTau_M-125_M-"},
        {"MMTT", "SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-"}
    };

    std::vector<TString> variations = {"", "_btagUp", "_btagDown", "_mistagUp", "_mistagDown", "_trkIsoUp", "_trkIsoDown"};

    TFile * file     = new TFile(dir+"/DoubleMuon_Run"+era+".root");
    TH2D * histOld    = (TH2D*)file->Get("InvMass2DH");

    // Fetch signal root files and open histograms
    auto histSignalOld = GetHistograms(dir, signalProcess, variations, mass);
    auto normFactors = CalculateNormFactors(dir, signalProcess, mass, era);

    TFile * fileCorr;
    if (era == "2016_preVFP" || era == "2016_postVFP"){
	fileCorr = new TFile("CorrCoefficients_data_2016.root");
    } else {
	fileCorr = new TFile(dir+"/CorrCoefficients_data_"+era+".root");
    } 
    TH2D * corrCoeff = (TH2D*)fileCorr->Get("corrCoeff");

    TFile * fileCorrCR;
    if (era == "2016_preVFP" || era == "2016_postVFP"){
        fileCorrCR = new TFile("CorrCoefficients_control_mc_2016.root");
    } else {
	fileCorrCR = new TFile(dir+"/CorrCoefficients_control_mc_"+era+".root");
    }
    TH2D * corrCoeffCR = (TH2D*)fileCorrCR->Get("corrCoeff");

    TFile * fileCorrSR;
    if (era == "2016_preVFP" || era == "2016_postVFP"){
        fileCorrSR = new TFile("CorrCoefficients_signal_mc_2016.root");
    } else {
	fileCorrSR = new TFile(dir+"/CorrCoefficients_signal_mc_"+era+".root");
    }
    TH2D * corrCoeffSR = (TH2D*)fileCorrSR->Get("corrCoeff");

    TH1D * hist1dN23Old  = (TH1D*)file->Get("InvMassN23H");
    TH1D * hist1dN45Old  = (TH1D*)file->Get("InvMassN45H");

    int nBinsNew = 6;
    double bins[7] = {0,1,2,3,4,5.2,20};

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

    // Rebin signal histograms
    auto rebinnedSigHistograms = RebinHistograms(histSignalOld, nBinsNew, bins);

    TH2D * histData = (TH2D*)TH2DtoTH2D(histOld,nBinsNew,bins,nBinsNew,bins,"_dataNew");

    double bkgNorm = histData->GetSumOfWeights();
    double sideBandNorm = 1/hist1d->GetSumOfWeights();

    std::cout << "Bkg  Norm = " << bkgNorm << std::endl;
    std::cout << "Signal expectations for BR(H->aa->4tau)=1.0" << std::endl;
    std::cout << "ggH  Norm = " << normFactors["GGH"]*rebinnedSigHistograms["histGGH"]->GetSumOfWeights() << std::endl;
    std::cout << "VBF  Norm = " << normFactors["VBF"]*rebinnedSigHistograms["histVBF"]->GetSumOfWeights() << std::endl; 
    std::cout << "VH  Norm = " << normFactors["VH"]*rebinnedSigHistograms["histVH"]->GetSumOfWeights() << std::endl;
    std::cout << "TTH  Norm = " << normFactors["TTH"]*rebinnedSigHistograms["histTTH"]->GetSumOfWeights() << std::endl;
    std::cout << "MMTT  Norm = " << normFactors["MMTT"]*rebinnedSigHistograms["histMMTT"]->GetSumOfWeights() << std::endl;
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
    UnrollHistograms(signalProcess, variations, rebinnedSigHistograms, normFactors, unrolledSigHistograms, nBinsNew);


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

    TString BaseName = "haa_"+era+"-13TeV_ma" + mass;
    TString rootFileName = BaseName+".root";
    TFile * fileInputs = new TFile(rootFileName,"recreate"); 
    if (Azimov)
        bkgd[0]->Write("data_obs");
    else 
        data->Write("data_obs");
  
    TString sysBkgd[3] = {"","_CMS_unc1d_"+era+"Up","_CMS_unc1d_"+era+"Down"}; 
    for (int i=0; i<3; ++i) {
        bkgd[i]->Write("bkgd"+sysBkgd[i]);
    } 
   
    bkgdCorrUp->Write("bkgd_CMS_uncCorrUp");
    bkgdCorrDown->Write("bkgd_CMS_uncCorrDown");

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
                newVar.ReplaceAll("_trkIso", "_CMS_trkiso_");
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

    if (era == "2016_preVFP" || era == "2016_postVFP") {
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

    auto jesUnc = CalculateJESUncertainties(dir, signalProcess, mass);
    auto Unc = CalculateAdditionalUncertainties(dir, signalProcess, mass);

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
    textFile << "process                  -4       -3       -2       -1        0       1" << std::endl;
    textFile << "rate     " 
	    << unrolledSigHistograms["mmtt"]->GetSumOfWeights() << "  " 
	    << unrolledSigHistograms["tth"]->GetSumOfWeights() << "  "
	    << unrolledSigHistograms["vh"]->GetSumOfWeights() << "  "
	    << unrolledSigHistograms["vbf"]->GetSumOfWeights() << "  " 
	    << unrolledSigHistograms["ggh"]->GetSumOfWeights() << "  " 
	    << bkgd[0]->GetSumOfWeights() << std::endl;
    textFile << "-----------------------------" << std::endl;
    textFile << "CMS_lumi_"+era+"               lnN   " << lumiUnc_era << "  " << lumiUnc_era << "   " << lumiUnc_era << "   " << lumiUnc_era << "   " << lumiUnc_era << "      -" << std::endl;
    textFile << "CMS_lumi_corr                  lnN   " << lumiUnc_corr << "  " << lumiUnc_corr << "  " << lumiUnc_corr << "  " << lumiUnc_corr << "  " << lumiUnc_corr << "     -" << std::endl;
    textFile << "CMS_lumi_1718_corr             lnN   " << lumiUnc_1718_corr << "  " << lumiUnc_1718_corr << "  " << lumiUnc_1718_corr << "  " << lumiUnc_1718_corr << "  " << lumiUnc_1718_corr << "     -" << std::endl;
    textFile << "CMS_eff_m                      lnN   1.03   1.03    1.03    1.03    1.03      -" << std::endl;
    textFile << "CMS_trkiso_"+ era+ "           shape   1.00   1.00    1.00    1.00    1.00      -" << std::endl;
    textFile << "CMS_unc1d_"+era+"              shape      -      -       -       -       -    1.00" << std::endl;
    textFile << "CMS_uncCorr                    shape      -      -       -       -       -    1.00" << std::endl;
    textFile << "CMS_btag_"+era+"               shape   1.00    1.00    1.00     1.00   1.00       -" << std::endl;
    textFile << "CMS_mistag_"+era+"             shape   1.00    1.00    1.00     1.00   1.00       -" << std::endl;
    textFile << "CMS_jes_"+era+"                lnN  "  << jesUnc["MMTT"] << " " << jesUnc["TTH"] << " " << jesUnc["VH"] << " " << jesUnc["VBF"] << " " << jesUnc["GGH"] << "     -" << std::endl;      
    textFile << "CMS_btag_corr                  lnN    "  << Unc["MMTT"]["btagCorrUp"] << " " << Unc["TTH"]["btagCorrUp"] << " " <<  Unc["VH"]["btagCorrUp"] << " " <<  Unc["VBF"]["btagCorrUp"] << " " <<  Unc["GGH"]["btagCorrUp"] <<  "   -" << std::endl;
    textFile << "CMS_btag_uncorr_"+era+"        lnN    "  << Unc["MMTT"]["btagUncorrUp"] << " " << Unc["TTH"]["btagUncorrUp"] << " " << Unc["VH"]["btagUncorrUp"] << " " << Unc["VBF"]["btagUncorrUp"] << " " << Unc["GGH"]["btagUncorrUp"] <<  "   -" << std::endl;
    textFile << "CMS_mistag_corr                lnN    "  << Unc["MMTT"]["mistagCorrUp"] << " " << Unc["TTH"]["mistagCorrUp"] << " " << Unc["VH"]["mistagCorrUp"] << " " <<  Unc["VBF"]["mistagCorrUp"] << " " << Unc["GGH"]["mistagCorrUp"] <<  "   -" << std::endl;
    textFile << "CMS_mistag_uncorr_"+era+"      lnN    "  << Unc["MMTT"]["mistagUncorrUp"] << " " << Unc["TTH"]["mistagUncorrUp"] << " " << Unc["VH"]["mistagUncorrUp"] << " " << Unc["VBF"]["mistagUncorrUp"] << " " << Unc["GGH"]["mistagUncorrUp"] << "   -" << std::endl; 
    textFile << "CMS_prefire_"+era+"            lnN     "  << Unc["MMTT"]["prefireUp"] << " " << Unc["TTH"]["prefireUp"] << " " << Unc["VH"]["prefireUp"] << " " << Unc["VBF"]["prefireUp"] << " " << Unc["GGH"]["prefireUp"] << "   -" << std::endl;

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

    textFile << "bkgNorm_"+era+"   rateParam  haa_"+era+"  bkgd  1  [0.5,1.5]" << std::endl;
    textFile << "* autoMCStats 25 1" << std::endl;
    textFile << std::endl;
    std::cout << "Datacards production completed for mass ma=" << mass << std::endl; 


}
