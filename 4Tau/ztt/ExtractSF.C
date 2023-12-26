std::vector<TString> ptbins = {"5to15","10to15","15to20","20toInf"};

void ExtractSF(TString era = "2018") {

  TString folder("/nfs/dust/cms/user/rasp/Run/MuTrk");

  std::cout << std::endl;
  std::cout << "Era : " << era << std::endl;
  std::cout << std::endl;
  
  for (auto ptbin : ptbins) {

    vector<TString> filenames;
    filenames.push_back(folder+"/fitDiagnostics_"+ptbin+"_"+era+".root");
    filenames.push_back(folder+"/fitDiagnostics_"+ptbin+"_"+era+"_simple.root");
    filenames.push_back(folder+"/fitDiagnostics_lowMT_"+ptbin+"_"+era+"_simple.root");
    double scale[3] = {1.,  1., 1.};
    double error[4] = {0.,  0., 0.};
    unsigned int counter = 0;
    for (auto filename : filenames) { 
  
      TFile * file = new TFile(filename);
      if (file->IsZombie()) { 
	std::cout << "file not found " << filename << std::endl; 
	return;
      }

      RooFitResult * fitres = (RooFitResult*)file->Get("fit_s");
      if (fitres==NULL) 
	continue;
      RooArgList list = fitres->floatParsFinal();

      int n = list.getSize();
      for (int iPar=0; iPar<n;++iPar) {
	RooRealVar* var = (RooRealVar*)list.at(iPar);
	//      std::cout << "nuisance : " << var << std::endl;
	if (var==NULL) continue;
	TString name = var->GetName();
	if (name=="r") {
	  scale[counter] = var->getVal();
	  error[counter] = var->getError();
	  break;
	}
      }
      counter++;
    }
    printf("%10s :  %5.3f +/- %5.3f | %5.3f +/- %5.3f | %5.3f +/- %5.3f |\n",
	   ptbin.Data(),scale[0],error[0],scale[1],error[1],scale[2],error[2]);

  }
}
