void ExtractSF( TString era = "2017" ) {

  std::vector<TString> ptbins = {"5to10","10to15","15to20","20toInf"};
  std::vector<TString> cones = {"cone0","cone0p15","cone0p30","cone0p45"};

  TString folder("/nfs/dust/cms/user/rasp/HtoAA/datacards/TrkID/PuppiMET");

  std::cout << std::endl;
  std::cout << "+----------------------------------------------+" << std::endl;
  std::cout << "|                   " << era << "                        |" << std::endl;
  // 
  std::cout << "+---------+------------------+-----------------+" << std::endl;
  std::cout << "|  pt bin |     robustFit    |   robustHesse   |" << std::endl;
  std::cout << "+---------+------------------+-----------------+" << std::endl;
  for (auto ptbin : ptbins) {

    vector<TString> filenames;
    filenames.push_back(folder+"/fitDiagnostics_"+ptbin+"_"+era+"_robustFit.root");
    filenames.push_back(folder+"/fitDiagnostics_"+ptbin+"_"+era+"_robustHesse.root");
    
    double scale[3] = {1., 1., 1.};
    double error[4] = {0., 0., 0.};
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
      file->Close();
      delete file;
      counter++;
    }
    printf("| %7s |  %5.3f +/- %5.3f | %5.3f +/- %5.3f |\n",
	   ptbin.Data(),scale[0],error[0],scale[1],error[1]);
    
  }
  std::cout << "+------------+------------------+--------------+" << std::endl;

  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "+-----------------------------------------------+" << std::endl;
  std::cout << "|                   " << era << "                        |" << std::endl;
  // 
  //           "  5to10 :  1.174 +/- 0.060 | 1.174 +/- 0.357 |"
  std::cout << "+----------+------------------+-----------------+" << std::endl;
  std::cout << "|   cone   |     robustFit    |   robustHesse   |" << std::endl;
  std::cout << "+----------+------------------+-----------------+" << std::endl;
  for (auto cone : cones) {
      
    vector<TString> filenames;
    filenames.push_back(folder+"/fitDiagnostics_"+cone+"_"+era+"_robustFit.root");
    filenames.push_back(folder+"/fitDiagnostics_"+cone+"_"+era+"_robustHesse.root");
    
    double scale[3] = {1., 1., 1.};
    double error[4] = {0., 0., 0.};
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
      file->Close();
      delete file;
      counter++;
    }
    printf("| %8s |  %5.3f +/- %5.3f | %5.3f +/- %5.3f |\n",
	   cone.Data(),scale[0],error[0],scale[1],error[1]);
  }
  std::cout << "+------------+------------------+---------------+" << std::endl;


}
