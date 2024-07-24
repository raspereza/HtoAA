#include "CreateCards.C"
void CreateAllCards(
                    bool Azimov = false, // create Asimov dataset
		    bool correlation = true, // use correlation coefficients C(i,j)
		    bool MassUncertPerEras = true, // decorrelate unc. in mass pdfs across eras
		    bool MassUncertPerBins = true, // decorrelate unc. in 1D mass pdfs across bins
		    TString bkgUnc = "0", // uncertainty in bkg norm, if set to 0, bkg norm is floated freely
		    bool multiProc = false // if set to true, bkg norm is treated by Combine as additional POI
		    ) {

  //  vector<TString> masses = {"4"};
  //  vector<TString> eras = {"2016_preVFP","2016_postVFP","2016"};
  vector<TString> masses = {"4","5","6","7","8","9","10","11","12","13","14","15"};
  vector<TString> eras = {"2016","2017","2018"};

  for (auto era : eras) {
    for (auto mass : masses) {
      CreateCards(mass,era,Azimov,correlation,MassUncertPerEras,MassUncertPerBins,bkgUnc,multiProc);
    }
  }

}
