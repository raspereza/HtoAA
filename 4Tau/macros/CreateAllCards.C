#include "CreateCards.C"
void CreateAllCards(bool Azimov = true,
		    bool correlation = true
		    ) {

  vector<TString> masses = {"4","5","6","7","8","9","10","11","12","13","14","15"};
  
  for (auto mass : masses)
    CreateCards(mass,Azimov,correlation);
  
}
