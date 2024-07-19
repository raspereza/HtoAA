// Systematic uncertainty (parton shower scale) 
// parameterized for convenience as a function of (mu,trk) mass
double sysUnc(double mass,
	      TString era,
	      bool isISR) {

  double alpha = 0.015;
  double beta = 12.2;
  if (isISR) {
    alpha = 0.014;
    beta = 21.6;
  }
  if (era=="2017") {
    alpha = 0.017;
    beta = 12.5;
    if (isISR) {
      alpha = 0.015;
      beta = 20.2;
    }
  }
  else if (era=="2018") {
    alpha = 0.017;
    beta = 12.7;
    if (isISR) {
      alpha = 0.016;
      beta = 20.8;
    }
  }
  return 1.0+alpha*TMath::Exp(mass/beta);    

}

void CheckSyst(TString era = "2018",
	       bool isISR = false) {

  TString PS = "FSR";
  if (isISR) PS = "ISR";
  std::cout << std::endl;
  std::cout << era << "  " << PS << std::endl;
  std::cout << std::endl;

  double mass[6] = {0.5, 1.5, 2.5, 3.5, 4.5, 9.0};
  for (int iB=0; iB<6; ++iB) {
    for (int jB=iB; jB<6; ++jB) {
      double mass1 = mass[iB];
      double mass2 = mass[jB];
      double corr1 = sysUnc(mass1,era,isISR);
      double corr2 = sysUnc(mass2,era,isISR);
      double corr = corr1 * corr2;
      printf("[%3.1f,%3.1f] -> %5.3f\n",mass1,mass2,corr);
    }
  }
  
}
