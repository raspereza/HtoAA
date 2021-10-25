#ifndef HTOAA_INTERFACE_FUNCTION_H
#define HTOAA_INTERFACE_FUNCTIO_H

// obsolete parameterization of muon trigger efficiencies
double MuLegEfficiency(double pt, double eta, double ptThres, double etaThres) {

  double absEta = fabs(eta);

  if (absEta>etaThres) return 0;

  double effEtaLt0p9 = 0.932;
  double effEta0p9to1p2 = 0.922;
  double effEtaGt1p2 = 0.950;

  double eff = 1.0;
  double ptThresLow = ptThres-1.0;
  double ptThresHigh = ptThres+1.0;
  if (ptThres>30) {
    ptThresLow = ptThres-3.0;
    ptThresHigh = ptThres+3.0;
    effEtaLt0p9 = 0.919;
    effEta0p9to1p2 = 0.841;
    effEtaGt1p2 = 0.866;
  }
  else if (ptThres>16) {
    effEtaLt0p9 =  0.931;
    effEta0p9to1p2 = 0.919;
    effEtaGt1p2 = 0.926;
  }

  if (pt<ptThresLow) {
    eff = 0;
  }
  else if (pt>=ptThresLow&&pt<ptThresHigh) {
    if (absEta<0.9) 
      eff = 0.5*effEtaLt0p9;
    else if (absEta>=0.9&&absEta<1.2)
      eff = 0.5*effEta0p9to1p2;
    else
      eff = 0.5*effEtaGt1p2;
  }
  else {
    if (absEta<0.9)
      eff = effEtaLt0p9;
    else if (absEta>=0.9&&absEta<1.2)
      eff = effEta0p9to1p2;
    else
      eff = effEtaGt1p2;
  }
  
  return eff;
  
}

bool looseId(
	     float energy,
	     float eta,
	     float energycorr,
	     float chargedhadronicenergy,
	     float neutralhadronicenergy,
	     float neutralemenergy,
	     float chargedemenergy,
	     float muonenergy,
	     unsigned int chargedmulti,
	     unsigned int neutralmulti
	     ) {
  bool looseJetID = false;
  energy *= energycorr; // uncorrected energy must be used                                                                                                             
  float chf = chargedhadronicenergy/energy;
  float nhf = neutralhadronicenergy/energy;
  float nem = neutralemenergy/energy;
  float elf = chargedemenergy/energy;
  float muf = muonenergy/energy;
  float chm = chargedmulti;
  float nm  = neutralmulti;
  float npr = chargedmulti + neutralmulti;

  if (fabs(eta)<=2.7)
    {
      looseJetID = (nhf < 0.99 && nem < 0.99 && npr > 1) && (fabs(eta)>2.4 || (chf>0 && chm > 0 && elf < 0.99));
    }
  else if (fabs(eta)<=3.0)
    {
      looseJetID = (nem < 0.9 && nm > 2);
    }
  else
    {
      looseJetID = nem < 0.9 && nm > 10;
    }

  return looseJetID;


}

float topPtWeight(float pt1,
                  float pt2,
		  bool run1) {
    
  float a = 0.0615;    // Run2 a parameter
  float b = -0.0005;  // Run2 b parameter

  if (run1) {
    if (pt1>400) pt1 = 400;
    if (pt2>400) pt2 = 400;
    a = 0.156;    // Run1 a parameter
    b = -0.00137;  // Run1 b parameter
  }
  float w1 = TMath::Exp(a+b*pt1);
  float w2 = TMath::Exp(a+b*pt2);
  
  return TMath::Sqrt(w1*w2);
    
}

double PtoEta(double Px, double Py, double Pz) {

  double P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
  double cosQ = Pz/P;
  double Q = TMath::ACos(cosQ);
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));
  return Eta;

}

double dPhiFrom2P(double Px1, double Py1,
                  double Px2, double Py2) {


  double prod = Px1*Px2 + Py1*Py2;
  double mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
  double mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);

  double cosDPhi = prod/(mod1*mod2);

  return TMath::ACos(cosDPhi);

}

double deltaR(double Eta1, double Phi1,
              double Eta2, double Phi2) {

  double Px1 = TMath::Cos(Phi1);
  double Py1 = TMath::Sin(Phi1);

  double Px2 = TMath::Cos(Phi2);
  double Py2 = TMath::Sin(Phi2);

  double dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);
  double dEta = Eta1 - Eta2;

  double dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

  return dR;

}

int binNumber(float x, int nbins, float * bins) {

  int binN = nbins-1;

  for (int iB=0; iB<nbins; ++iB) {
    if (x>=bins[iB]&&x<bins[iB+1]) {
      binN = iB;
      break;
    }
  }

  return binN;

}

void computeDzeta(float metX,  float metY,
                  float zetaX, float zetaY,
                  float pzetavis,
                  float & pzetamiss,
                  float & dzeta) {
    
    pzetamiss = metX*zetaX + metY*zetaY;
    dzeta = pzetamiss - 0.85*pzetavis;
    
}

const float PionMass = 0.13957;
const float MuMass = 0.105658367;

#endif
