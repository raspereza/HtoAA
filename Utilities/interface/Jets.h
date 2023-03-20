#ifndef Jets_h
#define Jets_h

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
bool tightJetID(float pfjet_e,
		float pfjet_eta,
		float pfjet_neutralhadronicenergy,
		float pfjet_neutralemenergy,
		float pfjet_muonenergy,
		float pfjet_chargedhadronicenergy,
		float pfjet_chargedemenergy,
		unsigned int pfjet_neutralmulti,
		unsigned int pfjet_chargedmulti,
		int era){

  bool tightJetID = false;
  float energy = pfjet_e;
  float eta = pfjet_eta;
  float NHF = pfjet_neutralhadronicenergy / energy;
  float NEMF = pfjet_neutralemenergy / energy;
  float NumConst = pfjet_chargedmulti + pfjet_neutralmulti;     
  float CHM = pfjet_chargedmulti;
  float MUF = pfjet_muonenergy / energy;
  float CHF = pfjet_chargedhadronicenergy / energy;
  float CEMF = pfjet_chargedemenergy / energy;
  float NumNeutralParticle  = pfjet_neutralmulti;	
  if (era == 2016){
    if (fabs(eta) <= 2.7)
      tightJetID = (NHF < 0.90 && NEMF < 0.90 && NumConst > 1) && ((abs(eta) <= 2.4 && CHF > 0 && CHM > 0 && CEMF < 0.99) || abs(eta) > 2.4);
    else if (fabs(eta) <= 3.0)
      tightJetID = NHF < 0.98 && NEMF > 0.01 && NumNeutralParticle > 2;
    else
      tightJetID = NEMF < 0.90 && NumNeutralParticle > 10;
  }
  else if (era == 2017){
    if (fabs(eta) <= 2.7)
      tightJetID = (NHF < 0.90 && NEMF < 0.90 && NumConst > 1) && ((abs(eta) <= 2.4 && CHF > 0 && CHM > 0) || abs(eta) > 2.4);
    else if (fabs(eta) <= 3.0)
      tightJetID = NEMF < 0.99 && NEMF > 0.02 && NumNeutralParticle > 2;
    else
      tightJetID = NEMF < 0.90 && NHF > 0.02 && NumNeutralParticle > 10;
  }
  else if (era == 2018){
    if (fabs(eta) <= 2.6)
      tightJetID = CEMF < 0.8 && CHM > 0 && CHF > 0 && NumConst > 1 && NEMF < 0.9 && MUF < 0.8 && NHF < 0.9;
    else if (fabs(eta) <= 2.7)
      tightJetID = CEMF < 0.8 && CHM > 0 && NEMF < 0.99 && MUF < 0.8 && NHF < 0.9;
    else if (fabs(eta) <= 3.0)
      tightJetID = NEMF > 0.02 && NEMF < 0.99 && NumNeutralParticle > 2;
    else
      tightJetID = NEMF < 0.90 && NHF > 0.2 && NumNeutralParticle > 10;
  }
  else
    {
      std::cout << "era is not 2016, 2017, 2018, exiting" << '\n';
      exit(-1);
    }
  return tightJetID;
}

bool jetPUID(float pt, float eta, float mva, TString wp="Tight"){
  // from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID#Working_points
  //4 Eta Categories 0-2.5 2.5-2.75 2.75-3.0 3.0-5.0
  //Tight Id
  size_t etaIndex=0;
  vector<float> Pt010 = { 0.69, -0.35, -0.26, -0.21};
  vector<float> Pt1020 = { 0.69, -0.35, -0.26, -0.21};
  vector<float> Pt2030 = { 0.69, -0.35, -0.26, -0.21};
  vector<float> Pt3050 = { 0.86, -0.10, -0.05, -0.01};
  if(wp=="Medium"){
    //Medium Id
    Pt010 = { 0.18, -0.55, -0.42, -0.36};
    Pt1020 = { 0.18, -0.55, -0.42, -0.36};
    Pt2030 = { 0.18, -0.55, -0.42, -0.36};
    Pt3050 = { 0.61, -0.35, -0.23, -0.17};
  }else if(wp=="Loose"){
    //Loose Id
    Pt010 = {-0.97, -0.68, -0.53, -0.47};
    Pt1020 = {-0.97, -0.68, -0.53, -0.47};
    Pt2030 = {-0.97, -0.68, -0.53, -0.47};
    Pt3050 = {-0.89, -0.52, -0.38, -0.30};
  }


  if(fabs(eta) < 2.5) etaIndex=0;
  else if(fabs(eta) < 2.75) etaIndex=1;
  else if(fabs(eta) < 3.0)  etaIndex=2;
  else if(fabs(eta) < 5.0)  etaIndex=3;

  float cut = -1.;

  if(pt<10) cut = Pt010[etaIndex];
  else if(pt<20) cut = Pt1020[etaIndex];
  else if(pt<30) cut = Pt2030[etaIndex];
  else if(pt<50) cut = Pt3050[etaIndex];

  bool passedPUID = (bool) (mva > cut);
  if (pt >= 50) passedPUID = true;
  if (eta >=5.0)passedPUID = false;

  return passedPUID;

}



#endif
