#!/bin/bash
for j in B-ver1 B-ver2 C D E F
do
    ./resubmitJobs.sh DoubleMuon_Run2016APV${j}
done

samples=(DYJetsToLL_M-10to50
DYJetsToLL_M-50
WW_13TeV-pythia8
WZ_13TeV-pythia8
ZZ_13TeV-pythia8
ST_t-channel_top
ST_t-channel_antitop
ST_tW_top
ST_tW_antitop
TTTo2L2Nu
TTToHadronic
TTToSemiLeptonic
WJetsToLNu
QCD_Pt-20To30_MuEnrichedPt5
QCD_Pt-30To50_MuEnrichedPt5
QCD_Pt-50To80_MuEnrichedPt5
QCD_Pt-80To120_MuEnrichedPt5
QCD_Pt-120To170_MuEnrichedPt5
QCD_Pt-170To300_MuEnrichedPt5
QCD_Pt-300To470_MuEnrichedPt5
QCD_Pt-470To600_MuEnrichedPt5
QCD_Pt-600To800_MuEnrichedPt5
QCD_Pt-800To1000_MuEnrichedPt5
QCD_Pt-1000_MuEnrichedPt5
)

j=0
while [ $j -lt ${#samples[@]} ] 
do
    ./resubmitJobs.sh ${samples[$j]}
    j=`expr $j + 1` 
done

for i in {4..21}
do
    ./resubmitJobs.sh SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-${i}	
    ./resubmitJobs.sh SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}
    ./resubmitJobs.sh SUSYVBFToHToAA_AToTauTau_M-125_M-${i}
    ./resubmitJobs.sh SUSYVH_HToAA_AToTauTau_M-125_M-${i}
    ./resubmitJobs.sh SUSYttH_HToAA_AToTauTau_M-125_M-${i}
done



   
