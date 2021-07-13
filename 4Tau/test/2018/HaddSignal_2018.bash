#!/bin/bash
for i in {4..15}
do
#    ./hadd.sh SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-${i}
#    ./hadd.sh SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}
#    ./hadd.sh SUSYVBFToHToAA_AToTauTau_M-125_M-${i}
    ./hadd.sh SUSYVH_HToAA_AToTauTau_M-125_M-${i}
#    ./hadd.sh SUSYttH_HToAA_AToTauTau_M-125_M-${i}
done
#./hadd.sh GluGluHToTauTau_M125
#./hadd.sh VBFHToTauTau_M125
#./hadd.sh WplusHToTauTau_M125
#./hadd.sh WminusHToTauTau_M125
#./hadd.sh ZHToTauTau_M125_13TeV


   
