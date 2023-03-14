#!/bin/bash
for i in {4..21}
do
    ./hadd.sh SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-${i}
    ./hadd.sh SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}
    ./hadd.sh SUSYVBFToHToAA_AToTauTau_M-125_M-${i}
    ./hadd.sh SUSYVH_HToAA_AToTauTau_M-125_M-${i}
    ./hadd.sh SUSYttH_HToAA_AToTauTau_M-125_M-${i}
done


   
