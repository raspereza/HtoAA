#!/bin/bash
for i in {4..21}
do
    ./hadd_mc.sh SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-${i}
    ./hadd_mc.sh SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}
    ./hadd_mc.sh SUSYVBFToHToAA_AToTauTau_M-125_M-${i}
    ./hadd_mc.sh SUSYVH_HToAA_AToTauTau_M-125_M-${i}
    ./hadd_mc.sh SUSYttH_HToAA_AToTauTau_M-125_M-${i}
done

   
