#!/bin/bash

echo creating file lists for HToAA_To4Tau samples
for i in {4..21}
do
    ls /pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2018/mc_2/SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}_pythia8/*root > SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}
done

