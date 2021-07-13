#!/bin/bash

for i in {4..15}
do
    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/SUSYVH_HToAA_AToTauTau_M-125_M-${i}/*.root > SUSYVH_HToAA_AToTauTau_M-125_M-${i}
done

