#!/bin/bash

echo submitting jobs for HToAA_AToMuMu_AToTauTau samples
for i in {4..21}
do
    ./HTC_submit_seq.sh analysis_macro_4tau analysisMacro_ggH_2mu2tau_2016APV.conf SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-${i} 5
done

echo submitting jobs for HToAA_To4Tau samples
for i in {4..21}
do
    ./HTC_submit_seq.sh analysis_macro_4tau analysisMacro_ggH_2016APV.conf SUSYGluGluToHToAA_AToTauTau_M-125_M-${i} 5
    ./HTC_submit_seq.sh analysis_macro_4tau analysisMacro_VBF_2016APV.conf SUSYVBFToHToAA_AToTauTau_M-125_M-${i} 5
    ./HTC_submit_seq.sh analysis_macro_4tau analysisMacro_VH_2016APV.conf SUSYVH_HToAA_AToTauTau_M-125_M-${i} 5
    
    ./HTC_submit_seq.sh analysis_macro_4tau analysisMacro_ttH_2016APV.conf SUSYttH_HToAA_AToTauTau_M-125_M-${i} 5
done

