#!/bin/bash

samples=(DYJetsToLL_M-10to50
DYJetsToLL_M-50
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
QCD_Pt-1000ToInf_MuEnrichedPt5
WW_13TeV-pythia8
WZ_13TeV-pythia8
ZZ_13TeV-pythia8
)

j=0
while [ $j -lt ${#samples[@]} ] 
do
    echo "Submitting jobs on sample " ${samples[$j]} 
    ./HTC_submit_seq.sh analysis_macro_4tau analysisMacro_mc_2017.conf ${samples[$j]} 100
    j=`expr $j + 1` 
done
