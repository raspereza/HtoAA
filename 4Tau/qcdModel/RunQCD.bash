#!/bin/bash

n=$#

if [[ $n -ne 1 ]]; then
    echo pass era as a single argument to the routine : ./RunQCD.bash [ERA]
    exit
fi

era=$1

samples=(QCD_Pt-15To20_MuEnrichedPt5
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
    echo "Submitting jobs on sample " ${samples[$j]} 
    ./HTC_submit_seq.sh mutrk_mass analysisMacro_mutrk_${era}.conf ${samples[$j]}
    j=`expr $j + 1` 
done

