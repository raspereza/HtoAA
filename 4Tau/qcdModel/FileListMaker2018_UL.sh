#!/bin/bash

names_QCD=(QCD_Pt-15To20_MuEnrichedPt5
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

if [ -f "parameters.txt" ]; then
    rm parameters.txt
fi
echo "CONFIGFILE,FILELIST" > parameters.txt

j=0
while [ $j -lt ${#names_QCD[@]} ] 
do
    echo "Creating file list for sample" ${names_QCD[$j]} 
    ls /pnfs/desy.de/cms/tier2/store/user/rasp/ntuples/H2AA/2018/mc/${names_QCD[$j]}/*root > ${names_QCD[$j]}
    ./split_filelist.sh analysisMacro_mutrk_2018.conf ${names_QCD[$j]} 100

    j=`expr $j + 1` 
done
