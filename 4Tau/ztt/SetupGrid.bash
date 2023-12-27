#!/bin/bash
#$1 - era (2016_preVFP, 2016_postVFP, 2017 and 2018)
n=$#

if [ $n -ne 1 ]; then
    echo Usage : ./SetupGrid.bash [ERA]
    echo ERA = 2016_preVFP, 2016_postVFP, 2017 and 2018
    exit
fi

ERA=$1

if [ ]
if [ $ERA != "2016_preVFP" ] && [ $ERA != "2016_postVFP" ] && [ $ERA != "2017" ] && [ $ERA != "2018" ]; then
    echo Unavailable option for era specified : $ERA
    echo Run script with one of available choices for era : 2016_preVFP, 2016_postVFP, 2017, 2018
    exit
fi
cp ${CMSSW_BASE}/src/HtoAA/4Tau/ztt/*.sh ./
cp ${CMSSW_BASE}/src/HtoAA/4Tau/ztt/gc_synch.conf ./
cp ${CMSSW_BASE}/src/HtoAA/4Tau/ztt/${ERA}/FileListMaker${ERA}.sh ./
cp ${CMSSW_BASE}/src/HtoAA/4Tau/ztt/${ERA}/analysisMacro_ztt.conf ./
./FileListMaker${ERA}.sh
echo "+++++++++++++++++++++++++++++++++++++++++++++"
echo "Don't forget to modify gc_synch.conf file!!!!"
echo "Then run grid control tool with this config. "
echo "+++++++++++++++++++++++++++++++++++++++++++++"
