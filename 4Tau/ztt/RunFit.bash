#!/bin/bash
era=$1
bin=$2
option=$3

folder=/nfs/dust/cms/user/rasp/HtoAA/datacards/TrkID/PuppiMET

eras=(2016 2017 2018)

options=(robustFit robustHesse)

n=$#

if [[ $n -ne 3 ]]; then
    echo Run scripts with 3 arguments: ./RunFit.bash [ERA] [BIN] [OPTION]
    echo ERA = 2016, 2017, 2018
    echo BIN = 5to10, 10to15, 5to15, 15to20, 20toInf, cone0, cone0p15, cone0p30, cone0p45
    echo OPTION = robustFit, robustHesse
    exit
fi

cd ${folder}

echo 
if [[ -f "zmm_${era}.txt" ]]; then
    echo datacards for ZMuMu control region : zmm_${era}.txt
else 
    echo datacards file zmm_${era}.txt does not exist
    exit
fi
if [[ -f "ztt_highMT_${bin}_${era}.txt" ]]; then
    echo datacards for highMT control region : ztt_highMT_${bin}_${era}.txt
else
    echo datacards file ztt_highMT_${bin}_${era}.txt foes not exist
    exit
fi
if [[ -f "ztt_lowMT_${bin}_${era}.txt" ]]; then
    echo datacards for highMT control region : ztt_lowMT_${bin}_${era}.txt
else
    echo datacards file ztt_lowMT_${bin}_${era}.txt foes not exist
    exit
fi
echo

# Combine cards
combineCards.py zmm_${era}.txt ztt_highMT_${bin}_${era}.txt ztt_lowMT_${bin}_${era}.txt > ztt_${bin}_${era}.txt
# Run fit
combineTool.py -M FitDiagnostics --${option} 1 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --cminDefaultMinimizerTolerance 0.05 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy 1 --rMin=0 --rMax=3 -m 91 ztt_${bin}_${era}.txt -n _${bin}_${era}_robustFit -v 3
