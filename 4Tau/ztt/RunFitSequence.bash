#!/bin/bash
era=$1

folder=/nfs/dust/cms/user/rasp/HtoAA/datacards/TrkID/PuppiMET
eras=(2016 2017 2018)
ptbins=(5to10 10to15 5to15 15to20 20toInf)
cones=(cone0 cone0p15 cone0p30 cone0p45)
n=$#

if [[ $n -ne 1 ]]; then
    echo pass era as a single argument to the routine : ./RunFit.bash [ERA]
    exit
fi

counter=0
for i in ${eras[@]} 
do
    if [[ $era -eq $i ]]; then
	counter=`expr $counter + 1`
    fi
done

if [[ $counter -eq 0 ]]; then
    echo unknown era $era
    echo available options : 2016, 2017, 2018
    exit
fi

cd ${folder}

# Loop over pt bins
for bin in ${ptbins[@]} 
do
    # Combine cards
    combineCards.py zmm_${era}.txt ztt_highMT_${bin}_${era}.txt ztt_lowMT_${bin}_${era}.txt > ztt_${bin}_${era}.txt
    # RobustFit option
    combineTool.py -M FitDiagnostics --robustFit 1 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --cminDefaultMinimizerTolerance 0.05 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy 0 --rMin=0 --rMax=2 -m 91 ztt_${bin}_${era}.txt -n _${bin}_${era}_robustFit 
    # RobustHesse option
    combineTool.py -M FitDiagnostics --robustHesse 1 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --cminDefaultMinimizerTolerance 0.05 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy 0 --rMin=0 --rMax=2 -m 91 ztt_${bin}_${era}.txt -n _${bin}_${era}_robustHesse 
done

##### Loop over cones ######
for bin in ${cones[@]} 
do
    # Combine cards 
    combineCards.py zmm_${era}.txt ztt_highMT_${bin}_${era}.txt ztt_lowMT_${bin}_${era}.txt > ztt_${bin}_${era}.txt
    # RobustFit option
    combineTool.py -M FitDiagnostics --robustFit 1 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --cminDefaultMinimizerTolerance 0.05 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy 0 --rMin=0 --rMax=2 -m 91 ztt_${bin}_${era}.txt -n _${bin}_${era}_robustFit 
    # RobustHesse option
    combineTool.py -M FitDiagnostics --robustHesse 1 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --cminDefaultMinimizerTolerance 0.05 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy 0 --rMin=0 --rMax=2 -m 91 ztt_${bin}_${era}.txt -n _${bin}_${era}_robustHesse 
done
rm higgsCombine*root

cd -
