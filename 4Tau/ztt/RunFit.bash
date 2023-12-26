#!/bin/bash
bin=$1
era=$2


folder=/nfs/dust/cms/user/rasp/Run/MuTrk

cd ${folder}

# advanced method
combineCards.py ztt_highMT_${bin}_${era}.txt ztt_lowMT_${bin}_${era}.txt > ztt_${bin}_${era}.txt

combineTool.py -M FitDiagnostics --robustFit 1 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --cminDefaultMinimizerTolerance 0.05 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy 1 --rMin=0 --rMax=2 -m 91 ztt_${bin}_${era}.txt  -n _${bin}_${era} -v 2

# simple method
combineCards.py ztt_highMT_${bin}_${era}_simple.txt ztt_lowMT_${bin}_${era}_simple.txt > ztt_${bin}_${era}_simple.txt

combineTool.py -M FitDiagnostics --robustFit 1 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --cminDefaultMinimizerTolerance 0.05 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy 1 --rMin=0 --rMax=2 -m 91 ztt_${bin}_${era}_simple.txt  -n _${bin}_${era}_simple -v 3

# fit only low MT region
combineTool.py -M FitDiagnostics --robustFit 1 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --cminDefaultMinimizerTolerance 0.05 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy 1 --rMin=0 --rMax=2 -m 91 ztt_lowMT_${bin}_${era}_simple.txt  -n _lowMT_${bin}_${era}_simple -v 3

cd -
