#!/bin/bash
sample=$1
mass=$2

folder=/nfs/dust/cms/user/rasp/Run/HtoAA/stat
combine -M FitDiagnostics --robustFit 1 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --cminDefaultMinimizerTolerance 0.05 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy=1 --rMin=-0.5 --rMax=0.5 -m ${mass} -d ${folder}/haa_${sample}_ma${mass}.root -v2
