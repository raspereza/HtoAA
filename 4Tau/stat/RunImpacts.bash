#!/bin/bash
# $1 - channel
ulimit -s unlimited
folder=/nfs/dust/cms/user/rasp/Run/HtoAA/stat

OUTDIR=impacts_${1}_${2}
if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR    
fi
cd $OUTDIR
rm *.root
combineTool.py -M Impacts -d ${folder}/haa_${1}_ma${2}.root -m ${2} --rMin=-1 --rMax=1 --robustFit 1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy=1 --doInitialFit 
#combineTool.py -M Impacts -d ${folder}/haa_${1}_ma${2}.root -m ${2} --rMin=-1 --rMax=1 --robustFit 1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy=1 --job-mode condor --sub-opts='+JobFlavour = "workday"' --merge 2 --doFits
combineTool.py -M Impacts -d ${folder}/haa_${1}_ma${2}.root -m ${2} --rMin=-1 --rMax=1 --robustFit 1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --cminDefaultMinimizerStrategy=1 --job-mode condor --sub-opts='+JobFlavour = "workday"' --merge 2 --doFits
cd -

