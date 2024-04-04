#!/bin/bash
################################################################
#    definition of parameters to steer running of GoF tests    #
################################################################
sample=$1        # options : 2016_postVFP-13TeV 2016_preVFP-13TeV 2017-13TeV 2018-13TeV Run2
ma=$2            # mass hypothesis
algo=saturated   # test-statistics, options saturated, KS, AD 
njobs=25         # number of jobs 
ntoys=40         # number of toys per job

npar=$#
if [ $npar -ne 2 ]; then
    echo 
    echo Execute script with two parameters:
    echo 2016_postVFP-13TeV, 2016_preVFP-13TeV, 2017-13TeV, 2018-13TeV, Run2
    echo ma
    echo Examples :
    echo ./RunGoF.bash datacards Run2 5
    echo ./RunGoF.bash datacards 2017-13TeV 6
    echo
    exit
fi


folder=$PWD
outdir=${folder}/GoF_${sample}_${ma} # output folder

if [ ! -d "$outdir" ]; then
    echo creating folder ${outdir}
    mkdir ${outdir}
    cd ${outdir}
else
    cd ${outdir}
    rm * # removing old stuff
fi

name=haa_${sample}_ma${ma}
combineTool.py -M GoodnessOfFit -d ${folder}/${name}.root -m ${ma} --algo ${algo} -n .obs
i=0
while [ $i -lt ${njobs} ]
do 
    random=$RANDOM
    echo running job $i with random seed $random
    combineTool.py -M GoodnessOfFit -d ${folder}/${name}.root --toysFreq -m ${ma} --algo ${algo} -n .exp -t ${ntoys} -s ${random} --job-mode condor --task-name exp.${random} --sub-opts='+JobFlavour = "workday"' 
    i=`expr $i + 1`
done
cd -
