#!/bin/bash
folder=/nfs/dust/cms/user/rasp/Run/HtoAA/stat
cd impacts_${1}_${2}
combineTool.py -M Impacts -d ${folder}/haa_${1}_ma${2}.root -m ${2} -o impacts_${1}_${2}.json
plotImpacts.py -i impacts_${1}_${2}.json -o impacts_${1}_${2} --blind
cd -
