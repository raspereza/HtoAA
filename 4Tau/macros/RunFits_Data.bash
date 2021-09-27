#!/bin/bash
# $1 : mass
combine -M FitDiagnostics --plots --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL  --robustFit 1 haa_2018-13TeV_ma${1}.txt --rMin=-10 --rMax=10 -m ${1} 
mv fitDiagnosticsTest.root fitDiagnostics_ma${1}.root

