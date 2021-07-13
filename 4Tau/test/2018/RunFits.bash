#!/bin/bash

for i in 4 7 10 15
do
    combine -M FitDiagnostics --plots --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL  --robustFit 1 haa-13TeV_ma${i}.txt -t -1 --rMin=-10 --rMax=10 -m ${i} 
    mv fitDiagnostics.root fitDiagnostics_ma${i}.root
done 
