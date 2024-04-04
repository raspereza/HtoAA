#!/bin/bash

cd GoF_${1}_${2}
mv higgsCombine.obs.GoodnessOfFit.mH${2}.root gof_obs.root
hadd gof_exp.root higgsCombine.exp.GoodnessOfFit.mH*.root
cd -
