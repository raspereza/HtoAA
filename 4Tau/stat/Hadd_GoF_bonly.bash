#!/bin/bash

cd GoF_${1}_bonly
mv higgsCombine.obs.GoodnessOfFit.mH4.root gof_obs.root
hadd gof_exp.root higgsCombine.exp.GoodnessOfFit.mH*.root
cd -
