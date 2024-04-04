#!/bin/bash
# $1 : common datacards name
#      options: haa_2016_preVFP-13TeV (2016_preVFP)
#		haa_2016_postVFP-13TeV (2016_postVFP)
#		haa_2017-13TeV (2017)
#               haa_2018-13TeV (2018)
#               haa-Run2 (2016+2018)  


rm limits_${1}.txt
cat > limits_${1}.txt <<EOF 
EOF

for i in {4..15}
#for i in 8
do
#    combine -M AsymptoticLimits --noFitAsimov --rMin=0 --rMax=1 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.01 ${1}_ma${i}.txt -t -1 -m ${i} 
    combine -M AsymptoticLimits --rMin=0 --rMax=1 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.01 haa_${1}_ma${i}.root -m ${i} 
    mv higgsCombineTest.AsymptoticLimits.mH${i}.root ${1}_limits_mH${i}.root
    cat >> limits_${1}.txt <<EOF1
${PWD}/${1}_limits_mH${i}.root
EOF1
done 
