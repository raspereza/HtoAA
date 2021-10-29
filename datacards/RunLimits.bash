#!/bin/bash
# $1 : common datacards name
#      options: haa-13TeV (2016)
#               haa_2018-13TeV (2018)
#               haa-Run2 (2016+2018)  


rm limits_${1}.txt
cat > limits_${1}.txt <<EOF 
EOF

for i in {4..15}
do
    combine -M AsymptoticLimits ${1}_ma${i}.txt -m ${i}
    mv higgsCombineTest.AsymptoticLimits.mH${i}.root ${1}_limits_mH${i}.root
    cat >> limits_${1}.txt <<EOF1
${PWD}/${1}_limits_mH${i}.root
EOF1
done 
