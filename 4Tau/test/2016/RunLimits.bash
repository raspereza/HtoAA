#!/bin/bash

rm limits.txt
cat > limits.txt <<EOF 
EOF

for i in {4..15}
do
    combine -M AsymptoticLimits haa_2016-13TeV_ma${i}.txt -m ${i}
    cat >> limits.txt <<EOF1
${PWD}/higgsCombineTest.AsymptoticLimits.mH${i}.root
EOF1
done 
