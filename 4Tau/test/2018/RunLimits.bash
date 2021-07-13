#!/bin/bash

for i in {4..15}
do
    combine -M AsymptoticLimits haa-13TeV_ma${i}.txt --noFitAsimov -m ${i}
done 
