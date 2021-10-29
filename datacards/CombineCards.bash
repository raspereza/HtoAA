#!/bin/bash

for mass in {4..15}
do
    combineCards.py haa-13TeV_ma${mass}.txt haa_2018-13TeV_ma${mass}.txt > haa-Run2_ma${mass}.txt
done
