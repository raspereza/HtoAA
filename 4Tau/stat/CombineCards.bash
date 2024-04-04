#!/bin/bash

for mass in {4..15}
#for mass in 10
do
    combineCards.py haa_2016-13TeV_ma${mass}.txt haa_2017-13TeV_ma${mass}.txt haa_2018-13TeV_ma${mass}.txt > haa_Run2_ma${mass}.txt
    combineTool.py -M T2W -o "haa_Run2_ma${mass}.root" -i haa_Run2_ma${mass}.txt
#    combineTool.py -M T2W -o "haa_2016_preVFP_ma${mass}.root" haa_2016_preVFP-13TeV_ma${mass}.txt
#    combineTool.py -M T2W -o "haa_2016_postVFP_ma${mass}.root" haa_2016_postVFP-13TeV_ma${mass}.txt
    combineTool.py -M T2W -o "haa_2016_ma${mass}.root" -i haa_2016-13TeV_ma${mass}.txt
    combineTool.py -M T2W -o "haa_2017_ma${mass}.root" -i haa_2017-13TeV_ma${mass}.txt
    combineTool.py -M T2W -o "haa_2018_ma${mass}.root" -i haa_2018-13TeV_ma${mass}.txt    
done
