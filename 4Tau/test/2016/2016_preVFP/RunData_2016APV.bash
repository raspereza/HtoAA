#!/bin/bash

for i in B-ver1 B-ver2 C D E F
do
    echo submitting jobs for sample DoubleMuon_Run2018${i}
    ./HTC_submit_seq.sh analysis_macro_4tau analysisMacro_2016APV.conf DoubleMuon_Run2016APV${i} 100
done
