#!/bin/bash

for i in B C D E F
do
    echo submitting jobs for sample DoubleMuon_Run2017${i}
    ./HTC_submit_seq.sh analysis_macro_4tau analysisMacro_2017.conf DoubleMuon_Run2017${i} 100
done
