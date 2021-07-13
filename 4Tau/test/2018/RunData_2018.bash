#!/bin/bash

for i in A B C D
do
    echo submitting jobs for sample DoubleMuon_Run2018${i}
    ./HTC_submit_seq.sh analysis_macro_4tau analysisMacro_2018.conf DoubleMuon_Run2018${i} 100
done
