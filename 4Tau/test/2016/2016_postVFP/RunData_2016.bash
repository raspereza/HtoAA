#!/bin/bash

for i in F G H
do
    echo submitting jobs for sample DoubleMuon_Run2016${i}
    ./HTC_submit_seq.sh analysis_macro_4tau analysisMacro_2016.conf DoubleMuon_Run2016${i} 100
done
