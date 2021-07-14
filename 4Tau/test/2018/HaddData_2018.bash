#!/bin/bash
for j in A B C D
do
    ./hadd.sh DoubleMuon_Run2018${j}
done
hadd DoubleMuon_Run2018.root DoubleMuon_Run2018A.root DoubleMuon_Run2018B.root DoubleMuon_Run2018C.root DoubleMuon_Run2018D.root
