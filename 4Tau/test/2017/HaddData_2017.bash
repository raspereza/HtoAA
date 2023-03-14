#!/bin/bash
for j in B C D E F
do
    ./hadd.sh DoubleMuon_Run2017${j}
done
hadd DoubleMuon_Run2017.root DoubleMuon_Run2017B.root DoubleMuon_Run2017C.root DoubleMuon_Run2017D.root DoubleMuon_Run2017E.root DoubleMuon_Run2017F.root 
