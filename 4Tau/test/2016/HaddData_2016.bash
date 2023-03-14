#!/bin/bash
for i in B-ver1 B-ver2 C D E F
do
    ./hadd.sh ./2016_preVFP/DoubleMuon_Run2016APV${i}
done    	
hadd DoubleMuon_Run2016_preVFP.root ./2016APV/DoubleMuon_Run2016APVB-ver1.root ./2016APV/DoubleMuon_Run2016APVB-ver2.root ./2016APV/DoubleMuon_Run2016APVC.root ./2016APV/DoubleMuon_Run2016APVD.root ./2016APV/DoubleMuon_Run2016APVE.root ./2016APV/DoubleMuon_Run2016APVF.root

for j in F G H
do
    ./hadd.sh ./2016_postVFP/DoubleMuon_Run2016${j}
done
hadd DoubleMuon_Run2016_postVFP.root ./2016/DoubleMuon_Run2016F.root ./2016/DoubleMuon_Run2016G.root ./2016/DoubleMuon_Run2016H.root 

hadd DoubleMuon_Run2016.root DoubleMuon_Run2016_preVFP.root DoubleMuon_Run2016_postVFP.root
