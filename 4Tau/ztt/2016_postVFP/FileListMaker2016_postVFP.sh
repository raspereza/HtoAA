#!/bin/bash

if [ -f "parameters.txt" ]; then
    rm parameters.txt
fi

echo "CONFIGFILE,FILELIST" > parameters.txt

#File lists for 2016 legacy data:
for i in F G H
do
    echo creating file list for data sample SingleMuon_Run2016${i}
    ls /pnfs/desy.de/cms/tier2/store/user/lsreelat/NTuples/2016/HtoAA/SingleMuon/SingleMuon-Run2016${i}-UL2016/*root > SingleMuon_Run2016${i}
    ./split_filelist.sh analysisMacro_ztt.conf SingleMuon_Run2016${i} 20

done

# TTSemileptonic sample should be split (long list for ls command)
echo "Creating file list for sample TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2016_postVFP/mc/TTToSemiLeptonic/*0.root > TTToSemiLeptonic
for index in {1..9}
do
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2016_postVFP/mc/TTToSemiLeptonic/*${index}.root >> TTToSemiLeptonic
done
./split_filelist.sh analysisMacro_ztt.conf TTToSemiLeptonic 20

# WJetsToLNu sample should be split (long list for ls command)
echo "Creating file list for sample WJetsToLNu"
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2016_postVFP/mc/WJetsToLNu/*0.root > WJetsToLNu
for index in {1..9}
do
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2016_postVFP/mc/WJetsToLNu/*${index}.root >> WJetsToLNu
done
./split_filelist.sh analysisMacro_ztt.conf WJetsToLNu 20

#File lists for VV background MC samples
samples_VV=(WW_TuneCP5_13TeV-pythia8_2016
WZ_TuneCP5_13TeV-pythia8_2016
ZZ_TuneCP5_13TeV-pythia8_2016
)

names_VV=(WW_13TeV-pythia8
WZ_13TeV-pythia8
ZZ_13TeV-pythia8
)

#File lists for background MC samples
samples=(DYJetsToLL_M-10to50
DYJetsToLL_M-50
ST_t-channel_antitop_4f
ST_t-channel_top_4f
ST_tW_top_5f_inclusiveDecays
ST_tW_antitop_5f_inclusiveDecays
TTTo2L2Nu
TTToHadronic
)

names=(DYJetsToLL_M-10to50
DYJetsToLL_M-50
ST_t-channel_top
ST_t-channel_antitop
ST_tW_top
ST_tW_antitop
TTTo2L2Nu
TTToHadronic
)

i=0
while [ $i -lt ${#samples[@]} ] 
do
    echo "Creating file list for sample" ${samples[$i]} 

    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2016_postVFP/mc/${samples[$i]}*/*root > ${names[$i]}
    ./split_filelist.sh analysisMacro_ztt.conf ${names[$i]} 20
      
    i=`expr $i + 1` 
done
echo "Creating file list for sample DYJetsToTT_M-50"
cp DYJetsToLL_M-50 DYJetsToTT_M-50
./split_filelist.sh analysisMacro_ztt.conf DYJetsToTT_M-50 20

k=0
while [ $k -lt ${#samples_VV[@]} ] 
do
    echo "Creating file list for sample" ${samples_VV[$k]} 

    ls /pnfs/desy.de/cms/tier2/store/user/lsreelat/NTuples/2016/HtoAA/VV_inclusive/${samples_VV[$k]}*/*root > ${names_VV[$k]}
    ./split_filelist.sh analysisMacro_ztt.conf ${names_VV[$k]} 20

    k=`expr $k + 1` 
done

