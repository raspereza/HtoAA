#!/bin/bash

if [ -f "parameters.txt" ]; then
    rm parameters.txt
fi

echo "CONFIGFILE,FILELIST" > parameters.txt

#File lists for 2018 legacy data:
for i in A B C D
do
    echo creating file list for data sample SingleMuon_Run2018${i}
    ls /pnfs/desy.de/cms/tier2/store/user/lsreelat/NTuples/2018/HtoAA/SingleMuon_v2/SingleMuon-Run2018${i}-UL2018/*0.root > SingleMuon_Run2018${i}
    for index in {1..9}
    do
	ls /pnfs/desy.de/cms/tier2/store/user/lsreelat/NTuples/2018/HtoAA/SingleMuon_v2/SingleMuon-Run2018${i}-UL2018/*${index}.root >> SingleMuon_Run2018${i}
    done
    ./split_filelist.sh analysisMacro_ztt.conf SingleMuon_Run2018${i} 20
done

# TTSemileptonic sample should be split (long list for ls command)
echo "Creating file list for sample TTToSemiLeptonic"
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2018/mc/TTToSemiLeptonic/*0.root > TTToSemiLeptonic
for index in {1..9}
do
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2018/mc/TTToSemiLeptonic/*${index}.root >> TTToSemiLeptonic
done
./split_filelist.sh analysisMacro_ztt.conf TTToSemiLeptonic 10

# WJetsToLNu sample should be split (long list for ls command)
echo "Creating file list for sample WJetsToLNu"
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2018/mc/WJetsToLNu/*0.root > WJetsToLNu
for index in {1..9}
do
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2018/mc/WJetsToLNu/*${index}.root >> WJetsToLNu
done
./split_filelist.sh analysisMacro_ztt.conf WJetsToLNu 20

#File lists for VV background MC samples
samples_VV=(WW_TuneCP5_13TeV-pythia8
WZ_TuneCP5_13TeV-pythia8
ZZ_TuneCP5_13TeV-pythia8
)

names_VV=(WW_13TeV-pythia8
WZ_13TeV-pythia8
ZZ_13TeV-pythia8
)

#File lists for background MC samples
samples=(DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8
DYJetsToLL_M-50
DYJetsToLL_M-50_amcatnlo
ST_t-channel_top_4f
ST_t-channel_antitop_4f
ST_tW_top_5f_inclusiveDecays
ST_tW_antitop_5f_inclusiveDecays
TTTo2L2Nu
TTToHadronic
)

names=(DYJetsToLL_M-10to50
DYJetsToLL_M-50
DYJetsToLL_M-50_amcatnlo
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
    echo "Creating file list for sample" ${names[$i]} 

    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2018/mc/${samples[$i]}/*root > ${names[$i]}
    ./split_filelist.sh analysisMacro_ztt.conf ${names[$i]} 20
      
    i=`expr $i + 1` 
done

echo "Creating file list for sample DYJetsToTT_M-50"
cp DYJetsToLL_M-50 DYJetsToTT_M-50
./split_filelist.sh analysisMacro_ztt.conf DYJetsToTT_M-50 20

echo "Creating file list for sample DYJetsToTT_M-50_amcatnlo"
cp DYJetsToLL_M-50_amcatnlo DYJetsToTT_M-50_amcatnlo
./split_filelist.sh analysisMacro_ztt.conf DYJetsToTT_M-50_amcatnlo 20

k=0
while [ $k -lt ${#samples_VV[@]} ] 
do
    echo "Creating file list for sample" ${names_VV[$k]} 

    ls /pnfs/desy.de/cms/tier2/store/user/lsreelat/NTuples/2018/HtoAA/VV_inclusive/${samples_VV[$k]}/*root > ${names_VV[$k]}
    ./split_filelist.sh analysisMacro_ztt.conf ${names_VV[$k]} 20
      
    k=`expr $k + 1` 
done

