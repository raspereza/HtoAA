#!/bin/bash

#File lists for 2018 legacy data:
for i in A B C D
do
    echo creating file list for data sample DoubleMuon_Run2018${i}
    ls /pnfs/desy.de/cms/tier2/store/user/rasp/ntuples/H2AA/2018/data/DoubleMuon_Run2018${i}/*root > DoubleMuon_Run2018${i}
done

#for i in A B C D
#do
#    echo creating file list for data sample SingleMuon_Run2018${i}
#    ls /pnfs/desy.de/cms/tier2/store/user/rasp/ntuples/H2AA/2018/data/SingleMuon_Run2018${i}/*root > SingleMuon_Run2018${i}
#done

#File lists for 2018 MC ntuples v2:
echo creating file lists for HToAA_AToMuMu_AToTauTau samples
for i in 3p6 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 
do
    ls /nfs/dust/cms/user/consuegs/ntuples/NMSSM_2018_v2/'SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-'$i''/*root > SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-$i
done

echo creating file lists for HToAA_To4Tau samples
for i in {4..21}
do
    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}/*root > SUSYGluGluToHToAA_AToTauTau_M-125_M-${i}
    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/SUSYVBFToHToAA_AToTauTau_M-125_M-${i}/*root > SUSYVBFToHToAA_AToTauTau_M-125_M-${i}
    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/SUSYVH_HToAA_AToTauTau_M-125_M-${i}/*root > SUSYVH_HToAA_AToTauTau_M-125_M-${i}
    ls /pnfs/desy.de/cms/tier2/store/user/sconsueg/ntuples/H2aa_4tau/2018/mc/SUSYttH_HToAA_AToTauTau_M-125_M-${i}/*root > SUSYttH_HToAA_AToTauTau_M-125_M-${i}
done

# TTSemileptonic sample should be splitted (long list for ls command)
echo "Creating file list for sample TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2018/mc/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*0.root > TTToSemiLeptonic
for index in {1..9}
do
    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2018/mc/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*${index}.root >> TTToSemiLeptonic
done

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
DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext
ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_erdON_13TeV-powheg-madspin-pythia8
ST_t-channel_top_4f_InclusiveDecays_TuneCP5_erdON_13TeV-powheg-madspin-pythia8
ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8
ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8
TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8
TTToHadronic_TuneCP5_13TeV-powheg-pythia8
WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
)

names=(DYJetsToLL_M-10to50
DYJetsToLL_M-50
ST_t-channel_top
ST_t-channel_antitop
ST_tW_top
ST_tW_antitop
TTTo2L2Nu
TTToHadronic
WJetsToLNu
)

names_QCD=(QCD_Pt-15To20_MuEnrichedPt5
QCD_Pt-20To30_MuEnrichedPt5  
QCD_Pt-30To50_MuEnrichedPt5
QCD_Pt-50To80_MuEnrichedPt5
QCD_Pt-80To120_MuEnrichedPt5
QCD_Pt-120To170_MuEnrichedPt5
QCD_Pt-170To300_MuEnrichedPt5
QCD_Pt-300To470_MuEnrichedPt5
QCD_Pt-470To600_MuEnrichedPt5
QCD_Pt-600To800_MuEnrichedPt5
QCD_Pt-800To1000_MuEnrichedPt5
QCD_Pt-1000_MuEnrichedPt5
)

i=0
while [ $i -lt ${#samples[@]} ] 
do
    echo "Creating file list for sample" ${samples[$i]} 

    ls /pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2018/mc/${samples[$i]}*/*root > ${names[$i]}
      
    i=`expr $i + 1` 
done

k=0
while [ $k -lt ${#samples_VV[@]} ] 
do
    echo "Creating file list for sample" ${samples_VV[$k]} 

    ls /pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/mc/${samples_VV[$k]}*/*root > ${names_VV[$k]}
      
    k=`expr $k + 1` 
done

j=0
while [ $j -lt ${#names_QCD[@]} ] 
do
    echo "Creating file list for sample" ${names_QCD[$j]} 
    ls /pnfs/desy.de/cms/tier2/store/user/rasp/ntuples/H2AA/2018/mc/${names_QCD[$j]}/*root > ${names_QCD[$j]}
      
    j=`expr $j + 1` 
done

