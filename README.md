# HtoAA

# Setup software

```
export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USER}/.cmsgit-cache
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_25
cd CMSSW_10_2_25/src
cmsenv
git clone https://github.com/CMS-HTT/LeptonEff-interface HTT-utilities
cd HTT-utilities/LepEffInterface/
rm -rf data
git clone https://github.com/CMS-HTT/LeptonEfficiencies.git data
cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/CorrectionsWorkspace HTT-utilities/CorrectionsWorkspace
cd ${CMSSW_BASE}/src/HTT-utilities/CorrectionsWorkspace
root -l -q CrystalBallEfficiency.cxx++
cd ${CMSSW_BASE}/src
git clone https://github.com/raspereza/HtoAA.git HtoAA
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.2.0
scramv1 b clean; scramv1 b -j 8
cd ${CMSSW_BASE}/src
scramv1 b -j 8
```

# Preparing filelists

First, you have to prepare filelists to run on different samples:
```
cd ${CMSSW_BASE}/src/HtoAA/4Tau/test/2018
./FileListMaker2018.sh
```

# Running analysis macro

The source code of the analysis macro is 
```
${CMSSW_BASE}/src/HtoAA/4Tau/bin/analysis_macro_4tau.cpp 
```

The macro is executed with two arguments: configuration file and the name of filelist.
```
analysis_macro_4tau [config] [filelist]
```

All config files and scripts to run analysis macro on the NAF batch system are found in the folder 
```
${CMSSW_BASE}/src/HtoAA/4Tau/test/2018
```

Running analysis macro on data (example)
```
analysis_macro_4tau analysisMacro_2018.conf DoubleMuon_Run2018B
```

Running on background MC samples (example)
```
analysis_macro_4tau analysisMacro_mc_2018.conf QCD_Pt-20to30_MuEnrichedPt5
```

Running on ggH H->aa->4tau signal samples (example)
```
analysis_macro_4tau analysisMacro_ggH_2018.conf SUSYGluGluToHToAA_AToTauTau_M-125_M-10
```

Running on VBF H->aa->4tau signal samples (example)
```
analysis_macro_4tau analysisMacro_VBF_2018.conf SUSYVBFToHToAA_AToTauTau_M-125_M-10
```

Running on VH H->aa->4tau signal samples (example)
```
analysis_macro_4tau analysisMacro_VH_2018.conf SUSYVH_HToAA_AToTauTau_M-125_M-10
```

Running on ttH H->aa->4tau signal samples (example)
```
analysis_macro_4tau analysisMacro_ttH_2018.conf SUSYttH_HToAA_AToTauTau_M-125_M-10
```

Running on ggH H->aa->2mu+2tau signal samples (example)
```
analysisMacro_ggH_2mu2tau_2018.conf SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-10
```

Note that different config files for different production modes of the signal contain respective names of the Higgs pT reweighting files.

Running analysis macro as given in these examples will produce output RooT file with histograms for further analysis and plotting. The output RooT file is named
```
[filelist].root
```