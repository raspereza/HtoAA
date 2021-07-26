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
analysis_macro_4tau analysisMacro_ggH_2mu2tau_2018.conf SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-10
```

Note that different config files for different production modes of the signal contain respective names of the Higgs pT reweighting files.

Running analysis macro as given in these examples will produce output RooT file with histograms for further analysis and plotting. The output RooT file is given the name
```
[filelist].root
```

# Running analysis macro on NAF batch

The analysis macro can be run on NAF batch system using the script HTC_submit_seq.sh. The script is executed with four arguments: the name of executable (in our case this is analysis_macro_4tau), the name of config file, the name of filelist and the number of input RooT files per job. Example of running the script on data sample DoubleMuon_Run2018B is given below
```
HTC_submit_seq.sh analysis_macro_4tau analysisMacro_2018.conf DoubleMuon_Run2018B 10
```
The script will split input filelist into several parts given specified value of the number of files per job (10), submit to NAF batch system one job for each part, create folder [filelist]_files (in example given above DoubleMuon_Run2018B_files) and put all output RooT files into this folder. You can monitor running of jobs on NAF system using utility
```
condor_q
```

Upon completion of jobs, make sure that all of them successfully finished. To check this execute the script resubmitJobs.sh with [filelist] as the only argument, e.g.
```
./resubmitJobs.sh DoubleMuon_Run2018B
```

This will check if all RooT files in the folder [filelist]_files are properly closed. If any RooT file is not closed properly, the script will resubmit the corresponding job.

Once all jobs successfully executed you can merge RooT files stored in the folder [filelist]_files by executing
```
./hadd.sh [filelist]
``` 
  
This will create RooT file named [filelist].root. 

# Useful scripts

```
RunData_2018.bash - submits jobs to NAF batch for all datasets
RunMC_2018.bash - submit jobs to NAF batch for all MC background samples
RunSignal_2018.bash - submit jobs to NAF batch for all signal samples
HaddData_2018.bash - merges RooT files for all datasets
HaddMC_2018.bash - merges RooT files for all background MC datasets
HaddSigal_2018.bash - merges RooT files for all signal MC datasets
resubmitAllJobs_2018.bash - examines all folders ([filelist]_files) for presence of missing, invalid and incomplete RooT files and resubmit corresponding jobs
```

# Study of trigger acceptance for the signal samples
 
One of important tasks is to evaluate the signal acceptance for various HLT paths. Presently, our primary ntuples contain information on the following interesting HLT paths:
```
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v
HLT_TripleMu_12_10_5_v
```
Unfortunately, no information is available for other interesting HLT paths:
```
HLT_TripleMu_10_5_5_DZ
HLT_TripleMu_5_3_3_Mass3p8to60_DZ
HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx
```
Information on this HLT paths will be added when processing UL samples.


In order to study signal acceptance for this HLT paths, modify config files.
Initial steering parameters
```
ApplyTriggerMatch = true
DiMuonTriggerName = HLT_Mu18_Mu9_SameSign_v
```
have to be changed "false" and the HLT path name to be studied. Example:
```
ApplyTriggerMatch = false
DiMuonTriggerName = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v
```

Run analysis macro with these modified config files on the signal samples. Once output RooT files are created you can compute signal acceptance. For that use histograms named 
```
histWeightsH - number of weighted generated events for a given sample
counter_FinalEventsH - number of weighted events selected into final sample
``` 
The acceptance can then be computed as
```
double acceptance = counter_FinalEventsH->GetSumOfWeights()/histWeightsH->GetSumOfWeights();
```

# Plotting variables at the generator level

The variables such as deltaR(mu,track) pT(track) can be plotted at the generator level for various masses of the pseudoscalar. The plotting macro is
```
PlotDistributions.C

#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"

//TH1D * deltaRMuonPionH = new TH1D("deltaRMuonPionH","",200,0,2);
//TH1D * pionPtH = new TH1D("pionPtH","",100,0,100);

void PlotDistributions(
TString histName = "deltaRMuonPionH", // histogram to plot
TString dir = "./", // folder where ntuples reside
float xLower = 0, // lower boundary of x-axis
float xUpper = 2, // upper boundary of y-axis
TString xtitle = "#delta R(#mu,track)", // title of x-axis
TString ytitle = "normalized to unity"  // title of y-axis
) 
{
...
}

```

Please note that macro rebins initial distribution (with very fine binning) and you will be requested to enter the new number of bins interactively after starting this macro.


