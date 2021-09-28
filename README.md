# HtoAA

# Setting up software

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


# Plotting variables at the generator level

The variables such as deltaR(mu,track) pT(track) can be plotted at the generator level for various masses of the pseudoscalar. The plotting macro is
```
${CMSSW_BASE}/src/HtoAA/4Tau/macros/PlotDistributions.C

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

# Plotting distributions at the reco level 

To plot various distributions saved as TH1 objects in the RooT files created by analysis executable, 
use the macro 
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/Plot.C

#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"

// these histograms are filled after applying same-sign dimuon selection
// no further cuts are applied
//
// muon kinematics ->
//TH1D * ptLeadingMuH = new TH1D("ptLeadingMuH","",50,0,100);
//TH1D * ptTrailingMuH = new TH1D("ptTrailingMuH","",50,0,100);
//TH1D * etaLeadingMuH = new TH1D("etaLeadingMuH","",50,-2.5,2.5);
//TH1D * etaTrailingMuH = new TH1D("etaTrailingMuH","",50,-2.5,2.5);
//TH1D * dimuonMassH = new TH1D("dimuonMassH","",500,0,500);
//
// number of tracks within dR<0.5 around muons
// we select events where there is only one track accompanies each of muons   
//TH1D * nTracksLeadingMuH = new TH1D("nTracksLeadingMuH","",21,-0.5,20.5);
//TH1D * nTracksTrailingMuH = new TH1D("nTracksTrailingMuH","",21,-0.5,20.5); 


void Plot(TString histName = "ptLeadingMuH", // histogram name
	  TString xtitle = "p_{T}^{#mu1} [GeV]", // title of x axis
	  TString ytitle = "Events / 2 GeV", // title of y axis
	  TString Mass = "10", // ma [GeV]
	  float xLower = 18, // lower boundary of x axis 
	  float xUpper = 100, // upper boundary of x axis
	  float yLower = 1000, // lower boundary of y axis (in case when plotting y axis in log scale)
	  int nBinsNew = 50, // new number of bins
	  bool logY = true, // log or linear scale of Y axis
	  bool drawLeg = true) { // draw or not legend
...
}

```
The macro produces file `$histName.png` in the working directory where the script is executed.
The list of histograms containing interesting distributions is given in the beginning of the macro. 
The distributions are obtained after selectring muon pairs of opposite charge. Further selection, in particular isolation requirement imposed on the muon-track pairs (each muon must be accompanied by only one track within dR cone of 0.5 around muon direction) significantly reduces statistics in the QCD multijet MC samples resulting in sparsely populated histograms with large bin-by-bin uncertainties. 
Other distributions filled in the process of signal selection can be found in the analysis macro 
```
$CMSSW_BASE/src/HtoAA/4Tau/bin/analysis_macro_4tau.cpp
```
Please fill free to edit this macro and add histograms of the distributions you are interested in.

# Creating datacards
To perform statistical analysis one needs to create datacards. This includes ascii file encoding uncertainty model and RooT files with the observed data, background and signal distributions. For more information on statistical tools used in the CMS analysis please consult [this documentation](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit). You are adviced also to familiarize with the analysis by reading [dedicated paper](https://arxiv.org/abs/1907.07235v2) on the CMS analysis of 2016 data or [latest presentation](https://indico.cern.ch/event/1053526/contributions/4447626/attachments/2279719/3873270/Haa_2018.pdf) on the strategy of the H->aa search with the 2018 dataset.  

At the last stage of the analysis, the signal is extracted from two-dimensional distribution of the muon-track pairs selected in the final sample. The purpose of the macro, producing datacards, is to create templates of the two-dimensional distributions for data, background and various signal processes as well as to define uncertainty model in the format compliant with the CMS statistical toolkit. 

The macro producing datacards is
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/CreateCards.C

#include "HtoH.h"

void CreateCards(TString mass="4", // mass of pseudoscalar boson
		 bool Azimov = true, // produce b-only Azimov dataset
		 bool correlation = true // apply correlation coefficients in background model
		 ) {

  TString dir = "/nfs/dust/cms/user/rasp/Run/Run2018/H2aa/";

....
}

```

It is steered with three input parameters: the probed pseudoscalar mass hypothesis (TString), and two boolean variables. First boolean (named Asimov) instructs code which type datacards should be produced. If parameter is set to true the background-only Asimov dataset is created, that is the data in each bin of distribution is replaced by background expectation. Otherwise real data are used to create datacards. It is advice to set this variable to true while analysis is kept blinded (we are not biasing ourselves by looking at data in the signal region). The second boolean, when set to true, instructs code to apply correlation mass correlation coefficients for the background model. Otherwise the masses of muon-track pairs in th background are assumed to be uncorrelated. 
The macro uses as inputs histograms in the RooT files produced by analysis 
macro analysis_macro_4tau.cpp. As background is estimated solely from data,
to produce datacards one needs only data and signal RooT files:
- DoubleMuon_Run2016.root 
- SUSYGluGluToHToAA_AToTauTau_M-125_M-${mass}.root
- SUSYVBFToHToAA_AToTauTau_M-125_M-${mass}.root
- SUSYVH_HToAA_AToTauTau_M-125_M-${mass}.root
- SUSYttH_HToAA_AToTauTau_M-125_M-${mass}.root
- SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-${mass}.root

Additionally one needs RooT files containing 2d histogram with mass correlation coefficients derived in the control (sideband) region from data and from the simulated QCD samples and in the signal region from the simulated QCD samples:  
- CorrCoefficients_data.root
- CorrCoefficients_control_mc.root
- CorrCoefficients_signal_mc.root

The first RooT file (correlation coefficients measured in control region in data) is produced by macro
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/CorrCoefficients.C
```
The correlation coefficients derived from simulation are presently inherited from the analysis of the 2016 data. The analysis code needs to be updated to rederive these coefficients for 2018 dataset.

All RooT files used in creation of datacards are located in the folder specified by parameter 
TString dir. By default this variable points to the directory `/nfs/dust/cms/user/rasp/Run/Run2018/H2aa`.
You should modify this parameter accordingly when using another working directory containing all necessary ingredients to produce datacards.  

The datacard producer creates two files named 
- haa_2018-13TeV_ma${mass}.root 
- haa_2018-13TeV_ma${mass}.txt 

These files are then used in the statistical analysis. The RooT file contains TH1 objects, representing 2d mass distributions unrolled in 1d array of bins. The following histograms are stored in this RooT file:
- data_obs: observed data (or background model if Azimov = true)
- bkg: background model
- ggh: gg->H signal with aa->4tau
- vbf: qqH signal with aa->4tau
- vh: ZH+WH signal with aa->4tau
- tth: ttH signal with aa->4tau
- mmtt: ggH signal with aa->2mu2tau (yield is rescaled to the inclusive cross section sigma(ggH)+sigma(qqH)+sigma(VH)+sigma(ttH)) 

The RooT file contains also up/down systematic variations of bkg model related to uncertainty in the estimate of 1D background shape in the distribution of muon-track invariant mass and uncertainty in mass correlation coefficients: 
- bkgd_CMS_unc1d_2018Up(Down)
- bkgd_CMS_uncCorr_2018Up(Down)

The ascii file haa_2018-13TeV_ma${mass}.txt defines the uncertainty model.
With the following line in the ascii file
```
bkgNorm_2018   rateParam  haa_2018  bkgd  1  [0.5,1.5]
```
we instruct the CMS statistical toolkit to perform fit with freely floating overall background normalization. 
You can produce datacards for all tested mass hypotheses (ma=4..15 GeV) at once by executing macro
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/CreateAllCards.C
```

# Running limits
The observed and expected limits are computed for all tested pseudoscalar mass hypotheses at once by running script
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/RunLimits.bash
```
The scripts call `combine` utility of the [CMS Higgs statistical package](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit) and creates RooT files named `higgsCombineTest.AsymptoticLimits.mH$mass$.root` (one file per tested mass hypothesis) where info on median/observed limits is stored. The script creates also filelist named `limits.txt`. The filelist is then used by the macro
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/PlotLimits.C

#include "HttStylesNew.cc"
#include "CMS_lumi.C"

void PlotLimits(char* fileList = "limits.txt",
		bool blindData = true) {

...
}

```

which plots expected median and observed upper limits on the signal strength (mu) defined as `mu=sigma(pp->H+X)*BR(H->aa->4tau)/sigma(pp->H+X,SM)`. Before unblinding the data the input parameter `bool blindData` is recommended to be set to true, which instructs the macro to hide the observed limits. If datacards are produced with data replaced by background-only Azimov dataset, the observed limits should be equal to the median expected ones.

# Running maximum likelihood fits
Maximum likelihood fits can be run either on data or on the Asimov dataset with predefined signal strength. 

To run fits on data execute the script
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/RunFits_Data.bash
```
It accepts one input parameter - tested mass hypothesis, e.g.
```
./RunFits_Data.bash 10
```
The script calls `combine` utility, computes the best-fit value of the signal strength and creates output RooT file named `fitDiagnostics_ma$mass.root`. This RooT file contains results of the fit (the best-fit mu value and fitted values of nuisance parameters), and the prefit and postfit distributions for signal and background. Please note if datacards are produced for background-only Asimov dataset, the fitted value mu and all nuisances would be zero as we the fit is performed on dataset which is obtained from background model. 

To run fits on the signal+background Asimov dataset execute the script
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/RunFits_Asimov.bash
```
You should pass to parameters to the script: pseudoscalar mass and the signal strength mu, e.g.
```
./RunFits_Asimov.bash 10 1
```
The script will create RooT file named `fitDiagnostics_ma$mass_mu$mu.root` containing the fit results and postfit shapes with uncertainties. The best-fit value of mu should be equal to signal strength injected into Asimov dataset (1 in the above example).

# Other useful macros

## Plotting unrolled 2D (m1,m2) signal distributions after final selection

The unrolled 2D (m1,m2) signal distributions are plotted by running the macro
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/PlotSignal.C

#include "CMS_lumi.C"
#include "HttStylesNew.cc"

void PlotSignal(TString mass = "4") {

  TString dir("/nfs/dust/cms/user/rasp/Run/Run2018/H2aa/");

...

}
```
The macro accepts one input parameter (TString mass): the mass of pseudoscalar.
Contributions from H->aa->4tau and H->aa->2mu2tau decays are plotted separately.
The macro uses the datacard RooT file `haa_2018-13TeV_ma$mass.root` 
in the folder defined by variable TString dir 
(the default folder is `/nfs/dust/cms/user/rasp/Run/Run2018/H2aa/`).

## Plotting 1D mass distribution of muon-track pairs after final selection

The 1D distribution of muon-track pairs is plotted with the macro
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/PlotMass1D.C

#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "HtoH.h"

void PlotMass1D(float xLower = 0, // lower boundary in x axis
		float xUpper = 12, // upper boundary in x axis
		TString xtitle = "m_{#mu,trk} [GeV]", // x-axis xtitle
		TString ytitle = "normalised to unity", // y-axis title
		bool drawLeg = true, // draw legend
		bool logY = true, // use log scale for y-axis
		bool blindData = true) { // blind data

  TString dir("/nfs/dust/cms/user/rasp/Run/Run2018/H2aa/");

...

}
```
If input parameter `blindData` is set to true, the data points are not shown for bins with m(mu,trak) > 6 GeV (these are bins which are most sensitive to the signal). 
The macro uses RooT files produced by the analysis executable `analysis_macro_4tau`. The directory, where the RooT files are stored, is defined by parameter `TString dir` (default location is `/nfs/dust/cms/user/rasp/Run/Run2018/H2aa/`).

## Plotting unrolled 2D (m1,m2) distribution (signal, background and data) after final selection

The unrolled 2D (m1,m2) distributions in data, background and signal samples are plotted using the macro
```
$CMSSW_BASE/src/HtoAA/4Tau/macros/PlotMass2D.C

#include "CMS_lumi.C"
#include "HttStylesNew.cc"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"
#include "boost/range/algorithm.hpp"
#include "boost/range/algorithm_ext.hpp"
#include "Plotting.h"
#include "Plotting_Style.h"

void PlotMass2D(
		bool prefit = true, // prefit (or postfit) distributions
		bool blindData = true, // blind data
		bool drawLeg = true, // draw legend
		bool logY = true // use log scale for
		) {
  
  TString dir("/nfs/dust/cms/user/rasp/Run/Run2018/H2aa/");
...

}

```

The macro uses datacard RooT files created for mass points ma=4,7,10,15 (`haa_2018-13TeV_ma$ma.root`). Additionally RooT file `fitDiagnostics_ma10.root` is required. This RooT file is created by executing command 
```
> RunFits_data.bash 10
``` 
All RooT files are supposed to be located in the directory specified by variable `TString dir` 
(default directory is `/nfs/dust/cms/user/rasp/Run/Run2018/H2aa/`).

# Task list

## Study of the signal acceptance for various HTL paths (done)
 
One of the important tasks is to evaluate the signal acceptance for various alternative triggers which can be potentially used in the analysis of 2017 data taken w/o same-sign dimuon HLT. Presently, our ntuples contain information on the following interesting HLT paths:
```
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v
HLT_TripleMu_12_10_5_v
```

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

Please compute acceptance for the studied HLT paths as a function of probed mass hypothesis in
the range [4..15] GeV.

## Optimization of cut on dR(muon,track)

The current version of analysis requires muon-track pair stemming from a->tautau decay to have dR(muon,track)<0.5. This cut may be sub-optimal for some tested pseudoscalar mass hypotheses. Therefore, it is desirable to study the sensitivity of the search for different mass hypothesis in dependence of upper cut on dR(mu,track). This parameter is specified in configuration file steering the running of analysis executable `analysis_macro_4tau`: 
```  
dRIsoMuon = 0.5
```

The workflow of the study would be
- create working directory to run analysis with a given value of `dRIsoMuon`;
- copy content of folder `$CMSSW_BASE/src/HtoAA/4Tau/test/2018` to the working directory;
- create filelists with the script `FileListMaker2018.sh`;
- modify appropriately the parameter `dRIsoMuon` in configuration files for data (`analysisMacro_2018.conf`) and signal samples (`analysisMacro_{ggH,VBF,VH,ttH,ggH_2mu2tau}_2018.conf`);
- run analysis executable `analysis_macro_4tau` on data and signal samples;
- create RooT file with (m1,m2) correlation coefficients by running macro `CorrCoefficients.C` in the working directory;
- create datacards; 
- compute expected limits for all mass hypothesis (ma=4..15 GeV).

Perform study for `dRIsoMuon=0.4,0.6,0.7`. Results for `dRIsoMuon=0.5` 
are already available in the folder `/nfs/dust/cms/user/rasp/Run/Run2018/H2aa`.

## Study correlation between muon+track momentum and dR(muon,track) in signal samples

It is well possible that the sensitivity of the search can be improved over entire range of probed mass hypotheses by using dynamic cut on dR(muon,track), e.g. by applying cut in dependence of pT(mu+track) or p(mu+track), where pT(p) is the transverse (full) momentum of the muon-track system.
To assess potential improvement of sensitivity one has to study the correlation between dR(muon,track) and pT(mu+track) or p(mu+track). The updated version of the analysis macro produces ntuple (named "muTrkTree") with the generator level information to facilitate such studies. Each entry of the ntuple stores kinematic properties of a muon-track pair resulting from the a->tau(mu)tau(1-prong) decay.
The variables stored in ntuple are:
```
genmu_P : generator muon momentum
genmu_Pt : generator muon transverse momentum
genmu_Eta : generator muon pseudorapidiry
gentrk_P : generator track momentum
gentrk_Pt : generator track transverse momentum
gentrk_Eta : generator track pseudorapidiry
genmutrk_P : generator momentum of muon-track system
genmutrk_Pt : generator transverse momentum of muon-track
genmutrk_Eta : generator pseudorapifity of muon-track system
genmutrk_DR : dR(mu,track)
```  

Using this information one can study correlations by plotting 2D distributions (genmutrk_DR,genmutrk_Pt) or (genmutrk_DR,genmutrk_P). Another way to study correlations is to plot genmutrk_DR in bins of genmutrk_Pt(genmutrk_P) (e.g. genmutrk_Pt=[10,15,20,25,100] GeV). 
One has to perform this study for several mass hypotheses, e.g. ma = 4, 7, 10, 12, 15 GeV. 
The results of these studies can be then used to derive an optimal formula for upper cut on dR(mu,track) as a function of pT(p) of the muon-track system.
