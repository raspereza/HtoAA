# Measurement of the track id and Isolation scale factors

## Setting up environment 

Selection of events for the measuremenr is implemented in macro [analysis_macro_ztt.cpp]()

Copy script SetupGrid.bash to the working directory where the code will be run.
Execute this script with one argument - ERA (ERA = 2016_preVFP, 2016_postVFP, 2017 or 2018), for example
```
./SetupGrid.bash 2018
```

The script copies into working directory all shell scripts necessary to run grid control from $CMSSW_BASE/src/HtoAA/4Tau/ztt. From $CMSSW_BASE/src/HtoAA/4Tau/ztt/$ERA it copies filelist maker FileListMaker$ERA.sh and configuration file analysisMacro_ztt.conf to run [analysis_macro_ztt.cpp](https://github.com/raspereza/HtoAA/blob/main/4Tau/bin/analysis_macro_ztt.cpp). Finally, it executes script FileListMaker$ERA.sh, which not only creates filelists, but also splits them into smaller sublists (each such sublist will be process by one grid-control job), and put them in the folder with the name of the sample and appended postfix `_files`. Job splitting per sample is done by the script [split_filelist.sh](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/split_filelist.sh).

The list of all subsamples to be processed is saved into the file named `parameters.txt`. This file also specifies configuration file to be used by the analysis macro (`analysisMacro_ztt.conf`). The same configuration file is used to process all different samples (data and MC) of a given era. The macro [analysis_macro_ztt.cpp](https://github.com/raspereza/HtoAA/blob/main/4Tau/bin/analysis_macro_ztt.cpp) identifies sample by the name of filelist and takes appropriating actions: selects good run-lumi sections from json for data, applies PU reweighting and muon ID, isolation and trigger scale factors for MC, performs Z mass/pT reweighting for simulated DY MC sample and selects same-sign muon-track pairs when processing RooT files contained in filelists with postfix `_SameSign`. Files `parameters.txt` and script [run_synchntuples.sh](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/run_synchntuples.sh) instruct grid-control what executable to call (`analysis_macro_ztt`), which files to process and what configration to use. The execution of the macro is controled with 
(`analysisMacro_ztt.conf`). 
 
The execution of grid-control is steered by configuration file [gc_synch.conf](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/gc_synch.conf) that need to be modified before running grid-control.
You have to specify correctly the working directory by changing [this parameter](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/gc_synch.conf#L11) and set the [full path to your CMSSW project area](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/gc_synch.conf#L37)

Once [gc_synch.conf](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/gc_synch.conf) file is appropriately modified you can proceed with running macro using grid-control. 

## Running grid-control and collecting output

It is recommended to run grid-control in the screen session. To run grid-control you need to have valid grid certificate registered at the CMS VO and in the dcms group. 

Proceed as follows
```
screen
cd $CMSSW_BASE/src
cmsenv 
setenv X509_USER_PROXY /afs/desy.de/user/${u}/${user}/public/k5-ca-proxy.pem
voms-proxy-init -voms cms:/cms/dcms -valid 96:00
```

Change to the working direcory where grid-control environment is setup and start grid-control
```
cd ${workdir}
/nfs/dust/cms/user/rasp/grid-control/go.py gc_synch.conf -iGc
```
This will run monitoring dashboard in the textual mode. It shows the progress with job
submission and running. You can leave screen session by typing `Ctrl A-D` and resume it
with the command `screen -r`.  

After all jobs are finished one collects output for all sample in one go with the script 
[add_samples.sh](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/add_samples.sh) run
in the working directory. This will merge output RooT files individualy for each sample
in loop. To merge output RooT files for one sample use script [hadd.sh](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/hadd.sh) with a name of the sample as an argument, for example
```
./hadd.sh DYJetsToLL_M-50
```

## Creation of datacards for the measurement of track ID and isolation scale factors 

Analysis macro [analysis_macro_ztt.cpp](https://github.com/raspereza/HtoAA/blob/main/4Tau/bin/analysis_macro_ztt.cpp) selects muon+track candidates, fulfilling set of loose requirements. Muon must meet nominal ID criteria and loose relative isolation criterion, relIso<0.4. The 1-prong tau candidate is required to meet all criteria of the nominal analysis. Both same-sign and opposite sign pairs of muon and 1-prong tau candidates are selected. The sample of same-sign pairs are used to estimate QCD background. No cut on transverse mass of muon and MET is applied in the selection. The macro `analysis_macro_ztt.cpp` fills tuple named `mutrkTree` with the [set of variables](https://github.com/raspereza/HtoAA/blob/main/4Tau/bin/analysis_macro_ztt.cpp#L441-L536) used to create datacards for the measurements. 

The datacards are created by macro (createCardsZtt.cpp)[https://github.com/raspereza/HtoAA/blob/main/4Tau/bin/createCardsZtt.cpp]. It is executed with two arguments, configuration file [Cards.conf](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/Cards.conf) and era (2016, 2017 or 2018; the code compines 2016_preVFP and 2016_post datasets in 2016 dataset). Example of running:
```
createCardsZtt Cards.conf 2017
``` 

The content of the default configuration file looks like this
```
# configuration file to create datacards
# for measurement of trk id/iso with Z->tautau
InputFolder = /nfs/dust/cms/user/rasp/Run/MuTrk/PuppiMET
DatacardsFolder = /nfs/dust/cms/user/rasp/HtoAA/datacards/TrkID/PuppiMET
PlotFolder = /nfs/dust/cms/user/rasp/HtoAA/datacards/TrkID/PuppiMET/plots
NBins = 14
XMin = 20
XMax = 160
NBinsZmm = 1
XMinZmm = 76
XMaxZmm = 106
CutLowMT = 30
CutHighMT = 60
DeltaPhiCut = 2.0
TrkPtCutForCones = 10.0
UseZMuMu = true
FloatFakes = true 
```

The meaning of the configuration parameters is explained below
* `InputFolder` : the folder where data and MC RooT files with tuple are stored
* `DatacardsFolder` : the folder where datacards for measurements are saved
* `PlotFolder` : the folder where control plots are saved
* `NBins` : number of bins in the fitted distributions of the muon-track invariant mass
* `XMin` : lower boundary of the fitted distributions of the muon-track invariant mass
* `XMax` : upper boundary of the fitted	distributions of the muon-track invariant mass
* `NBinsZmm` : number of bins in the distribution of dimuon invariant mass in the Z->mumu control region : since we are interested primarily in normalization of the DY process, it is recommeded to set this parameter to 1
* `XMinZmm` : lower boundary of the dimuon invariant mass distribution in the Z->mumu control region
* `XMaxZmm` : upper boundary of the dimuon invariant mass distribution in the Z->mumu control region 
* `CutLowMT` : upper cut on mT(muon,MET), defining `lowMT` measurement region enriched in Z->tautau events
* `CutHighMT` : lower cut on mT(muon,MET), defining `highMT` control region t fake background
* `DeltaPhiCut` : lower cut on deltaPhi(muon,track)
* `TrkPtCutForCones` : lower cut on track pT in measurements for random isolation cones
* `UseZMuMu` : if true, Z->mumu control region is included in the fit and rate parameter, controlling normalisation of the DY process, is introduced; otherwise Z->mumu control region is removed from the measurement.
* `FloatFakes` : if true, normalisation of jet->tau fake background is floated freely in the fit, otherwise log-normal uncertainty of 30% is assigned to this background.   

Before running the code, don't forget to properly set paths to the folder where datacards will be stored (argument `DatacardsFolder`) and where plots will be saved (argument `PlotFolder`).

Measurement of the track ID/Iso scale factors is performed in three (two) regions:
* `Z->mumu` control region (excluded if `UseZMuMu = false`)
* `lowMT` region (measurement region enriched in Z->tautau signal) and
* `highMT` region (control region enriched in jet->tau fakes).

The macro loops over measurement bins and creates datacards for `lowMT` and `highMT` regions. 

Measurement bins include:
* bins of track pT : 5to10, 10to15, 5to15, 15to20, 20toInf (bin 5to15 is introduced as we might consider merging two ranges 5to10 and 10to15 into one);
* bins of isolation axis orientation : cone0, cone0p15, cone0p30, cone0p45.

Additionally, datacards are created for the Z->mumu control region (if `UseZMuMu = true`).

QCD multijet background is extrapolated from `lowMT` (`highMT`) sideband region with same-sign muon-track pairs. Extrapolation factors from same-sign region to the opposite sign region are determined with sample of events, where muon is loosely isolated (0.2<relIso<0.4). RooT files with shape and ascii files with cards for statistical analysis are saved in the folder specified by argument `DatacardsFolder` under the names:
* `zmm_$era.txt(root)` : datacards and shapes in the Z->mumu control region;
* `ztt_lowMT_$bin_$era.txt(root)` : datacards and shapes in the `lowMT` measurement region;
* `ztt_highMT_$bin_$era.txt(root)` : datacards and shapes in the `highMT` control region;
Here `$era = 2016, 2017, 2018` and `$bin` is one of aforementioned measurement bins.

The produced datacards are used to measure track ID/Iso scale factors.
Bash script [RunFit.bash](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/RunFit.bash) illustrates execution of fits using [Higgs Combination toolkit](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit). Script takes three arguments as an input:
```
./RunFit.bash $era $bin $option
```
where
* `$era = 2016, 2017, 2018` : data-taking period;
* `$bin` : one of the measurement bins (see above);
* `$option = robustFit or robustHesse` : fit option of combine utility; the second option provides more reliable estimate of the fitted parameters. 

The output of the fit is saved in the same folder, where datacards are stored, under the name
```
fitDiagnostics_$bin_$era_$option.root
``` 
The output file contains results of the fit, including computed value of track ID/Iso scale factor (nuisance parameter named `r`), as well as prefit and postfit shapes which can be used to produce nice looking plots. 
