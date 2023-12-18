# Running analysis_macro_ztt using grid control

## Setting up environment 
Copy script SetupGrid.bash to the working directory where the code will be run.
Execute this script with one argument - ERA (ERA = 2016_preVFP, 2016_postVFP, 2017 or 2018), for example```
./SetupGrid.bash 2018
```
The script copies into working directory all shell scripts necessary to run grid control from ($CMSSW_BASE/src/HtoAA/4Tau/ztt). From `$CMSSW_BASE/src/HtoAA/4Tau/ztt/$ERA` it copies filelist maker `FileListMaker${ERA}.sh and configuration file `analysisMacro_ztt.conf` to run [analysis_macro_ztt.cpp](https://github.com/raspereza/HtoAA/blob/main/4Tau/bin/analysis_macro_ztt.cpp). Finally, it executes script `FileListMaker${ERA}.sh`, which not only creates filelists, but also splits them into smaller sublists (each such sublist will be process by one grid-control job), and put them in the folder with the name of the sample and appended postfix `_files`. Job splitting per sample is done by the script [split_filelist.sh](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/split_filelist.sh).

For each sample a duplicate filelist with the name containing postfix `_SameSign` is created. This is needed for selection of same-sign muon-track pairs. The list of all subsamples to be processed is saved into the file named `parameter.txt`. This file also specifies configuration file to be used by the analysis macro (`analysisMacro_ztt.con`). The same configuration file is used to process all different samples (data and MC) of a given era. The macro [analysis_macro_ztt.cpp](https://github.com/raspereza/HtoAA/blob/main/4Tau/bin/analysis_macro_ztt.cpp) identifies sample by the name of filelist and takes appropriating actions: selects good run-lumi sections from json for data, applies PU reweighting and muon ID, isolation and trigger scale factors for MC, performs Z mass/pT reweighting for simulated DY MC sample and selects same-sign muon-track pairs when processing RooT files contained in filelists with postfix `_SameSign`. Files `parameter.txt` and script (run_synchntuples.sh)[https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/run_synchntuples.sh] instruct grid-control what executable to call (`analysis_macro_ztt`), which files to process and what configration
 
The execution of grid-control is steered by configuration file (gc_synch.conf)[https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/gc_synch.conf] that need to be modified before running grid-control.
You have to specify correctly the working directory by changing (this parameter)[https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/gc_synch.conf#L11] 
```
[storage]
; please modify according to your working directory
se path = /nfs/dust/cms/user/rasp/Run/
```
and set the (full path to your CMSSW project area)[https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/gc_synch.conf#L37]
```
; please modify acoording to your project area
project area      = /nfs/dust/cms/user/rasp/CMSSW/CMSSW_10_2_25
```

Once (gc_synch.conf)[https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/gc_synch.conf] file is appropriately modified you can proceed with running macro using grid-control 

## Running grid-control and collecting output

It is recommended to run grid-control in the screen session. To run grid-control you need to have valid grid certificate registered at the CMS VO and in the dcms group. 

Proceed as follows
```
screen
cd $CMSSW_BASE/src
screen
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
(add_samples.sh)[https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/add_samples.sh] run
in the working directory. This will merge output RooT files individualy for each sample
in loop. To merge output RooT files for one sample use script [hadd.sh](https://github.com/raspereza/HtoAA/blob/main/4Tau/ztt/hadd.sh) with a name of the sample as an argument, for example
```
./hadd.sh DYJetsToLL_M-50
```
 
