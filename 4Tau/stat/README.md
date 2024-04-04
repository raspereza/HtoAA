# Statistical inference in the HtoAA analysis

## Before running the code
It is suggested to run statistical code in a separate folder.
```
mkdir /nfs/dust/cms/user/username/HtoAA
cd /nfs/dust/cms/user/username/HtoAA
cp $CMSSW_BASE/src/HtoAA/4Tau/stat/* ./
```

Make sure that the outputs of [analysis_macro_4tau.cpp](https://github.com/raspereza/HtoAA/blob/main/4Tau/bin/analysis_macro_4tau.cpp) macro are located in separate subfolders under the same base directory
```
${base_directory}/${year}
```

where `${year} = {2016_preVFP, 2016_postVFP, 2017, 2018}`.

The script [CreateCards.C](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/CreateAllCards.C) will look for inputs in these folders. In this script specify correct name of the [base directory](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/CreateCards.C#L314).

## Creating datacards and workspaces

Datacards are created with RooT macro [CreateAllCards.C](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/CreateAllCards.C) which takes four input arguments.
```
#include "CreateCards.C"
void CreateAllCards(
bool Azimov = false, // create Asimov dataset
bool correlation = true, // use correlation coefficients C(i,j)
bool MassUncertPerEras = true, // decorrelate unc. in mass pdfs across eras
bool MassUncertPerBins = true // decorrelate unc. in 1D mass pdfs across bins
) {
```

Input parameters are: 
* `Azimov (bool)` : if true, background-only Asimov dataset is created, otherwise real data are used. 
* `correlation (bool)` : if true, background model is built with mass correlation coefficients according to the formula `b(i,j)=C(i,j)f1D(i)f1D(j)`, otherwise possible correlations are ignored `b(i,j)=f1D(i)f1D(j)`. 
* `(bool) MassUncertPerEra` : if true, uncertainies in C(i,j) and f1D(i) are decorrelated across data taking periods, otherwise uncertainties are correlated across data taking periods. Recommended option is `true`.
* `(bool) MassUncertPerBins` : if true, uncertainties in f1D(i) are decorrelated across bins, e.i. variation of f1D in a given bin is controlled by a separate nuisance parameter. If false, variations in different bins are floated in fits in a correlated way and controlled by common nuisance parameter. Recommended option is `true`. 

Datacards of different data taking periods are combined by executing script [CombineCards.bash](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/CombineCards.bash). Additionally, this script creates workspaces which are used as inputs to the `combine` utility of the [Higgs combination package](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit). Creation of workspaces may take a while. Be patient. At the end you will have one workspace per data taking period and mass point. They are named `haa_${year}_ma${mass}.root` for individual data taking periods and `haa_Run2_ma${mass}.root` for Run 2 combination.

## Running and plotting limits

Limits are run with the script [RunLimits.bash](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/RunLimits.bash). It takes as an argument data taking period (or `Run2` for combination), e.g.
```
./RunLimits.bash 2018
```

or
```
./RunLimits.bash Run2
```

The script produces RooT files (one RooT file per mass point) with observed and expected 95% CL upper limits, as well as 68% and 95% CL intervals around expected limits. The filelist of RooT files is created and named `limits_Run2.txt` or `limits_${year}.txt`, depending on the argument you passed to the [RunLimits.bash](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/RunLimits.bash) script. This filelist is then used by macro [PlotLimits.C](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/PlotLimits.C) 
```
void PlotLimits(
char* fileList = "limits_Run2.txt",
bool blindData = false) {
```
  
## Running Goodness-of-Fit test   

Goodness-of-fit test is performed by executing macro [RunGoF.bash](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/RunGoF.bash) with two parameters : data taking period (or `Run2`) and mass point:
```
./RunGoF.bash Run2 10
```

The script will create folder `GoF_Run2_10` and put all output files there. First, observed value of the test-statistics (saturated likelihood), is computed and stored in the RooT file. Second, jobs to generate toys of test-statistics are submitted to condor. One can configure the number of jobs and number of toys per job by editing the script [RunGoF.bash](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/RunGoF.bash):
```
#!/bin/bash
################################################################
#    definition of parameters to steer running of GoF tests    #
################################################################
sample=$1        # options : 2016 2017 2018 Run2
ma=$2            # mass hypothesis
algo=saturated   # test-statistics, options : saturated, KS, AD 
njobs=25         # number of jobs 
ntoys=40         # number of toys per job

``` 

Wait until jobs are finished and collect the results with script [Hadd_GoF.bash](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/Hadd_GoF.bash). You have to pass as arguments year (or `Run2`) and mass point, for instance
```
./Hadd_GoF.bash Run2 10   
```

Afterwards you can plot results of GoF test and compute p-value, quantifying consistency of your data with the model. Usually, in the search for low yield signal, p-value of GoF test is driven by the quality of background modeling, and one expect comparable values of GoF test for different mass hypothesis. The plotting of GoF test is done with RooT macro [Compatibility.C](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/Compatibility.C)
```
void Compatibility(
TString folder = "GoF_2018_8", // folder with RooT files (output of GoF test)
TString Algo = "saturated", // algorithm
TString legend = "H(125)#rightarrowaa#rightarrow4#tau (2018)",
int bins = 50 // number of bins in the histogram of toys
) {
```

where the first argument is the folder where results of GoF test are stored.

## Running fit and plotting prefit and postfit distributions

The maximum-likelihood fit is done by executing the script [RunFit_Data.bash](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/RunFit_Data.bash) with two input arguments: data taking period (or `Run2`) and mass point:
```
./RunFit_Data.bash Run2 10
```

The script will produce RooT file fitDiagnosticsTest.root containing results of the fit as well as prefit shapes of the background and signal model, postfit shape of the background model after background-only fit ignoring signal, and postfit shapes of background and given signal model after background+signal fit. This RooT file is used as an input to the macro [PlotMass2D.C](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/PlotMass2D.C) which plots prefit and postfit unrolled (m1,m2) distribution:

```
void PlotMass2D(
TString era = "2016", // 2016, 2017, 2018 or Run2
bool prefit = false, // prefit (or postfit) distributions
bool blindData = false, // blind data
bool drawLeg = true, // draw legend
bool logY = true // use log scale for Y axis
) 
```

## Running impacts and pulls

Pulls and impacts of nuisance parameters are computed with the script 
[RunImpacts.bash](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/RunImpacts.bash), which requires two input parameters: data taking period (or `Run2`) and mass point, for exampl
```
./RunImpacts.bash Run2 10
```

or
```
./RunImpacts.bash 2017 4
```

The routine performs initial fit of data and computes signal strength and its uncertainty and then for each nuisance parameter likelihood scan is performed. The task is parallelized by submitting jobs to condor, with each job processing two nuisances. The output is stored in the folder `impacts_${era}_${mass}`. Once jobs are finished, you can create pdf file with plot showing pulls and impacts of nuisance parameters with the script [PlotImpacts.bash](https://github.com/raspereza/HtoAA/blob/main/4Tau/stat/PlotImpacts.bash). And you have to specify data taking period (or `Run2`) and mass point for which you ran impacts. With examples above one should execute
```
./PlotImpacts.bash Run2 10
```

or
```
./RunImpacts.bash 2017 4
```

The pdf file with impact plot `impacts_${era}_${mass}.pdf` will appear in folder the `impacts_${era}_${mass}`.    

## Useful links

* [Higgs combination package](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit)
* [CombineHarvester toolkit](https://cms-analysis.github.io/CombineHarvester/index.html) 




 