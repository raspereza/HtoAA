# Statistical inference in the HtoAA analysis

## Before running the code
It is suggested to run statistical code in the separate folder.
```
mkdir /nfs/dust/cms/user/username/HtoAA
cd /nfs/dust/cms/user/username/HtoAA
cp /nfs/dust/cms/user/rasp/CMSSW/CMSSW_10_2_25/src/HtoAA/4Tau/stat/* ./
```

Make sure that the outputs of [analysis_macro_4tau.cpp]() macro are located in separate subfolders
under the same base directory
```
$base_directory/$year
```

where `$year = {2016_postVFP, 2016_postVFP, }

The script will look for outputs in these folders. Specify correct name of the base directory in the RooT macro [CreateCards.C]()

