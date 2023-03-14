#!/bin/csh
rm $1.root 
hadd $1.root ./2016_postVFP/$1_files/*.root ./2016_preVFP/$1_files/*.root

