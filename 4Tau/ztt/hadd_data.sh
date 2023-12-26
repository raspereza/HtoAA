#!/bin/bash
if [ -f "Data.root" ]; then
    rm Data.root
fi
hadd Data.root SingleMuon*root


