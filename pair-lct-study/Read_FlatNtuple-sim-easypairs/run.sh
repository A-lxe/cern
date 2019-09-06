#!/bin/bash

if [ $# -lt 2 ]; then
  echo "Not enough arguments!"
  exit 1
fi

reco_cut=$1
dataset=$2
max_events=$3

if [ "$4" != "" ]; then
  filename="$4"
else
  filename="dataset-${dataset}_events-${max_events}.root"
fi

# HtCondor can't copy directory structure, so here's a workaround for local runs
if [ -e ./src ]; then
  cd src
  filename=../$filename
fi

echo "Running cmsenv"
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_9_2_8
eval $(/cvmfs/cms.cern.ch/common/scramv1 runtime -sh)
cd ~-

echo "Running Read_FlatNtuple.cpp on dataset $dataset with $max_events events"

root -b -q -l "Read_FlatNtuple.cpp+(\"$filename\", $dataset, false, $reco_cut, $max_events)"

echo "Outputting results to $filename"
