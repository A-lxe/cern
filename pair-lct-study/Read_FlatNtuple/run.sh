#!/bin/bash

if [ $# -lt 2 ]; then
  echo "Not enough arguments!"
  exit 1
fi

dataset=$1
max_events=$2

if [ "$3" != "" ]; then
  filename="$3"
else
  filename="dataset-${1}_events-${2}.root"
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

root -b -q -l "Read_FlatNtuple.cpp+(\"$filename\", $1, false, $2)"

echo "Outputting results to $filename"
