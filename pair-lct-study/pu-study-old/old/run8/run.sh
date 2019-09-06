#!/bin/bash

if [ $# ] <2; then
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

echo "Running Read_FlatNtuple.cpp on dataset $dataset with $max_events events"
echo "Outputting results to $filename"

eval $(scramv1 runtime -sh)
root -b -q -l "Read_FlatNtuple.cpp(\"$filename\", $1, false, $2)"
echo "Outputting results to $filename"
