#!/bin/sh

max_events=3000000

if [ $# -lt 1 ]; then
  config=all
else
  config=$1
fi

echo "Running config $config with $max_events events"

case "$config" in
"all")
  tmux-sync.sh "Read_FlatNTuples" \
    "./run.sh 'false' 3 $max_events; read" \
    "sleep 30; ./run.sh 'false' 4 $max_events; read" \
    "sleep 30; ./run.sh 'false' 5 $max_events; read" \
    "sleep 30; ./run.sh 'true' 4 $max_events ./dataset-4_events-${max_events}_reco-cut.root; read"
  ;;
"0")
  ./run.sh 'false' 0 $max_events
  ;;
"1")
  ./run.sh 'false' 1 $max_events
  ;;
"2")
  ./run.sh 'false' 2 $max_events
  ;;
"3")
  ./run.sh 'false' 3 $max_events
  ;;
"4")
  ./run.sh 'false' 4 $max_events
  ;;
"5")
  ./run.sh 'false' 5 $max_events
  ;;
"6")
  ./run.sh 'true' 4 $max_events ./dataset-4_events-${max_events}_reco-cut.root
  ;;
esac
