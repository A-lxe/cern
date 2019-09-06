kinit
eval `scramv1 runtime -sh`
root -b -q -l "Read_FlatNtuple.cpp(\"output_dataset-$1.root\", $1, false, $2)"