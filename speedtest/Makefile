.PHONY: cpp python

cpp:
	time root -b -q -l "cpp/NTupleTest.cpp+(\"NTuple_SingleMuon_FlatNtuple_Run_306092_2018_03_02_SingleMu.root\", 2000000, \"cpp.root\")"

python:
	time python python/NTupleTest.py "NTuple_SingleMuon_FlatNtuple_Run_306092_2018_03_02_SingleMu.root" 2000000 "python.root"
