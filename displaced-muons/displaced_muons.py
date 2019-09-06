#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Analysis code for studying the CSC L1T DAQ data rate through EMTF LCTs.

See
https://github.com/chadfreer/cmssw/tree/LCT-Matched-Plotter/EMTFAnalyzer/NTupleMaker
for a list of values available in the NTuples used.

Author: Alex Aubuchon
"""
import ROOT
from plots import get_plots, post_fill
from helpers import possible_tracks


# Constants -------------------------------------------------------------------

MAX_EVT = 1
OUTFILE = "output.root"
FILES = [
    "NTuple_SingleMuon_FlatNtuple_Run_306092_2018_03_02_SingleMu.root",
]


def run():

    # Set ROOT options --------------------------------------------------------

    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    # Prepare input/output ----------------------------------------------------

    output_file = ROOT.TFile(OUTFILE, "RECREATE")
    p = get_plots(output_file)  # Load plots from ./plots.py

    chain = ROOT.TChain("FlatNtupleData/tree")
    for f in FILES:
        chain.Add(f)

    # Event loop --------------------------------------------------------------

    counter = 0
    # Technically event is the same TChain pointer as chain, but iterating
    # updates the index
    for event in chain:
        if counter % 1000 == 0:
            print("Processed {0} events".format(counter))
        if counter > MAX_EVT:
            break
        counter += 1

        lcts = {(endcap, station): [] for endcap in [+1,-1] for station in [1,2,3,4]}
        for lct in range(event.nHits):
            stat = event.hit_station[lct]
            endcap = event.hit_endcap[lct]
            lcts[(endcap,stat)].append(lct)
        
        print(list(possible_tracks(lcts)))

        for trk in range(event.nTracks):
            pass

    post_fill(output_file, p)


if __name__ == "__main__":
    run()
