#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import ROOT


def run(in_file, max_events, out_file):

    # Set ROOT options --------------------------------------------------------

    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    # Prepare input/output ----------------------------------------------------

    # Create output file and load empty histograms
    output_file = ROOT.TFile(out_file, "RECREATE")
    lctHist = ROOT.TH1D("LCTTheta", "LCT Theta", 100, -50, 50)
    trkHist = ROOT.TH1D("TrackPt", "Track pT", 100, 0, 100)

    # Load ntuple files
    chain = ROOT.TChain("FlatNtupleData/tree")
    chain.Add(in_file)

    # Event loop --------------------------------------------------------------

    # Technically event is the same TChain pointer as chain, but iterating
    # updates the event index. Enumerate just lets counter hold the number of
    # events processed
    for counter, event in enumerate(chain):
        if counter % 10000 == 0:
            print("Processed {0} events".format(counter))
        if counter > max_events:
            break

        for lct in range(event.nHits):
            if not event.hit_isCSC[lct]:
                continue
            lctHist.Fill(event.hit_theta[lct])

        for trk in range(event.nTracks):
            trkHist.Fill(event.trk_pt[trk])


if __name__ == "__main__":
    run(sys.argv[1], int(sys.argv[2]), sys.argv[3])
