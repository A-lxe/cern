#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"

#include "Read_FlatNtuple.h" // List of input branches and functions to return values

void NTupleTest(const char* inputFile, int maxEvents, const char* outFile)
{

  TFile* out = new TFile(outFile, "RECREATE");
  TH1D* lctHist = new TH1D("LCTTheta", "LCT Theta", 100, -50, 50);
  TH1D* trkHist = new TH1D("TrackPt", "Track pT", 100, 0, 100);

  TChain* in_chain = new TChain("FlatNtupleData/tree");
  in_chain->Add(inputFile);

  InitializeMaps();
  SetBranchAddresses(in_chain);

  int nEvents = in_chain->GetEntries();

  for (int iEvt = 0; iEvt < nEvents; iEvt++) {
    if (iEvt > maxEvents)
      break;
    if ((iEvt % 10000) == 0) {
      std::cout << "Looking at event " << iEvt << " out of " << nEvents
                << std::endl;
    }

    in_chain->GetEntry(iEvt);

    for (int i = 0; i < I("nHits"); i++) {
      if (I("hit_isCSC", i) != 1)
        continue;
      lctHist->Fill(F("hit_theta", i));
    }

    for (int trk = 0; trk < I("nTracks"); trk++) {
      trkHist->Fill(F("trk_pt", trk));
    }
  }

  delete in_chain;
  out->Write();
}
