
////////////////////////////////////////////////////////////
///        Simple macro to read EMTF Flat NTuples        ///
///              Andrew Brinkerhoff 29.09.17             ///
///                                                      ///
///   TChain can be used to read multiple files.         ///
///   Format: interface/FlatNtupleBranches.h             ///
////////////////////////////////////////////////////////////

#include <list>
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "Read_FlatNtuple.h"  // List of input branches and functions to return values

void fillOccupancyHit(TH2D *plot, int hit);

const int MAX_FILES = 10;  // Max number of files to process
const int MAX_EVT = 2000000;  // Max number of events to process
const int PRT_EVT = 1000;  // Print every N events
const bool verbose =
    false;  // Print information about the event and RECO and L1T muons

///////////////////////////////////////////////
/////  Prepare the plots and their labels  ////
///////////////////////////////////////////////
TFile *f = new TFile("PURateTrainRPC.root", "RECREATE");

/////Alex's PU rate study/////
/////      :)START(:      ////
//////////////////////////////
TH1D *EventPairCount = new TH1D("EventPairCount", "EventPairCount", 10, 0, 10);

TH1D *SingleMuPtLoWPair =
    new TH1D("SingleMuonPtLoWPair", "SingleMuonPtLoWPair", 22, 1, 23);
TH1D *SingleMuPtHiWPair =
    new TH1D("SingleMuonPtHiWPair", "SingleMuonPtHiWPair", 78, 22, 101);
TH1D *SingleMuPtAllWPair =
    new TH1D("SingleMuonPtAllWPair", "SingleMuonPtAllWPair", 100, 1, 101);

TH1D *SingleMuPtLoNoPair =
    new TH1D("SingleMuonPtLoNoPair", "SingleMuonPtLoNoPair", 22, 1, 23);
TH1D *SingleMuPtHiNoPair =
    new TH1D("SingleMuonPtHiNoPair", "SingleMuonPtHiNoPair", 78, 22, 101);
TH1D *SingleMuPtAllNoPair =
    new TH1D("SingleMuonPtAllNoPair", "SingleMuonPtAllNoPair", 100, 1, 101);

TH2D *SingleMuOccupancy =
    new TH2D("SingleMuOccupancy", "SingleMuOccupancy", 54, 1, 55, 12, -6, 6);
TH2D *SingleMuPairOccupancy = new TH2D(
    "SingleMuPairOccupancy", "SingleMuPairOccupancy", 54, 1, 55, 12, -6, 6);
TH2D *SingleMuPairInTrack = new TH2D(
    "SingleMuPairInTrack", "SingleMuPairInTrack", 54, 1, 55, 12, -6, 6);

//////////////////////////////////////////////
/////  Main function to read the NTuples  ////
//////////////////////////////////////////////

void Read_FlatNtuple() {
  // BookHistos();

  gROOT->SetStyle("Plain");
  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);

  // List of input files
  std::vector<TString> in_file_names;
  TString store = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/";
  TString in_dir = "ntuples/HADD/";
  TString file_name;

  for (int i = 1; i < MAX_FILES + 1; i++) {
    // file_name.Form("%s/%s/NTuple_ZeroBias8b4e_FlatNtuple_Run_302674_2017_09_30.root",
    // store.Data(), in_dir.Data());
    // file_name.Form("%s/%s/NTuple_ZeroBias8b4e_FlatNtuple_Skim_Run_302674_2017_09_30.root",
    // store.Data(), in_dir.Data());
    // file_name.Form("%s/%s/NTuple_ZeroBiasIsolatedBunches_FlatNtuple_Run_302674_2017_09_30.root",
    // store.Data(), in_dir.Data());
    // file_name.Form("%s/%s/NTuple_ZeroBiasIsolatedBunches_FlatNtuple_Skim_Run_302674_2017_09_30.root",
    // store.Data(), in_dir.Data());
    file_name.Form(
        "%s/%s/"
        "NTuple_ZeroBiasNominalTrains_FlatNtuple_Run_302674_2017_09_"
        "30.root",
        store.Data(), in_dir.Data());
    // file_name.Form("%s/%s/NTuple_ZeroBiasNominalTrains_FlatNtuple_Skim_Run_302674_2017_09_30.root",
    // store.Data(), in_dir.Data());
    std::cout << "Adding file " << file_name.Data() << std::endl;
    in_file_names.push_back(file_name.Data());
  }

  // Open all input files
  for (int i = 0; i < in_file_names.size(); i++) {
    if (!gSystem->AccessPathName(in_file_names.at(i)))
      file_tmp = TFile::Open(in_file_names.at(i));  // Check if file exists
    if (!file_tmp) {
      std::cout << "ERROR: could not open data file " << in_file_names.at(i)
                << std::endl;
      return;
    }
  }

  // Add tree from the input files to the TChain
  TChain *in_chain = new TChain("ntuple/tree");
  for (int i = 0; i < in_file_names.size(); i++) {
    in_chain->Add(in_file_names.at(i));
  }

  InitializeMaps();
  SetBranchAddresses(in_chain);

  // if (verbose) in_chain->GetListOfBranches()->Print();

  std::cout << "\n******* About to loop over the events *******" << std::endl;
  int nEvents = in_chain->GetEntries();
  for (int iEvt = 0; iEvt < nEvents; iEvt++) {
    if (iEvt > MAX_EVT) break;
    if ((iEvt % PRT_EVT) == 0) {
      std::cout << "\n*************************************" << std::endl;
      std::cout << "Looking at event " << iEvt << " out of " << nEvents
                << std::endl;
      std::cout << "*************************************" << std::endl;
    }

    in_chain->GetEntry(iEvt);

    // From Read_FlatNtuple.h, use 'I("branch_name")' to get an integer branch
    // value, 'F("branch_name") to get a float
    if (verbose)
      std::cout << "\nRun = " << I("evt_run") << ", LS = " << I("evt_LS")
                << ", event = " << I("evt_event") << std::endl;

    // Print info for emulated EMTF hits
    if (verbose)
      std::cout << "\n"
                << I("nHits") << " emulated EMTF hits in the event"
                << std::endl;

    // Alex PU Study: Find events with pair LCTs with same chamber, BX, and phi,
    // but differnet Theta
    int pairCount = 0;
    std::list<std::pair<int, int>> pairs;
    for (int i = 0; i < I("nHits"); i++) {
      if (I("hit_isCSC", i) != 1 || I("hit_neighbor", i) == 1) continue;
      // cscid + cscid_offset, endcap * (station - 0.5)
      SingleMuOccupancy->Fill(
          I("hit_CSC_ID", i) + (I("hit_sector", i) - 1) * 9,
          I("hit_endcap", i) * (I("hit_station", i) +
                                (I("hit_subsector", i) == 1 ? -0.5 : 0.5)));
      for (int j = i + 1; j < I("nHits"); j++) {
        if (I("hit_isCSC", j) != 1 || I("hit_neighbor", j) == 1) continue;
        // What about hits with same theta? Do I need to account for sector etc?
        bool isPair = I("hit_endcap", i) == I("hit_endcap", j) &&
                      I("hit_station", i) == I("hit_station", j) &&
                      I("hit_ring", i) == I("hit_ring", j) &&
                      I("hit_chamber", i) == I("hit_chamber", j) &&
                      I("hit_BX", i) == I("hit_BX", j) &&
                      I("hit_phi_int", i) == I("hit_phi_int", j) &&
                      I("hit_theta_int", i) != I("hit_BX", j);
        if (isPair) {
          SingleMuPairOccupancy->Fill(
              I("hit_CSC_ID", i) + (I("hit_sector", i) - 1) * 9,
              I("hit_endcap", i) * (I("hit_station", i) +
                                    (I("hit_subsector", i) == 1 ? -0.5 : 0.5)));
          pairCount += 1;
          std::pair<int, int> pair(i, j);
          pairs.push_back(pair);
        }
      }
    }
    EventPairCount->Fill(pairCount, 1);

    // Loop over unpacked tracks
    for (int trk = 0; trk < I("nUnpTracks"); trk++) {
      bool trackUsesPair = false;

      // Loop over hits in this tracks
      for (int trkHitIdx = 0; trkHitIdx < I("unp_trk_found_hits", trk);
           trkHitIdx++) {
        int iHit = I("unp_trk_iHit", trk, trkHitIdx);

        for (auto pair : pairs) {
          if (iHit == pair.first || iHit == pair.second) {
            trackUsesPair = true;
            if(I("unp_trk_nHits", trk) >= 3)
              fillOccupancyHit(SingleMuPairInTrack, iHit);
          }
        }
      }

      if (/*I("unp_trk_nRPC",i)==0 &&*/ I("unp_trk_nHits", trk) >= 3) {
        if (trackUsesPair) {
          SingleMuPtAllWPair->Fill(F("unp_trk_pt", trk), 1);
          F("unp_trk_pt", trk) < 22
              ? SingleMuPtLoWPair->Fill(F("unp_trk_pt", trk), 1)
              : SingleMuPtHiWPair->Fill(F("unp_trk_pt", trk), 1);
        } else {
          SingleMuPtAllNoPair->Fill(F("unp_trk_pt", trk), 1);
          F("unp_trk_pt", trk) < 22
              ? SingleMuPtLoNoPair->Fill(F("unp_trk_pt", trk), 1)
              : SingleMuPtHiNoPair->Fill(F("unp_trk_pt", trk), 1);
        }
      }
    }

  }  // End loop: for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Finished looping over the events *******"
            << std::endl;

  delete in_chain;

  std::cout << "\nDone with Read_FlatNtuple(). Exiting.\n" << std::endl;
  f->Write();
}  // End function: void Read_FlatNtuple()

void fillOccupancyHit(TH2D *plot, int hit) {
  plot->Fill(
      I("hit_CSC_ID", hit) + (I("hit_sector", hit) - 1) * 9,
      I("hit_endcap", hit) * (I("hit_station", hit) +
                              (I("hit_subsector", hit) == 1 ? -0.5 : 0.5)));
}
