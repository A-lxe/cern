
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
bool isPhiPair(int i, int j);
bool isThetaPair(int i, int j);

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
TH1D *EventPhiPairCount =
    new TH1D("EventPhiPairCount", "EventPhiPairCount", 10, 0, 10);
TH1D *EventPhiPairTrackBreakdown = new TH1D(
    "EventPhiPairTrackBreakdown", "EventPhiPairTrackBreakdown", 4, 0, 4);
TH1D *EventThetaPairCount =
    new TH1D("EventThetaPairCount", "EventThetaPairCount", 10, 0, 10);
TH1D *EventThetaPairTrackBreakdown = new TH1D(
    "EventThetaPairTrackBreakdown", "EventThetaPairTrackBreakdown", 4, 0, 4);
TH1D *EventQuadCount = new TH1D("EventQuadCount", "EventQuadCount", 10, 0, 10);
TH1D *EventQuadTrackBreakdown =
    new TH1D("EventQuadTrackBreakdown", "EventQuadTrackBreakdown", 4, 0, 4);

TH1D *SingleMuPtLoWPhiPair =
    new TH1D("SingleMuonPtLoWPhiPair", "SingleMuonPtLoWPhiPair", 22, 1, 23);
TH1D *SingleMuPtHiWPhiPair =
    new TH1D("SingleMuonPtHiWPhiPair", "SingleMuonPtHiWPhiPair", 78, 22, 101);
TH1D *SingleMuPtAllWPhiPair =
    new TH1D("SingleMuonPtAllWPhiPair", "SingleMuonPtAllWPhiPair", 100, 1, 101);
TH1D *SingleMuPtLoWThetaPair =
    new TH1D("SingleMuonPtLoWThetaPair", "SingleMuonPtLoWThetaPair", 22, 1, 23);
TH1D *SingleMuPtHiWThetaPair = new TH1D(
    "SingleMuonPtHiWThetaPair", "SingleMuonPtHiWThetaPair", 78, 22, 101);
TH1D *SingleMuPtAllWThetaPair = new TH1D(
    "SingleMuonPtAllWThetaPair", "SingleMuonPtAllWThetaPair", 100, 1, 101);
TH1D *SingleMuPtLoWQuad =
    new TH1D("SingleMuonPtLoWQuad", "SingleMuonPtLoWQuad", 22, 1, 23);
TH1D *SingleMuPtHiWQuad =
    new TH1D("SingleMuonPtHiWQuad", "SingleMuonPtHiWQuad", 78, 22, 101);
TH1D *SingleMuPtAllWQuad =
    new TH1D("SingleMuonPtAllWQuad", "SingleMuonPtAllWQuad", 100, 1, 101);

TH1D *SingleMuPtLoNoPair =
    new TH1D("SingleMuonPtLoNoPair", "SingleMuonPtLoNoPair", 22, 1, 23);
TH1D *SingleMuPtHiNoPair =
    new TH1D("SingleMuonPtHiNoPair", "SingleMuonPtHiNoPair", 78, 22, 101);
TH1D *SingleMuPtAllNoPair =
    new TH1D("SingleMuonPtAllNoPair", "SingleMuonPtAllNoPair", 100, 1, 101);

TH1D *SingleMuPhiPairsInTracks =
    new TH1D("SingleMuPhiPairsInTracks", "SingleMuPhiPairsInTracks", 10, 0, 10);
TH1D *SingleMuThetaPairsInTracks = new TH1D(
    "SingleMuThetaPairsInTracks", "SingleMuThetaPairsInTracks", 10, 0, 10);
TH1D *SingleMuQuadInTracks =
    new TH1D("SingleMuQuadInTracks", "SingleMuQuadInTracks", 10, 0, 10);

TH2D *SingleMuOccupancy =
    new TH2D("ChamberOccupancy", "ChamberOccupancy", 54, 1, 55, 12, -6, 6);
TH2D *SingleMuPhiPairOccupancy = new TH2D(
    "PhiPairChamberOccupancy", "PhiPairChamberOccupancy", 54, 1, 55, 12, -6, 6);
TH2D *SingleMuPhiPairInTrack = new TH2D(
    "SingleMuTrackPhiPair", "SingleMuTrackPhiPair", 54, 1, 55, 12, -6, 6);
TH2D *SingleMuThetaPairOccupancy =
    new TH2D("ThetaPairChamberOccupancy", "ThetaPairChamberOccupancy", 54, 1,
             55, 12, -6, 6);
TH2D *SingleMuThetaPairInTrack = new TH2D(
    "SingleMuTrackThetaPair", "SingleMuTrackThetaPair", 54, 1, 55, 12, -6, 6);
TH2D *SingleMuQuadOccupancy = new TH2D(
    "QuadChamberOccupancy", "QuadChamberOccupancy", 54, 1, 55, 12, -6, 6);
TH2D *SingleMuQuadInTrack =
    new TH2D("SingleMuTrackQuad", "SingleMuTrackQuad", 54, 1, 55, 12, -6, 6);

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
     file_name.Form("%s/%s/NTuple_ZeroBias8b4e_FlatNtuple_Run_302674_2017_09_30.root",
     store.Data(), in_dir.Data());
    // file_name.Form("%s/%s/NTuple_ZeroBias8b4e_FlatNtuple_Skim_Run_302674_2017_09_30.root",
    // store.Data(), in_dir.Data());
    // file_name.Form("%s/%s/NTuple_ZeroBiasIsolatedBunches_FlatNtuple_Run_302674_2017_09_30.root",
    // store.Data(), in_dir.Data());
    // file_name.Form("%s/%s/NTuple_ZeroBiasIsolatedBunches_FlatNtuple_Skim_Run_302674_2017_09_30.root",
    // store.Data(), in_dir.Data());
  //  file_name.Form(
  //      "%s/%s/"
  //      "NTuple_ZeroBiasNominalTrains_FlatNtuple_Run_302674_2017_09_"
  //      "30.root",
  //      store.Data(), in_dir.Data());
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
    int phiPairCount = 0;
    int thetaPairCount = 0;
    int quadCount = 0;
    // <HitA <HitB, #TrksWHitA>>
    std::map<int, std::pair<int, int>> phiPairs;
    std::map<int, std::pair<int, int>> thetaPairs;
    std::map<int, std::tuple<int, int, int, int, int> *> quads;

    for (int i = 0; i < I("nHits"); i++) {
      if (I("hit_isCSC", i) != 1 || I("hit_neighbor", i) == 1) continue;
      fillOccupancyHit(SingleMuOccupancy, i);
      for (int j = i + 1; j < I("nHits"); j++) {
        if (I("hit_isCSC", j) != 1 || I("hit_neighbor", j) == 1) continue;
        // What about hits with same theta? Do I need to account for sector etc?
        bool isPhiP = isPhiPair(i, j);
        bool isThetaP = isThetaPair(i, j);
        if (isPhiP) {
          fillOccupancyHit(SingleMuPhiPairOccupancy, i);
          phiPairCount += 1;
          std::pair<int, int> p1(j, 0);
          std::pair<int, int> p2(i, 0);
          phiPairs[i] = p1;
          phiPairs[j] = p2;
        }
        if (isThetaP) {
          fillOccupancyHit(SingleMuThetaPairOccupancy, i);
          thetaPairCount += 1;
          std::pair<int, int> p1(j, 0);
          std::pair<int, int> p2(i, 0);
          thetaPairs[i] = p1;
          thetaPairs[j] = p2;
        }
      }
    }
    EventPhiPairCount->Fill(phiPairCount, 1);
    EventThetaPairCount->Fill(thetaPairCount, 1);

    for (auto const & phiPair : phiPairs) {
      int phiHitA = phiPair.first;
      int phiHitB = phiPair.second.first;
      if (quads.find(phiHitA) == quads.end() &&
          thetaPairs.find(phiHitA) != thetaPairs.end() &&
          thetaPairs.find(phiHitB) != thetaPairs.end() &&
          phiPairs.find(thetaPairs[phiHitA].first) != phiPairs.end() &&
          phiPairs.find(thetaPairs[phiHitB].first) != phiPairs.end() &&
          phiPairs[thetaPairs[phiHitB].first].first ==
              thetaPairs[phiHitB].first) {
        quadCount += 1;
        int hit3 = thetaPairs[phiHitA].first;
        int hit4 = thetaPairs[phiHitB].first;

        std::tuple<int,int,int,int,int> quadTuple = {phiHitA, phiHitB, hit3, hit4, 0};
        quads[phiHitA] = &quadTuple;
        quads[phiHitB] = &quadTuple;
        quads[hit3] = &quadTuple;
        quads[hit4] = &quadTuple;
      }
    }
    EventQuadCount->Fill(quadCount, 1);

    // Loop over unpacked tracks
    for (int trk = 0; trk < I("nUnpTracks"); trk++) {
      bool trackUsesPhiPair = false;
      bool trackUsesThetaPair = false;
      bool trackUsesQuad = false;

      if (/*I("unp_trk_nRPC",i)==0 &&*/ I("unp_trk_nHits", trk) >= 3) {
        // Loop over hits in this tracks
        for (int trkHitIdx = 0; trkHitIdx < I("unp_trk_found_hits", trk);
             trkHitIdx++) {
          int iHit = I("unp_trk_iHit", trk, trkHitIdx);

          auto phiP = phiPairs.find(iHit);
          if (phiP != phiPairs.end()) {
            phiPairs[iHit].second += 1;
            trackUsesPhiPair = true;
            fillOccupancyHit(SingleMuPhiPairInTrack, iHit);
          }
          auto thetaP = thetaPairs.find(iHit);
          if (thetaP != thetaPairs.end()) {
            thetaPairs[iHit].second += 1;
            trackUsesThetaPair = true;
            fillOccupancyHit(SingleMuThetaPairInTrack, iHit);
          }
          auto quadP = quads.find(iHit);
          if (quadP != quads.end()) {
            std::get<4>(*quads[iHit]) += 1;
            trackUsesQuad = true;
            fillOccupancyHit(SingleMuQuadInTrack, iHit);
          }
        }

        if (trackUsesPhiPair) {
          SingleMuPtAllWPhiPair->Fill(F("unp_trk_pt", trk), 1);
          F("unp_trk_pt", trk) < 22
              ? SingleMuPtLoWPhiPair->Fill(F("unp_trk_pt", trk), 1)
              : SingleMuPtHiWPhiPair->Fill(F("unp_trk_pt", trk), 1);
        }
        if (trackUsesThetaPair) {
          SingleMuPtAllWThetaPair->Fill(F("unp_trk_pt", trk), 1);
          F("unp_trk_pt", trk) < 22
              ? SingleMuPtLoWThetaPair->Fill(F("unp_trk_pt", trk), 1)
              : SingleMuPtHiWThetaPair->Fill(F("unp_trk_pt", trk), 1);
        }
        if (trackUsesQuad) {
          SingleMuPtAllWQuad->Fill(F("unp_trk_pt", trk), 1);
          F("unp_trk_pt", trk) < 22
              ? SingleMuPtLoWQuad->Fill(F("unp_trk_pt", trk), 1)
              : SingleMuPtHiWQuad->Fill(F("unp_trk_pt", trk), 1);
        }

        if (!(trackUsesPhiPair || trackUsesThetaPair || trackUsesQuad)) {
          SingleMuPtAllNoPair->Fill(F("unp_trk_pt", trk), 1);
          F("unp_trk_pt", trk) < 22
              ? SingleMuPtLoNoPair->Fill(F("unp_trk_pt", trk), 1)
              : SingleMuPtHiNoPair->Fill(F("unp_trk_pt", trk), 1);
        }

        if (F("unp_trk_pt", trk) < 22) {
          EventPhiPairTrackBreakdown->Fill(phiPairCount > 0 ? 2 : 0, 1);
          EventThetaPairTrackBreakdown->Fill(thetaPairCount > 0 ? 2 : 0, 1);
          EventQuadTrackBreakdown->Fill(quadCount > 0 ? 2 : 0, 1);

        } else {
          EventPhiPairTrackBreakdown->Fill(phiPairCount > 0 ? 3 : 1, 1);
          EventThetaPairTrackBreakdown->Fill(thetaPairCount > 0 ? 3 : 1, 1);
          EventQuadTrackBreakdown->Fill(quadCount > 0 ? 3 : 1, 1);
        }
      }
    }

    // Fill in the track counts histogram
    for (std::map<int, std::pair<int, int>>::iterator it = phiPairs.begin();
         it != phiPairs.end(); ++it) {
      int tracks = it->second.second + phiPairs[it->second.first].second;
      SingleMuPhiPairsInTracks->Fill(tracks, 0.5);
    }
    for (std::map<int, std::pair<int, int>>::iterator it = thetaPairs.begin();
         it != thetaPairs.end(); ++it) {
      int tracks = it->second.second + thetaPairs[it->second.first].second;
      SingleMuThetaPairsInTracks->Fill(tracks, 0.5);
    }
    for (auto const & quad : quads) {
      SingleMuQuadInTracks->Fill(std::get<4>(*quad.second), 1);
    }
  }  // End loop: for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Finished looping over the events *******"
            << std::endl;

  delete in_chain;

  std::cout << "\nDone with Read_FlatNtuple(). Exiting.\n" << std::endl;
  f->Write();
}  // End function: void Read_FlatNtuple()

// cscid + cscid_offset, endcap * (station - 0.5)
void fillOccupancyHit(TH2D *plot, int hit) {
  plot->Fill(
      I("hit_CSC_ID", hit) + (I("hit_sector", hit) - 1) * 9,
      I("hit_endcap", hit) * (I("hit_station", hit) +
                              (I("hit_subsector", hit) == 1 ? -0.5 : 0.5)));
}

bool isPhiPair(int i, int j) {
  return I("hit_endcap", i) == I("hit_endcap", j) &&
         I("hit_station", i) == I("hit_station", j) &&
         I("hit_ring", i) == I("hit_ring", j) &&
         I("hit_chamber", i) == I("hit_chamber", j) &&
         I("hit_BX", i) == I("hit_BX", j) &&
         I("hit_phi_int", i) == I("hit_phi_int", j) &&
         I("hit_theta_int", i) != I("hit_theta_int", j);
}

bool isThetaPair(int i, int j) {
  return I("hit_endcap", i) == I("hit_endcap", j) &&
         I("hit_station", i) == I("hit_station", j) &&
         I("hit_ring", i) == I("hit_ring", j) &&
         I("hit_chamber", i) == I("hit_chamber", j) &&
         I("hit_BX", i) == I("hit_BX", j) &&
         I("hit_theta_int", i) == I("hit_theta_int", j) &&
         I("hit_phi_int", i) != I("hit_phi_int", j);
}
