
////////////////////////////////////////////////////////////
///        Simple macro to read EMTF Flat NTuples        ///
///              Andrew Brinkerhoff 29.09.17             ///
///                                                      ///
///   TChain can be used to read multiple files.         ///
///   Format: interface/FlatNtupleBranches.h             ///
////////////////////////////////////////////////////////////

#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <iostream>
#include <list>

#include "Read_FlatNtuple.h" // List of input branches and functions to return values

void fillOccupancyHit(TH2D *plot, int hit, int weight);
void fillOccupancyHit(TH2D *plot, int hit);
bool isPhiPair(int i, int j);
bool isThetaPair(int i, int j);

const int MAX_FILES = 10;  // Max number of files to process
const int MAX_EVT = 1000;  // Max number of events to process
const int PRT_EVT = 10000; // Print every N events
const bool verbose =
    false; // Print information about the event and RECO and L1T muons

///////////////////////////////////////////////
/////  Prepare the plots and their labels  ////
///////////////////////////////////////////////

//////////////////////////////////////////////
/////  Main function to read the NTuples  ////
//////////////////////////////////////////////

void Read_FlatNtuple(const char *outFile = "output.root", int dataset = 0,
                     bool skim = false, int maxEvents = MAX_EVT,
                     int maxFiles = MAX_FILES) {
  TFile *f = new TFile(outFile, "RECREATE");
#include "plots.h"

  // BookHistos();
  gROOT->SetStyle("Plain");
  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);

  // List of input files
  std::vector<TString> in_file_names;
  TString store = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/";
  TString in_dir = "ntuples/HADD/";
  TString file_name;

  for (int i = 1; i < maxFiles + 1; i++) {
    switch (dataset) {
    case 0:
      file_name.Form(
          "%s/%s/NTuple_ZeroBias8b4e_FlatNtuple%sRun_302674_2017_09_30.root",
          store.Data(), in_dir.Data(), skim ? "_Skim_" : "_");
      break;
    case 1:
      file_name.Form(
          "%s/%s/"
          "NTuple_ZeroBiasIsolatedBunches_FlatNtuple%sRun_302674_2017_09_30."
          "root",
          store.Data(), in_dir.Data(), skim ? "_Skim_" : "_");
      break;
    case 2:
      file_name.Form(
          "%s/%s/"
          "NTuple_ZeroBiasNominalTrains_FlatNtuple%sRun_302674_2017_09_30."
          "root",
          store.Data(), in_dir.Data(), skim ? "_Skim_" : "_");
      break;
    }
    std::cout << "Adding file " << file_name.Data() << std::endl;
    in_file_names.push_back(file_name.Data());
  }

  // Open all input files
  for (int i = 0; i < in_file_names.size(); i++) {
    if (!gSystem->AccessPathName(in_file_names.at(i)))
      file_tmp = TFile::Open(in_file_names.at(i)); // Check if file exists
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
    if (iEvt > maxEvents)
      break;
    if ((iEvt % PRT_EVT) == 0) {
      std::cout << "Looking at event " << iEvt << " out of " << nEvents
                << std::endl;
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

    // Alex PU Study: Find events with pair LCTs with same chamber, BX, and
    // phi, but differnet Theta
    int phiPairCount = 0;
    int thetaPairCount = 0;
    int quadCount = 0;
    // <HitA <HitB, #TrksWHitA>>
    std::map<int, std::pair<int, int>> phiPairs;
    std::map<int, std::pair<int, int>> thetaPairs;
    std::map<int, std::tuple<int, int, int, int, int> *> quads;

    for (int i = 0; i < I("nHits"); i++) {
      if (I("hit_isCSC", i) != 1)
        continue;
      fillOccupancyHit(SingleMuOccupancy, i);

      AllLctStatVRing->Fill(I("hit_ring", i), I("hit_station", i));
      AllLctStat->Fill(I("hit_station", i));
      AllLctQual->Fill(I("hit_quality", i));
      AllLctTheta->Fill(F("hit_theta", i));

      for (int j = i + 1; j < I("nHits"); j++) {
        if (I("hit_isCSC", j) != 1)
          continue;
        // What about hits with same theta? Do I need to account for sector
        // etc?
        bool isPhiP = isPhiPair(i, j);
        bool isThetaP = isThetaPair(i, j);
        if (isPhiP) {
          fillOccupancyHit(SingleMuPhiPairOccupancy, i, 2);
          phiPairCount += 1;
          std::pair<int, int> p1(j, 0);
          std::pair<int, int> p2(i, 0);
          phiPairs[i] = p1;
          phiPairs[j] = p2;

          PhiPairDTheta->Fill(
              abs(F("hit_theta", i) - F("hit_theta", j)));
          PhiPairDWire->Fill(abs(I("hit_wire", i) - I("hit_wire", j)));

          PhiPairLctStatVRing->Fill(I("hit_ring", i), I("hit_station", i));
          PhiPairLctStatVRing->Fill(I("hit_ring", i), I("hit_station", i));

          PhiPairLctStat->Fill(I("hit_station", i));
          PhiPairLctStat->Fill(I("hit_station", i));

          PhiPairLctQual->Fill(I("hit_quality", i));
          PhiPairLctQual->Fill(I("hit_quality", j));
          PhiPairLctTheta->Fill(F("hit_theta", i));
          PhiPairLctTheta->Fill(F("hit_theta", j));
        }
        if (isThetaP) {
          fillOccupancyHit(SingleMuThetaPairOccupancy, i, 2);
          thetaPairCount += 1;
          std::pair<int, int> p1(j, 0);
          std::pair<int, int> p2(i, 0);
          thetaPairs[i] = p1;
          thetaPairs[j] = p2;

          ThetaPairDPhi->Fill(abs(I("hit_phi_int", i) - I("hit_phi_int", j)));
          ThetaPairDStrip->Fill(abs(I("hit_strip", i) - I("hit_strip", j)));

          ThetaPairLctStatVRing->Fill(I("hit_ring", i), I("hit_station", i));
          ThetaPairLctStatVRing->Fill(I("hit_ring", i), I("hit_station", i));

          ThetaPairLctStat->Fill(I("hit_station", i));
          ThetaPairLctStat->Fill(I("hit_station", i));

          ThetaPairLctQual->Fill(I("hit_quality", i));
          ThetaPairLctQual->Fill(I("hit_quality", j));
          ThetaPairLctTheta->Fill(F("hit_theta", i));
          ThetaPairLctTheta->Fill(F("hit_theta", j));
        }
      }
    }
    EventPhiPairCount->Fill(phiPairCount);
    EventThetaPairCount->Fill(thetaPairCount);

    for (auto const &phiPair : phiPairs) {
      int phiHitA = phiPair.first;
      int phiHitB = phiPair.second.first;
      if (quads.find(phiHitA) == quads.end() &&
          thetaPairs.find(phiHitA) != thetaPairs.end() &&
          thetaPairs.find(phiHitB) != thetaPairs.end() &&
          phiPairs.find(thetaPairs[phiHitA].first) != phiPairs.end() &&
          phiPairs.find(thetaPairs[phiHitB].first) != phiPairs.end() &&
          phiPairs[thetaPairs[phiHitA].first].first ==
              thetaPairs[phiHitB].first) {
        quadCount += 1;
        int hit3 = thetaPairs[phiHitA].first;
        int hit4 = thetaPairs[phiHitB].first;

        std::tuple<int, int, int, int, int> quadTuple{phiHitA, phiHitB, hit3,
                                                      hit4, 0};
        quads[phiHitA] = &quadTuple;
        quads[phiHitB] = &quadTuple;
        quads[hit3] = &quadTuple;
        quads[hit4] = &quadTuple;
        fillOccupancyHit(SingleMuQuadOccupancy, phiHitA, 4);

        QuadLctStatVRing->Fill(I("hit_ring", phiHitA),
                               I("hit_station", phiHitA), 4);
        QuadLctStat->Fill(I("hit_station", phiHitA));
        QuadLctStat->Fill(I("hit_station", phiHitA));
        QuadLctStat->Fill(I("hit_station", phiHitA));
        QuadLctStat->Fill(I("hit_station", phiHitA));

        QuadLctQual->Fill(I("hit_quality", phiHitA));
        QuadLctQual->Fill(I("hit_quality", phiHitB));
        QuadLctQual->Fill(I("hit_quality", hit3));
        QuadLctQual->Fill(I("hit_quality", hit4));
        QuadLctTheta->Fill(F("hit_theta", phiHitA));
        QuadLctTheta->Fill(F("hit_theta", phiHitB));
        QuadLctTheta->Fill(F("hit_theta", hit3));
        QuadLctTheta->Fill(F("hit_theta", hit4));
      }
    }
    EventQuadCount->Fill(quadCount);

    // Loop over unpacked tracks
    for (int trk = 0; trk < I("nUnpTracks"); trk++) {
      bool trackUsesPhiPair = false;
      bool trackUsesThetaPair = false;
      bool trackUsesQuad = false;

      if (I("unp_trk_nHits", trk) >= 3) {
        // Loop over hits in this tracks
        for (int trkHitIdx = 0; trkHitIdx < I("unp_trk_found_hits", trk);
             trkHitIdx++) {
          int iHit = I("unp_trk_iHit", trk, trkHitIdx);
          if (I("hit_isCSC", iHit) != 1) {
            continue;
          }
          int ring = I("hit_ring", iHit);
          int stat = I("hit_station", iHit);
          int qual = I("hit_quality", iHit);
          int theta = F("hit_theta", iHit);

          AllLctStatVRingAllPt->Fill(ring, stat);
          AllLctStatAllPt->Fill(stat);
          AllLctQualAllPt->Fill(qual);
          AllLctThetaAllPt->Fill(theta);

          F("unp_trk_pt", trk) < 22 ? AllLctStatVRingLowPt->Fill(ring, stat)
                                    : AllLctStatVRingHighPt->Fill(ring, stat);
          F("unp_trk_pt", trk) < 22 ? AllLctStatLowPt->Fill(stat)
                                    : AllLctStatHighPt->Fill(stat);
          F("unp_trk_pt", trk) < 22 ? AllLctQualLowPt->Fill(qual)
                                    : AllLctQualHighPt->Fill(qual);
          F("unp_trk_pt", trk) < 22 ? AllLctThetaLowPt->Fill(theta)
                                    : AllLctThetaHighPt->Fill(theta);

          auto phiP = phiPairs.find(iHit);
          if (phiP != phiPairs.end()) {
            phiPairs[iHit].second += 1;
            trackUsesPhiPair = true;
            fillOccupancyHit(SingleMuPhiPairInTrack, iHit);

            int dTheta = abs(F("hit_theta", iHit) -
                             F("hit_theta", phiPairs[iHit].first));
            int dWire =
                abs(I("hit_wire", iHit) - I("hit_wire", phiPairs[iHit].first));

            PhiPairTrackAllPtDTheta->Fill(dTheta);
            PhiPairTrackAllPtDWire->Fill(dWire);
            F("unp_trk_pt", trk) < 22 ? PhiPairTrackLoPtDTheta->Fill(dTheta)
                                      : PhiPairTrackHiPtDTheta->Fill(dTheta);
            F("unp_trk_pt", trk) < 22 ? PhiPairTrackLoPtDWire->Fill(dWire)
                                      : PhiPairTrackHiPtDWire->Fill(dWire);

            PhiPairLctStatVRingAllPt->Fill(ring, stat);
            PhiPairLctStatAllPt->Fill(stat);
            PhiPairLctQualAllPt->Fill(qual);
            PhiPairLctThetaAllPt->Fill(theta);

            F("unp_trk_pt", trk) < 22
                ? PhiPairLctStatVRingLowPt->Fill(ring, stat)
                : PhiPairLctStatVRingHighPt->Fill(ring, stat);
            F("unp_trk_pt", trk) < 22 ? PhiPairLctStatLowPt->Fill(stat)
                                      : PhiPairLctStatHighPt->Fill(stat);
            F("unp_trk_pt", trk) < 22 ? PhiPairLctQualLowPt->Fill(qual)
                                      : PhiPairLctQualHighPt->Fill(qual);
            F("unp_trk_pt", trk) < 22 ? PhiPairLctThetaLowPt->Fill(theta)
                                      : PhiPairLctThetaHighPt->Fill(theta);
          }
          auto thetaP = thetaPairs.find(iHit);
          if (thetaP != thetaPairs.end()) {
            thetaPairs[iHit].second += 1;
            trackUsesThetaPair = true;
            fillOccupancyHit(SingleMuThetaPairInTrack, iHit);

            int dPhi = abs(I("hit_phi_int", iHit) -
                           I("hit_phi_int", thetaPairs[iHit].first));
            int dStrip = abs(I("hit_strip", iHit) -
                             I("hit_strip", thetaPairs[iHit].first));

            ThetaPairTrackAllPtDPhi->Fill(dPhi);
            ThetaPairTrackAllPtDStrip->Fill(dStrip);
            F("unp_trk_pt", trk) < 22 ? ThetaPairTrackLoPtDPhi->Fill(dPhi)
                                      : ThetaPairTrackHiPtDPhi->Fill(dPhi);
            F("unp_trk_pt", trk) < 22 ? ThetaPairTrackLoPtDStrip->Fill(dStrip)
                                      : ThetaPairTrackHiPtDStrip->Fill(dStrip);

            ThetaPairLctStatVRingAllPt->Fill(ring, stat);
            ThetaPairLctStatAllPt->Fill(stat);
            ThetaPairLctQualAllPt->Fill(qual);
            ThetaPairLctThetaAllPt->Fill(theta);

            F("unp_trk_pt", trk) < 22
                ? ThetaPairLctStatVRingLowPt->Fill(ring, stat)
                : ThetaPairLctStatVRingHighPt->Fill(ring, stat);
            F("unp_trk_pt", trk) < 22 ? ThetaPairLctStatLowPt->Fill(stat)
                                      : ThetaPairLctStatHighPt->Fill(stat);
            F("unp_trk_pt", trk) < 22 ? ThetaPairLctQualLowPt->Fill(qual)
                                      : ThetaPairLctQualHighPt->Fill(qual);
            F("unp_trk_pt", trk) < 22 ? ThetaPairLctThetaLowPt->Fill(theta)
                                      : ThetaPairLctThetaHighPt->Fill(theta);
          }
          auto quadP = quads.find(iHit);
          if (quadP != quads.end()) {
            // std::get<4>(*quads[iHit]) += 1;
            trackUsesQuad = true;
            fillOccupancyHit(SingleMuQuadInTrack, iHit);

            QuadLctStatVRingAllPt->Fill(ring, stat);
            QuadLctStatAllPt->Fill(stat);
            QuadLctQualAllPt->Fill(qual);
            QuadLctThetaAllPt->Fill(theta);

            F("unp_trk_pt", trk) < 22
                ? QuadLctStatVRingLowPt->Fill(ring, stat)
                : QuadLctStatVRingHighPt->Fill(ring, stat);
            F("unp_trk_pt", trk) < 22 ? QuadLctStatLowPt->Fill(stat)
                                      : QuadLctStatHighPt->Fill(stat);
            F("unp_trk_pt", trk) < 22 ? QuadLctQualLowPt->Fill(qual)
                                      : QuadLctQualHighPt->Fill(qual);
            F("unp_trk_pt", trk) < 22 ? QuadLctThetaLowPt->Fill(theta)
                                      : QuadLctThetaHighPt->Fill(theta);
          }
        }

        SingleMuPtAll->Fill(F("unp_trk_pt", trk));

        if (trackUsesPhiPair) {
          SingleMuPtAllWPhiPair->Fill(F("unp_trk_pt", trk));
          F("unp_trk_pt", trk) < 22
              ? SingleMuPtLoWPhiPair->Fill(F("unp_trk_pt", trk))
              : SingleMuPtHiWPhiPair->Fill(F("unp_trk_pt", trk));
        }
        if (trackUsesThetaPair) {
          SingleMuPtAllWThetaPair->Fill(F("unp_trk_pt", trk));
          F("unp_trk_pt", trk) < 22
              ? SingleMuPtLoWThetaPair->Fill(F("unp_trk_pt", trk))
              : SingleMuPtHiWThetaPair->Fill(F("unp_trk_pt", trk));
        }
        if (trackUsesQuad) {
          SingleMuPtAllWQuad->Fill(F("unp_trk_pt", trk));
          F("unp_trk_pt", trk) < 22
              ? SingleMuPtLoWQuad->Fill(F("unp_trk_pt", trk))
              : SingleMuPtHiWQuad->Fill(F("unp_trk_pt", trk));
        }

        if (!(trackUsesPhiPair || trackUsesThetaPair || trackUsesQuad)) {
          SingleMuPtAllNoPair->Fill(F("unp_trk_pt", trk));
          F("unp_trk_pt", trk) < 22
              ? SingleMuPtLoNoPair->Fill(F("unp_trk_pt", trk))
              : SingleMuPtHiNoPair->Fill(F("unp_trk_pt", trk));
        }

        if (F("unp_trk_pt", trk) < 22) {
          EventPhiPairTrackBreakdown->Fill(phiPairCount > 0 ? 2 : 0);
          EventThetaPairTrackBreakdown->Fill(thetaPairCount > 0 ? 2 : 0);
          EventQuadTrackBreakdown->Fill(quadCount > 0 ? 2 : 0);

        } else {
          EventPhiPairTrackBreakdown->Fill(phiPairCount > 0 ? 3 : 1);
          EventThetaPairTrackBreakdown->Fill(thetaPairCount > 0 ? 3 : 1);
          EventQuadTrackBreakdown->Fill(quadCount > 0 ? 3 : 1);
        }
      }
    }
  } // End loop: for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Finished looping over the events *******"
            << std::endl;

  delete in_chain;

  std::cout << "\nDone with Read_FlatNtuple(). Exiting.\n" << std::endl;
  f->Write();
} // End function: void Read_FlatNtuple()

const std::map<std::pair<int, int>, int> histIndexCSC = {
    {{1, 4}, 9}, {{1, 1}, 8}, {{1, 2}, 7}, {{1, 3}, 6}, {{2, 1}, 5},
    {{2, 2}, 4}, {{3, 1}, 3}, {{3, 2}, 2}, {{4, 1}, 1}, {{4, 2}, 0}};
// cscid + cscid_offset, endcap * (station - 0.5)
void fillOccupancyHit(TH2D *plot, int hit) { fillOccupancyHit(plot, hit, 1); }

void fillOccupancyHit(TH2D *plot, int hit, int weight) {
  int station = I("hit_station", hit);
  int ring = I("hit_ring", hit);
  int chamber = I("hit_chamber", hit);
  int endcap = I("hit_endcap", hit);
  int sector = I("hit_sector", hit);

  int chamber_bin_index;
  if (station > 1 && (ring % 2) == 1) {
    chamber_bin_index = (chamber * 2) + ((chamber + 1) / 3);
  } else {
    chamber_bin_index = chamber + ((chamber + 3) / 6);
  }

  if (I("hit_isCSC", hit) != 1) {
    std::cout << "Hit isn't CSC: " << station << " " << ring << std::endl;
  }
  if (histIndexCSC.find({station, ring}) == histIndexCSC.end()) {
    std::cout << "Strange station, ring: " << station << " " << ring
              << std::endl;
    return;
  }

  int hist_index = histIndexCSC.at({station, ring});
  if (endcap > 0)
    hist_index = 19 - hist_index;

  if (I("hit_neighbor", hit) == 1) {
    plot->Fill(sector * 7 - 4, hist_index, weight);
  } else {
    plot->Fill(chamber_bin_index, hist_index, weight);
    if (station > 1 && (ring % 2) == 1) {
      plot->Fill(chamber_bin_index - 1, hist_index, weight);
    }
  }
}

bool isPhiPair(int i, int j) {
  return I("hit_endcap", i) == I("hit_endcap", j) &&
         I("hit_station", i) == I("hit_station", j) &&
         I("hit_ring", i) == I("hit_ring", j) &&
         I("hit_chamber", i) == I("hit_chamber", j) &&
         I("hit_BX", i) == I("hit_BX", j) &&
         I("hit_phi_int", i) == I("hit_phi_int", j) &&
         I("hit_theta_int", i) != I("hit_theta_int", j) &&
         I("hit_neighbor", i) == I("hit_neighbor", j);
}

bool isThetaPair(int i, int j) {
  return I("hit_endcap", i) == I("hit_endcap", j) &&
         I("hit_station", i) == I("hit_station", j) &&
         I("hit_ring", i) == I("hit_ring", j) &&
         I("hit_chamber", i) == I("hit_chamber", j) &&
         I("hit_BX", i) == I("hit_BX", j) &&
         I("hit_theta_int", i) == I("hit_theta_int", j) &&
         I("hit_phi_int", i) != I("hit_phi_int", j) &&
         I("hit_neighbor", i) == I("hit_neighbor", j);
}
