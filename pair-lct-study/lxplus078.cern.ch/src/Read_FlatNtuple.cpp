
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
    case 3:
      file_name.Form(
          "%s/%s/"
          "NTuple_ZeroBias2_FlatNtuple_Run_306091_2018_02_24_ZB2.root",
          store.Data(), in_dir.Data());
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
  TChain *in_chain = new TChain("FlatNtupleData/tree");
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
                << I("nSimHits") << " emulated EMTF hits in the event"
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

    for (int i = 0; i < I("nSimHits"); i++) {
      if (I("sim_hit_isCSC", i) != 1)
        continue;
      fillOccupancyHit(SingleMuOccupancy, i);

      int station = I("sim_hit_station", i);
      int ring = I("sim_hit_ring", i);
      int qual = I("sim_hit_quality", i);
      int alctQual = I("sim_hit_alct_quality", i);
      int clctQual = I("sim_hit_clct_quality", i);
      float theta = F("sim_hit_theta", i);

      AllLctStatVRing->Fill(ring, station);
      AllLctStat->Fill(station);
      AllLctQual->Fill(qual);
      AllLctTheta->Fill(theta);

      AllAlctQuality->Fill(alctQual);
      AllClctQuality->Fill(clctQual);

      AllLctAlctQualityByStat[station - 1]->Fill(alctQual);
      AllLctClctQualityByStat[station - 1]->Fill(clctQual);
      AllLctAlctVsClctByStat[station - 1]->Fill(clctQual, alctQual);
      AllLctAlctQualVsStat->Fill(station, alctQual);
      AllLctClctQualVsStat->Fill(station, clctQual);

      for (int j = i + 1; j < I("nSimHits"); j++) {
        if (I("sim_hit_isCSC", j) != 1)
          continue;

        float thetaJ = F("sim_hit_theta", j);
        bool isPhiP = isPhiPair(i, j);
        bool isThetaP = isThetaPair(i, j);
        if (isPhiP) {
          fillOccupancyHit(SingleMuPhiPairOccupancy, i, 2);
          phiPairCount += 1;
          std::pair<int, int> p1(j, 0);
          std::pair<int, int> p2(i, 0);
          phiPairs[i] = p1;
          phiPairs[j] = p2;

          PhiPairDTheta->Fill(abs(theta - thetaJ));
          PhiPairDWire->Fill(abs(I("sim_hit_wire", i) - I("sim_hit_wire", j)));

          PhiPairLctStatVRing->Fill(ring, station);
          PhiPairLctStatVRing->Fill(ring, station);

          PhiPairLctStat->Fill(station);
          PhiPairLctStat->Fill(station);

          PhiPairLctQual->Fill(qual);
          PhiPairLctQual->Fill(I("sim_hit_quality", j));

          PhiPairAlctQuality->Fill(alctQual);
          PhiPairAlctQuality->Fill(I("sim_hit_alct_quality", j));
          PhiPairClctQuality->Fill(clctQual);
          PhiPairClctQuality->Fill(I("sim_hit_clct_quality", j));

          PhiPairLctTheta->Fill(theta);
          PhiPairLctTheta->Fill(thetaJ);

          PhiPairAlctQualityByStat[station - 1]->Fill(alctQual);
          PhiPairAlctQualityByStat[station - 1]->Fill(
              I("sim_hit_alct_quality", j));
          PhiPairClctQualityByStat[station - 1]->Fill(clctQual);
          PhiPairClctQualityByStat[station - 1]->Fill(
              I("sim_hit_alct_quality", j));
          PhiPairAlctVsClctByStat[station - 1]->Fill(clctQual, alctQual);
          PhiPairAlctVsClctByStat[station - 1]->Fill(
              I("sim_hit_clct_quality", j), I("sim_hit_alct_quality", j));
          PhiPairAlctQualVsStat->Fill(station, alctQual);
          PhiPairAlctQualVsStat->Fill(station, I("sim_hit_alct_quality", j));
          PhiPairClctQualVsStat->Fill(station, clctQual);
          PhiPairClctQualVsStat->Fill(station, I("sim_hit_clct_quality", j));

          PhiPairQual1VsQual2ByStat[station - 1]->Fill(I("sim_hit_quality", i),
                                                       I("sim_hit_quality", j));
          PhiPairAlct1VsAlct2ByStat[station - 1]->Fill(
              I("sim_hit_alct_quality", i), I("sim_hit_alct_quality", j));
          PhiPairClct1VsClct2ByStat[station - 1]->Fill(
              I("sim_hit_clct_quality", i), I("sim_hit_clct_quality", j));
        }
        if (isThetaP) {
          fillOccupancyHit(SingleMuThetaPairOccupancy, i, 2);
          thetaPairCount += 1;
          std::pair<int, int> p1(j, 0);
          std::pair<int, int> p2(i, 0);
          thetaPairs[i] = p1;
          thetaPairs[j] = p2;

          ThetaPairDPhi->Fill(
              abs(I("sim_hit_phi_int", i) - I("sim_hit_phi_int", j)));
          ThetaPairDStrip->Fill(
              abs(I("sim_hit_strip", i) - I("sim_hit_strip", j)));

          ThetaPairLctStatVRing->Fill(ring, station);
          ThetaPairLctStatVRing->Fill(ring, station);

          ThetaPairLctStat->Fill(station);
          ThetaPairLctStat->Fill(station);

          ThetaPairLctQual->Fill(qual);
          ThetaPairLctQual->Fill(I("sim_hit_quality", j));

          ThetaPairAlctQuality->Fill(alctQual);
          ThetaPairAlctQuality->Fill(I("sim_hit_alct_quality", j));
          ThetaPairClctQuality->Fill(clctQual);
          ThetaPairClctQuality->Fill(I("sim_hit_clct_quality", j));

          ThetaPairLctTheta->Fill(theta);
          ThetaPairLctTheta->Fill(thetaJ);

          ThetaPairAlctQualityByStat[station - 1]->Fill(alctQual);
          ThetaPairAlctQualityByStat[station - 1]->Fill(
              I("sim_hit_alct_quality", j));
          ThetaPairClctQualityByStat[station - 1]->Fill(clctQual);
          ThetaPairClctQualityByStat[station - 1]->Fill(
              I("sim_hit_alct_quality", j));
          ThetaPairAlctVsClctByStat[station - 1]->Fill(clctQual, alctQual);
          ThetaPairAlctVsClctByStat[station - 1]->Fill(
              I("sim_hit_clct_quality", j), I("sim_hit_alct_quality", j));
          ThetaPairAlctQualVsStat->Fill(station, alctQual);
          ThetaPairAlctQualVsStat->Fill(station, I("sim_hit_alct_quality", j));
          ThetaPairClctQualVsStat->Fill(station, clctQual);
          ThetaPairClctQualVsStat->Fill(station, I("sim_hit_clct_quality", j));

          ThetaPairQual1VsQual2ByStat[station - 1]->Fill(
              I("sim_hit_quality", i), I("sim_hit_quality", j));
          ThetaPairAlct1VsAlct2ByStat[station - 1]->Fill(
              I("sim_hit_alct_quality", i), I("sim_hit_alct_quality", j));
          ThetaPairClct1VsClct2ByStat[station - 1]->Fill(
              I("sim_hit_clct_quality", i), I("sim_hit_clct_quality", j));
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

        int station = I("sim_hit_station", phiHitA);

        QuadLctStatVRing->Fill(I("sim_hit_ring", phiHitA),
                               I("sim_hit_station", phiHitA), 4);
        QuadLctStat->Fill(I("sim_hit_station", phiHitA));
        QuadLctStat->Fill(I("sim_hit_station", phiHitA));
        QuadLctStat->Fill(I("sim_hit_station", phiHitA));
        QuadLctStat->Fill(I("sim_hit_station", phiHitA));

        QuadLctQual->Fill(I("sim_hit_quality", phiHitA));
        QuadLctQual->Fill(I("sim_hit_quality", phiHitB));
        QuadLctQual->Fill(I("sim_hit_quality", hit3));
        QuadLctQual->Fill(I("sim_hit_quality", hit4));

        QuadAlctQuality->Fill(I("sim_hit_alct_quality", phiHitA));
        QuadAlctQuality->Fill(I("sim_hit_alct_quality", phiHitB));
        QuadAlctQuality->Fill(I("sim_hit_alct_quality", hit3));
        QuadAlctQuality->Fill(I("sim_hit_alct_quality", hit4));

        QuadClctQuality->Fill(I("sim_hit_clct_quality", phiHitA));
        QuadClctQuality->Fill(I("sim_hit_clct_quality", phiHitB));
        QuadClctQuality->Fill(I("sim_hit_clct_quality", hit3));
        QuadClctQuality->Fill(I("sim_hit_clct_quality", hit4));

        QuadLctTheta->Fill(I("sim_hit_theta_int", phiHitA));
        QuadLctTheta->Fill(I("sim_hit_theta_int", phiHitB));
        QuadLctTheta->Fill(I("sim_hit_theta_int", hit3));
        QuadLctTheta->Fill(I("sim_hit_theta_int", hit4));

        QuadAlctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", phiHitA));
        QuadAlctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", phiHitB));
        QuadAlctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", hit3));
        QuadAlctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", hit4));

        QuadClctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", phiHitA));
        QuadClctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", phiHitB));
        QuadClctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", hit3));
        QuadClctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", hit4));

        QuadAlctVsClctByStat[station - 1]->Fill(
            I("sim_hit_clct_quality", phiHitA),
            I("sim_hit_alct_quality", phiHitA));
        QuadAlctVsClctByStat[station - 1]->Fill(
            I("sim_hit_clct_quality", phiHitB),
            I("sim_hit_alct_quality", phiHitB));
        QuadAlctVsClctByStat[station - 1]->Fill(
            I("sim_hit_clct_quality", hit3), I("sim_hit_alct_quality", hit3));
        QuadAlctVsClctByStat[station - 1]->Fill(
            I("sim_hit_clct_quality", hit4), I("sim_hit_alct_quality", hit4));

        QuadAlctQualVsStat->Fill(station, I("sim_hit_alct_quality", phiHitA));
        QuadAlctQualVsStat->Fill(station, I("sim_hit_alct_quality", phiHitB));
        QuadAlctQualVsStat->Fill(station, I("sim_hit_alct_quality", hit3));
        QuadAlctQualVsStat->Fill(station, I("sim_hit_alct_quality", hit4));

        QuadClctQualVsStat->Fill(station, I("sim_hit_clct_quality", phiHitA));
        QuadClctQualVsStat->Fill(station, I("sim_hit_clct_quality", phiHitB));
        QuadClctQualVsStat->Fill(station, I("sim_hit_clct_quality", hit3));
        QuadClctQualVsStat->Fill(station, I("sim_hit_clct_quality", hit4));
      }
    }
    EventQuadCount->Fill(quadCount);

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
  int station = I("sim_hit_station", hit);
  int ring = I("sim_hit_ring", hit);
  int chamber = I("sim_hit_chamber", hit);
  int endcap = I("sim_hit_endcap", hit);
  int sector = I("sim_hit_sector", hit);

  int chamber_bin_index;
  if (station > 1 && (ring % 2) == 1) {
    chamber_bin_index = (chamber * 2) + ((chamber + 1) / 3);
  } else {
    chamber_bin_index = chamber + ((chamber + 3) / 6);
  }

  if (I("sim_hit_isCSC", hit) != 1) {
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

  if (I("sim_hit_neighbor", hit) == 1) {
    plot->Fill(sector * 7 - 4, hist_index, weight);
  } else {
    plot->Fill(chamber_bin_index, hist_index, weight);
    if (station > 1 && (ring % 2) == 1) {
      plot->Fill(chamber_bin_index - 1, hist_index, weight);
    }
  }
}

bool isPhiPair(int i, int j) {
  return I("sim_hit_endcap", i) == I("sim_hit_endcap", j) &&
         I("sim_hit_station", i) == I("sim_hit_station", j) &&
         I("sim_hit_ring", i) == I("sim_hit_ring", j) &&
         I("sim_hit_chamber", i) == I("sim_hit_chamber", j) &&
         I("sim_hit_BX", i) == I("sim_hit_BX", j) &&
         I("sim_hit_phi_int", i) == I("sim_hit_phi_int", j) &&
         I("sim_hit_theta_int", i) != I("sim_hit_theta_int", j) &&
         I("sim_hit_neighbor", i) == I("sim_hit_neighbor", j);
}

bool isThetaPair(int i, int j) {
  return I("sim_hit_endcap", i) == I("sim_hit_endcap", j) &&
         I("sim_hit_station", i) == I("sim_hit_station", j) &&
         I("sim_hit_ring", i) == I("sim_hit_ring", j) &&
         I("sim_hit_chamber", i) == I("sim_hit_chamber", j) &&
         I("sim_hit_BX", i) == I("sim_hit_BX", j) &&
         I("sim_hit_theta_int", i) == I("sim_hit_theta_int", j) &&
         I("sim_hit_phi_int", i) != I("sim_hit_phi_int", j) &&
         I("sim_hit_neighbor", i) == I("sim_hit_neighbor", j);
}
