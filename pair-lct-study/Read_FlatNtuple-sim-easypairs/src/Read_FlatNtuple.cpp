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

#include "Helpers.cpp"

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
                     bool skim = false, bool recoTrackCut = false,
                     int maxEvents = MAX_EVT, int maxFiles = MAX_FILES) {
  TFile *f = new TFile(outFile, "RECREATE");
#include "Plots.h"

  // BookHistos();
  gROOT->SetStyle("Plain");

  auto in_chain = prepareInput(dataset, skim);
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

    std::map<int, int> simHitsInTrackRecoMatch;

    if (recoTrackCut) {
      for (int trk = 0; trk < I("nTracks"); trk++) {
        for (int hitIdx = 0; hitIdx < I("trk_nHits", trk); hitIdx++) {
          int hit = I("trk_iHit", trk, hitIdx);
          if (hit < 0)
            continue;
          int simHit = I("hit_match_iSimHit", hit);
          if (simHit < 0)
            continue;
          bool simHitIsCSC = I("sim_hit_isCSC", simHit) == 1;
          int simHitBx = I("sim_hit_BX", simHit);
          int trkBx = I("trk_BX", trk);
          int trkMode = I("trk_mode", trk);
          bool recoMatch = I("trk_dR_match_unique", trk) == 1;
          if (!recoMatch)
            continue;
          int recoMuon = I("trk_dR_match_iReco", trk);
          float recoPt = F("reco_pt", recoMuon);
          int recoIdMed = I("reco_ID_medium", recoMuon);
          int recoIdStat = I("reco_ID_station", recoMuon);

          if (simHit >= 0 && simHitIsCSC && simHitBx == 0) {
            if (trkBx == 0 && trkMode > 10 && trkMode < 16) {
              if (recoMatch && recoPt >= 10 && recoPt <= 20) {
                if (recoIdMed == 1 && recoIdStat == 1) {
                  simHitsInTrackRecoMatch[simHit] = trk;
                }
              }
            }
          }
        }
      }
    }

    for (int i = 0; i < I("nSimHits"); i++) {
      if (I("sim_hit_isCSC", i) != 1)
        continue;

      if (recoTrackCut &&
          simHitsInTrackRecoMatch.find(i) == simHitsInTrackRecoMatch.end()) {
        continue;
      }
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
        //if (recoTrackCut &&
            //simHitsInTrackRecoMatch.find(i) == simHitsInTrackRecoMatch.end()) {
          //continue;
        //}

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
              I("sim_hit_clct_quality", j));
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
              I("sim_hit_alct_quality", j), I("sim_hit_alct_quality", i));
          PhiPairClct1VsClct2ByStat[station - 1]->Fill(
              I("sim_hit_clct_quality", j), I("sim_hit_clct_quality", i));
        }
        if (isThetaP) {
          fillOccupancyHit(SingleMuThetaPairOccupancy, i, 2);
          thetaPairCount += 1;
          std::pair<int, int> p1(j, 0);
          std::pair<int, int> p2(i, 0);
          thetaPairs[i] = p1;
          thetaPairs[j] = p2;

          ThetaPairDPhi->Fill(abs(F("sim_hit_phi", i) - F("sim_hit_phi", j)));
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
              I("sim_hit_clct_quality", j));
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
              I("sim_hit_alct_quality", j), I("sim_hit_alct_quality", i));
          ThetaPairClct1VsClct2ByStat[station - 1]->Fill(
              I("sim_hit_clct_quality", j), I("sim_hit_clct_quality", i));
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

        QuadLctTheta->Fill(F("sim_hit_theta", phiHitA));
        QuadLctTheta->Fill(F("sim_hit_theta", phiHitB));
        QuadLctTheta->Fill(F("sim_hit_theta", hit3));
        QuadLctTheta->Fill(F("sim_hit_theta", hit4));

        QuadAlctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", phiHitA));
        QuadAlctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", phiHitB));
        QuadAlctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", hit3));
        QuadAlctQualityByStat[station - 1]->Fill(
            I("sim_hit_alct_quality", hit4));

        QuadClctQualityByStat[station - 1]->Fill(
            I("sim_hit_clct_quality", phiHitA));
        QuadClctQualityByStat[station - 1]->Fill(
            I("sim_hit_clct_quality", phiHitB));
        QuadClctQualityByStat[station - 1]->Fill(
            I("sim_hit_clct_quality", hit3));
        QuadClctQualityByStat[station - 1]->Fill(
            I("sim_hit_clct_quality", hit4));

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

    for (int trk = 0; trk < I("nTracks"); trk++) {
      bool trackUsesPhiPair = false;
      bool trackUsesThetaPair = false;
      bool trackUsesQuad = false;

      bool recoMatch = I("trk_dR_match_unique", trk) == 1;
      if (recoTrackCut && !recoMatch) continue;

        if (I("trk_nHits", trk) >= 3) {
          // Loop over hits in this tracks
          for (int trkHitIdx = 0; trkHitIdx < I("trk_nHits", trk);
               trkHitIdx++) {
            int trkHit = I("trk_iHit", trk, trkHitIdx);
            if (trkHit < 0)
              continue;
            try {
              int iHit = I("hit_match_iSimHit", trkHit);
              if (iHit < 0)
                continue;
              if (recoTrackCut && simHitsInTrackRecoMatch.find(iHit) ==
                                      simHitsInTrackRecoMatch.end())
                continue;
              if (I("sim_hit_isCSC", iHit) != 1) {
                continue;
              }
              int ring = I("sim_hit_ring", iHit);
              int stat = I("sim_hit_station", iHit);
              int qual = I("sim_hit_quality", iHit);
              float theta = F("sim_hit_theta", iHit);

              AllLctStatVRingAllPt->Fill(ring, stat);
              AllLctStatAllPt->Fill(stat);
              AllLctQualAllPt->Fill(qual);
              AllLctThetaAllPt->Fill(theta);

              F("trk_pt", trk) < 22 ? AllLctStatVRingLowPt->Fill(ring, stat)
                                    : AllLctStatVRingHighPt->Fill(ring, stat);
              F("trk_pt", trk) < 22 ? AllLctStatLowPt->Fill(stat)
                                    : AllLctStatHighPt->Fill(stat);
              F("trk_pt", trk) < 22 ? AllLctQualLowPt->Fill(qual)
                                    : AllLctQualHighPt->Fill(qual);
              F("trk_pt", trk) < 22 ? AllLctThetaLowPt->Fill(theta)
                                    : AllLctThetaHighPt->Fill(theta);

              auto phiP = phiPairs.find(iHit);
              if (phiP != phiPairs.end()) {
                phiPairs[iHit].second += 1;
                trackUsesPhiPair = true;
                fillOccupancyHit(SingleMuPhiPairInTrack, iHit);

                float dTheta = abs(F("sim_hit_theta", iHit) -
                                   F("sim_hit_theta", phiPairs[iHit].first));
                int dWire = abs(I("sim_hit_wire", iHit) -
                                I("sim_hit_wire", phiPairs[iHit].first));

                PhiPairTrackAllPtDTheta->Fill(dTheta);
                PhiPairTrackAllPtDWire->Fill(dWire);
                F("trk_pt", trk) < 22 ? PhiPairTrackLoPtDTheta->Fill(dTheta)
                                      : PhiPairTrackHiPtDTheta->Fill(dTheta);
                F("trk_pt", trk) < 22 ? PhiPairTrackLoPtDWire->Fill(dWire)
                                      : PhiPairTrackHiPtDWire->Fill(dWire);

                PhiPairLctStatVRingAllPt->Fill(ring, stat);
                PhiPairLctStatAllPt->Fill(stat);
                PhiPairLctQualAllPt->Fill(qual);
                PhiPairLctThetaAllPt->Fill(theta);

                F("trk_pt", trk) < 22
                    ? PhiPairLctStatVRingLowPt->Fill(ring, stat)
                    : PhiPairLctStatVRingHighPt->Fill(ring, stat);
                F("trk_pt", trk) < 22 ? PhiPairLctStatLowPt->Fill(stat)
                                      : PhiPairLctStatHighPt->Fill(stat);
                F("trk_pt", trk) < 22 ? PhiPairLctQualLowPt->Fill(qual)
                                      : PhiPairLctQualHighPt->Fill(qual);
                F("trk_pt", trk) < 22 ? PhiPairLctThetaLowPt->Fill(theta)
                                      : PhiPairLctThetaHighPt->Fill(theta);
              }
              auto thetaP = thetaPairs.find(iHit);
              if (thetaP != thetaPairs.end()) {
                thetaPairs[iHit].second += 1;
                trackUsesThetaPair = true;
                fillOccupancyHit(SingleMuThetaPairInTrack, iHit);

                float dPhi = abs(F("sim_hit_phi", iHit) -
                                 F("sim_hit_phi", thetaPairs[iHit].first));
                int dStrip = abs(I("sim_hit_strip", iHit) -
                                 I("sim_hit_strip", thetaPairs[iHit].first));

                ThetaPairTrackAllPtDPhi->Fill(dPhi);
                ThetaPairTrackAllPtDStrip->Fill(dStrip);
                F("trk_pt", trk) < 22 ? ThetaPairTrackLoPtDPhi->Fill(dPhi)
                                      : ThetaPairTrackHiPtDPhi->Fill(dPhi);
                F("trk_pt", trk) < 22 ? ThetaPairTrackLoPtDStrip->Fill(dStrip)
                                      : ThetaPairTrackHiPtDStrip->Fill(dStrip);

                ThetaPairLctStatVRingAllPt->Fill(ring, stat);
                ThetaPairLctStatAllPt->Fill(stat);
                ThetaPairLctQualAllPt->Fill(qual);
                ThetaPairLctThetaAllPt->Fill(theta);

                F("trk_pt", trk) < 22
                    ? ThetaPairLctStatVRingLowPt->Fill(ring, stat)
                    : ThetaPairLctStatVRingHighPt->Fill(ring, stat);
                F("trk_pt", trk) < 22 ? ThetaPairLctStatLowPt->Fill(stat)
                                      : ThetaPairLctStatHighPt->Fill(stat);
                F("trk_pt", trk) < 22 ? ThetaPairLctQualLowPt->Fill(qual)
                                      : ThetaPairLctQualHighPt->Fill(qual);
                F("trk_pt", trk) < 22 ? ThetaPairLctThetaLowPt->Fill(theta)
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

                F("trk_pt", trk) < 22
                    ? QuadLctStatVRingLowPt->Fill(ring, stat)
                    : QuadLctStatVRingHighPt->Fill(ring, stat);
                F("trk_pt", trk) < 22 ? QuadLctStatLowPt->Fill(stat)
                                      : QuadLctStatHighPt->Fill(stat);
                F("trk_pt", trk) < 22 ? QuadLctQualLowPt->Fill(qual)
                                      : QuadLctQualHighPt->Fill(qual);
                F("trk_pt", trk) < 22 ? QuadLctThetaLowPt->Fill(theta)
                                      : QuadLctThetaHighPt->Fill(theta);
              }
            } catch (const std::invalid_argument &e) {
              continue;
            }
          }

          SingleMuPtAll->Fill(F("trk_pt", trk));

          if (trackUsesPhiPair) {
            SingleMuPtAllWPhiPair->Fill(F("trk_pt", trk));
            F("trk_pt", trk) < 22
                ? SingleMuPtLoWPhiPair->Fill(F("trk_pt", trk))
                : SingleMuPtHiWPhiPair->Fill(F("trk_pt", trk));
          }
          if (trackUsesThetaPair) {
            SingleMuPtAllWThetaPair->Fill(F("trk_pt", trk));
            F("trk_pt", trk) < 22
                ? SingleMuPtLoWThetaPair->Fill(F("trk_pt", trk))
                : SingleMuPtHiWThetaPair->Fill(F("trk_pt", trk));
          }
          if (trackUsesQuad) {
            SingleMuPtAllWQuad->Fill(F("trk_pt", trk));
            F("trk_pt", trk) < 22 ? SingleMuPtLoWQuad->Fill(F("trk_pt", trk))
                                  : SingleMuPtHiWQuad->Fill(F("trk_pt", trk));
          }

          if (!(trackUsesPhiPair || trackUsesThetaPair || trackUsesQuad)) {
            SingleMuPtAllNoPair->Fill(F("trk_pt", trk));
            F("trk_pt", trk) < 22 ? SingleMuPtLoNoPair->Fill(F("trk_pt", trk))
                                  : SingleMuPtHiNoPair->Fill(F("trk_pt", trk));
          }

          if (F("trk_pt", trk) < 22) {
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
