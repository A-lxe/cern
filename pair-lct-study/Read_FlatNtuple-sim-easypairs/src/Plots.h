/////Alex's PU rate study/////
/////      :)START(:      ////
//////////////////////////////

// Event Breakdowns
f->mkdir("Event Breakdowns")->cd();

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

// Track Cuts
f->mkdir("Track Cuts")->cd();

TH1D *SingleMuPtAll =
    new TH1D("SingleMuonPtAll", "Single Muon Tracks", 100, 1, 101);

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

// Occupancy
f->mkdir("Occupancies")->cd();

TH2D *SingleMuOccupancy =
    new TH2D("ChamberOccupancy", "ChamberOccupancy", 42, 1, 43, 20, 0, 20);
TH2D *SingleMuPhiPairOccupancy = new TH2D(
    "PhiPairChamberOccupancy", "PhiPairChamberOccupancy", 42, 1, 43, 20, 0, 20);
TH2D *SingleMuPhiPairInTrack = new TH2D(
    "SingleMuTrackPhiPair", "SingleMuTrackPhiPair", 42, 1, 43, 20, 0, 20);
TH2D *SingleMuThetaPairOccupancy =
    new TH2D("ThetaPairChamberOccupancy", "ThetaPairChamberOccupancy", 42, 1,
             43, 20, 0, 20);
TH2D *SingleMuThetaPairInTrack = new TH2D(
    "SingleMuTrackThetaPair", "SingleMuTrackThetaPair", 42, 1, 43, 20, 0, 20);
TH2D *SingleMuQuadOccupancy = new TH2D(
    "QuadChamberOccupancy", "QuadChamberOccupancy", 42, 1, 43, 20, 0, 20);
TH2D *SingleMuQuadInTrack =
    new TH2D("SingleMuTrackQuad", "SingleMuTrackQuad", 42, 1, 43, 20, 0, 20);

TH2D *occupancyPlots[7] = {SingleMuOccupancy,        SingleMuPhiPairOccupancy,
                           SingleMuPhiPairInTrack,   SingleMuThetaPairOccupancy,
                           SingleMuThetaPairInTrack, SingleMuQuadOccupancy,
                           SingleMuQuadInTrack};
const std::array<std::string, 10> suffix_label{
    {"4/2", "4/1", "3/2", "3/1", " 2/2", "2/1", "1/3", "1/2", "1/1b", "1/1a"}};

for (auto op : occupancyPlots) {
  op->GetXaxis()->SetTitle("10#circ Chamber (N=neighbor)");
  int count = 0;
  for (int xbin = 1; xbin < 43; ++xbin) {
    op->GetXaxis()->SetBinLabel(xbin, std::to_string(xbin - count).c_str());
    if (xbin == 2 || xbin == 9 || xbin == 16 || xbin == 23 || xbin == 30 ||
        xbin == 37) {
      ++xbin;
      ++count;
      op->GetXaxis()->SetBinLabel(xbin, "N");
    }
  }
  for (int ybin = 1; ybin <= 10; ++ybin) {
    op->GetYaxis()->SetBinLabel(ybin, ("ME-" + suffix_label[ybin - 1]).c_str());
    op->GetYaxis()->SetBinLabel(21 - ybin,
                                ("ME+" + suffix_label[ybin - 1]).c_str());
  }
  // op->GetXaxis()->SetCanExtend(false); // Needed to stop multi-thread summing
}

// LCT Plots
TDirectory *lctStatDir = f->mkdir("LCT Stats");
lctStatDir->cd();

TH2D *AllLctStatVRing =
    new TH2D("AllLctStatVRing", "All LCTs Station vs Ring", 4, 1, 5, 4, 1, 5);
TH2D *PhiPairLctStatVRing = new TH2D(
    "PhiPairLctStatVRing", "Phi Pair LCTs Station vs Ring", 4, 1, 5, 4, 1, 5);
TH2D *ThetaPairLctStatVRing =
    new TH2D("ThetaPairLctStatVRing", "Theta Pair LCTs Station vs Ring", 4, 1,
             5, 4, 1, 5);
TH2D *QuadLctStatVRing =
    new TH2D("QuadLctStatVRing", "Quad LCTs Station vs Ring", 4, 1, 5, 4, 1, 5);

TH1D *AllLctStat = new TH1D("AllLctStat", "All LCTs Station", 4, 1, 5);
TH1D *PhiPairLctStat =
    new TH1D("PhiPairLctStat", "Phi Pair LCTs Station", 4, 1, 5);
TH1D *ThetaPairLctStat =
    new TH1D("ThetaPairLctStat", "Theta Pair LCTs Station", 4, 1, 5);
TH1D *QuadLctStat = new TH1D("QuadLctStat", "Quad LCTs Station", 4, 1, 5);

TH1D *AllLctQual = new TH1D("AllLctQual", "All LCTs Quality", 16, 0, 16);
TH1D *PhiPairLctQual =
    new TH1D("PhiPairLctQual", "Phi Pair LCTs Quality", 16, 0, 16);
TH1D *ThetaPairLctQual =
    new TH1D("ThetaPairLctQual", "Theta Pair LCTs Quality", 16, 0, 16);
TH1D *QuadLctQual = new TH1D("QuadLctQual", "Quad LCTs Quality", 16, 0, 16);

TH1D *AllLctTheta = new TH1D("AllLctTheta", "All LCTs Theta", 100, 0, 50);
TH1D *PhiPairLctTheta =
    new TH1D("PhiPairLctTheta", "Phi Pair LCTs Theta", 100, 0, 50);
TH1D *ThetaPairLctTheta =
    new TH1D("ThetaPairLctTheta", "Theta Pair LCTs Theta", 100, 0, 50);
TH1D *QuadLctTheta = new TH1D("QuadLctTheta", "Quad LCTs Theta", 100, 0, 50);

// All pT
lctStatDir->mkdir("All pT")->cd();

TH2D *AllLctStatVRingAllPt =
    new TH2D("AllLctStatVRingAllPt", "All LCTs Station vs Ring All pT", 4, 1, 5,
             4, 1, 5);
TH2D *PhiPairLctStatVRingAllPt =
    new TH2D("PhiPairLctStatVRingAllPt", "Phi Pair LCTs Station vs Ring All pT",
             4, 1, 5, 4, 1, 5);
TH2D *ThetaPairLctStatVRingAllPt =
    new TH2D("ThetaPairLctStatVRingAllPt",
             "Theta Pair LCTs Station vs Ring All pT", 4, 1, 5, 4, 1, 5);
TH2D *QuadLctStatVRingAllPt =
    new TH2D("QuadLctStatVRingAllPt", "Quad LCTs Station vs Ring All pT", 4, 1,
             5, 4, 1, 5);

TH1D *AllLctStatAllPt =
    new TH1D("AllLctStatAllPt", "All LCTs Station All pT", 4, 1, 5);
TH1D *PhiPairLctStatAllPt =
    new TH1D("PhiPairLctStatAllPt", "Phi Pair LCTs Station All pT", 4, 1, 5);
TH1D *ThetaPairLctStatAllPt = new TH1D(
    "ThetaPairLctStatAllPt", "Theta Pair LCTs Station All pT", 4, 1, 5);
TH1D *QuadLctStatAllPt =
    new TH1D("QuadLctStatAllPt", "Quad LCTs Station All pT", 4, 1, 5);

TH1D *AllLctQualAllPt =
    new TH1D("AllLctQualAllPt", "All LCTs Quality All pT", 16, 0, 16);
TH1D *PhiPairLctQualAllPt =
    new TH1D("PhiPairLctQualAllPt", "Phi Pair LCTs Quality All pT", 16, 0, 16);
TH1D *ThetaPairLctQualAllPt = new TH1D(
    "ThetaPairLctQualAllPt", "Theta Pair LCTs Quality All pT", 16, 0, 16);
TH1D *QuadLctQualAllPt =
    new TH1D("QuadLctQualAllPt", "Quad LCTs Quality All pT", 16, 0, 16);

TH1D *AllLctThetaAllPt =
    new TH1D("AllLctThetaAllPt", "All LCTs Theta All pT", 100, 0, 50);
TH1D *PhiPairLctThetaAllPt =
    new TH1D("PhiPairLctThetaAllPt", "Phi Pair LCTs Theta All pT", 100, 0, 50);
TH1D *ThetaPairLctThetaAllPt = new TH1D(
    "ThetaPairLctThetaAllPt", "Theta Pair LCTs Theta All pT", 100, 0, 50);
TH1D *QuadLctThetaAllPt =
    new TH1D("QuadLctThetaAllPt", "Quad LCTs Theta All pT", 100, 0, 50);

// Low pT
lctStatDir->mkdir("Low pT")->cd();

TH2D *AllLctStatVRingLowPt =
    new TH2D("AllLctStatVRingLowPt", "All LCTs Station vs Ring Low pT", 4, 1, 5,
             4, 1, 5);
TH2D *PhiPairLctStatVRingLowPt =
    new TH2D("PhiPairLctStatVRingLowPt", "Phi Pair LCTs Station vs Ring Low pT",
             4, 1, 5, 4, 1, 5);
TH2D *ThetaPairLctStatVRingLowPt =
    new TH2D("ThetaPairLctStatVRingLowPt",
             "Theta Pair LCTs Station vs Ring Low pT", 4, 1, 5, 4, 1, 5);
TH2D *QuadLctStatVRingLowPt =
    new TH2D("QuadLctStatVRingLowPt", "Quad LCTs Station vs Ring Low pT", 4, 1,
             5, 4, 1, 5);

TH1D *AllLctStatLowPt =
    new TH1D("AllLctStatLowPt", "All LCTs Station Low pT", 4, 1, 5);
TH1D *PhiPairLctStatLowPt =
    new TH1D("PhiPairLctStatLowPt", "Phi Pair LCTs Station Low pT", 4, 1, 5);
TH1D *ThetaPairLctStatLowPt = new TH1D(
    "ThetaPairLctStatLowPt", "Theta Pair LCTs Station Low pT", 4, 1, 5);
TH1D *QuadLctStatLowPt =
    new TH1D("QuadLctStatLowPt", "Quad LCTs Station Low pT", 4, 1, 5);

TH1D *AllLctQualLowPt =
    new TH1D("AllLctQualLowPt", "All LCTs Quality Low pT", 16, 0, 16);
TH1D *PhiPairLctQualLowPt =
    new TH1D("PhiPairLctQualLowPt", "Phi Pair LCTs Quality Low pT", 16, 0, 16);
TH1D *ThetaPairLctQualLowPt = new TH1D(
    "ThetaPairLctQualLowPt", "Theta Pair LCTs Quality Low pT", 16, 0, 16);
TH1D *QuadLctQualLowPt =
    new TH1D("QuadLctQualLowPt", "Quad LCTs Quality Low pT", 16, 0, 16);

TH1D *AllLctThetaLowPt =
    new TH1D("AllLctThetaLowPt", "All LCTs Theta Low pT", 100, 0, 50);
TH1D *PhiPairLctThetaLowPt =
    new TH1D("PhiPairLctThetaLowPt", "Phi Pair LCTs Theta Low pT", 100, 0, 50);
TH1D *ThetaPairLctThetaLowPt = new TH1D(
    "ThetaPairLctThetaLowPt", "Theta Pair LCTs Theta Low pT", 100, 0, 50);
TH1D *QuadLctThetaLowPt =
    new TH1D("QuadLctThetaLowPt", "Quad LCTs Theta Low pT", 100, 0, 50);

// High pT
lctStatDir->mkdir("High pT")->cd();

TH2D *AllLctStatVRingHighPt =
    new TH2D("AllLctStatVRingHighPt", "All LCTs Station vs Ring High pT", 4, 1,
             5, 4, 1, 5);
TH2D *PhiPairLctStatVRingHighPt =
    new TH2D("PhiPairLctStatVRingHighPt",
             "Phi Pair LCTs Station vs Ring High pT", 4, 1, 5, 4, 1, 5);
TH2D *ThetaPairLctStatVRingHighPt =
    new TH2D("ThetaPairLctStatVRingHighPt",
             "Theta Pair LCTs Station vs Ring High pT", 4, 1, 5, 4, 1, 5);
TH2D *QuadLctStatVRingHighPt =
    new TH2D("QuadLctStatVRingHighPt", "Quad LCTs Station vs Ring High pT", 4,
             1, 5, 4, 1, 5);

TH1D *AllLctStatHighPt =
    new TH1D("AllLctStatHighPt", "All LCTs Station High pT", 4, 1, 5);
TH1D *PhiPairLctStatHighPt =
    new TH1D("PhiPairLctStatHighPt", "Phi Pair LCTs Station High pT", 4, 1, 5);
TH1D *ThetaPairLctStatHighPt = new TH1D(
    "ThetaPairLctStatHighPt", "Theta Pair LCTs Station High pT", 4, 1, 5);
TH1D *QuadLctStatHighPt =
    new TH1D("QuadLctStatHighPt", "Quad LCTs Station High pT", 4, 1, 5);

TH1D *AllLctQualHighPt =
    new TH1D("AllLctQualHighPt", "All LCTs Quality High pT", 16, 0, 16);
TH1D *PhiPairLctQualHighPt = new TH1D(
    "PhiPairLctQualHighPt", "Phi Pair LCTs Quality High pT", 16, 0, 16);
TH1D *ThetaPairLctQualHighPt = new TH1D(
    "ThetaPairLctQualHighPt", "Theta Pair LCTs Quality High pT", 16, 0, 16);
TH1D *QuadLctQualHighPt =
    new TH1D("QuadLctQualHighPt", "Quad LCTs Quality High pT", 16, 0, 16);

TH1D *AllLctThetaHighPt =
    new TH1D("AllLctThetaHighPt", "All LCTs Theta High pT", 100, 0, 50);
TH1D *PhiPairLctThetaHighPt = new TH1D(
    "PhiPairLctThetaHighPt", "Phi Pair LCTs Theta High pT", 100, 0, 50);
TH1D *ThetaPairLctThetaHighPt = new TH1D(
    "ThetaPairLctThetaHighPt", "Theta Pair LCTs Theta High pT", 100, 0, 50);
TH1D *QuadLctThetaHighPt =
    new TH1D("QuadLctThetaHighPt", "Quad LCTs Theta High pT", 100, 0, 50);

// Pair Deltas
f->mkdir("Pair Deltas")->cd();

TH1D *PhiPairDTheta =
    new TH1D("PhiPairDTheta", "Phi Pair LCT dTheta", 40, 0, 20);
TH1D *PhiPairDWire =
    new TH1D("PhiPairDWire", "Phi Pair LCT dWire", 120, 0, 120);
TH1D *ThetaPairDPhi =
    new TH1D("ThetaPairDPhi", "Theta Pair LCT dPhi", 60, 0, 30);
TH1D *ThetaPairDStrip =
    new TH1D("ThetaPairDStrip", "Theta Pair LCT dStrip", 200, 0, 200);

TH1D *PhiPairTrackLoPtDTheta =
    new TH1D("PhiPairTrackLoPtDTheta", "Phi Pair LCT in Track (pT < 22) dTheta",
             40, 0, 20);
TH1D *PhiPairTrackLoPtDWire =
    new TH1D("PhiPairTrackLoPtDWire", "Phi Pair LCT in Track (pT < 22) dWire",
             120, 0, 120);
TH1D *ThetaPairTrackLoPtDPhi =
    new TH1D("ThetaPairTrackLoPtDPhi", "Theta Pair LCT in Track (pT < 22) dPhi",
             60, 0, 30);
TH1D *ThetaPairTrackLoPtDStrip =
    new TH1D("ThetaPairTrackLoPtDStrip",
             "Theta Pair LCT in Track (pT < 22) dStrip", 200, 0, 200);

TH1D *PhiPairTrackHiPtDTheta =
    new TH1D("PhiPairTrackHiPtDTheta",
             "Phi Pair LCT in Track (pT >= 22) dTheta", 40, 0, 20);
TH1D *PhiPairTrackHiPtDWire =
    new TH1D("PhiPairTrackHiPtDWire", "Phi Pair LCT in Track (pT >= 22) dWire",
             120, 0, 120);
TH1D *ThetaPairTrackHiPtDPhi =
    new TH1D("ThetaPairTrackHiPtDPhi",
             "Theta Pair LCT in Track (pT >= 22) dPhi", 60, 0, 30);
TH1D *ThetaPairTrackHiPtDStrip =
    new TH1D("ThetaPairTrackHiPtDStrip",
             "Theta Pair LCT in Track (pT >= 22) dStrip", 200, 0, 200);

TH1D *PhiPairTrackAllPtDTheta = new TH1D(
    "PhiPairTrackAllPtDTheta", "Phi Pair LCT in Track dTheta", 40, 0, 20);
TH1D *PhiPairTrackAllPtDWire = new TH1D(
    "PhiPairTrackAllPtDWire", "Phi Pair LCT in Track dWire", 120, 0, 120);
TH1D *ThetaPairTrackAllPtDPhi = new TH1D(
    "ThetaPairTrackAllPtDPhi", "Theta Pair LCT in Track dPhi", 60, 0, 30);
TH1D *ThetaPairTrackAllPtDStrip = new TH1D(
    "ThetaPairTrackAllPtDStrip", "Theta Pair LCT in Track dStrip", 200, 0, 200);

THStack *DWireStack = new THStack("dWireStack", "PhiPairDWire");
// h1->SetFillColor(kRed);
DWireStack->Add(PhiPairDWire);
DWireStack->Add(PhiPairTrackAllPtDWire);
DWireStack->Add(PhiPairTrackLoPtDWire);
DWireStack->Add(PhiPairTrackHiPtDWire);
DWireStack->Draw("nostack");

// LCT Quality
auto qualDir = f -> mkdir("LCT Quality");
qualDir->cd();

qualDir->mkdir("All Lct")->cd();

TH1D *AllAlctQuality =
    new TH1D("AllAlctQual", "All Alcts (in LCT) Quality", 16, 0, 16);
TH1D *PhiPairAlctQuality = new TH1D(
    "PhiPairAlctQual", "Phi Pair Alcts (in Phi Pair LCT) Quality", 16, 0, 16);
TH1D *ThetaPairAlctQuality =
    new TH1D("ThetaPairAlctQual",
             "Theta Pair Alcts (in Theta Pair LCT) Quality", 16, 0, 16);
TH1D *QuadAlctQuality =
    new TH1D("QuadAlctQual", "Quad Alcts (in LCT) Quality", 16, 0, 16);

TH1D *AllClctQuality =
    new TH1D("AllClctQual", "All Clcts (in LCT) Quality", 16, 0, 16);
TH1D *PhiPairClctQuality = new TH1D(
    "PhiPairClctQual", "Phi Pair Clcts (in Phi Pair LCT) Quality", 16, 0, 16);
TH1D *ThetaPairClctQuality =
    new TH1D("ThetaPairClctQual",
             "Theta Pair Clcts (in Theta Pair LCT) Quality", 16, 0, 16);
TH1D *QuadClctQuality =
    new TH1D("QuadClctQual", "Quad Clcts (in LCT) Quality", 16, 0, 16);

TH1D *AllLctAlctQualityByStat[4];
TH1D *AllLctClctQualityByStat[4];
TH2D *AllLctAlctVsClctByStat[4];
TH2D *AllLctAlctQualVsStat =
    new TH2D("AllLctAlctQualVsStat", "All Lcts - Alct qual vs Station", 4, 1, 5,
             16, 0, 16);
TH2D *AllLctClctQualVsStat =
    new TH2D("AllLctClctQualVsStat", "All Lcts - Clct qual vs Station", 4, 1, 5,
             16, 0, 16);
for (int s = 1; s < 5; s++) {
  AllLctAlctQualityByStat[s - 1] =
      new TH1D(("AllLctAlctQualStat" + std::to_string(s)).c_str(),
               ("All Lcts - Alct Qual - Station " + std::to_string(s)).c_str(),
               16, 0, 16);
  AllLctClctQualityByStat[s - 1] =
      new TH1D(("AllLctClctQualStat" + std::to_string(s)).c_str(),
               ("All Lcts - Clct Qual - Station " + std::to_string(s)).c_str(),
               16, 0, 16);
  AllLctAlctVsClctByStat[s - 1] = new TH2D(
      ("AllLctAlctVsClctStat" + std::to_string(s)).c_str(),
      ("All Lcts - Alct qual vs Clct qual - Station " + std::to_string(s))
          .c_str(),
      16, 0, 16, 16, 0, 16);
}

qualDir->mkdir("Phi Pair Lct")->cd();

TH1D *PhiPairAlctQualityByStat[4];
TH1D *PhiPairClctQualityByStat[4];
TH2D *PhiPairAlctVsClctByStat[4];
TH2D *PhiPairQual1VsQual2ByStat[4];
TH2D *PhiPairAlct1VsAlct2ByStat[4];
TH2D *PhiPairClct1VsClct2ByStat[4];
TH2D *PhiPairAlctQualVsStat =
    new TH2D("PhiPairAlctQualVsStat", "Phi Pair Lcts - Alct qual vs Station", 4,
             1, 5, 16, 0, 16);
TH2D *PhiPairClctQualVsStat =
    new TH2D("PhiPairClctQualVsStat", "Phi Pair Lcts - Clct qual vs Station", 4,
             1, 5, 16, 0, 16);
for (int s = 1; s < 5; s++) {
  PhiPairAlctQualityByStat[s - 1] = new TH1D(
      ("PhiPairAlctQualStat" + std::to_string(s)).c_str(),
      ("Phi Pair Lcts - Alct Qual - Station " + std::to_string(s)).c_str(), 16,
      0, 16);
  PhiPairClctQualityByStat[s - 1] = new TH1D(
      ("PhiPairClctQualStat" + std::to_string(s)).c_str(),
      ("Phi Pair Lcts - Clct Qual - Station " + std::to_string(s)).c_str(), 16,
      0, 16);
  PhiPairAlctVsClctByStat[s - 1] = new TH2D(
      ("PhiPairAlctVsClctStat" + std::to_string(s)).c_str(),
      ("Phi Pair Lcts - Alct qual vs Clct qual - Station " + std::to_string(s))
          .c_str(),
      16, 0, 16, 16, 0, 16);
  PhiPairQual1VsQual2ByStat[s - 1] =
      new TH2D(("PhiPairQual1VsQual2Stat" + std::to_string(s)).c_str(),
               ("Phi Pair Lcts - Qual 1 qual vs Qual 2 qual - Station " +
                std::to_string(s))
                   .c_str(),
               16, 0, 16, 16, 0, 16);
  PhiPairAlct1VsAlct2ByStat[s - 1] =
      new TH2D(("PhiPairAlct1VsAlct2Stat" + std::to_string(s)).c_str(),
               ("Phi Pair Lcts - Alct 1 qual vs Alct 2 qual - Station " +
                std::to_string(s))
                   .c_str(),
               16, 0, 16, 16, 0, 16);
  PhiPairClct1VsClct2ByStat[s - 1] =
      new TH2D(("PhiPairClct1VsClct2Stat" + std::to_string(s)).c_str(),
               ("Phi Pair Lcts - Clct 1 qual vs Clct 2 qual - Station " +
                std::to_string(s))
                   .c_str(),
               16, 0, 16, 16, 0, 16);
}

qualDir->mkdir("Theta Pair Lct")->cd();

TH1D *ThetaPairAlctQualityByStat[4];
TH1D *ThetaPairClctQualityByStat[4];
TH2D *ThetaPairAlctVsClctByStat[4];
TH2D *ThetaPairQual1VsQual2ByStat[4];
TH2D *ThetaPairAlct1VsAlct2ByStat[4];
TH2D *ThetaPairClct1VsClct2ByStat[4];
TH2D *ThetaPairAlctQualVsStat =
    new TH2D("ThetaPairAlctQualVsStat",
             "Theta Pair Lcts - Alct qual vs Station", 4, 1, 5, 16, 0, 16);
TH2D *ThetaPairClctQualVsStat =
    new TH2D("ThetaPairClctQualVsStat",
             "Theta Pair Lcts - Clct qual vs Station", 4, 1, 5, 16, 0, 16);
for (int s = 1; s < 5; s++) {
  ThetaPairAlctQualityByStat[s - 1] = new TH1D(
      ("ThetaPairAlctQualStat" + std::to_string(s)).c_str(),
      ("Theta Pair Lcts - Alct Qual - Station " + std::to_string(s)).c_str(),
      16, 0, 16);
  ThetaPairClctQualityByStat[s - 1] = new TH1D(
      ("ThetaPairClctQualStat" + std::to_string(s)).c_str(),
      ("Theta Pair Lcts - Clct Qual - Station " + std::to_string(s)).c_str(),
      16, 0, 16);
  ThetaPairAlctVsClctByStat[s - 1] =
      new TH2D(("ThetaPairAlctVsClctStat" + std::to_string(s)).c_str(),
               ("Theta Pair Lcts - Alct qual vs Clct qual - Station " +
                std::to_string(s))
                   .c_str(),
               16, 0, 16, 16, 0, 16);
  ThetaPairQual1VsQual2ByStat[s - 1] =
      new TH2D(("ThetaPairQual1VsQual2Stat" + std::to_string(s)).c_str(),
               ("Theta Pair Lcts - Qual 1 qual vs Qual 2 qual - Station " +
                std::to_string(s))
                   .c_str(),
               16, 0, 16, 16, 0, 16);
  ThetaPairAlct1VsAlct2ByStat[s - 1] =
      new TH2D(("ThetaPairAlct1VsAlct2Stat" + std::to_string(s)).c_str(),
               ("Theta Pair Lcts - Alct 1 qual vs Alct 2 qual - Station " +
                std::to_string(s))
                   .c_str(),
               16, 0, 16, 16, 0, 16);
  ThetaPairClct1VsClct2ByStat[s - 1] =
      new TH2D(("ThetaPairClct1VsClct2Stat" + std::to_string(s)).c_str(),
               ("Theta Pair Lcts - Clct 1 qual vs Clct 2 qual - Station " +
                std::to_string(s))
                   .c_str(),
               16, 0, 16, 16, 0, 16);
}

qualDir->mkdir("Quad Lct")->cd();

TH1D *QuadAlctQualityByStat[4];
TH1D *QuadClctQualityByStat[4];
TH2D *QuadAlctVsClctByStat[4];
TH2D *QuadAlctQualVsStat =
    new TH2D("QuadAlctQualVsStat", "Quad Lcts - Alct qual vs Station", 4, 1, 5,
             16, 0, 16);
TH2D *QuadClctQualVsStat =
    new TH2D("QuadClctQualVsStat", "Quad Lcts - Clct qual vs Station", 4, 1, 5,
             16, 0, 16);
for (int s = 1; s < 5; s++) {
  QuadAlctQualityByStat[s - 1] =
      new TH1D(("QuadAlctQualStat" + std::to_string(s)).c_str(),
               ("Quad Lcts - Alct Qual - Station " + std::to_string(s)).c_str(),
               16, 0, 16);
  QuadClctQualityByStat[s - 1] =
      new TH1D(("QuadClctQualStat" + std::to_string(s)).c_str(),
               ("Quad Lcts - Clct Qual - Station " + std::to_string(s)).c_str(),
               16, 0, 16);
  QuadAlctVsClctByStat[s - 1] = new TH2D(
      ("QuadAlctVsClctStat" + std::to_string(s)).c_str(),
      ("Quad Lcts - Alct qual vs Clct qual - Station " + std::to_string(s))
          .c_str(),
      16, 0, 16, 16, 0, 16);
}
