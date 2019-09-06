/////Alex's PU rate study/////
/////      :)START(:      ////
//////////////////////////////

// Event Breakdowns
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

// Pair Deltas
TH1D *PhiPairDTheta =
    new TH1D("PhiPairDTheta", "Phi Pair LCT dTheta", 100, 0, 100);
TH1D *PhiPairDWire =
    new TH1D("PhiPairDWire", "Phi Pair LCT dWire", 120, 0, 120);
TH1D *ThetaPairDPhi =
    new TH1D("ThetaPairDPhi", "Theta Pair LCT dPhi", 100, 0, 100);
TH1D *ThetaPairDStrip =
    new TH1D("ThetaPairDStrip", "Theta Pair LCT dStrip", 200, 0, 200);

TH1D *PhiPairTrackLoPtDTheta =
    new TH1D("PhiPairTrackLoPtDTheta", "Phi Pair LCT in Track (pT < 22) dTheta", 100, 0, 100);
TH1D *PhiPairTrackLoPtDWire =
    new TH1D("PhiPairTrackLoPtDWire", "Phi Pair LCT in Track (pT < 22) dWire", 120, 0, 120);
TH1D *ThetaPairTrackLoPtDPhi =
    new TH1D("ThetaPairTrackLoPtDPhi", "Theta Pair LCT in Track (pT < 22) dPhi", 100, 0, 100);
TH1D *ThetaPairTrackLoPtDStrip =
    new TH1D("ThetaPairTrackLoPtDStrip", "Theta Pair LCT in Track (pT < 22) dStrip", 200, 0, 200);

TH1D *PhiPairTrackHiPtDTheta =
    new TH1D("PhiPairTrackHiPtDTheta", "Phi Pair LCT in Track (pT >= 22) dTheta", 100, 0, 100);
TH1D *PhiPairTrackHiPtDWire =
    new TH1D("PhiPairTrackHiPtDWire", "Phi Pair LCT in Track (pT >= 22) dWire", 120, 0, 120);
TH1D *ThetaPairTrackHiPtDPhi =
    new TH1D("ThetaPairTrackHiPtDPhi", "Theta Pair LCT in Track (pT >= 22) dPhi", 100, 0, 100);
TH1D *ThetaPairTrackHiPtDStrip =
    new TH1D("ThetaPairTrackHiPtDStrip", "Theta Pair LCT in Track (pT >= 22) dStrip", 200, 0, 200);

TH1D *PhiPairTrackAllPtDTheta =
    new TH1D("PhiPairTrackAllPtDTheta", "Phi Pair LCT in Track dTheta", 100, 0, 100);
TH1D *PhiPairTrackAllPtDWire =
    new TH1D("PhiPairTrackAllPtDWire", "Phi Pair LCT in Track dWire", 120, 0, 120);
TH1D *ThetaPairTrackAllPtDPhi =
    new TH1D("ThetaPairTrackAllPtDPhi", "Theta Pair LCT in Track dPhi", 100, 0, 100);
TH1D *ThetaPairTrackAllPtDStrip =
    new TH1D("ThetaPairTrackAllPtDStrip", "Theta Pair LCT in Track dStrip", 200, 0, 200);





