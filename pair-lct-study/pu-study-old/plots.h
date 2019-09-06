using std::map;
using std::string;

string[] typeKeys = {"All", "SamePhiPair", "SameThetaPair", "Quad"};
string[] energyKeys = {"All", "Low-pT", "High-pT"};
map<string, map<string, map<string, TH1D*>>> hists1D;
map<string, map<string, map<string, TH2D*>>> hists2D;

for (auto type : typeKeys) {
    for (auto energy : energyKeys) {
        auto slot1D = hists1D[type][energy];
        auto slot2D = hists1D[type][energy];
        hists
    }
}



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