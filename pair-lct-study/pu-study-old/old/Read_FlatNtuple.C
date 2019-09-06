
////////////////////////////////////////////////////////////
///        Simple macro to read EMTF Flat NTuples        ///
///              Andrew Brinkerhoff 29.09.17             ///
///                                                      ///
///   TChain can be used to read multiple files.         ///
///   Format: interface/FlatNtupleBranches.h             ///
////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TStyle.h"
#include "TMath.h"

#include "Read_FlatNtuple.h" // List of input branches and functions to return values

const int MAX_FILES = 10;   // Max number of files to process
const int MAX_EVT   = 5000000;   // Max number of events to process
const int PRT_EVT   =  1000;   // Print every N events
const bool verbose  = false; // Print information about the event and RECO and L1T muons

///////////////////////////////////////////////
/////  Prepare the plots and their labels  ////
///////////////////////////////////////////////
TFile *f = new TFile("PURateTrainRPC.root","RECREATE");
//Chad's LCT plots//
//    :)START(:   //
////////////////////
TH1D *LCTBX = new TH1D("LCT BX","LCT BX",7,-3.5,3.5);
TH2D *ChamberBXPos42 = new TH2D("Chamber vs BX ME+4/2","Chamber vs BX ME+4/2",36,1,37,7,-3.5,3.5);
TH2D *ChamberBXPos32 = new TH2D("Chamber vs BX ME+3/2","Chamber vs BX ME+3/2",36,1,37,7,-3.5,3.5);
TH2D *ChamberBXNeg42 = new TH2D("Chamber vs BX ME-4/2","Chamber vs BX ME-4/2",36,1,37,7,-3.5,3.5);
TH2D *ChamberBXNeg32 = new TH2D("Chamber vs BX ME-3/2","Chamber vs BX ME-3/2",36,1,37,7,-3.5,3.5);
TH2D *EtaBXPos42 = new TH2D("BX vs Eta ME+4/2","BX vs Eta ME+4/2",7,-3.5,3.5,40,1.1,1.8);
TH2D *EtaBXPos32 = new TH2D("BX vs Eta ME+3/2","BX vs Eta ME+3/2",7,-3.5,3.5,40,1.1,1.8);
TH2D *EtaBXNeg42 = new TH2D("BX vs Eta ME-4/2","BX vs Eta ME-4/2",7,-3.5,3.5,40,-1.8,-1.1);
TH2D *EtaBXNeg32 = new TH2D("BX vs Eta ME-3/2","BX vs Eta ME-3/2",7,-3.5,3.5,40,-1.8,-1.1);
TH1D *Pos42 = new TH1D("+4/2 occupancy","chamber occupancy",36,1,37);
TH1D *Pos41 = new TH1D("+4/1 occupancy","chamber occupancy",18,1,19);
TH1D *Pos32 = new TH1D("+3/2 occupancy","chamber occupancy",36,1,37);
TH1D *Pos31 = new TH1D("+3/1 occupancy","chamber occupancy",18,1,19);
TH1D *Pos22 = new TH1D("+2/2 occupancy","chamber occupancy",36,1,37);
TH1D *Neg42 = new TH1D("-4/2 occupancy","chamber occupancy",36,1,37);
TH1D *Neg41 = new TH1D("-4/1 occupancy","chamber occupancy",18,1,19);
TH1D *Neg32 = new TH1D("-3/2 occupancy","chamber occupancy",36,1,37);
TH1D *Neg31 = new TH1D("-3/1 occupancy","chamber occupancy",18,1,19);
TH1D *Neg22 = new TH1D("-2/2 occupancy","chamber occupancy",36,1,37);
//Chad's Track plots//
//    :)START(:     //
//////////////////////
TH1D *TrackLCTBX = new TH1D("Track LCT BX","Track LCT BX",7,-3.5,3.5);
TH1D *TrackLCTBX42 = new TH1D("Track LCT BX 4/2","Track LCT BX",7,-3.5,3.5);
TH1D *TrackLCTBX42HiPt = new TH1D("Track LCT BX 4/2 with hi Pt","Track LCT BX with High Pt",7,-3.5,3.5);
TH1D *TrackBX = new TH1D("Track BX", "Track BX",7,-3.5,3.5);
TH1D *TrackNo42BX = new TH1D("Track no 42 BX", "Track with no 4/2 BX",7,-3.5,3.5);
TH1D *TrackBXHiPt = new TH1D("Track BX with Hi Pt", "Track BX with High Pt",7,-3.5,3.5);
TH1D *TrackBXNo42HiPt = new TH1D("Track no 42 BX with Hi Pt", "Track BX with no 4/2 and High Pt",7,-3.5,3.5);
TH1D *HiPtPhi = new TH1D("High Pt Phi","High Pt Phi",100,-3.15,3.15);
TH2D *BXvsPhi = new TH2D("BX vs Phi","BX vs Phi",7,-3.5,3.5,100,-3.15,3.15);
TH2D *EtavsPhi = new TH2D("Eta vs Phi","Eta vs Phi",50,-2.5,2.5,100,-3.15,3.15);
TH2D *HiPt_EtavsPhi = new TH2D("Hi Pt Eta vs Phi","Hi Pt Eta vs Phi",50,-2.5,2.5,100,-3.15,3.15);
TH2D *HierPt_EtavsPhi = new TH2D("Hier Pt Eta vs Phi","Hier Pt Eta vs Phi",50,-2.5,2.5,100,-3.15,3.15);
TH2D *TrackPtvsPhiPos42 = new TH2D("Pt vs Phi positive endcap","Pt vs Phi",250,1,251,100,-3.15,3.15);
TH2D *TrackPtvsPhiNeg42 = new TH2D("Pt vs Phi negative endcap","Pt vs Phi",250,1,251,100,-3.15,3.15);
//Chad's Timing plots////
//      :)START(:      //
/////////////////////////
//TDirectory *Timing = f->mkdir("Time_Plots");
//Time_Plots->cd();    // make the "Time_Plots" directory the current directory
TH1D *TrackdBX = new TH1D("Track dBX","Track dBX",7,0,6);
TH1D *TrackdPhi = new TH1D("Track dPhi","Track dPhi",1001,0,1000);
TH1D *TrackdTheta = new TH1D("Track dTheta","Track dTheta",51,0,50);
TH1D *TrackPT = new TH1D("Track PT", "Track PT",100,1,101);
TH1D *TrackPTdBX = new TH1D("Track PT w/ dBX", "Track PT w/ dBX",100,1,101);
TH1D *TrackPTdPhi = new TH1D("Track PT w/ dPhi", "Track PT w/ dPhi",100,1,101);
TH1D *TrackPTdTheta = new TH1D("Track PT w/ dTheta", "Track PT w/ dTheta",100,1,101);

TH2D *ME1vME2BX = new TH2D("ME1 vs ME2 BX","ME1 vs ME2 BX",7,-3.5,3.5,7,-3.5,3.5);
TH2D *ME1vME3BX = new TH2D("ME1 vs ME3 BX","ME1 vs ME3 BX",7,-3.5,3.5,7,-3.5,3.5);
TH2D *ME1vME4BX = new TH2D("ME1 vs ME4 BX","ME1 vs ME4 BX",7,-3.5,3.5,7,-3.5,3.5);
TH2D *ME2vME3BX = new TH2D("ME2 vs ME3 BX","ME2 vs ME3 BX",7,-3.5,3.5,7,-3.5,3.5);
TH2D *ME2vME4BX = new TH2D("ME2 vs ME4 BX","ME2 vs ME4 BX",7,-3.5,3.5,7,-3.5,3.5);
TH2D *ME3vME4BX = new TH2D("ME3 vs ME4 BX","ME3 vs ME4 BX",7,-3.5,3.5,7,-3.5,3.5);
TH2D *ME1vME2Phi = new TH2D("ME1 vs ME2 Phi","ME1 vs ME2 Phi",1500,500,2001,1500,500,2001);
TH2D *ME1vME3Phi = new TH2D("ME1 vs ME3 Phi","ME1 vs ME3 Phi",1500,500,2001,1500,500,2001);
TH2D *ME1vME4Phi = new TH2D("ME1 vs ME4 Phi","ME1 vs ME4 Phi",1500,500,2001,1500,500,2001);
TH2D *ME2vME3Phi = new TH2D("ME2 vs ME3 Phi","ME2 vs ME3 Phi",1500,500,2001,1500,500,2001);
TH2D *ME2vME4Phi = new TH2D("ME2 vs ME4 Phi","ME2 vs ME4 Phi",1500,500,2001,1500,500,2001);
TH2D *ME3vME4Phi = new TH2D("ME3 vs ME4 Phi","ME3 vs ME4 Phi",1500,500,2001,1500,500,2001);
TH2D *ME1vME2Theta = new TH2D("ME1 vs ME2 Theta","ME1 vs ME2 Theta",20,1,101,20,1,101);
TH2D *ME1vME3Theta = new TH2D("ME1 vs ME3 Theta","ME1 vs ME3 Theta",20,1,101,20,1,101);
TH2D *ME1vME4Theta = new TH2D("ME1 vs ME4 Theta","ME1 vs ME4 Theta",20,1,101,20,1,101);
TH2D *ME2vME3Theta = new TH2D("ME2 vs ME3 Theta","ME2 vs ME3 Theta",20,1,101,20,1,101);
TH2D *ME2vME4Theta = new TH2D("ME2 vs ME4 Theta","ME2 vs ME4 Theta",20,1,101,20,1,101);
TH2D *ME3vME4Theta = new TH2D("ME3 vs ME4 Theta","ME3 vs ME4 Theta",20,1,101,20,1,101);
int ME1BX;
int ME2BX;
int ME3BX;
int ME4BX;
int ME1Phi;
int ME2Phi;
int ME3Phi;
int ME4Phi;
int ME1Theta;
int ME2Theta;
int ME3Theta;
int ME4Theta;
int dThetaFlag3;
int dThetaFlag4;
//Chad's dTheta plots////
////      :)START(:      //
///////////////////////////
TH2D *Station3dThetaCheck = new TH2D("3 Station dTheta Check","3 Station dTheta Check",18,-1,16,10,-1,8);
TH2D *Station4dThetaCheck = new TH2D("4 Station dTheta Check","4 Station dTheta Check",18,-1,16,10,-1,8);
int station1ThetaInt;
int station2ThetaInt;
int station3ThetaInt;
int station4ThetaInt;
int dTheta1x;
int dThetaxy;
int dTheta41x;
int dTheta4xy;

/////Chad's PU rate study/////
/////      :)START(:      ////
//////////////////////////////
TH1D *SingleMuPtLo = new TH1D("SingleMuonPtLo","SingleMuonPtLo",22,1,23);
TH1D *SingleMuPtLodBX1 = new TH1D("SingleMuonPtLodBX1","SingleMuonPtLodBX1",22,1,23);
TH1D *SingleMuPtLodBX1ME11 = new TH1D("SingleMuonPtLodBX1ME11","SingleMuonPtLodBX1ME11",22,1,23);

TH1D *SingleMuPtHi = new TH1D("SingleMuonPtHi","SingleMuonPtHi",78,22,101);
TH1D *SingleMuPtHidBX1 = new TH1D("SingleMuonPtHidBX1","SingleMuonPtHidBX1",78,22,101);
TH1D *SingleMuPtHidBX1ME11 = new TH1D("SingleMuonPtHidBX1ME11","SingleMuonPtHidBX1ME11",78,22,101);

TH1D *SingleMuPtAll = new TH1D("SingleMuonPtAll","SingleMuonPtAll",100,1,101);
TH1D *SingleMuPtAlldBX1 = new TH1D("SingleMuonPtAlldBX1","SingleMuonPtAlldBX1",100,1,101);
TH1D *SingleMuPtAlldBX1ME11 = new TH1D("SingleMuonPtAlldBX1ME11","SingleMuonPtAlldBX1ME11",100,1,101);

TH1D *SingleMuPtLoPattern = new TH1D("SingleMuonPtLoPattern","SingleMuonPtLoPattern",22,1,23);
TH1D *SingleMuPtHiPattern = new TH1D("SingleMuonPtHiPattern","SingleMuonPtHiPattern",78,22,101);
TH1D *SingleMuPtAllPattern = new TH1D("SingleMuonPtAllPattern","SingleMuonPtAllPattern",100,1,101);

/////Alex's PU rate study/////
/////      :)START(:      ////
//////////////////////////////
TH1D *SingleMuPtLoNoPair = new TH1D("SingleMuonPtLoNoPair","SingleMuonPtLoNoPair",22,1,23);
TH1D *SingleMuPtHiNoPair = new TH1D("SingleMuonPtHiNoPair","SingleMuonPtHiNoPair",78,22,101);
TH1D *SingleMuPtAllNoPair = new TH1D("SingleMuonPtAllNoPair","SingleMuonPtAllNoPair",100,1,101);

//////////////////////////////////////////////
/////  Main function to read the NTuples  ////
//////////////////////////////////////////////

void Read_FlatNtuple() {
   //BookHistos();
 
  gROOT->SetStyle("Plain");
 // Initialize empty file to access each file in the list
  TFile *file_tmp(0);

  // List of input files
  std::vector<TString> in_file_names;
  TString store = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/";
  TString in_dir = "ntuples/HADD/";
  TString file_name;

  //file_name.Form("%s/%s/NTuple_ZeroBiasIsolatedBunches_FlatNtuple_Skim_Run_302674_2017_09_30.root", store.Data(), in_dir.Data());
  //std::cout << "Adding file " << file_name.Data() << std::endl;
  //in_file_names.push_back(file_name.Data());
  
   for (int i = 1; i < MAX_FILES+1; i++) {
     //file_name.Form("%s/%s/NTuple_ZeroBias8b4e_FlatNtuple_Run_302674_2017_09_30.root", store.Data(), in_dir.Data());
     //file_name.Form("%s/%s/NTuple_ZeroBias8b4e_FlatNtuple_Skim_Run_302674_2017_09_30.root", store.Data(), in_dir.Data());
     //file_name.Form("%s/%s/NTuple_ZeroBiasIsolatedBunches_FlatNtuple_Run_302674_2017_09_30.root", store.Data(), in_dir.Data());
     //file_name.Form("%s/%s/NTuple_ZeroBiasIsolatedBunches_FlatNtuple_Skim_Run_302674_2017_09_30.root", store.Data(), in_dir.Data());
     file_name.Form("%s/%s/NTuple_ZeroBiasNominalTrains_FlatNtuple_Run_302674_2017_09_30.root", store.Data(), in_dir.Data());
     //file_name.Form("%s/%s/NTuple_ZeroBiasNominalTrains_FlatNtuple_Skim_Run_302674_2017_09_30.root", store.Data(), in_dir.Data());
     std::cout << "Adding file " << file_name.Data() << std::endl;
     in_file_names.push_back(file_name.Data());
   }

  // Open all input files
  for (int i = 0; i < in_file_names.size(); i++) {
    if ( !gSystem->AccessPathName(in_file_names.at(i)) )
      file_tmp = TFile::Open( in_file_names.at(i) ); // Check if file exists
    if (!file_tmp) {
      std::cout << "ERROR: could not open data file " << in_file_names.at(i) << std::endl;
      return;
    }
  }

  // Add tree from the input files to the TChain
  TChain *in_chain= new TChain("ntuple/tree");
  for (int i = 0; i < in_file_names.size(); i++) {
    in_chain->Add( in_file_names.at(i) );
  }

  InitializeMaps();
  SetBranchAddresses(in_chain);

  // if (verbose) in_chain->GetListOfBranches()->Print();
  
  std::cout << "\n******* About to loop over the events *******" << std::endl;
  int nEvents = in_chain->GetEntries();
  for (int iEvt = 0; iEvt < nEvents; iEvt++) {
    if (iEvt > MAX_EVT) break;
    if ( (iEvt % PRT_EVT) == 0 ) {
      std::cout << "\n*************************************" << std::endl;
      std::cout << "Looking at event " << iEvt << " out of " << nEvents << std::endl;
      std::cout << "*************************************" << std::endl;
    }

    in_chain->GetEntry(iEvt);
    
    // From Read_FlatNtuple.h, use 'I("branch_name")' to get an integer branch value, 'F("branch_name") to get a float
    if (verbose) std::cout << "\nRun = " << I("evt_run") << ", LS = " << I("evt_LS") << ", event = " << I("evt_event") << std::endl;

    // Print info for emulated EMTF hits
    if (verbose) std::cout << "\n" << I("nHits") << " emulated EMTF hits in the event" << std::endl;
    for (int i = 0; i < I("nHits"); i++) {

      if        (I("hit_isCSC", i) == 1) {
    	if (verbose) std::cout << " * CSC LCT with BX = " << I("hit_BX", i) << ", endcap = " << I("hit_endcap", i)
    			       << ", station = " << I("hit_station", i) << ", ring = " << I("hit_ring", i) << std::endl; 

        /////////////Chad filling LCT plots//////////
        //                :)START(:                //
        /////////////////////////////////////////////
        LCTBX->Fill(I("hit_BX", i),1);
        LCTBX->GetXaxis()->SetTitle("BX");
        LCTBX->GetYaxis()->SetTitle("Entries");
        if (I("hit_endcap",i)>0 && I("hit_neighbor",i)==0){
           if(I("hit_station",i)==4 && I("hit_ring",i)==2) Pos42->Fill(I("hit_chamber",i),1);
           if(I("hit_station",i)==4 && I("hit_ring",i)==1) Pos41->Fill(I("hit_chamber",i),1);
           if(I("hit_station",i)==3 && I("hit_ring",i)==2) Pos32->Fill(I("hit_chamber",i),1);
           if(I("hit_station",i)==3 && I("hit_ring",i)==1) Pos31->Fill(I("hit_chamber",i),1);
           if(I("hit_station",i)==2 && I("hit_ring",i)==2) Pos22->Fill(I("hit_chamber",i),1);
           if(I("hit_station",i)==4 && I("hit_ring",i)==2) ChamberBXPos42->Fill(I("hit_chamber",i),I("hit_BX",i),1);
           if(I("hit_station",i)==3 && I("hit_ring",i)==2) ChamberBXPos32->Fill(I("hit_chamber",i),I("hit_BX",i),1);
           if(I("hit_station",i)==4 && I("hit_ring",i)==2) EtaBXPos42->Fill(I("hit_BX",i),F("hit_eta",i),1);
           if(I("hit_station",i)==3 && I("hit_ring",i)==2) EtaBXPos32->Fill(I("hit_BX",i),F("hit_eta",i),1);
        }
        if (I("hit_endcap",i)<0 && I("hit_neighbor",i)==0){
           if(I("hit_station",i)==4 && I("hit_ring",i)==2) Neg42->Fill(I("hit_chamber",i),1);
           if(I("hit_station",i)==4 && I("hit_ring",i)==1) Neg41->Fill(I("hit_chamber",i),1);
           if(I("hit_station",i)==3 && I("hit_ring",i)==2) Neg32->Fill(I("hit_chamber",i),1);
           if(I("hit_station",i)==3 && I("hit_ring",i)==1) Neg31->Fill(I("hit_chamber",i),1);
           if(I("hit_station",i)==2 && I("hit_ring",i)==2) Neg22->Fill(I("hit_chamber",i),1);
           if(I("hit_station",i)==4 && I("hit_ring",i)==2) ChamberBXNeg42->Fill(I("hit_chamber",i),I("hit_BX",i),1);
           if(I("hit_station",i)==3 && I("hit_ring",i)==2) ChamberBXNeg32->Fill(I("hit_chamber",i),I("hit_BX",i),1);
           if(I("hit_station",i)==4 && I("hit_ring",i)==2) EtaBXNeg42->Fill(I("hit_BX",i),F("hit_eta",i),1);
           if(I("hit_station",i)==3 && I("hit_ring",i)==2) EtaBXNeg32->Fill(I("hit_BX",i),F("hit_eta",i),1);
        }
        //Set up the X bin labels
        Pos42->GetXaxis()->SetTitle("Chamber");
        Pos41->GetXaxis()->SetTitle("Chamber");
        Pos32->GetXaxis()->SetTitle("Chamber");
        Pos31->GetXaxis()->SetTitle("Chamber");
        Pos22->GetXaxis()->SetTitle("Chamber");
        Neg42->GetXaxis()->SetTitle("Chamber");
        Neg41->GetXaxis()->SetTitle("Chamber");
        Neg32->GetXaxis()->SetTitle("Chamber");
        Neg31->GetXaxis()->SetTitle("Chamber");
        Neg22->GetXaxis()->SetTitle("Chamber");
        Pos42->GetYaxis()->SetTitle("Entries");
        Pos41->GetYaxis()->SetTitle("Entries");
        Pos32->GetYaxis()->SetTitle("Entries");
        Pos31->GetYaxis()->SetTitle("Entries");
        Pos22->GetYaxis()->SetTitle("Entries");
        Neg42->GetYaxis()->SetTitle("Entries");
        Neg41->GetYaxis()->SetTitle("Entries");
        Neg32->GetYaxis()->SetTitle("Entries");
        Neg31->GetYaxis()->SetTitle("Entries");
        Neg22->GetYaxis()->SetTitle("Entries");
        ChamberBXPos42->GetXaxis()->SetTitle("Chamber");
        ChamberBXPos32->GetXaxis()->SetTitle("Chamber");
        ChamberBXNeg42->GetXaxis()->SetTitle("Chamber");
        ChamberBXNeg32->GetXaxis()->SetTitle("Chamber");
        EtaBXPos42->GetYaxis()->SetTitle("#eta");
        EtaBXPos32->GetYaxis()->SetTitle("#eta");
        EtaBXNeg42->GetYaxis()->SetTitle("#eta");
        EtaBXNeg32->GetYaxis()->SetTitle("#eta");
        ChamberBXPos42->GetYaxis()->SetTitle("BX");
        ChamberBXPos32->GetYaxis()->SetTitle("BX");
        ChamberBXNeg42->GetYaxis()->SetTitle("BX");
        ChamberBXNeg32->GetYaxis()->SetTitle("BX");
        EtaBXPos42->GetXaxis()->SetTitle("BX");
        EtaBXPos32->GetXaxis()->SetTitle("BX");
        EtaBXNeg42->GetXaxis()->SetTitle("BX");
        EtaBXNeg32->GetXaxis()->SetTitle("BX");
        /////////////Chad filling LCT plots//////////
        //                :(END):                  //
        /////////////////////////////////////////////
      } else if (I("hit_isRPC", i) == 1) {
    	if (verbose) std::cout << " * RPC hit with BX = " << I("hit_BX", i) << ", endcap = " << I("hit_endcap", i)
    			       << ", station = " << I("hit_station", i) << ", ring = " << I("hit_ring", i) << std::endl; 
      }
    }
    
    // Print info for emulated EMTF tracks
    if (verbose) std::cout << "\n" << I("nTracks") << " emulated EMTF tracks in the event" << std::endl;
    if ( I("nTracks") != VF("trk_pt")->size() )  // From Read_FlatNtuple.h, can use 'VF("branch_name")' to get a vector of floats
      throw std::out_of_range("nTracks does not equal size of vector of tracks");
    for (int i = 0; i < I("nTracks"); i++) {
      // // From Read_FlatNtuple.h, can use '("branch_name", i)' to get the ith element of a vector of ints
      if (verbose) std::cout << " * Mode " << I("trk_mode", i) << " track with BX = " << I("trk_BX", i) 
    			     << ", pT = " << F("trk_pt", i) << ", eta = " << F("trk_eta", i) << ", phi = " << F("trk_phi", i)
			     << ", maximum dPhi_int among hits = " << I("trk_dPhi_int", i) << std::endl;
    }
    
    // Print info for unpacked EMTF tracks
    if (verbose) std::cout << "\n" << I("nUnpTracks") << " unpacked EMTF tracks in the event" << std::endl;
    for (int i = 0; i < I("nUnpTracks"); i++) {
      
      if (verbose) std::cout << " * Mode " << I("unp_trk_mode", i) << " track with BX = " << I("unp_trk_BX", i) 
    			     << ", pT = " << F("unp_trk_pt", i) << ", eta = " << F("unp_trk_eta", i) << ", phi = " << F("unp_trk_phi", i)
			     << ", maximum dPhi_int among hits = " << I("unp_trk_dPhi_int", i) << std::endl;

      /////////////Chad fills Time plots///////////
      //                :)START(:                //
      /////////////////////////////////////////////
      if(F("unp_trk_pt", i)>25){
         TrackdBX->Fill(I("unp_trk_dBX",i),1);
         TrackdPhi->Fill(I("unp_trk_dPhi_int",i),1);
         TrackdTheta->Fill(I("unp_trk_dTheta_int",i),1);       
      }
      TrackPT->Fill(F("unp_trk_pt", i),1);
      if(I("unp_trk_dBX",i)>1) TrackPTdBX->Fill(F("unp_trk_pt", i),1);     
      if(I("unp_trk_dPhi_int",i)>300) TrackPTdPhi->Fill(F("unp_trk_pt", i),1);
      if(I("unp_trk_dTheta_int",i)>10) TrackPTdTheta->Fill(F("unp_trk_pt", i),1);

      ME1BX = 100;
      ME2BX = 100;
      ME3BX = 100;
      ME4BX = 100;

        station1ThetaInt=-900;
        station2ThetaInt=-900;
        station3ThetaInt=-900;
        station4ThetaInt=-900;
        dTheta1x = 0;
        dThetaxy = 0;
        dTheta41x = 0;
        dTheta4xy = 0;
        dThetaFlag3=-1;
        dThetaFlag4=-1;


      for (int j = 0; j < I("unp_trk_nHits", i); j++) {
	if (I("unp_trk_found_hits", i) != 1) continue;  // For a very small fraction of unpacked tracks, can't find all hits (mostly BX != 0)
	int iHit = I("unp_trk_iHit", i, j);  // Access the index of each hit in the track
	if        (I("hit_isCSC", iHit) == 1) {
        /////////////Chad fills Track plots//////////
        //                :)START(:                //
        /////////////////////////////////////////////
		TrackLCTBX->Fill(I("hit_BX",iHit),1);
                	TrackLCTBX->GetXaxis()->SetTitle("BX");
                	TrackLCTBX->GetYaxis()->SetTitle("Entries");
                if (I("hit_station",iHit)==4 && I("hit_ring",iHit)==2)TrackLCTBX42->Fill(I("hit_BX",iHit),1);
                	TrackLCTBX42->GetXaxis()->SetTitle("BX");
                	TrackLCTBX42->GetYaxis()->SetTitle("Entries");
                if (I("hit_station",iHit)==4 && I("hit_ring",iHit)==2 && F("unp_trk_pt", i)>25)TrackLCTBX42HiPt->Fill(I("hit_BX",iHit),1);
                        TrackLCTBX42HiPt->GetXaxis()->SetTitle("BX");
                        TrackLCTBX42HiPt->GetYaxis()->SetTitle("Entries");
                TrackBX->Fill(I("unp_trk_BX",i),1);
                        TrackBX->GetXaxis()->SetTitle("Track BX");
                        TrackBX->GetYaxis()->SetTitle("Entries");
                if (I("hit_station",iHit)<4 || I("hit_ring",iHit)==1)TrackNo42BX->Fill(I("unp_trk_BX",i),1);
                        TrackNo42BX->GetXaxis()->SetTitle("Track BX");
                        TrackNo42BX->GetYaxis()->SetTitle("Entries");
                if (F("unp_trk_pt", i)>25)TrackBXHiPt->Fill(I("unp_trk_BX",i),1);
                        TrackBXHiPt->GetXaxis()->SetTitle("Track BX");
                        TrackBXHiPt->GetYaxis()->SetTitle("Entries");
                if (I("hit_station",iHit)<4 || I("hit_ring",iHit)==1){ 
			if (F("unp_trk_pt", i)>25){
				TrackBXNo42HiPt->Fill(I("unp_trk_BX",i),1);
			}
		}
                        TrackBXNo42HiPt->GetXaxis()->SetTitle("Track BX");
                        TrackBXNo42HiPt->GetYaxis()->SetTitle("Entries");
                BXvsPhi->Fill(I("hit_BX",iHit),F("unp_trk_phi",i)*3.14159/180,1);
                	BXvsPhi->GetXaxis()->SetTitle("LCT BX");
                	BXvsPhi->GetYaxis()->SetTitle("Track #phi");
	        if (F("unp_trk_pt", i)>25 && I("unp_trk_nNeighbor", i)<3) HiPtPhi->Fill(F("unp_trk_phi",i)*3.14159/180,1);
                	HiPtPhi->GetXaxis()->SetTitle("Track #phi");
                	HiPtPhi->GetYaxis()->SetTitle("Entries");
                if (I("hit_station",iHit)==4 && I("hit_ring",iHit)==2) EtavsPhi->Fill(F("unp_trk_eta", i),F("unp_trk_phi",i)*3.14159/180,1);
                	EtavsPhi->GetXaxis()->SetTitle("Track #eta");
                	EtavsPhi->GetYaxis()->SetTitle("Track #phi");
                if (I("hit_station",iHit)==4 && I("hit_ring",iHit)==2 && F("unp_trk_pt", i)>25) HiPt_EtavsPhi->Fill(F("unp_trk_eta", i),F("unp_trk_phi",i)*3.14159/180,1);
                	HiPt_EtavsPhi->GetXaxis()->SetTitle("Track #eta");
                	HiPt_EtavsPhi->GetYaxis()->SetTitle("Track #phi");
                if (I("hit_station",iHit)==4 && I("hit_ring",iHit)==2 && F("unp_trk_pt", i)>50) HierPt_EtavsPhi->Fill(F("unp_trk_eta", i),F("unp_trk_phi",i)*3.14159/180,1);
                	HierPt_EtavsPhi->GetXaxis()->SetTitle("Track #eta");
                	HierPt_EtavsPhi->GetYaxis()->SetTitle("Track #phi");
                if (I("hit_station",iHit)==4 && I("hit_ring",iHit)==2 && I("hit_endcap",iHit)>0) TrackPtvsPhiPos42->Fill(F("unp_trk_pt", i),F("unp_trk_phi",i)*3.14159/180,1);
                        TrackPtvsPhiPos42->GetXaxis()->SetTitle("Track Pt");
                        TrackPtvsPhiPos42->GetYaxis()->SetTitle("Track #phi"); 
                if (I("hit_station",iHit)==4 && I("hit_ring",iHit)==2 && I("hit_endcap",iHit)<0) TrackPtvsPhiNeg42->Fill(F("unp_trk_pt", i),F("unp_trk_phi",i)*3.14159/180,1);
                        TrackPtvsPhiNeg42->GetXaxis()->SetTitle("Track Pt");
                        TrackPtvsPhiNeg42->GetYaxis()->SetTitle("Track #phi");
         /////////////Chad fills Track plots//////////
         //                 :(END):                 //
         /////////////////////////////////////////////
	

        ////////////Chad fills dTheta plots//////////
        //                :)START(:                //
        /////////////////////////////////////////////
	if (F("unp_trk_pt", i)>22 && fabs(I("unp_trk_BX", i))<2 && I("unp_trk_nRPC",i)==0 && I("unp_trk_found_hits",i)==1 && I("unp_trk_nHits",i)==3){
		if (I("hit_station",iHit)==1) station1ThetaInt=I("hit_theta_int",iHit);
		if (I("hit_station",iHit)==2) station2ThetaInt=I("hit_theta_int",iHit);
		if (I("hit_station",iHit)==3) station3ThetaInt=I("hit_theta_int",iHit);
		if (I("hit_station",iHit)==4) station4ThetaInt=I("hit_theta_int",iHit);
		dThetaFlag3=1;
		if (station1ThetaInt!=-900 && station2ThetaInt!=-900) dTheta1x=fabs(station1ThetaInt-station2ThetaInt);
		if (station1ThetaInt!=-900 && station3ThetaInt!=-900 && fabs(station1ThetaInt-station3ThetaInt)>dTheta1x) dTheta1x=fabs(station1ThetaInt-station3ThetaInt);
                if (station1ThetaInt!=-900 && station4ThetaInt!=-900 && fabs(station1ThetaInt-station4ThetaInt)>dTheta1x) dTheta1x=fabs(station1ThetaInt-station4ThetaInt);
                if (station2ThetaInt!=-900 && station3ThetaInt!=-900) dThetaxy=fabs(station2ThetaInt-station3ThetaInt);
		if (station2ThetaInt!=-900 && station4ThetaInt!=-900 && fabs(station2ThetaInt-station4ThetaInt)>dThetaxy) dThetaxy=fabs(station2ThetaInt-station4ThetaInt);
                if (station3ThetaInt!=-900 && station4ThetaInt!=-900 && fabs(station3ThetaInt-station4ThetaInt)>dThetaxy) dThetaxy=fabs(station3ThetaInt-station4ThetaInt);
	}
        if (F("unp_trk_pt", i)>22 && fabs(I("unp_trk_BX", i))<2 && I("unp_trk_nRPC",i)==0 && I("unp_trk_found_hits",i)==1 && I("unp_trk_nHits",i)==4){
                if (I("hit_station",iHit)==1) station1ThetaInt=I("hit_theta_int",iHit);
                if (I("hit_station",iHit)==2) station2ThetaInt=I("hit_theta_int",iHit);
                if (I("hit_station",iHit)==3) station3ThetaInt=I("hit_theta_int",iHit);
                if (I("hit_station",iHit)==4) station4ThetaInt=I("hit_theta_int",iHit);
		dThetaFlag4=1;
                if (station1ThetaInt!=-900 && station2ThetaInt!=-900) dTheta41x=fabs(station1ThetaInt-station2ThetaInt);
                if (station1ThetaInt!=-900 && station3ThetaInt!=-900 && fabs(station1ThetaInt-station3ThetaInt)>dTheta41x) dTheta41x=fabs(station1ThetaInt-station3ThetaInt);
                if (station1ThetaInt!=-900 && station4ThetaInt!=-900 && fabs(station1ThetaInt-station4ThetaInt)>dTheta41x) dTheta41x=fabs(station1ThetaInt-station4ThetaInt);
                if (station2ThetaInt!=-900 && station3ThetaInt!=-900) dTheta4xy=fabs(station2ThetaInt-station3ThetaInt);
                if (station2ThetaInt!=-900 && station4ThetaInt!=-900 && fabs(station2ThetaInt-station4ThetaInt)>dTheta4xy) dTheta4xy=fabs(station2ThetaInt-station4ThetaInt);
                if (station3ThetaInt!=-900 && station4ThetaInt!=-900 && fabs(station3ThetaInt-station4ThetaInt)>dTheta4xy) dTheta4xy=fabs(station3ThetaInt-station4ThetaInt);
        }
         ////////////Chad fills dTheta plots//////////
         //                 :(END):                 //
         /////////////////////////////////////////////


         /////////////Chad fills Time plots///////////
         //                :)START(:                //
         /////////////////////////////////////////////
                if (I("unp_trk_nHits", i)==4 /*&& I("unp_trk_endcap",i)>0*/){
                 if (I("hit_station",iHit)==1 /*&& F("unp_trk_pt", i)>25 */&& I("unp_trk_dBX",i) >1 && I("unp_trk_dPhi_int",i)>300 && I("unp_trk_dTheta_int",i)> 10){
                     ME1BX = I("hit_BX",iHit);
                     ME1Phi = I("hit_phi_int", iHit);
		     ME1Theta = I("hit_theta_int", iHit);
                 } else if (I("hit_station",iHit)==2 /*&& F("unp_trk_pt", i)>25 */&& I("unp_trk_dBX",i) >1 && I("unp_trk_dPhi_int",i)>300 && I("unp_trk_dTheta_int",i)> 10){
                     ME2BX = I("hit_BX",iHit);
                     ME2Phi = I("hit_phi_int", iHit);
                     ME2Theta = I("hit_theta_int", iHit);
                 } else if (I("hit_station",iHit)==3 /*&& F("unp_trk_pt", i)>25 */&& I("unp_trk_dBX",i) >1 && I("unp_trk_dPhi_int",i)>300 && I("unp_trk_dTheta_int",i)> 10){
                     ME3BX = I("hit_BX",iHit);
                     ME3Phi = I("hit_phi_int", iHit);
                     ME3Theta = I("hit_theta_int", iHit);
                 } else if (I("hit_station",iHit)==4 /*&& F("unp_trk_pt", i)>25 */&& I("unp_trk_dBX",i) >1 && I("unp_trk_dPhi_int",i)>300 && I("unp_trk_dTheta_int",i)> 10){
                     ME4BX = I("hit_BX",iHit);
                     ME4Phi = I("hit_phi_int", iHit);
                     ME4Theta = I("hit_theta_int", iHit);
                 }
                }

         //////////////Chad PU rate study/////////////
         //                :)START(:                //
         /////////////////////////////////////////////
         if(/*I("unp_trk_nRPC",i)==0 &&*/ I("unp_trk_nHits",i)>=3 && F("unp_trk_pt", i)<22){
            SingleMuPtLo->Fill(F("unp_trk_pt", i),1);      
            if(I("hit_station", iHit)==1 && I("hit_ring", iHit)==1 && I("hit_pattern",iHit)<6){
               } else {
               SingleMuPtLoPattern->Fill(F("unp_trk_pt", i),1);
            }
            if(I("unp_trk_dBX",i)<=1){ 
               SingleMuPtLodBX1->Fill(F("unp_trk_pt", i),1); 
               if(I("hit_station", iHit)==1 && I("hit_ring", iHit)==1 && I("hit_BX", iHit)!=I("unp_trk_BX",i)){
               } else { 
                  SingleMuPtLodBX1ME11->Fill(F("unp_trk_pt", i),1);
               } 
            }
         }
         if(/*I("unp_trk_nRPC",i)==0 &&*/ I("unp_trk_nHits",i)>=3 && F("unp_trk_pt", i)>=22){
            SingleMuPtHi->Fill(F("unp_trk_pt", i),1);
            if(I("hit_station", iHit)==1 && I("hit_ring", iHit)==1 && I("hit_pattern",iHit)<6){
               } else {
               SingleMuPtHiPattern->Fill(F("unp_trk_pt", i),1);
            }
            if(I("unp_trk_dBX",i)<=1){ 
               SingleMuPtHidBX1->Fill(F("unp_trk_pt", i),1);
               if(I("hit_station", iHit)==1 && I("hit_ring", iHit)==1 && I("hit_BX", iHit)!=I("unp_trk_BX",i)){
                  } else {
                  SingleMuPtHidBX1ME11->Fill(F("unp_trk_pt", i),1);
               }
            }
         }
         if(/*I("unp_trk_nRPC",i)==0 &&*/ I("unp_trk_nHits",i)>=3){
            SingleMuPtAll->Fill(F("unp_trk_pt", i),1);
            if(I("hit_station", iHit)==1 && I("hit_ring", iHit)==1 && I("hit_pattern",iHit)<6){
               } else {
               SingleMuPtAllPattern->Fill(F("unp_trk_pt", i),1);
            }
            if(I("unp_trk_dBX",i)<=1){
               SingleMuPtAlldBX1->Fill(F("unp_trk_pt", i),1);
               if(I("hit_station", iHit)==1 && I("hit_ring", iHit)==1 && I("hit_BX", iHit)!=I("unp_trk_BX",i)){
                  } else {
                  SingleMuPtAlldBX1ME11->Fill(F("unp_trk_pt", i),1);
               }
            }
         }
         //////////////Chad PU rate study/////////////
         //                 :(END):                 //
         /////////////////////////////////////////////



	  if (verbose) std::cout << "  - CSC LCT with BX = " << I("hit_BX", iHit) << ", endcap = " << I("hit_endcap", iHit)
				 << ", station = " << I("hit_station", iHit) << ", ring = " << I("hit_ring", iHit) << std::endl; 
	} else if (I("hit_isRPC", iHit) == 1) {
	  if (verbose) std::cout << "  - RPC hit with BX = " << I("hit_BX", iHit) << ", endcap = " << I("hit_endcap", iHit)
				 << ", station = " << I("hit_station", iHit) << ", ring = " << I("hit_ring", iHit) << std::endl; 
	}
      }//End of iHit for loop
	

      if (dThetaFlag3!=-1){
      if (dTheta1x>=16 || dThetaxy>=8){ 
         dTheta1x=-1;
         dThetaxy=-1;
      }
	Station3dThetaCheck->Fill(dTheta1x,dThetaxy); 
      }
      if (dThetaFlag4!=-1){  
      if (dTheta41x>=16 || dTheta4xy>=8){
	 dTheta41x=-1;
      	 dTheta4xy=-1;
      }
	Station4dThetaCheck->Fill(dTheta41x,dTheta4xy);
      }
      Station3dThetaCheck->GetXaxis()->SetTitle("dTheta(1-X)");
      Station3dThetaCheck->GetYaxis()->SetTitle("dTheta(X-Y)");

      Station4dThetaCheck->GetXaxis()->SetTitle("dTheta(1-X)");
      Station4dThetaCheck->GetYaxis()->SetTitle("dTheta(X-Y)");

      if (ME1BX != 100 && ME2BX != 100){
         ME1vME2BX->Fill(ME1BX,ME2BX);
         ME1vME2Phi->Fill(ME1Phi,ME2Phi);
         ME1vME2Theta->Fill(ME1Theta,ME2Theta);
      } 
      if (ME1BX != 100 && ME3BX != 100){
         ME1vME3BX->Fill(ME1BX,ME3BX);
         ME1vME3Phi->Fill(ME1Phi,ME3Phi);
         ME1vME3Theta->Fill(ME1Theta,ME3Theta);
      }
      if (ME1BX != 100 && ME4BX != 100){
         ME1vME4BX->Fill(ME1BX,ME4BX);
         ME1vME4Phi->Fill(ME1Phi,ME4Phi);
         ME1vME4Theta->Fill(ME1Theta,ME4Theta);
      }
      if (ME2BX != 100 && ME3BX != 100){
         ME2vME3BX->Fill(ME2BX,ME3BX);
         ME2vME3Phi->Fill(ME2Phi,ME3Phi);
         ME2vME3Theta->Fill(ME2Theta,ME3Theta);
      }
      if (ME2BX != 100 && ME4BX != 100){
         ME2vME4BX->Fill(ME2BX,ME4BX);
         ME2vME4Phi->Fill(ME2Phi,ME4Phi);
         ME2vME4Theta->Fill(ME2Theta,ME4Theta);
      }
      if (ME3BX != 100 && ME4BX != 100){
         ME3vME4BX->Fill(ME3BX,ME4BX);    
         ME3vME4Phi->Fill(ME3Phi,ME4Phi);
         ME3vME4Theta->Fill(ME3Theta,ME4Theta);
      }
    }
    
  } // End loop: for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Finished looping over the events *******" << std::endl;

  delete in_chain;

  std::cout << "\nDone with Read_FlatNtuple(). Exiting.\n" << std::endl;
  f->Write();
} // End function: void Read_FlatNtuple()

