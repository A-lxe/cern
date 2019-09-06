
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

const std::map<std::pair<int, int>, int> histIndexCSC = {
    {{1, 4}, 9}, {{1, 1}, 8}, {{1, 2}, 7}, {{1, 3}, 6}, {{2, 1}, 5},
    {{2, 2}, 4}, {{3, 1}, 3}, {{3, 2}, 2}, {{4, 1}, 1}, {{4, 2}, 0}};
// cscid + cscid_offset, endcap * (station - 0.5)

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

void fillOccupancyHit(TH2D *plot, int hit) { fillOccupancyHit(plot, hit, 1); }

TChain *prepareInput(int dataset, bool skim) {
  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);

  // List of input files
  std::vector<TString> in_file_names;
  TString store = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/";
  TString in_dir = "ntuples/HADD/";
  TString file_name;
  TString treeName;

  switch (dataset) {
  case 0:
    file_name.Form(
        "%s/%s/NTuple_ZeroBias8b4e_FlatNtuple%sRun_302674_2017_09_30.root",
        store.Data(), in_dir.Data(), skim ? "_Skim_" : "_");
    std::cout << "Adding file " << file_name.Data() << std::endl;
    in_file_names.push_back(file_name.Data());
    treeName = "ntuple/tree";
    break;
  case 1:
    file_name.Form(
        "%s/%s/"
        "NTuple_ZeroBiasIsolatedBunches_FlatNtuple%sRun_302674_2017_09_30."
        "root",
        store.Data(), in_dir.Data(), skim ? "_Skim_" : "_");
    std::cout << "Adding file " << file_name.Data() << std::endl;
    in_file_names.push_back(file_name.Data());
    treeName = "ntuple/tree";
    break;
  case 2:
    file_name.Form(
        "%s/%s/"
        "NTuple_ZeroBiasNominalTrains_FlatNtuple%sRun_302674_2017_09_30."
        "root",
        store.Data(), in_dir.Data(), skim ? "_Skim_" : "_");
    std::cout << "Adding file " << file_name.Data() << std::endl;
    in_file_names.push_back(file_name.Data());
    treeName = "ntuple/tree";
    break;
  case 3:
    for (int i = 2; i < 5; i++) {
      file_name.Form(
          "%s/%s/"
          "NTuple_ZeroBias%i_FlatNtuple_Run_306091_2018_02_24_ZB%i.root",
          store.Data(), in_dir.Data(), i, i);
      std::cout << "Adding file " << file_name.Data() << std::endl;
      in_file_names.push_back(file_name.Data());
    }
    treeName = "FlatNtupleData/tree";
    break;
  case 4:
    file_name.Form(
        "%s/%s/"
        "NTuple_SingleMuon_FlatNtuple_Run_306092_2018_03_02_SingleMu.root",
        store.Data(), in_dir.Data());
    std::cout << "Adding file " << file_name.Data() << std::endl;
    in_file_names.push_back(file_name.Data());
    treeName = "FlatNtupleData/tree";
    break;
  case 5:
    for (int i = 1; i < 5; i++) {
      file_name.Form(
          "%s/%s/"
          "NTuple_ZeroBias%i_FlatNtuple_Run_306091_2018_03_02_ZB%i.root",
          store.Data(), in_dir.Data(), i, i);
      std::cout << "Adding file " << file_name.Data() << std::endl;
      in_file_names.push_back(file_name.Data());
    }
    treeName = "FlatNtupleData/tree";
    break;
  }

  // Open all input files
  for (int i = 0; i < in_file_names.size(); i++) {
    if (!gSystem->AccessPathName(in_file_names.at(i)))
      file_tmp = TFile::Open(in_file_names.at(i)); // Check if file exists
    if (!file_tmp) {
      std::cout << "ERROR: could not open data file " << in_file_names.at(i)
                << std::endl;
      std::exit(1);
    }
  }

  // Add tree from the input files to the TChain
  TChain *in_chain = new TChain(treeName);
  for (int i = 0; i < in_file_names.size(); i++) {
    in_chain->Add(in_file_names.at(i));
  }
  return in_chain;
}
