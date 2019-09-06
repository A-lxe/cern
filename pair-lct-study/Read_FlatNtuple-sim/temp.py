import ROOT

d = f.Get("LCT Quality")
c = ROOT.TCanvas("AlctVsClct", "Alct qual Vs Clct qual by Station", 1200, 1200)
c.Divide(4, 4)

padIter = iter(range(1, 17))
for stat in [1, 2, 3, 4]:
    for lctType in [("All Lct", "AllLct"), ("Phi Pair Lct", "PhiPair"), ("Theta Pair Lct", "ThetaPair"), ("Quad Lct", "Quad")]:
        subDir = d.Get(lctType[0])
        pad = c.cd(next(padIter))
        pad.SetRightMargin(0.05)
        hist = subDir.Get("{0}AlctVsClctStat{1}".format(lctType[1], stat))
        hist.Draw("colz")
        hist.GetXaxis().SetRange(5,7)
        hist.GetYaxis().SetRange(2,4)
        hist.GetXaxis().SetTitle("Clct Quality")
        hist.GetYaxis().SetTitle("Alct Quality")
        hist.SetTitleSize(1)
        hist.SetStats(False)
        pad.Draw()

c.Draw()


zerobias_8b4e = {
    "idx": 0,
    "source": "NTuple_ZeroBias8b4e_FlatNtuple_Run_302674_2017_09_30.root",
    "name": "Run 302674 ZeroBias 8b4e"
}
zerobias_isolated = {
    "idx": 1,
    "source": "NTuple_ZeroBiasIsolatedBunches_FlatNtuple_Run_302674_2017_09_30.root",
    "name": "Run 302674 ZeroBias Isolated Bunches"
}
zerobias_nominal = {
    "idx": 2,
    "source": "NTuple_ZeroBiasNominalTrains_FlatNtuple_Run_302674_2017_09_30.root",
    "name": "Run 302674 ZeroBias Nominal Trains"
}
zerobias_306091 = {
    "idx": 3,
    "source": "NTuple_ZeroBias2_FlatNtuple_Run_306091_2018_02_24_ZB2.root",
    "name": "Run 306091 2018/02/24 ZeroBias"
}
singlemu = {
    "idx": 4,
    "source": "NTuple_SingleMuon_FlatNtuple_Run_306092_2018_03_02_SingleMu.root",
    "name": "Run 306092 SingleMu"
}
zerobias_306091_2 = {
    "idx": 5,
    "source": "NTuple_ZeroBias[#]_FlatNtuple_Run_306091_2018_03_02_ZB[#].root",
    "name": "Run 306091 2018/03/02 ZeroBias"
}
singlemu_recomatch = {
    "idx": 5,
    "source": "NTuple_SingleMuon_FlatNtuple_Run_306092_2018_03_02_SingleMu.root"
    "name": "Run 306092 SingleMu Reco Muon Match"
}


ds = zerobias_306091
event_count = 2000000  # This is just used to select the file, not control #events
base_dir = "/eos/cms/store/user/abrinke1/EMTF/Emulator/ntuples/HADD/"

fname = "./dataset-{0}_events-{1}.root".format(ds.["idx"], event_count)
source_file = base_dir + ds.["source"]
print("The dataset being used is {0} events of {1} from {2}".format(
    event_count, ds.["name"], source_file))

f = ROOT.TFile(fname)

useLogY = True
