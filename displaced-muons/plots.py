import ROOT


def get_plots(outfile):
    plots = {}

    plots['EventsByLcts'] = ROOT.TH1D(
        'EventsByLcts', 'Events by LCT Count', 100, 0, 100)
    plots['trkPt'] = ROOT.TH1D('trkPt', 'Track pT', 100, 0, 100)
    plots['trkNLcts'] = ROOT.TH1D('trkNLcts', 'Track # LCTs', 5, 0, 5)

    outfile.mkdir("XbyLS").cd()
    for obj in ['Events', 'LCTs', 'Chambers', 'Tracks']:
        plots[obj + 'ByLS'] = ROOT.TH1D(obj + 'ByLS',
                                        obj + ' by LS', 1000, 0, 1000)
        plots[obj + 'ByLS'].GetXaxis().SetTitle('Lumi Section #')

    outfile.mkdir("XbyPU").cd()
    for obj in ['Events', 'LCTs', 'Chambers', 'Tracks']:
        plots[obj + 'ByPU'] = ROOT.TH1D(obj +
                                        'ByPU', obj + ' by PU', 100, 0, 100)
        plots[obj + 'ByPU'].GetXaxis().SetTitle('Lumi Section Avg Pileup')

    outfile.mkdir("XbyDelLumi").cd()
    for obj in ['Events', 'LCTs', 'Chambers', 'Tracks']:
        plots[obj + 'ByDelLumi'] = ROOT.TH1D(obj + 'ByDelLumi',
                                             obj + ' by Delivered Lumi', 20000, 0, 20000)
        plots[obj + 'ByDelLumi'].GetXaxis().SetTitle('Lumi Section Delivered Lumi')

    return plots


def post_fill(outfile, plots):

    # Derivative plots --------------------------------------------------------

    outfile.cd()
    plots['LCTsOEvtsByLS'] = plots['LCTsByLS'].Clone('LCTsOEvtsByLS')
    plots['LCTsOEvtsByLS'].SetTitle('LCTs/Event by Lumi Section')
    plots['LCTsOEvtsByLS'].Divide(plots['EventsByLS'])

    plots['LCTsOEvtsByPU'] = plots['LCTsByPU'].Clone('LCTsOEvtsByPU')
    plots['LCTsOEvtsByPU'].SetTitle('LCTs/Event by PU')
    plots['LCTsOEvtsByPU'].Divide(plots['EventsByPU'])

    plots['LCTsOEvtsByDelLumi'] = plots['LCTsByDelLumi'].Clone(
        'LCTsOEvtsByDelLumi')
    plots['LCTsOEvtsByDelLumi'].SetTitle('LCTs/Event by Delivered Lumi')
    plots['LCTsOEvtsByDelLumi'].Divide(plots['EventsByDelLumi'])

    # Make fits ---------------------------------------------------------------

    ROOT.gStyle.SetOptFit(1)
    for p in plots:
        plots[p].Fit("pol1")

    # Cleanup -----------------------------------------------------------------

    for p in plots:
        plots[p].Draw()
    outfile.Write()
