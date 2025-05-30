import sys
from collections import OrderedDict
import ROOT
sys.path.append("./CMSPLOTS")
from myFunction import DrawHistos


print("Starting")

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(4)


def makePlots(elefile, pionfile, suffix=""):
    chains = OrderedDict()
    chains['ele'] = ROOT.TChain("tree")
    chains['pion'] = ROOT.TChain("tree")

    loops = [('ele', elefile), ('pion', pionfile)]

    rdfs = OrderedDict()

    for part, filename in loops:
        nfiles = 0
        with open(filename) as f:
            print(f"Reading {filename}")
            elefiles = f.readlines()
            for line in elefiles:
                line = line.strip()
                if line.startswith("#"):
                    continue
                nfiles += 1
                print(f"{part} " + line)
                chains[part].Add(line)

                if nfiles > 100:
                    break
        rdfs[part] = ROOT.RDataFrame(chains[part])

    rdfs_new = OrderedDict()
    for part in rdfs.keys():
        rdfs_new[part] = rdfs[part].Define( "eTotaltruth", "eLeaktruth + eCalotruth + eWorldtruth + eInvisible") \
            .Define("eTotalGeant", "eLeaktruth + eCalotruth + eWorldtruth") \
            .Define("truthhit_r", f"sqrt((truthhit_x-2.5)*(truthhit_x-2.5) + (truthhit_y+2.5)*(truthhit_y+2.5))") \
            .Define("nOPsCer_Cer_perGeV", "nOPsCer_Cer / eCentruth") \
            .Define("nOPsCer_Sci_perGeV", "nOPsCer_Sci / eScintruth") \
            .Define("eScinSummed", "Sum(truthhit_edep * (truthhit_calotype == 2))") \
            .Define("eCentrSummed", "Sum(truthhit_edep * (truthhit_calotype == 3))")

    rdfs = rdfs_new

    nEvts = OrderedDict()
    nEvts['ele'] = rdfs['ele'].Count().GetValue()
    nEvts['pion'] = rdfs['pion'].Count().GetValue()
    print("Number of events for electrons: ", nEvts['ele'])
    print("Number of events for pions: ", nEvts['pion'])

    rdfs['ele'] = rdfs['ele'].Define("eweight", f"truthhit_edep/ {nEvts['ele']}")
    rdfs['pion'] = rdfs['pion'].Define("eweight", f"truthhit_edep/ {nEvts['pion']}")

    histos = OrderedDict()
    figures = ['eLeaktruth', 'eCalotruth', 'eTotaltruth', 'eTotalGeant',
               'eRodtruth', 'eCentruth', 'eScintruth',
               'truthhit_x', 'truthhit_y', 'truthhit_z', 'truthhit_r',
               'time', 'time_zoomed',
               "nOPsCer_Cer", "nOPsCer_Sci", "nOPsCer_Cer_perGeV", "nOPsCer_Sci_perGeV",
               "nOPsCer_Cer_vs_nOPsCer_Sci", "nOPsCer_Cer_vs_eCentruth", "nOPsCer_Sci_vs_eScintruth",
               "nOPsCer_Cer_vs_eCentrSummed", "nOPsCer_Sci_vs_eScinSummed",
               "eScinSummed_vs_eScintruth",
               "nOPsCer_Cer_vs_eScintruth",
               'truthhit_x_vs_truthhit_y', 'truthhit_x_vs_truthhit_z', 'truthhit_r_vs_truthhit_z',
               'time_vs_truthhit_z', 'time_vs_truthhit_r']
    #evtlist = [1, 10, 20, 30, 40]
    evtlist = []
    for i in evtlist:
        figures.append(f"event_{i}_truthhit_x_vs_truthhit_y")
        figures.append(f"event_{i}_truthhit_x_vs_truthhit_z")
        figures.append(f"event_{i}_truthhit_r_vs_truthhit_z")
        figures.append(f"event_{i}_time_vs_truthhit_z")
        figures.append(f"event_{i}_time_vs_truthhit_r")

    for fig in figures:
        histos[fig] = OrderedDict()

    values = OrderedDict()
    values['eCentruth'] = 1.2
    values['eScintruth'] = 0.8
    values['nOPsCer_Cer'] = 1.3e5
    values['nOPsCer_Sci'] = 1.0e5
    values['nOPsCer_Cer_perGeV'] = values['nOPsCer_Cer'] / values['eCentruth']
    values['nOPsCer_Sci_perGeV'] = values['nOPsCer_Sci'] / values['eScintruth']

    for part, rdf in rdfs.items():
        suffix_h = suffix + "_" + part

        histos['eLeaktruth'][part] = rdf.Histo1D(
            ("eLeaktruth" + suffix_h, "eLeaktruth", 50, 0, 20.0), "eLeaktruth")
        histos['eCalotruth'][part] = rdf.Histo1D(
            ("eCalotruth" + suffix_h, "eCalotruth", 50, 0., 102.0), "eCalotruth")
        histos['eTotaltruth'][part] = rdf.Histo1D(
            ("eTotaltruth" + suffix_h, "eTotaltruth", 50, 90., 110.0), "eTotaltruth")
        histos['eTotalGeant'][part] = rdf.Histo1D(
            ("eTotalGeant" + suffix_h, "eTotalGeant", 50, 80., 110.0), "eTotalGeant")
        histos['eRodtruth'][part] = rdf.Histo1D(
            ("eRodtruth" + suffix_h, "eRodtruth", 50, 50, 102.0), "eRodtruth")
        histos['eCentruth'][part] = rdf.Histo1D(
            ("eCentruth" + suffix_h, "eCentruth", 50, 0, values["eCentruth"]), "eCentruth")
        histos['eScintruth'][part] = rdf.Histo1D(
            ("eScintruth" + suffix_h, "eScintruth", 50, 0, values['eScintruth']), "eScintruth")

        # energy weighted
        histos['truthhit_x'][part] = rdf.Histo1D(
            ("truthhit_x" + suffix_h, "truthhit_x", 50, -20, 20), "truthhit_x", "eweight")
        histos['truthhit_y'][part] = rdf.Histo1D(
            ("truthhit_y" + suffix_h, "truthhit_y", 50, -20, 20), "truthhit_y", "eweight")
        histos['truthhit_z'][part] = rdf.Histo1D(
            ("truthhit_z" + suffix_h, "truthhit_z", 100, -100, 100), "truthhit_z", "eweight")
        histos['truthhit_r'][part] = rdf.Histo1D(
            ("truthhit_r" + suffix_h, "truthhit_r", 50, 0, 30), "truthhit_r", "eweight")

        histos['time'][part] = rdf.Histo1D(
            ("time" + suffix_h, "time", 50, 0, 50), "truthhit_globaltime", "eweight")
        histos['time_zoomed'][part] = rdf.Histo1D(
            ("time_zoomed" + suffix_h, "time_zoomed", 50, 0, 10), "truthhit_globaltime", "eweight")

        histos['nOPsCer_Cer'][part] = rdf.Histo1D(
            ("nOPsCer_Cer" + suffix_h, "nOPsCer_Cer", 50, 0, values["nOPsCer_Cer"]), "nOPsCer_Cer", "eweight")
        histos['nOPsCer_Sci'][part] = rdf.Histo1D(
            ("nOPsCer_Sci" + suffix_h, "nOPsCer_Sci", 50, 0, values["eScintruth"]), "nOPsCer_Sci", "eweight")
        histos['nOPsCer_Cer_perGeV'][part] = rdf.Histo1D(
            ("nOPsCer_Cer_perGeV" + suffix_h, "nOPsCer_Cer_perGeV", 50, 0, values["nOPsCer_Cer_perGeV"]), "nOPsCer_Cer_perGeV", "eweight")
        histos['nOPsCer_Sci_perGeV'][part] = rdf.Histo1D(
            ("nOPsCer_Sci_perGeV" + suffix_h, "nOPsCer_Sci_perGeV", 50, 0, values["nOPsCer_Sci_perGeV"]), "nOPsCer_Sci_perGeV", "eweight")
        histos['nOPsCer_Cer_vs_nOPsCer_Sci'][part] = rdf.Histo2D(
            ("nOPsCer_Cer_vs_nOPsCer_Sci" + suffix_h, "nOPsCer_Cer_vs_nOPsCer_Sci", 50, 0, values["nOPsCer_Cer"], 50, 0, values["nOPsCer_Sci"]), "nOPsCer_Cer", "nOPsCer_Sci", "eweight")
        histos['nOPsCer_Cer_vs_eCentruth'][part] = rdf.Histo2D(
            ("nOPsCer_Cer_vs_eCentruth" + suffix_h, "nOPsCer_Cer_vs_eCentruth", 50, 0, values["nOPsCer_Cer"], 50, 0, values["eCentruth"]), "nOPsCer_Cer", "eCentruth", "eweight")
        histos['nOPsCer_Sci_vs_eScintruth'][part] = rdf.Histo2D(
            ("nOPsCer_Sci_vs_eScintruth" + suffix_h, "nOPsCer_Sci_vs_eScintruth", 50, 0, values["nOPsCer_Sci"], 50, 0, values['eScintruth']), "nOPsCer_Sci", "eScintruth", "eweight")
        histos['nOPsCer_Cer_vs_eCentrSummed'][part] = rdf.Histo2D(
            ("nOPsCer_Cer_vs_eCentrSummed" + suffix_h, "nOPsCer_Cer_vs_eCentrSummed", 50, 0, values["nOPsCer_Cer"], 50, 0, values["eCentruth"]), "nOPsCer_Cer", "eCentrSummed", "eweight")
        histos['nOPsCer_Sci_vs_eScinSummed'][part] = rdf.Histo2D(
            ("nOPsCer_Sci_vs_eScinSummed" + suffix_h, "nOPsCer_Sci_vs_eScinSummed", 50, 0, values["nOPsCer_Sci"], 50, 0, values["eScintruth"]), "nOPsCer_Sci", "eScinSummed", "eweight")
        histos['eScinSummed_vs_eScintruth'][part] = rdf.Histo2D(
            ("eScinSummed_vs_eScintruth" + suffix_h, "eScinSummed_vs_eScintruth", 50, 0, values["eCentruth"], 50, 0, values["eScintruth"]), "eScinSummed", "eScintruth", "eweight")
        histos['nOPsCer_Cer_vs_eScintruth'][part] = rdf.Histo2D(
            ("nOPsCer_Cer_vs_eScintruth" + suffix_h, "nOPsCer_Cer_vs_eScintruth", 50, 0, values["nOPsCer_Cer"], 50, 0, values["eScintruth"]), "nOPsCer_Cer", "eScintruth", "eweight")


        histos['truthhit_x_vs_truthhit_y'][part] = rdf.Histo2D(
            ("truthhit_x_vs_truthhit_y" + suffix_h, "truthhit_x_vs_truthhit_y", 50, -20, 20, 50, -20, 20), "truthhit_x", "truthhit_y", "eweight")
        histos['truthhit_x_vs_truthhit_z'][part] = rdf.Histo2D(
            ("truthhit_x_vs_truthhit_z" + suffix_h, "truthhit_x_vs_truthhit_z", 50, -20, 20, 100, -100, 100), "truthhit_x", "truthhit_z", "eweight")
        histos['truthhit_r_vs_truthhit_z'][part] = rdf.Histo2D(
            ("truthhit_r_vs_truthhit_z" + suffix_h, "truthhit_r_vs_truthhit_z", 50, 0, 30, 100, -100, 100), "truthhit_r", "truthhit_z", "eweight")

        histos['time_vs_truthhit_z'][part] = rdf.Histo2D(
            ("time_vs_truthhit_z" + suffix_h, "time_vs_truthhit_z", 50, 0, 20, 100, -100, 100), "truthhit_globaltime", "truthhit_z", "eweight")
        histos['time_vs_truthhit_r'][part] = rdf.Histo2D(
            ("time_vs_truthhit_r" + suffix_h, "time_vs_truthhit_r", 50, 0, 20, 50, 0, 30), "truthhit_globaltime", "truthhit_r", "eweight")

        # some event displays
        for i in evtlist:
            rdf_event = rdf.Filter(f"rdfentry_ == {i}")
            histos[f"event_{i}_truthhit_x_vs_truthhit_y"][part] = rdf_event.Histo2D(
                (f"event_{i}_truthhit_x_vs_truthhit_y" + suffix_h, f"event_{i}_truthhit_x_vs_truthhit_y", 50, -20, 20, 50, -20, 20), "truthhit_x", "truthhit_y", "truthhit_edep")
            histos[f"event_{i}_truthhit_x_vs_truthhit_z"][part] = rdf_event.Histo2D(
                (f"event_{i}_truthhit_x_vs_truthhit_z" + suffix_h, f"event_{i}_truthhit_x_vs_truthhit_z", 50, -20, 20, 100, -100, 100), "truthhit_x", "truthhit_z", "truthhit_edep")
            histos[f"event_{i}_truthhit_r_vs_truthhit_z"][part] = rdf_event.Histo2D(
                (f"event_{i}_truthhit_r_vs_truthhit_z" + suffix_h, f"event_{i}_truthhit_r_vs_truthhit_z", 50, 0, 30, 100, -100, 100), "truthhit_r", "truthhit_z", "truthhit_edep")

            histos[f"event_{i}_time_vs_truthhit_z"][part] = rdf_event.Histo2D(
                (f"event_{i}_time_vs_truthhit_z" + suffix_h, f"event_{i}_time_vs_truthhit_z", 50, 0, 20, 100, -100, 100), "truthhit_globaltime", "truthhit_z", "truthhit_edep")
            histos[f"event_{i}_time_vs_truthhit_r"][part] = rdf_event.Histo2D(
                (f"event_{i}_time_vs_truthhit_r" + suffix_h, f"event_{i}_time_vs_truthhit_r", 50, 0, 20, 50, 0, 30), "truthhit_globaltime", "truthhit_r", "truthhit_edep")

    colormaps = {
        'ele': 2,
        'pion': 3,
        'neu': 4
    }


    def GetColors(ene_fracs):
        colors = []
        for str in ene_fracs.keys():
            part, _ = str.split("_")
            color = colormaps[part]
            colors.append(color)
        return colors


    args = {
        'dology': True,
        'donormalize': True,
        'mycolors': [colormaps[part] for part in rdfs.keys()],
        "MCOnly": True,
        'addOverflow': True,
        'addUnderflow': True
    }

    print("Drawing")

    DrawHistos(list(histos['eLeaktruth'].values()), list(histos['eLeaktruth'].keys(
    )), 0, 20, "Leakage Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eLeaktruth_" + suffix, **args)
    DrawHistos(list(histos['eCalotruth'].values()), list(histos['eCalotruth'].keys(
    )), 0., 102, "Calo Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eCalotruth_" + suffix, **args)
    DrawHistos(list(histos['eTotaltruth'].values()), list(histos['eTotaltruth'].keys(
    )), 80., 110, "Total Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eTotaltruth_" + suffix, **args)
    DrawHistos(list(histos['eTotalGeant'].values()), list(histos['eTotalGeant'].keys(
    )), 80., 110, "Total 'Visible' Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eTotalGeant_" + suffix, **args)
    DrawHistos(list(histos['eRodtruth'].values()), list(histos['eRodtruth'].keys(
    )), 50, 102, "Rod Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eRodtruth_" + suffix, **args)
    DrawHistos(list(histos['eCentruth'].values()), list(histos['eCentruth'].keys(
    )), 0, values["eCentruth"], "CFiber Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eCentruth_" + suffix, **args)
    DrawHistos(list(histos['eScintruth'].values()), list(histos['eScintruth'].keys(
    )), 0, values['eScintruth'], "SFiber Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eScintruth_" + suffix, **args)


    args['donormalize'] = False
    DrawHistos(list(histos['truthhit_x'].values()), list(histos['truthhit_x'].keys(
    )), -20, 20, "x [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "truthhit_x_" + suffix, **args)
    DrawHistos(list(histos['truthhit_y'].values()), list(histos['truthhit_y'].keys(
    )), -20, 20, "y [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "truthhit_y_" + suffix, **args)
    DrawHistos(list(histos['truthhit_z'].values()), list(histos['truthhit_z'].keys(
    )), -100, 100, "z [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "truthhit_z_" + suffix, **args)
    DrawHistos(list(histos['truthhit_r'].values()), list(histos['truthhit_r'].keys(
    )), 0, 30, "r [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "truthhit_r_" + suffix, **args)

    DrawHistos(list(histos['time'].values()), list(histos['time'].keys(
    )), 0, 100, "Time [ns]", 1e-3, 1e2, "Deposited Energy [GeV]", "time_" + suffix, **args)
    DrawHistos(list(histos['time_zoomed'].values()), list(histos['time_zoomed'].keys(
    )), 0, 10, "Time [ns]", 1e-3, 1e2, "Deposited Energy [GeV]", "time_zoomed_" + suffix, **args)

    DrawHistos(list(histos['nOPsCer_Cer'].values()), list(histos['nOPsCer_Cer'].keys(
    )), 0, values["nOPsCer_Cer"], "Number of OPs in Cerenkov", 1e-3, 1e2, "Fraction of events", "nOPsCer_Cer_" + suffix, **args)
    DrawHistos(list(histos['nOPsCer_Sci'].values()), list(histos['nOPsCer_Sci'].keys(
    )), 0, values["nOPsCer_Sci"], "Number of OPs in Scintillator", 1e-3, 1e2, "Fraction of events", "nOPsCer_Sci_" + suffix, **args)
    DrawHistos(list(histos['nOPsCer_Cer_perGeV'].values()), list(histos['nOPsCer_Cer_perGeV'].keys(
    )), 0, values["nOPsCer_Cer_perGeV"], "Number of OPs in Cerenkov per GeV", 1e-3, 1e2, "Fraction of events", "nOPsCer_Cer_perGeV_" + suffix, **args)
    DrawHistos(list(histos['nOPsCer_Sci_perGeV'].values()), list(histos['nOPsCer_Sci_perGeV'].keys(
    )), 0, values["nOPsCer_Sci_perGeV"], "Number of OPs in Scintillator per GeV", 1e-3, 1e2, "Fraction of events", "nOPsCer_Sci_perGeV_" + suffix, **args)


    # 2D plots
    args['dology'] = False
    args['drawoptions'] = "colz"
    args['dologz'] = True
    args['zmax'] = 1e-1
    args['zmin'] = 1e-4
    args['doth2'] = True
    args['addOverflow'] = False
    args['addUnderflow'] = False
    for part in rdfs.keys():
        DrawHistos([histos['truthhit_x_vs_truthhit_y'][part]], [], -20, 20,
                   "x [cm]", -20, 20, "y [cm]", f"truthhit_x_vs_truthhit_y_{part}_{suffix}", **args)
        DrawHistos([histos['truthhit_x_vs_truthhit_z'][part]], [], -20, 20,
                   "x [cm]", -100, 100, "z [cm]", f"truthhit_x_vs_truthhit_z_{part}_{suffix}", **args)
        DrawHistos([histos['truthhit_r_vs_truthhit_z'][part]], [], 0, 30,
                   "r [cm]", -100, 100, "z [cm]", f"truthhit_r_vs_truthhit_z_{part}_{suffix}", **args)

        DrawHistos([histos['time_vs_truthhit_z'][part]], [], 0, 20,
                   "Time [ns]", -100, 100, "z [cm]", f"time_vs_truthhit_z_{part}_{suffix}", **args)
        DrawHistos([histos['time_vs_truthhit_r'][part]], [], 0, 20,
                   "Time [ns]", 0, 30, "r [cm]", f"time_vs_truthhit_r_{part}_{suffix}", **args)

        DrawHistos([histos['nOPsCer_Cer_vs_nOPsCer_Sci'][part]], [], 0, values["nOPsCer_Cer"],
                   "Number of OPs in Cerenkov", 0, values["nOPsCer_Sci"], "Number of OPs in Scintillator", f"nOPsCer_Cer_vs_nOPsCer_Sci_{part}_{suffix}", **args)
        DrawHistos([histos['nOPsCer_Cer_vs_eCentruth'][part]], [], 0, values["nOPsCer_Cer"],
                   "Number of OPs in Cerenkov", 0, values["eCentruth"], "CFiber Energy [GeV]", f"nOPsCer_Cer_vs_eCentruth_{part}_{suffix}", **args)  
        DrawHistos([histos['nOPsCer_Sci_vs_eScintruth'][part]], [], 0, values["nOPsCer_Sci"],
                   "Number of OPs in Scintillator", 0, values["eScintruth"], "SFiber Energy [GeV]", f"nOPsCer_Sci_vs_eScintruth_{part}_{suffix}", **args)
        DrawHistos([histos['nOPsCer_Cer_vs_eCentrSummed'][part]], [], 0, values["nOPsCer_Cer"],
                    "Number of OPs in Cerenkov", 0, values["eCentruth"], "CFiber Energy Summed [GeV]", f"nOPsCer_Cer_vs_eCentrSummed_{part}_{suffix}", **args)
        DrawHistos([histos['nOPsCer_Sci_vs_eScinSummed'][part]], [], 0, values["nOPsCer_Sci"],
                    "Number of OPs in Scintillator", 0, values["eScintruth"], "SFiber Energy Summed [GeV]", f"nOPsCer_Sci_vs_eScinSummed_{part}_{suffix}", **args)
        DrawHistos([histos['eScinSummed_vs_eScintruth'][part]], [], 0, values["eScintruth"],
                    "SFiber Energy Summed [GeV]", 0, values["eScintruth"], "SFiber Energy [GeV]", f"eScinSummed_vs_eScintruth_{part}_{suffix}", **args)
        DrawHistos([histos['nOPsCer_Cer_vs_eScintruth'][part]], [], 0, values["nOPsCer_Cer"],
                    "Number of OPs in Cerenkov", 0, values["eScintruth"], "SFiber Energy [GeV]", f"nOPsCer_Cer_vs_eScintruth_{part}_{suffix}", **args)

        # event displays
        for i in evtlist:
            DrawHistos([histos[f"event_{i}_truthhit_x_vs_truthhit_y"][part]], [
            ], -20, 20, "x [cm]", -20, 20, "y [cm]", f"event_{i}_truthhit_x_vs_truthhit_y_{part}_{suffix}", **args)
            DrawHistos([histos[f"event_{i}_truthhit_x_vs_truthhit_z"][part]], [
            ], -20, 20, "x [cm]", -100, 100, "z [cm]", f"event_{i}_truthhit_x_vs_truthhit_z_{part}_{suffix}", **args)
            DrawHistos([histos[f"event_{i}_truthhit_r_vs_truthhit_z"][part]], [
            ], 0, 30, "r [cm]", -100, 100, "z [cm]", f"event_{i}_truthhit_r_vs_truthhit_z_{part}_{suffix}", **args)

            DrawHistos([histos[f"event_{i}_time_vs_truthhit_z"][part]], [
            ], 0, 20, "Time [ns]", -100, 100, "z [cm]", f"event_{i}_time_vs_truthhit_z_{part}_{suffix}", **args)
            DrawHistos([histos[f"event_{i}_time_vs_truthhit_r"][part]], [
            ], 0, 20, "Time [ns]", 0, 30, "r [cm]", f"event_{i}_time_vs_truthhit_r_{part}_{suffix}", **args)
            
if __name__ == "__main__":
    elefile = "inputs/electrons_10GeV.txt"
    pionfile = "inputs/pions_10GeV.txt"
    makePlots(elefile, pionfile, suffix="_10GeV")
    elefile = "inputs/electrons_20GeV.txt"
    pionfile = "inputs/pions_20GeV.txt"
    makePlots(elefile, pionfile, suffix="_20GeV")
    elefile = "inputs/electrons_30GeV.txt"
    pionfile = "inputs/pions_30GeV.txt"
    makePlots(elefile, pionfile, suffix="_30GeV")

print("Done")


