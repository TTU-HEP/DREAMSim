import sys
from collections import OrderedDict
import ROOT
sys.path.append("./CMSPLOTS")
from myFunction import DrawHistos


print("Starting")

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(4)

elefile = "inputs/electrons.txt"
pionfile = "inputs/pions.txt"
# neufile = "inputs/neutrons.txt"

chains = OrderedDict()
chains['ele'] = ROOT.TChain("tree")
chains['pion'] = ROOT.TChain("tree")
# chains['neu'] = ROOT.TChain("tree")

# loops = [('ele', elefile), ('pion', pionfile), ('neu', neufile)]
loops = [('ele', elefile), ('pion', pionfile)]
#loops = [('pion', pionfile)]

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

nEvts = OrderedDict()
nEvts['ele'] = rdfs['ele'].Count().GetValue()
nEvts['pion'] = rdfs['pion'].Count().GetValue()
# nEvts['neu'] = rdfs['neu'].Count().GetValue()
print("Number of events for electrons: ", nEvts['ele'])
print("Number of events for pions: ", nEvts['pion'])
# print("Number of events for neutrons: ", nEvts['neu'])

rdfs['ele'] = rdfs['ele'].Define("eventweight", f"1.0 / {nEvts['ele']}")
rdfs['pion'] = rdfs['pion'].Define("eventweight", f"1.0 / {nEvts['pion']}")

rdfs['ele'] = rdfs['ele'].Define("eweight", f"truthhit_edep/ {nEvts['ele']}")
rdfs['pion'] = rdfs['pion'].Define(
    "eweight", f"truthhit_edep/ {nEvts['pion']}")
# rdfs['neu'] = rdfs['neu'].Define("eweight", f"truthhit_edep/ {nEvts['neu']}")


rdfs_new = OrderedDict()
for part in rdfs.keys():
    rdfs_new[part] = rdfs[part].Define("eTotaltruth", "eLeaktruth + eCalotruth + eWorldtruth + eInvisible") \
        .Define("eTotalGeant", "eLeaktruth + eCalotruth + eWorldtruth") \
        .Define("truthhit_r", f"sqrt((truthhit_x-2.5)*(truthhit_x-2.5) + (truthhit_y+2.5)*(truthhit_y+2.5))")

    # optics features
    rdfs_new[part] = rdfs_new[part].Define("OP_passEnd", "OP_pos_final_z > 99.0") \
        .Define("OP_passEnd_isCoreC", "OP_passEnd && OP_isCoreC") \
        .Define("OP_passEnd_isCoreS", "OP_passEnd && OP_isCoreS") \
        .Define("OP_passEnd_normalized", "OP_passEnd * eventweight") \
        .Define("OP_passEnd_isCoreC_normalized", "OP_passEnd_isCoreC * eventweight") \
        .Define("OP_passEnd_isCoreS_normalized", "OP_passEnd_isCoreS * eventweight") \
        .Define("OP_isCoreC_normalized", "OP_isCoreC * eventweight") \
        .Define("OP_isCoreS_normalized", "OP_isCoreS * eventweight")
        
    rdfs_new[part] = rdfs_new[part].Define("OP_energy_final", "sqrt(OP_mom_final_z*OP_mom_final_z+OP_mom_final_y*OP_mom_final_y+OP_mom_final_x*OP_mom_final_x)") \
        .Define("OP_wavelength_final", "1239.8e-9 / OP_energy_final") \
        .Define("OP_energy_produced", "sqrt(OP_mom_produced_z*OP_mom_produced_z+OP_mom_produced_y*OP_mom_produced_y+OP_mom_produced_x*OP_mom_produced_x)") \
        .Define("OP_wavelength_produced", "1239.8e-9 / OP_energy_produced")
        
    rdfs_new[part] = rdfs_new[part].Define("nOPs_produced", "Sum(OP_time_produced > 0.)") \
        .Define("nOPs_produced_isCoreC", "Sum(OP_time_produced > 0. && OP_isCoreC)") \
        .Define("nOPs_produced_isCoreS", "Sum(OP_time_produced > 0. && OP_isCoreS)") \
        .Define("nOPs_passEnd", "Sum(OP_passEnd)") \
        .Define("nOPs_passEnd_isCoreC", "Sum(OP_passEnd_isCoreC)") \
        .Define("nOPs_passEnd_isCoreS", "Sum(OP_passEnd_isCoreS)") \
        .Define("capEff", "float(nOPs_passEnd) / nOPs_produced") \
        .Define("capEff_isCoreC", "float(nOPs_passEnd_isCoreC) / nOPs_produced_isCoreC") \
        .Define("capEff_isCoreS", "float(nOPs_passEnd_isCoreS) / nOPs_produced_isCoreS")
        
    # energy deposit in the "center of the calorimeter"
    rdfs_new[part] = rdfs_new[part].Define("eDep_center", "Sum(truthhit_edep * (truthhit_layerNumber == 40 && truthhit_rodNumber == 45))")
            
rdfs = rdfs_new

histos = OrderedDict()
figures = ['eLeaktruth', 'eCalotruth', 'eTotaltruth', 'eTotalGeant',
           'eRodtruth', 'eCentruth', 'eScintruth',
           'truthhit_x', 'truthhit_y', 'truthhit_z', 'truthhit_r',
           'time', 'time_zoomed',
           'truthhit_x_vs_truthhit_y', 'truthhit_x_vs_truthhit_z', 'truthhit_r_vs_truthhit_z',
           'time_vs_truthhit_z', 'time_vs_truthhit_r',
           "OP_time_produced", "OP_time_final", "OP_time_final_vs_OP_pos_produced_z",
           "eDep_center", "eDep_center_vs_nOPs_produced", "eDep_center_vs_nOPs_passEnd",
           "nOPs_produced", "nOPs_passEnd", "capEff",
           "OP_wavelength_produced", "OP_wavelength_final"]


evtlist = [1, 3, 5, 10, 15]
for i in evtlist:
    figures.append(f"event_{i}_truthhit_x_vs_truthhit_y")
    figures.append(f"event_{i}_truthhit_x_vs_truthhit_z")
    figures.append(f"event_{i}_truthhit_r_vs_truthhit_z")
    figures.append(f"event_{i}_time_vs_truthhit_z")
    figures.append(f"event_{i}_time_vs_truthhit_r")

for fig in figures:
    histos[fig] = OrderedDict()

for part, rdf in rdfs.items():
    suffix = "_" + part
    
    histos['nOPs_produced'][part] = rdf.Histo1D(
        ("nOPs_produced" + suffix, "nOPs_produced", 50, 0, 6e4), "nOPs_produced")
    histos['nOPs_passEnd'][part] = rdf.Histo1D(
        ("nOPs_passEnd" + suffix, "nOPs_passEnd", 50, 0, 2e3), "nOPs_passEnd")
    histos['capEff'][part] = rdf.Histo1D(
        ("capEff" + suffix, "capEff", 50, 0, 0.1), "capEff")
    histos['eDep_center'][part] = rdf.Histo1D(
        ("eDep_center" + suffix, "eDep_center", 50, 0, 0.5), "eDep_center")
    histos['eDep_center_vs_nOPs_produced'][part] = rdf.Histo2D(
        ("eDep_center_vs_nOPs_produced" + suffix, "eDep_center_vs_nOPs_produced", 50, 0, 0.5, 50, 0, 6e4), "eDep_center", "nOPs_produced")
    histos['eDep_center_vs_nOPs_passEnd'][part] = rdf.Histo2D(
        ("eDep_center_vs_nOPs_passEnd" + suffix, "eDep_center_vs_nOPs_passEnd", 50, 0, 0.5, 50, 0, 2e3), "eDep_center", "nOPs_passEnd")

    histos['eLeaktruth'][part] = rdf.Histo1D(
        ("eLeaktruth" + suffix, "eLeaktruth", 50, 0, 20.0), "eLeaktruth")
    histos['eCalotruth'][part] = rdf.Histo1D(
        ("eCalotruth" + suffix, "eCalotruth", 50, 0., 102.0), "eCalotruth")
    histos['eTotaltruth'][part] = rdf.Histo1D(
        ("eTotaltruth" + suffix, "eTotaltruth", 50, 90., 110.0), "eTotaltruth")
    histos['eTotalGeant'][part] = rdf.Histo1D(
        ("eTotalGeant" + suffix, "eTotalGeant", 50, 80., 110.0), "eTotalGeant")
    histos['eRodtruth'][part] = rdf.Histo1D(
        ("eRodtruth" + suffix, "eRodtruth", 50, 50, 102.0), "eRodtruth")
    histos['eCentruth'][part] = rdf.Histo1D(
        ("eCentruth" + suffix, "eCentruth", 50, 1.0, 5.0), "eCentruth")
    histos['eScintruth'][part] = rdf.Histo1D(
        ("eScintruth" + suffix, "eScintruth", 50, 1.0, 5.0), "eScintruth")

    # energy weighted
    histos['truthhit_x'][part] = rdf.Histo1D(
        ("truthhit_x" + suffix, "truthhit_x", 50, -20, 20), "truthhit_x", "eweight")
    histos['truthhit_y'][part] = rdf.Histo1D(
        ("truthhit_y" + suffix, "truthhit_y", 50, -20, 20), "truthhit_y", "eweight")
    histos['truthhit_z'][part] = rdf.Histo1D(
        ("truthhit_z" + suffix, "truthhit_z", 100, -100, 100), "truthhit_z", "eweight")
    histos['truthhit_r'][part] = rdf.Histo1D(
        ("truthhit_r" + suffix, "truthhit_r", 50, 0, 30), "truthhit_r", "eweight")

    histos['time'][part] = rdf.Histo1D(
        ("time" + suffix, "time", 50, 0, 50), "truthhit_globaltime", "eweight")
    histos['time_zoomed'][part] = rdf.Histo1D(
        ("time_zoomed" + suffix, "time_zoomed", 50, 0, 10), "truthhit_globaltime", "eweight")
    
    histos['truthhit_x_vs_truthhit_y'][part] = rdf.Histo2D(
        ("truthhit_x_vs_truthhit_y" + suffix, "truthhit_x_vs_truthhit_y", 50, -20, 20, 50, -20, 20), "truthhit_x", "truthhit_y", "eweight")
    histos['truthhit_x_vs_truthhit_z'][part] = rdf.Histo2D(
        ("truthhit_x_vs_truthhit_z" + suffix, "truthhit_x_vs_truthhit_z", 50, -20, 20, 100, -100, 100), "truthhit_x", "truthhit_z", "eweight")
    histos['truthhit_r_vs_truthhit_z'][part] = rdf.Histo2D(
        ("truthhit_r_vs_truthhit_z" + suffix, "truthhit_r_vs_truthhit_z", 50, 0, 30, 100, -100, 100), "truthhit_r", "truthhit_z", "eweight")

    histos['truthhit_r_vs_truthhit_z'][part] = rdf.Histo2D(
        ("truthhit_r_vs_truthhit_z" + suffix, "truthhit_r_vs_truthhit_z", 50, 0, 30, 100, -100, 100), "truthhit_r", "truthhit_z", "eweight")

    histos['time_vs_truthhit_z'][part] = rdf.Histo2D(
        ("time_vs_truthhit_z" + suffix, "time_vs_truthhit_z", 50, 0, 20, 100, -100, 100), "truthhit_globaltime", "truthhit_z", "eweight")
    histos['time_vs_truthhit_r'][part] = rdf.Histo2D(
        ("time_vs_truthhit_r" + suffix, "time_vs_truthhit_r", 50, 0, 20, 50, 0, 30), "truthhit_globaltime", "truthhit_r", "eweight")

    # some event displays
    for i in evtlist:
        rdf_event = rdf.Filter(f"rdfentry_ == {i}")
        histos[f"event_{i}_truthhit_x_vs_truthhit_y"][part] = rdf_event.Histo2D(
            (f"event_{i}_truthhit_x_vs_truthhit_y" + suffix, f"event_{i}_truthhit_x_vs_truthhit_y", 50, -20, 20, 50, -20, 20), "truthhit_x", "truthhit_y", "truthhit_edep")
        histos[f"event_{i}_truthhit_x_vs_truthhit_z"][part] = rdf_event.Histo2D(
            (f"event_{i}_truthhit_x_vs_truthhit_z" + suffix, f"event_{i}_truthhit_x_vs_truthhit_z", 50, -20, 20, 100, -100, 100), "truthhit_x", "truthhit_z", "truthhit_edep")
        histos[f"event_{i}_truthhit_r_vs_truthhit_z"][part] = rdf_event.Histo2D(
            (f"event_{i}_truthhit_r_vs_truthhit_z" + suffix, f"event_{i}_truthhit_r_vs_truthhit_z", 50, 0, 30, 100, -100, 100), "truthhit_r", "truthhit_z", "truthhit_edep")

        histos[f"event_{i}_time_vs_truthhit_z"][part] = rdf_event.Histo2D(
            (f"event_{i}_time_vs_truthhit_z" + suffix, f"event_{i}_time_vs_truthhit_z", 50, 0, 20, 100, -100, 100), "truthhit_globaltime", "truthhit_z", "truthhit_edep")
        histos[f"event_{i}_time_vs_truthhit_r"][part] = rdf_event.Histo2D(
            (f"event_{i}_time_vs_truthhit_r" + suffix, f"event_{i}_time_vs_truthhit_r", 50, 0, 20, 50, 0, 30), "truthhit_globaltime", "truthhit_r", "truthhit_edep")

    for fib in ["isCoreC", "isCoreS"]:

        suffix = "_" + part + "_" + fib
        oppart = part + "_" + fib

        histos['OP_time_produced'][oppart] = rdf.Histo1D(
            ("OP_time_produced" + suffix, "OP_time_produced", 50, 0, 5), "OP_time_produced", f"OP_passEnd_{fib}_normalized")
        histos['OP_time_final'][oppart] = rdf.Histo1D(
            ("OP_time_final" + suffix, "OP_time_final", 50, 5, 17), "OP_time_final", f"OP_passEnd_{fib}_normalized")
        histos["OP_time_final_vs_OP_pos_produced_z"][oppart] = rdf.Histo2D(
            ("OP_time_final_vs_OP_pos_produced_z" + suffix, "OP_time_final_vs_OP_pos_produced_z", 50, 6, 17, 50, -100, 100), "OP_time_final", "OP_pos_produced_z", f"OP_passEnd_{fib}_normalized")
        
        histos['OP_wavelength_produced'][oppart] = rdf.Histo1D(
            ("OP_wavelength_produced" + suffix, "OP_wavelength_produced", 100, 300, 700), "OP_wavelength_produced", f"OP_{fib}_normalized")
        histos['OP_wavelength_final'][oppart] = rdf.Histo1D(
            ("OP_wavelength_final" + suffix, "OP_wavelength_final", 100, 300, 700), "OP_wavelength_final", f"OP_passEnd_{fib}_normalized")


colormaps = {
    'ele_isCoreS': 2,
    'ele_isCoreC': 6,
    'pion_isCoreS': 3,
    'pion_isCoreC': 4,
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
    'mycolors': [colormaps[part] for part in histos["OP_time_produced"].keys()],
    "MCOnly": True,
    'addOverflow': True,
    'addUnderflow': True
}

print("Drawing")

DrawHistos(list(histos['eLeaktruth'].values()), list(histos['eLeaktruth'].keys(
)), 0, 20, "Leakage Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eLeaktruth", **args)
DrawHistos(list(histos['eCalotruth'].values()), list(histos['eCalotruth'].keys(
)), 0., 102, "Calo Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eCalotruth", **args)
DrawHistos(list(histos['eTotaltruth'].values()), list(histos['eTotaltruth'].keys(
)), 80., 110, "Total Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eTotaltruth", **args)
DrawHistos(list(histos['eTotalGeant'].values()), list(histos['eTotalGeant'].keys(
)), 80., 110, "Total 'Visible' Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eTotalGeant", **args)
DrawHistos(list(histos['eRodtruth'].values()), list(histos['eRodtruth'].keys(
)), 50, 102, "Rod Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eRodtruth", **args)
DrawHistos(list(histos['eCentruth'].values()), list(histos['eCentruth'].keys(
)), 1, 5, "CFiber Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eCentruth", **args)
DrawHistos(list(histos['eScintruth'].values()), list(histos['eScintruth'].keys(
)), 1, 5, "SFiber Energy [GeV]", 1e-3, 1e2, "Fraction of events", "eScintruth", **args)

DrawHistos(list(histos['nOPs_produced'].values()), list(histos['nOPs_produced'].keys(
)), 0, 6e4, "Number of OPs produced", 1e-3, 1e2, "Fraction of events", "nOPs_produced", **args)
DrawHistos(list(histos['nOPs_passEnd'].values()), list(histos['nOPs_passEnd'].keys(
)), 0, 2e3, "Number of OPs reaching the end", 1e-3, 1e2, "Fraction of events", "nOPs_passEnd", **args)
DrawHistos(list(histos['capEff'].values()), list(histos['capEff'].keys(
)), 0, 0.1, "Capture Efficiency", 1e-3, 1e2, "Fraction of events", "capEff", **args)
DrawHistos(list(histos['eDep_center'].values()), list(histos['eDep_center'].keys(
)), 0, 0.5, "Energy deposited in the center of the calorimeter [GeV]", 1e-3, 1e2, "Fraction of events", "eDep_center", **args)


args['donormalize'] = False
DrawHistos(list(histos['truthhit_x'].values()), list(histos['truthhit_x'].keys(
)), -20, 20, "x [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "truthhit_x", **args)
DrawHistos(list(histos['truthhit_y'].values()), list(histos['truthhit_y'].keys(
)), -20, 20, "y [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "truthhit_y", **args)
DrawHistos(list(histos['truthhit_z'].values()), list(histos['truthhit_z'].keys(
)), -100, 100, "z [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "truthhit_z", **args)
DrawHistos(list(histos['truthhit_r'].values()), list(histos['truthhit_r'].keys(
)), 0, 30, "r [cm]", 1e-3, 1e2, "Deposited Energy [GeV]", "truthhit_r", **args)

DrawHistos(list(histos['time'].values()), list(histos['time'].keys(
)), 0, 100, "Time [ns]", 1e-3, 1e2, "Deposited Energy [GeV]", "time", **args)
DrawHistos(list(histos['time_zoomed'].values()), list(histos['time_zoomed'].keys(
)), 0, 10, "Time [ns]", 1e-3, 1e2, "Deposited Energy [GeV]", "time_zoomed", **args)
DrawHistos(list(histos['OP_time_produced'].values()), list(histos['OP_time_produced'].keys(
)), 0, 5, "Time [ns]", 1e1, 1e8, "Deposited Energy [GeV]", "OP_time_produced", **args)
DrawHistos(list(histos['OP_time_final'].values()), list(histos['OP_time_final'].keys(
)), 5, 17, "Time [ns]", 1e1, 1e8, "Deposited Energy [GeV]", "OP_time_final", **args)
DrawHistos(list(histos['OP_wavelength_produced'].values()), list(histos['OP_wavelength_produced'].keys(
)), 300, 700, "Wavelength [nm]", 1e-1, 1e5, "# Photons", "OP_wavelength_produced", **args)
DrawHistos(list(histos['OP_wavelength_final'].values()), list(histos['OP_wavelength_final'].keys(
)), 300, 700, "Wavelength [nm]", 1e-1, 1e5, "# Photons", "OP_wavelength_final", **args)


# 2D plots
args['dology'] = False
args['drawoptions'] = "colz"
args['dologz'] = True
args['zmax'] = 1e0
args['zmin'] = 1e-6
args['doth2'] = True
args['addOverflow'] = False
args['addUnderflow'] = False
for part in rdfs.keys():
    DrawHistos([histos['truthhit_x_vs_truthhit_y'][part]], [], -20, 20,
               "x [cm]", -20, 20, "y [cm]", f"truthhit_x_vs_truthhit_y_{part}", **args)
    DrawHistos([histos['truthhit_x_vs_truthhit_z'][part]], [], -20, 20,
               "x [cm]", -100, 100, "z [cm]", f"truthhit_x_vs_truthhit_z_{part}", **args)
    DrawHistos([histos['truthhit_r_vs_truthhit_z'][part]], [], 0, 30,
               "r [cm]", -100, 100, "z [cm]", f"truthhit_r_vs_truthhit_z_{part}", **args)

    DrawHistos([histos['time_vs_truthhit_z'][part]], [], 0, 20,
               "Global time [ns]", -100, 100, "Production z [cm]", f"time_vs_truthhit_z_{part}", **args)
    DrawHistos([histos['time_vs_truthhit_r'][part]], [], 0, 20,
               "Global time [ns]", 0, 30, "Production r [cm]", f"time_vs_truthhit_r_{part}", **args)
    

    # event displays
    for i in evtlist:
        DrawHistos([histos[f"event_{i}_truthhit_x_vs_truthhit_y"][part]], [
        ], -20, 20, "x [cm]", -20, 20, "y [cm]", f"event_{i}_truthhit_x_vs_truthhit_y_{part}", **args)
        DrawHistos([histos[f"event_{i}_truthhit_x_vs_truthhit_z"][part]], [
        ], -20, 20, "x [cm]", -100, 100, "z [cm]", f"event_{i}_truthhit_x_vs_truthhit_z_{part}", **args)
        DrawHistos([histos[f"event_{i}_truthhit_r_vs_truthhit_z"][part]], [
        ], 0, 30, "r [cm]", -100, 100, "z [cm]", f"event_{i}_truthhit_r_vs_truthhit_z_{part}", **args)

        DrawHistos([histos[f"event_{i}_time_vs_truthhit_z"][part]], [
        ], 0, 20, "Time [ns]", -100, 100, "z [cm]", f"event_{i}_time_vs_truthhit_z_{part}", **args)
        DrawHistos([histos[f"event_{i}_time_vs_truthhit_r"][part]], [
        ], 0, 20, "Time [ns]", 0, 30, "r [cm]", f"event_{i}_time_vs_truthhit_r_{part}", **args)
        

args['zmax'] = 1e3
args['zmin'] = 1e-2
for part in rdfs.keys():
    DrawHistos([histos['eDep_center_vs_nOPs_produced'][part]], [], 0, 0.5,
               "E dep in the center hole [GeV]", 0, 6e4, "Number of OPs produced", f"eDep_center_vs_nOPs_produced_{part}", **args)
    DrawHistos([histos['eDep_center_vs_nOPs_passEnd'][part]], [], 0, 0.5,
                "E dep in the center hole [GeV]", 0, 2e3, "Number of OPs reaching the end", f"eDep_center_vs_nOPs_passEnd_{part}", **args)

args['zmax'] = 1e2
args['zmin'] = 1e-4
for oppart in histos["OP_time_produced"].keys():
    DrawHistos([histos["OP_time_final_vs_OP_pos_produced_z"][oppart]], [], 7, 14,
               "Arrival time [ns]", -100, 100, "Production z [cm]", f"OP_time_final_vs_OP_pos_produced_z_{oppart}", **args)
    



print("Done")
