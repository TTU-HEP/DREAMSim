"""
makeHitsZ.py — energy deposition and Cherenkov photon z-profile plots

Usage examples
--------------
# Run from simulation files (default)
python makeHitsZ.py --ele inputs/electrons_40GeV.txt --pion inputs/pions_40GeV.txt \
                    --energy 40 --suffix 40GeV

# Add neutrons
python makeHitsZ.py --ele inputs/electrons_40GeV.txt --pion inputs/pions_40GeV.txt \
                    --neu inputs/neutrons.txt --energy 40 --suffix 40GeV

# Replot from a previously saved ROOT file (fast, no re-reading simulation)
python makeHitsZ.py --from-root hitsZ_40GeV.root --energy 40 --suffix 40GeV_replot
"""

import argparse
import sys
import os
from collections import OrderedDict

import ROOT
sys.path.append("./CMSPLOTS")
from myFunction import DrawHistos

from utils import (
    build_rdfs,
    ensure_dir,
    save_histos_to_root,
    load_histos_from_root,
    load_calibration,
    energy_axis_ranges,
    get_colors,
)

# ---------------------------------------------------------------------------
# Detector / beam constants — update if geometry changes
# ---------------------------------------------------------------------------
DETECTOR_HALF_Z  = 125.0   # cm  (half of 250 cm detector length)
BEAM_X           =   2.5   # cm
BEAM_Y           =  -2.5   # cm
CENTER_ROD_MIN   =  35
CENTER_ROD_MAX   =  51
CENTER_LAYER_MIN =  40
CENTER_LAYER_MAX =  52

# Histogram names that are saved to / loaded from the ROOT cache file
SAVE_HNAMES = ["truthhit_z", "truthhit_z_Cer", "truthhit_z_Cer_center"]


# ---------------------------------------------------------------------------
# Step 1 — build RDataFrames and add derived columns
# ---------------------------------------------------------------------------

def build_rdfs_with_defines(particle_files, nEvts=None, max_files=100):
    """
    Read input files, build RDataFrames, and attach all derived columns.
    Returns (rdfs, chains, nEvts).
    """
    rdfs, chains = build_rdfs(particle_files, max_files=max_files)

    # event counts (needed for per-event weights)
    nEvts = {p: float(rdfs[p].Count().GetValue()) for p in rdfs}
    for p, n in nEvts.items():
        print(f"  {p}: {int(n)} events")

    for part in rdfs:
        n = nEvts[part]
        rdfs[part] = (
            rdfs[part]
            .Define("eTotaltruth",
                    "eLeaktruth + eCalotruth + eWorldtruth + eInvisible")
            .Define("eTotalGeant",
                    "eLeaktruth + eCalotruth + eWorldtruth")
            .Define("truthhit_r",
                    f"sqrt((truthhit_x-{BEAM_X})*(truthhit_x-{BEAM_X})"
                    f"+(truthhit_y-({BEAM_Y}))*(truthhit_y-({BEAM_Y})))")
            .Define("nOPsCer_Cer_perGeV", "nOPsCer_Cer / eCentruth")
            .Define("nOPsCer_Sci_perGeV", "nOPsCer_Sci / eScintruth")
            .Define("eScinSummed",
                    "Sum(truthhit_edep * (truthhit_calotype == 2))")
            .Define("eCentrSummed",
                    "Sum(truthhit_edep * (truthhit_calotype == 3))")
            .Define("isCenterHit",
                    f"(truthhit_rodNumber   > {CENTER_ROD_MIN}   &&"
                    f" truthhit_rodNumber   < {CENTER_ROD_MAX}   &&"
                    f" truthhit_layerNumber > {CENTER_LAYER_MIN} &&"
                    f" truthhit_layerNumber < {CENTER_LAYER_MAX})")
            .Define("nCerPho",            "truthhit_ncercap * isCenterHit")
            .Define("eweight",            f"truthhit_edep / {n}")
            .Define("evt_weight",         f"1.0 / {n}")
            .Define("nCerPho_perEvt",     f"truthhit_ncercap / {n}")
            .Define("nCerPho_perEvt_center", f"nCerPho / {n}")
        )
        # define skew and Meridional selections
        rdfs[part] = (
            rdfs[part]
            .Define("is_skew", f"(OP_mom_produced_z > 0 && OP_isCoreC && OP_isCaptured_s && OP_isAttenuated_s == 0) / {n}")
            .Define("is_meridional", f"(OP_mom_produced_z > 0 && OP_isCoreC && OP_isCaptured_m && OP_isAttenuated_m == 0) / {n}")
            .Define("is_center_op", f"OP_productionRod > {CENTER_ROD_MIN} && OP_productionRod < {CENTER_ROD_MAX} && OP_productionLayer > {CENTER_LAYER_MIN} && OP_productionLayer < {CENTER_LAYER_MAX}")
            .Define("is_skew_center", "is_skew * is_center_op")
            .Define("is_meridional_center", "is_meridional * is_center_op")
            .Define("angle_skew", "ROOT::VecOps::cos(OP_captureAngle_s)") 
            .Define("angle_meridional", "ROOT::VecOps::cos(OP_captureAngle_m)")
        )

    return rdfs, chains, nEvts


# ---------------------------------------------------------------------------
# Step 2 — book histograms
# ---------------------------------------------------------------------------

def book_histos(rdfs, book, suffix):
    """
    Book all histograms against the RDataFrames.
    Returns a nested OrderedDict {hname: {label: lazy histogram proxy}}.
    All ROOT histogram names are globally unique via the suffix.
    """
    histos = OrderedDict()

    def h(name):
        histos[name] = OrderedDict()

    for name in [
        "eLeaktruth", "eCalotruth", "eTotaltruth", "eTotalGeant",
        "eRodtruth", "eCentruth", "eScintruth",
        "nOPsCer_Cer", "nOPsCer_Sci", "nOPsCer_Cer_perGeV", "nOPsCer_Sci_perGeV",
        "truthhit_x", "truthhit_y", "truthhit_z", "truthhit_r",
        "truthhit_z_Cer", "truthhit_z_Cer_center",
        "angle_skew", "angle_meridional",
        "time_skew", "time_meridional",
        # 2D
        "truthhit_x_vs_truthhit_y", "truthhit_x_vs_truthhit_z",
        "truthhit_r_vs_truthhit_z",
        "time_vs_truthhit_z", "time_vs_truthhit_r",
        "time_skew_vs_truthhit_z", "time_meridional_vs_truthhit_z",
        "time_skew_vs_z_prof", "time_meridional_vs_z_prof",

    ]:
        h(name)

    z_lo, z_hi = book["z"]
    t_min, t_max = book["time"]

    for part, rdf in rdfs.items():
        tag = f"{suffix}_{part}"   # unique per (run, particle)

        # --- per-event scalars ---
        histos["eLeaktruth"][part] = rdf.Histo1D(
            (f"eLeaktruth_{tag}",    "eLeaktruth",   50, *book["eLeak"]),    "eLeaktruth")
        histos["eCalotruth"][part] = rdf.Histo1D(
            (f"eCalotruth_{tag}",    "eCalotruth",   50, *book["eCalo"]),    "eCalotruth")
        histos["eTotaltruth"][part] = rdf.Histo1D(
            (f"eTotaltruth_{tag}",   "eTotaltruth",  50, *book["eTotal"]),   "eTotaltruth")
        histos["eTotalGeant"][part] = rdf.Histo1D(
            (f"eTotalGeant_{tag}",   "eTotalGeant",  50, *book["eTotalGeant"]), "eTotalGeant")
        histos["eRodtruth"][part] = rdf.Histo1D(
            (f"eRodtruth_{tag}",     "eRodtruth",    50, *book["eRod"]),     "eRodtruth")
        histos["eCentruth"][part] = rdf.Histo1D(
            (f"eCentruth_{tag}",     "eCentruth",    50, *book["eCentruth"]),"eCentruth")
        histos["eScintruth"][part] = rdf.Histo1D(
            (f"eScintruth_{tag}",    "eScintruth",   50, *book["eScintruth"]),"eScintruth")

        # --- energy-weighted hit profiles ---
        histos["truthhit_x"][part] = rdf.Histo1D(
            (f"truthhit_x_{tag}", "truthhit_x", 50, -20, 20), "truthhit_x", "eweight")
        histos["truthhit_y"][part] = rdf.Histo1D(
            (f"truthhit_y_{tag}", "truthhit_y", 50, -20, 20), "truthhit_y", "eweight")
        histos["truthhit_z"][part] = rdf.Histo1D(
            (f"truthhit_z_{tag}", "truthhit_z", 200, z_lo, z_hi), "truthhit_z", "eweight")
        histos["truthhit_r"][part] = rdf.Histo1D(
            (f"truthhit_r_{tag}", "truthhit_r", 50, 0, 30), "truthhit_r", "eweight")
        histos["truthhit_z_Cer"][part] = rdf.Histo1D(
            (f"truthhit_z_Cer_{tag}", "truthhit_z_Cer", 200, z_lo, z_hi),
            "truthhit_z", "nCerPho_perEvt")
        histos["truthhit_z_Cer_center"][part] = rdf.Histo1D(
            (f"truthhit_z_Cer_center_{tag}", "truthhit_z_Cer_center", 200, z_lo, z_hi),
            "truthhit_z", "nCerPho_perEvt_center")
        histos["angle_skew"][part] = rdf.Histo1D(
            (f"angle_skew_{tag}", "angle_skew", 50, 0, 1), "angle_skew", "is_skew_center")
        histos["angle_meridional"][part] = rdf.Histo1D(
            (f"angle_meridional_{tag}", "angle_meridional", 50, 0, 1), "angle_meridional", "is_meridional_center")
        histos["time_skew"][part] = rdf.Histo1D(
            (f"time_skew_{tag}", "time_skew", 50, t_min, t_max), "OP_analyticalArrivalTime_s", "is_skew_center")
        histos["time_meridional"][part] = rdf.Histo1D(
            (f"time_meridional_{tag}", "time_meridional", 50, t_min, t_max), "OP_analyticalArrivalTime_m", "is_meridional_center")

        # --- 2D spatial ---
        histos["truthhit_x_vs_truthhit_y"][part] = rdf.Histo2D(
            (f"truthhit_x_vs_truthhit_y_{tag}", "", 50, -20, 20, 50, -20, 20),
            "truthhit_x", "truthhit_y", "eweight")
        histos["truthhit_x_vs_truthhit_z"][part] = rdf.Histo2D(
            (f"truthhit_x_vs_truthhit_z_{tag}", "", 50, -20, 20, 100, z_lo, z_hi),
            "truthhit_x", "truthhit_z", "eweight")
        histos["truthhit_r_vs_truthhit_z"][part] = rdf.Histo2D(
            (f"truthhit_r_vs_truthhit_z_{tag}", "", 50, 0, 30, 100, z_lo, z_hi),
            "truthhit_r", "truthhit_z", "eweight")
        
        histos["time_skew_vs_truthhit_z"][part] = rdf.Histo2D(
            (f"time_skew_vs_truthhit_z_{tag}", "", 50, 0, t_max, 100, z_lo, z_hi),
            "OP_analyticalArrivalTime_s", "OP_pos_produced_z", "is_skew_center")
        histos["time_meridional_vs_truthhit_z"][part] = rdf.Histo2D(
            (f"time_meridional_vs_truthhit_z_{tag}", "", 50, 0, t_max, 100, z_lo, z_hi),
            "OP_analyticalArrivalTime_m", "OP_pos_produced_z", "is_meridional_center")
        
        # --- profile: mean analytical arrival time vs production z (RMS error bars) ---
        histos["time_skew_vs_z_prof"][part] = rdf.Profile1D(
            ROOT.RDF.TProfile1DModel(
                f"time_skew_vs_z_prof_{tag}", "", 100, z_lo, z_hi, "s"),
            "OP_pos_produced_z", "OP_analyticalArrivalTime_s", "is_skew_center")
        histos["time_meridional_vs_z_prof"][part] = rdf.Profile1D(
            ROOT.RDF.TProfile1DModel(
                f"time_meridional_vs_z_prof_{tag}", "", 100, z_lo, z_hi, "s"),
            "OP_pos_produced_z", "OP_analyticalArrivalTime_m", "is_meridional_center")
        
    return histos


# ---------------------------------------------------------------------------
# Step 3 — draw histograms
# ---------------------------------------------------------------------------

def draw_histos(histos, draw, suffix, outdir, labels, calib=None):
    """
    Produce all output plots.

    Parameters
    ----------
    histos  : nested dict {hname: {label: TH1/TH2}}
    draw    : draw-range dict from energy_axis_ranges()
    suffix  : string appended to every output filename
    outdir  : directory for output PNG/PDF files
    labels  : ordered list of particle labels present in histos
    calib   : optional calibration dict (from load_calibration)
    """
    ensure_dir(outdir)
    colors = get_colors(labels)
    z_lo, z_hi = draw["z"]
    t_min, t_max       = draw["time"]

    def vals(n): return [histos[n][l] for l in labels if l in histos[n]]
    def keys(n): return [l           for l in labels if l in histos[n]]

    args1d = dict(dology=True, donormalize=True, mycolors=colors,
                  MCOnly=True, addOverflow=True, addUnderflow=True,
                  outdir=outdir)
    
    def d1(hname, xlo, xhi, xlabel, outname):
        DrawHistos(vals(hname), keys(hname),
                   xlo, xhi, xlabel, 1e-3, 1e2, "Fraction of events",
                   outname + "_" + suffix, **args1d)

    d1("eLeaktruth",        *draw["eLeak"],            "Leakage Energy [GeV]",        "eLeaktruth")
    d1("eCalotruth",        *draw["eCalo"],             "Calo Energy [GeV]",           "eCalotruth")
    d1("eTotaltruth",       *draw["eTotal"],            "Total Energy [GeV]",          "eTotaltruth")
    d1("eTotalGeant",       *draw["eTotalGeant"],       "Total Visible Energy [GeV]",  "eTotalGeant")
    d1("eRodtruth",         *draw["eRod"],              "Rod Energy [GeV]",            "eRodtruth")
    d1("eCentruth",         *draw["eCentruth"],         "C-Fiber Energy [GeV]",        "eCentruth")
    d1("eScintruth",        *draw["eScintruth"],        "S-Fiber Energy [GeV]",        "eScintruth")
    d1("nOPsCer_Cer",       *draw["nOPsCer_Cer"],       "N Cer OPs (C-fiber)",         "nOPsCer_Cer")
    d1("nOPsCer_Sci",       *draw["nOPsCer_Sci"],       "N Cer OPs (S-fiber)",         "nOPsCer_Sci")
    d1("nOPsCer_Cer_perGeV",*draw["nOPsCer_Cer_perGeV"],"N Cer OPs / GeV (C-fiber)",  "nOPsCer_Cer_perGeV")
    d1("nOPsCer_Sci_perGeV",*draw["nOPsCer_Sci_perGeV"],"N Cer OPs / GeV (S-fiber)",  "nOPsCer_Sci_perGeV")

    args_edep = {**args1d, "donormalize": False}
    
    args_op = {**args1d, "donormalize": False, "drawoptions": ["hist,C"]*5}
    args_op = {**args1d, "donormalize": False}

    DrawHistos(vals("truthhit_x"), keys("truthhit_x"),
               -20, 20, "x [cm]", 1e-3, 1e2, "Deposited Energy [GeV]",
               f"truthhit_x_{suffix}", **args_edep)
    DrawHistos(vals("truthhit_y"), keys("truthhit_y"),
               -20, 20, "y [cm]", 1e-3, 1e2, "Deposited Energy [GeV]",
               f"truthhit_y_{suffix}", **args_edep)
    DrawHistos(vals("truthhit_z"), keys("truthhit_z"),
               z_lo, z_hi, "z [cm]", 1e-3, 1, "Deposited Energy [GeV]",
               f"truthhit_z_{suffix}", **args_edep)
    DrawHistos(vals("truthhit_r"), keys("truthhit_r"),
               0, 30, "r [cm]", 1e-3, 1e2, "Deposited Energy [GeV]",
               f"truthhit_r_{suffix}", **args_edep)
    DrawHistos(vals("truthhit_z_Cer"), keys("truthhit_z_Cer"),
               z_lo, z_hi, "z [cm]", 1e0, 5e3, "# C Photons",
               f"truthhit_z_Cer_{suffix}", **args_edep)
    DrawHistos(vals("truthhit_z_Cer_center"), keys("truthhit_z_Cer_center"),
               z_lo, z_hi, "z [cm]", 1e0, 5e3, "# C Photons (Central)",
               f"truthhit_z_Cer_center_{suffix}", **args_edep)
    DrawHistos(vals("angle_skew"), keys("angle_skew"),
               0, 1, "Cosine of Capture Angle (Skew)", 1e0, 1e6, "Counts",
               f"angle_skew_{suffix}", **args_edep)
    DrawHistos(vals("angle_meridional"), keys("angle_meridional"),
               0, 1, "Cosine of Capture Angle (Meridional)", 1e0, 1e6, "Counts",
               f"angle_meridional_{suffix}", **args_edep)
    DrawHistos(vals("time_skew"), keys("time_skew"),
               t_min, t_max, "Analytical Arrival Time (Skew) [ns]", 1, 1e4, "Counts",
               f"time_skew_{suffix}", **args_op)
    DrawHistos(vals("time_meridional"), keys("time_meridional"),
               t_min, t_max, "Analytical Arrival Time (Meridional) [ns]", 1, 1e4, "Counts",
               f"time_meridional_{suffix}", **args_op)

    # --- 2D plots ---
    args2d = {**args_edep, "dology": False, "drawoptions": "colz",
              "dologz": True, "zmax": 1e2, "zmin": 1e-4,
              "doth2": True, "addOverflow": False, "addUnderflow": False}
    
    args2d_op = args2d.copy()
    args2d_op["zmax"] = 2e3
    args2d_op["zmin"] = 1e0

    for label in labels:
        if label not in histos["truthhit_x_vs_truthhit_y"]:
            continue
        DrawHistos([histos["truthhit_x_vs_truthhit_y"][label]], [],
                   -20, 20, "x [cm]", -20, 20, "y [cm]",
                   f"truthhit_x_vs_y_{label}_{suffix}", **args2d)
        DrawHistos([histos["truthhit_x_vs_truthhit_z"][label]], [],
                   -20, 20, "x [cm]", z_lo, z_hi, "z [cm]",
                   f"truthhit_x_vs_z_{label}_{suffix}", **args2d)
        DrawHistos([histos["truthhit_r_vs_truthhit_z"][label]], [],
                   0, 30, "r [cm]", z_lo, z_hi, "z [cm]",
                   f"truthhit_r_vs_z_{label}_{suffix}", **args2d)
        DrawHistos([histos["time_skew_vs_truthhit_z"][label]], [],
                   t_min, t_max, "Time [ns] (skew)", z_lo, z_hi, "z [cm]",
                   f"time_skew_vs_z_{label}_{suffix}", **args2d_op)
        DrawHistos([histos["time_meridional_vs_truthhit_z"][label]], [],
                   t_min, t_max, "Time [ns] (Meridional)", z_lo, z_hi, "z [cm]",
                   f"time_meridional_vs_z_{label}_{suffix}", **args2d_op)
        
    # ── Profile overlays: mean arrival time vs z, RMS error bars ─────────────
    args_prof = dict(
        dology=False, donormalize=False, MCOnly=True,
        addOverflow=False, addUnderflow=False,
    )
    markers_solid = [20, 21, 22, 23]
    markers_open  = [24, 25, 26, 27]

    prof_skew = vals("time_skew_vs_z_prof")
    prof_meri = vals("time_meridional_vs_z_prof")
    n = len(labels)

    # skew only: all particles overlaid
    if prof_skew:
        DrawHistos(
            prof_skew, keys("time_skew_vs_z_prof"),
            z_lo, z_hi, "Production z [cm]",
            t_min, t_max + 5, "Mean arrival time [ns] (skew)",
            f"time_skew_vs_z_prof_{suffix}",
            mycolors=get_colors(keys("time_skew_vs_z_prof")),
            drawoptions=["E1"] * len(prof_skew),
            markerstyles=markers_solid[:len(prof_skew)],
            outdir=outdir, **args_prof)

    # meridional only: all particles overlaid
    if prof_meri:
        DrawHistos(
            prof_meri, keys("time_meridional_vs_z_prof"),
            z_lo, z_hi, "Production z [cm]",
            t_min, t_max, "Mean arrival time [ns] (meridional)",
            f"time_meridional_vs_z_prof_{suffix}",
            mycolors=get_colors(keys("time_meridional_vs_z_prof")),
            drawoptions=["E1"] * len(prof_meri),
            markerstyles=markers_solid[:len(prof_meri)],
            outdir=outdir, **args_prof)



# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def run(particle_files, suffix, beam_energy_gev,
        from_root=None, calib_json=None, outdir=None, max_files=100):
    """
    Parameters
    ----------
    particle_files : dict {label: file_list_path}
    suffix         : str, appended to all output names
    beam_energy_gev: float
    from_root      : str or None — if set, reload histos from this ROOT file
                     instead of re-reading simulation
    calib_json     : str or None — path to JSON calibration sidecar
    outdir         : str or None — output directory (default: plots_<suffix>)
    """
    if outdir is None:
        outdir = f"plots_{suffix}"
    ensure_dir(outdir)

    calib  = load_calibration(calib_json) if calib_json else {}
    labels = list(particle_files.keys())
    book, draw = energy_axis_ranges(beam_energy_gev, DETECTOR_HALF_Z)

    root_cache = os.path.join(outdir, f"hitsZ_{suffix}.root")

    if from_root:
        # ----- fast replot mode -----
        print(f"\nLoading histograms from {from_root}")
        all_hnames = [
            "eLeaktruth", "eCalotruth", "eTotaltruth", "eTotalGeant",
            "eRodtruth", "eCentruth", "eScintruth",
            "nOPsCer_Cer", "nOPsCer_Sci", "nOPsCer_Cer_perGeV", "nOPsCer_Sci_perGeV",
            "truthhit_x", "truthhit_y", "truthhit_z", "truthhit_r",
            "truthhit_z_Cer", "truthhit_z_Cer_center",
            "angle_skew", "angle_meridional",
            "time_skew", "time_meridional",
            "truthhit_x_vs_truthhit_y", "truthhit_x_vs_truthhit_z",
            "truthhit_r_vs_truthhit_z",
            "time_skew_vs_truthhit_z", "time_meridional_vs_truthhit_z"
        ]
        histos = load_histos_from_root(from_root, all_hnames, labels)
    else:
        # ----- full run mode -----
        ROOT.ROOT.EnableImplicitMT(10)
        print("\nBuilding RDataFrames...")
        rdfs, chains, nEvts = build_rdfs_with_defines(particle_files, max_files=max_files)

        print("\nBooking histograms...")
        histos = book_histos(rdfs, book, suffix)

        print("\nSaving histograms to ROOT cache...")
        save_histos_to_root(histos, root_cache)

    print("\nDrawing plots...")
    draw_histos(histos, draw, suffix, outdir, labels, calib=calib)
    print(f"\nDone. Plots in: {outdir}/")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Make energy deposition and Cherenkov z-profile plots.")
    p.add_argument("--ele",    help="Text file listing electron ROOT files")
    p.add_argument("--pion",   help="Text file listing pion ROOT files")
    p.add_argument("--neu",    help="Text file listing neutron ROOT files (optional)")
    p.add_argument("--mu",     help="Text file listing muon ROOT files (optional)")
    p.add_argument("--energy", type=float, default=100.0,
                   help="Nominal beam energy in GeV (default: 100)")
    p.add_argument("--suffix", default="",
                   help="Label appended to all output file names")
    p.add_argument("--outdir", default=None,
                   help="Output directory (default: plots_<suffix>)")
    p.add_argument("--from-root", default=None, metavar="FILE",
                   help="Skip simulation reading; reload histos from this ROOT file")
    p.add_argument("--calib",  default=None, metavar="FILE",
                   help="JSON calibration sidecar file (optional)")
    p.add_argument("--max-files", type=int, default=100, metavar="N",
                   help="Max ROOT files to read per particle type (default: 100; use small N for quick tests)")
    return p.parse_args()


if __name__ == "__main__":
    import ROOT
    ROOT.gROOT.SetBatch(True)

    args = parse_args()

    # Build particle_files dict from whatever was provided
    particle_files = OrderedDict()
    for label, path in [("ele", args.ele), ("pion", args.pion),
                        ("neu", args.neu), ("mu", args.mu)]:
        if path:
            particle_files[label] = path

    if not particle_files and not args.from_root:
        print("ERROR: provide at least one of --ele / --pion / --neu / --mu, "
              "or use --from-root.")
        sys.exit(1)

    # If replotting from ROOT we still need labels — infer from what was given,
    # or fall back to ele+pion as the historical default.
    if args.from_root and not particle_files:
        particle_files = OrderedDict([("ele", ""), ("pion", "")])

    run(particle_files,
        suffix         = args.suffix or f"{int(args.energy)}GeV",
        beam_energy_gev= args.energy,
        from_root      = args.from_root,
        calib_json     = args.calib,
        outdir         = args.outdir,
        max_files      = args.max_files)

    print("Done")