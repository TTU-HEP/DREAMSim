"""
makeHitsZ.py — energy deposition and Cherenkov photon z-profile plots

Produces plots for scintillation (S-fiber) and two Cherenkov fiber types
(C-Plastic, C-Quartz) overlaid on shared canvases where relevant.

Output plots
------------
  - Energy depositions: eCalo, eTotal, eRod, ePlatruth, eQuatruth, eScintruth, ...
  - Hit distributions vs x / y / z / r
  - Cherenkov photon z-profiles: Plastic and Quartz overlaid
  - Capture-angle distributions: skew and meridional modes
  - Analytical arrival-time distributions and 2D (time vs z) histograms
  - Profile plots (mean arrival time vs z) with pol1 fits and fit-parameter
    text boxes; fit ranges: ele [-100, -50] cm, pion [-100, 10] cm

Color coding
------------
  Color encodes particle type (by label index):
    index 0 — Plastic: kRed,    Quartz: kPink+6
    index 1 — Plastic: kBlue,   Quartz: kAzure+6
    index 2 — Plastic: kGreen+2, Quartz: kTeal+3
  Legend height is auto-sized to 0.05 * n_entries;
  two-column layout is used automatically when entries > 6.

Suffix conventions
------------------
  If "rotate" appears in the suffix, the meridional time axis is widened
  to [0, 16] ns instead of the default [8, 15] ns.

Usage examples
--------------
# Single particle type
python makeHitsZ.py --ele inputs/electrons_40GeV.txt --energy 40 --suffix ele_40GeV

# Overlay electrons and pions (Plastic/Quartz get distinct color families)
python makeHitsZ.py \\
    --ele  inputs/electrons_40GeV.txt \\
    --pion inputs/pions_40GeV.txt \\
    --energy 40 --suffix elePion_40GeV

# Position scan (each label gets its own color)
python makeHitsZ.py \\
    --input ele_center inputs/electrons_center.txt \\
    --input ele_x5y0   inputs/electrons_x5y0.txt  \\
    --input ele_x10y0  inputs/electrons_x10y0.txt \\
    --energy 40 --suffix pos_scan_40GeV

# Rotated geometry (widens meridional time range)
python makeHitsZ.py --ele inputs/electrons_rotated.txt --energy 40 --suffix ele_40GeV_rotate

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
    lpos,
    make_fit_pave,
    PLA_COLORS, QUA_COLORS,
    PLA_MARKERS, QUA_MARKERS,
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

hnames =  [
        "eLeaktruth", "eCalotruth", "eTotaltruth", "eTotalGeant",
        "eRodtruth", "ePlatruth", "eQuatruth", "eScintruth",
        "nOPsCer_Pla", "nOPsCer_Qua", "nOPsCer_Sci",
        "nOPsCer_Pla_perGeV", "nOPsCer_Qua_perGeV", "nOPsCer_Sci_perGeV",
        "truthhit_x", "truthhit_y", "truthhit_z", "truthhit_r",
        "truthhit_z_Pla", "truthhit_z_Pla_center",
        "truthhit_z_Qua", "truthhit_z_Qua_center",
        "angle_skew_Pla", "angle_skew_Qua",
        "angle_meridional_Pla", "angle_meridional_Qua",
        "time_skew_Pla", "time_skew_Qua",
        "time_meridional_Pla", "time_meridional_Qua",
        # 2D
        "truthhit_x_vs_truthhit_y", "truthhit_x_vs_truthhit_z", "truthhit_y_vs_truthhit_z",
        "truthhit_r_vs_truthhit_z",
        "time_vs_truthhit_z", "time_vs_truthhit_r",
        "time_skew_vs_truthhit_z_Pla", "time_skew_vs_truthhit_z_Qua",
        "time_meridional_vs_truthhit_z_Pla", "time_meridional_vs_truthhit_z_Qua",
        "time_skew_vs_z_prof_Pla", "time_skew_vs_z_prof_Qua",
        "time_meridional_vs_z_prof_Pla", "time_meridional_vs_z_prof_Qua",
    ]

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
            .Define("nOPsCer_Pla_perGeV", "nOPsCer_Pla / ePlatruth")
            .Define("nOPsCer_Qua_perGeV", "nOPsCer_Qua / eQuatruth")
            .Define("nOPsCer_Sci_perGeV", "nOPsCer_Sci / eScintruth")
            .Define("eScinSummed",
                    "Sum(truthhit_edep * (truthhit_calotype == 2))")
            .Define("eCentrSummed",
                    "Sum(truthhit_edep * (truthhit_calotype == 3))")
            .Define("eQentrSummed",
                    "Sum(truthhit_edep * (truthhit_calotype == 4))")
            .Define("isCenterHit",
                    f"(truthhit_rodNumber   > {CENTER_ROD_MIN}   &&"
                    f" truthhit_rodNumber   < {CENTER_ROD_MAX}   &&"
                    f" truthhit_layerNumber > {CENTER_LAYER_MIN} &&"
                    f" truthhit_layerNumber < {CENTER_LAYER_MAX})")
            .Define("nCerPho_Pla",        "truthhit_ncercap * (truthhit_calotype == 3) * isCenterHit")
            .Define("nCerPho_Qua",        "truthhit_ncercap * (truthhit_calotype == 4) * isCenterHit")
            .Define("eweight",            f"truthhit_edep / {n}")
            .Define("evt_weight",         f"1.0 / {n}")
            .Define("nCerPho_Pla_perEvt", f"truthhit_ncercap * (truthhit_calotype == 3) / {n}")
            .Define("nCerPho_Qua_perEvt", f"truthhit_ncercap * (truthhit_calotype == 4) / {n}")
            .Define("nCerPho_Pla_perEvt_center", f"nCerPho_Pla / {n}")
            .Define("nCerPho_Qua_perEvt_center", f"nCerPho_Qua / {n}")
        )
        # define skew and meridional selections, split by fiber type (Plastic / Quartz)
        rdfs[part] = (
            rdfs[part]
            .Define("is_skew_Pla",      f"(OP_mom_produced_z > 0 && OP_isCoreC && OP_isCaptured_s && OP_isAttenuated_s == 0) / {n}")
            .Define("is_skew_Qua",      f"(OP_mom_produced_z > 0 && OP_isCoreQ && OP_isCaptured_s && OP_isAttenuated_s == 0) / {n}")
            .Define("is_meridional_Pla",f"(OP_mom_produced_z > 0 && OP_isCoreC && OP_isCaptured_m && OP_isAttenuated_m == 0) / {n}")
            .Define("is_meridional_Qua",f"(OP_mom_produced_z > 0 && OP_isCoreQ && OP_isCaptured_m && OP_isAttenuated_m == 0) / {n}")
            .Define("is_center_op",     f"OP_productionRod > {CENTER_ROD_MIN} && OP_productionRod < {CENTER_ROD_MAX} && OP_productionLayer > {CENTER_LAYER_MIN} && OP_productionLayer < {CENTER_LAYER_MAX}")
            .Define("is_skew_Pla_center",       "is_skew_Pla       * is_center_op")
            .Define("is_skew_Qua_center",       "is_skew_Qua       * is_center_op")
            .Define("is_meridional_Pla_center", "is_meridional_Pla * is_center_op")
            .Define("is_meridional_Qua_center", "is_meridional_Qua * is_center_op")
            .Define("angle_skew",       "ROOT::VecOps::cos(OP_captureAngle_s)")
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
        
    for name in hnames:
        h(name)

    z_lo, z_hi = book["z"]
    t_min, t_max = book["time"]

    for part, rdf in rdfs.items():
        tag = f"{suffix}_{part}"   # unique per (run, label)

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
        histos["ePlatruth"][part] = rdf.Histo1D(
            (f"ePlatruth_{tag}",     "ePlatruth",    50, *book["ePlatruth"]),"ePlatruth")
        histos["eQuatruth"][part] = rdf.Histo1D(
            (f"eQuatruth_{tag}",     "eQuatruth",    50, *book["ePlatruth"]),"eQuatruth")
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
        # Cherenkov z-profiles — Plastic and Quartz separately
        histos["truthhit_z_Pla"][part] = rdf.Histo1D(
            (f"truthhit_z_Pla_{tag}", "truthhit_z_Pla", 200, z_lo, z_hi),
            "truthhit_z", "nCerPho_Pla_perEvt")
        histos["truthhit_z_Pla_center"][part] = rdf.Histo1D(
            (f"truthhit_z_Pla_center_{tag}", "truthhit_z_Pla_center", 200, z_lo, z_hi),
            "truthhit_z", "nCerPho_Pla_perEvt_center")
        histos["truthhit_z_Qua"][part] = rdf.Histo1D(
            (f"truthhit_z_Qua_{tag}", "truthhit_z_Qua", 200, z_lo, z_hi),
            "truthhit_z", "nCerPho_Qua_perEvt")
        histos["truthhit_z_Qua_center"][part] = rdf.Histo1D(
            (f"truthhit_z_Qua_center_{tag}", "truthhit_z_Qua_center", 200, z_lo, z_hi),
            "truthhit_z", "nCerPho_Qua_perEvt_center")
        # Capture angle — Plastic and Quartz
        histos["angle_skew_Pla"][part] = rdf.Histo1D(
            (f"angle_skew_Pla_{tag}", "angle_skew_Pla", 50, 0, 1), "angle_skew", "is_skew_Pla_center")
        histos["angle_skew_Qua"][part] = rdf.Histo1D(
            (f"angle_skew_Qua_{tag}", "angle_skew_Qua", 50, 0, 1), "angle_skew", "is_skew_Qua_center")
        histos["angle_meridional_Pla"][part] = rdf.Histo1D(
            (f"angle_meridional_Pla_{tag}", "angle_meridional_Pla", 50, 0, 1), "angle_meridional", "is_meridional_Pla_center")
        histos["angle_meridional_Qua"][part] = rdf.Histo1D(
            (f"angle_meridional_Qua_{tag}", "angle_meridional_Qua", 50, 0, 1), "angle_meridional", "is_meridional_Qua_center")
        # Arrival time — Plastic and Quartz
        histos["time_skew_Pla"][part] = rdf.Histo1D(
            (f"time_skew_Pla_{tag}", "time_skew_Pla", 500, t_min, t_max), "OP_analyticalArrivalTime_s", "is_skew_Pla_center")
        histos["time_skew_Qua"][part] = rdf.Histo1D(
            (f"time_skew_Qua_{tag}", "time_skew_Qua", 500, t_min, t_max), "OP_analyticalArrivalTime_s", "is_skew_Qua_center")
        histos["time_meridional_Pla"][part] = rdf.Histo1D(
            (f"time_meridional_Pla_{tag}", "time_meridional_Pla", 500, t_min, t_max), "OP_analyticalArrivalTime_m", "is_meridional_Pla_center")
        histos["time_meridional_Qua"][part] = rdf.Histo1D(
            (f"time_meridional_Qua_{tag}", "time_meridional_Qua", 500, t_min, t_max), "OP_analyticalArrivalTime_m", "is_meridional_Qua_center")

        # --- 2D spatial ---
        # Convention: a_vs_b → a on y-axis, b on x-axis
        histos["truthhit_x_vs_truthhit_y"][part] = rdf.Histo2D(
            (f"truthhit_x_vs_truthhit_y_{tag}", "", 50, -20, 20, 60, -30, 30),
            "truthhit_y", "truthhit_x", "eweight")
        histos["truthhit_x_vs_truthhit_z"][part] = rdf.Histo2D(
            (f"truthhit_x_vs_truthhit_z_{tag}", "", 100, z_lo, z_hi, 60, -30, 30),
            "truthhit_z", "truthhit_x", "eweight")
        histos["truthhit_y_vs_truthhit_z"][part] = rdf.Histo2D(
            (f"truthhit_y_vs_truthhit_z_{tag}", "", 100, z_lo, z_hi, 60, -30, 30),
            "truthhit_z", "truthhit_y", "eweight")
        histos["truthhit_r_vs_truthhit_z"][part] = rdf.Histo2D(
            (f"truthhit_r_vs_truthhit_z_{tag}", "", 100, z_lo, z_hi, 50, 0, 40),
            "truthhit_z", "truthhit_r", "eweight")

        # 2D time vs z — Plastic and Quartz (time on y-axis, z on x-axis)
        histos["time_skew_vs_truthhit_z_Pla"][part] = rdf.Histo2D(
            (f"time_skew_vs_truthhit_z_Pla_{tag}", "", 100, z_lo, z_hi, 500, 0, t_max),
            "OP_pos_produced_z", "OP_analyticalArrivalTime_s", "is_skew_Pla_center")
        histos["time_skew_vs_truthhit_z_Qua"][part] = rdf.Histo2D(
            (f"time_skew_vs_truthhit_z_Qua_{tag}", "", 100, z_lo, z_hi, 500, 0, t_max),
            "OP_pos_produced_z", "OP_analyticalArrivalTime_s", "is_skew_Qua_center")
        t_mer_lo, t_mer_hi = book["time_meridional"]
        histos["time_meridional_vs_truthhit_z_Pla"][part] = rdf.Histo2D(
            (f"time_meridional_vs_truthhit_z_Pla_{tag}", "", 100, z_lo, z_hi, 500, t_mer_lo, t_mer_hi),
            "OP_pos_produced_z", "OP_analyticalArrivalTime_m", "is_meridional_Pla_center")
        histos["time_meridional_vs_truthhit_z_Qua"][part] = rdf.Histo2D(
            (f"time_meridional_vs_truthhit_z_Qua_{tag}", "", 100, z_lo, z_hi, 500, t_mer_lo, t_mer_hi),
            "OP_pos_produced_z", "OP_analyticalArrivalTime_m", "is_meridional_Qua_center")

        ## --- profile: mean analytical arrival time vs production z — Plastic and Quartz ---
        histos["time_skew_vs_z_prof_Pla"][part] = rdf.Profile1D(
            ROOT.RDF.TProfile1DModel(
                f"time_skew_vs_z_prof_Pla_{tag}", "", 100, z_lo, z_hi, "s"),
            "OP_pos_produced_z", "OP_analyticalArrivalTime_s", "is_skew_Pla_center")
        histos["time_skew_vs_z_prof_Qua"][part] = rdf.Profile1D(
            ROOT.RDF.TProfile1DModel(
                f"time_skew_vs_z_prof_Qua_{tag}", "", 100, z_lo, z_hi, "s"),
            "OP_pos_produced_z", "OP_analyticalArrivalTime_s", "is_skew_Qua_center")
        histos["time_meridional_vs_z_prof_Pla"][part] = rdf.Profile1D(
            ROOT.RDF.TProfile1DModel(
                f"time_meridional_vs_z_prof_Pla_{tag}", "", 100, z_lo, z_hi, "s"),
            "OP_pos_produced_z", "OP_analyticalArrivalTime_m", "is_meridional_Pla_center")
        histos["time_meridional_vs_z_prof_Qua"][part] = rdf.Profile1D(
            ROOT.RDF.TProfile1DModel(
                f"time_meridional_vs_z_prof_Qua_{tag}", "", 100, z_lo, z_hi, "s"),
            "OP_pos_produced_z", "OP_analyticalArrivalTime_m", "is_meridional_Qua_center")
        
    return histos


# ---------------------------------------------------------------------------
# Step 3 — draw histograms
# ---------------------------------------------------------------------------

def draw_histos(histos, draw, suffix, outdir, labels, calib=None, fit_group=None):
    """
    Produce all output plots.

    Parameters
    ----------
    histos  : nested dict {hname: {label: TH1/TH2}}
    draw    : draw-range dict from energy_axis_ranges()
    suffix  : string appended to every output filename
    outdir  : directory for output PNG/PDF files
    labels  : ordered list of labels present in histos
    calib   : optional calibration dict (from load_calibration)
    """
    ensure_dir(outdir)
    colors = get_colors(labels)
    z_lo, z_hi = draw["z"]
    t_min, t_max = draw["time"]
    t_range = {"skew": draw["time_skew"], "meridional": draw["time_meridional"]}
    if "rotate" in suffix:
        t_range["meridional"] = (0.0, 16.0)

    def vals(n): return [histos[n][l] for l in labels if l in histos[n]]
    def keys(n): return [l           for l in labels if l in histos[n]]

    # For Cherenkov plots: overlay Plastic and Quartz on the same canvas,
    # labelling as "<label> (Plastic)" and "<label> (Quartz)".
    # Colors are paired per particle: Plastic and Quartz of the same particle
    # share a color family; different particles get different color families.
    # This avoids duplicate colors when multiple particle types are plotted.
    def cer_vals(n_pla, n_qua):
        v = []
        for l in labels:
            if l in histos[n_pla]: v.append(histos[n_pla][l])
            if l in histos[n_qua]: v.append(histos[n_qua][l])
        return v
    def cer_keys(n_pla, n_qua):
        k = []
        for l in labels:
            if l in histos[n_pla]: k.append(f"{l} (Plastic)")
            if l in histos[n_qua]: k.append(f"{l} (Quartz)")
        return k
    def cer_colors(n_pla, n_qua):
        c = []
        for i, l in enumerate(labels):
            if l in histos[n_pla]: c.append(PLA_COLORS[i % len(PLA_COLORS)])
            if l in histos[n_qua]: c.append(QUA_COLORS[i % len(QUA_COLORS)])
        return c
    def cer_markers(n_pla, n_qua):
        m = []
        for i, l in enumerate(labels):
            if l in histos[n_pla]: m.append(PLA_MARKERS[i % len(PLA_MARKERS)])
            if l in histos[n_qua]: m.append(QUA_MARKERS[i % len(QUA_MARKERS)])
        return m

    args1d = dict(dology=True, donormalize=True, mycolors=colors,
                  MCOnly=True, addOverflow=True, addUnderflow=True,
                  outdir=outdir)

    def d1(hname, xlo, xhi, xlabel, outname):
        k = keys(hname)
        DrawHistos(vals(hname), k,
                   xlo, xhi, xlabel, 1e-3, 1e2, "Fraction of events",
                   outname + "_" + suffix, **lpos(k), **args1d)

    d1("eLeaktruth",        *draw["eLeak"],            "Leakage Energy [GeV]",        "eLeaktruth")
    d1("eCalotruth",        *draw["eCalo"],             "Calo Energy [GeV]",           "eCalotruth")
    d1("eTotaltruth",       *draw["eTotal"],            "Total Energy [GeV]",          "eTotaltruth")
    d1("eTotalGeant",       *draw["eTotalGeant"],       "Total Visible Energy [GeV]",  "eTotalGeant")
    d1("eRodtruth",         *draw["eRod"],              "Rod Energy [GeV]",            "eRodtruth")
    d1("ePlatruth",         *draw["ePlatruth"],         "C-Fiber Energy, Plastic [GeV]", "ePlatruth")
    d1("eQuatruth",         *draw["ePlatruth"],         "C-Fiber Energy, Quartz [GeV]",  "eQuatruth")
    d1("eScintruth",        *draw["eScintruth"],        "S-Fiber Energy [GeV]",          "eScintruth")

    args_edep = {**args1d, "donormalize": False}
    
    args_op = {**args1d, "donormalize": False}

    DrawHistos(vals("truthhit_x"), keys("truthhit_x"),
               -20, 20, "x [cm]", 1e-3, 1e2, "Deposited Energy [GeV]",
               f"truthhit_x_{suffix}", **lpos(keys("truthhit_x")), **args_edep)
    DrawHistos(vals("truthhit_y"), keys("truthhit_y"),
               -20, 20, "y [cm]", 1e-3, 1e2, "Deposited Energy [GeV]",
               f"truthhit_y_{suffix}", **lpos(keys("truthhit_y")), **args_edep)
    DrawHistos(vals("truthhit_z"), keys("truthhit_z"),
               z_lo, z_hi, "z [cm]", 1e-3, 1, "Deposited Energy [GeV]",
               f"truthhit_z_{suffix}", **lpos(keys("truthhit_z")), **args_edep)
    DrawHistos(vals("truthhit_r"), keys("truthhit_r"),
               0, 30, "r [cm]", 1e-3, 1e2, "Deposited Energy [GeV]",
               f"truthhit_r_{suffix}", **lpos(keys("truthhit_r")), **args_edep)
    # Cherenkov z-profiles: Plastic and Quartz overlaid
    _ck = cer_keys("truthhit_z_Pla", "truthhit_z_Qua")
    DrawHistos(cer_vals("truthhit_z_Pla", "truthhit_z_Qua"), _ck,
               z_lo, z_hi, "z [cm]", 1e0, 5e3, "# C Photons / event",
               f"truthhit_z_Cer_{suffix}", **lpos(_ck),
               **{**args_edep, "mycolors": cer_colors("truthhit_z_Pla", "truthhit_z_Qua")})
    _ck = cer_keys("truthhit_z_Pla_center", "truthhit_z_Qua_center")
    DrawHistos(cer_vals("truthhit_z_Pla_center", "truthhit_z_Qua_center"), _ck,
               z_lo, z_hi, "z [cm]", 1e0, 5e3, "# C Photons / event (Central)",
               f"truthhit_z_Cer_center_{suffix}", **lpos(_ck),
               **{**args_edep, "mycolors": cer_colors("truthhit_z_Pla_center", "truthhit_z_Qua_center")})
    # Capture angle: Plastic and Quartz overlaid
    _ck = cer_keys("angle_skew_Pla", "angle_skew_Qua")
    DrawHistos(cer_vals("angle_skew_Pla", "angle_skew_Qua"), _ck,
               0, 1, "Cosine of Capture Angle (Skew)", 1e0, 1e6, "Counts",
               f"angle_skew_{suffix}", **lpos(_ck),
               **{**args_edep, "mycolors": cer_colors("angle_skew_Pla", "angle_skew_Qua")})
    _ck = cer_keys("angle_meridional_Pla", "angle_meridional_Qua")
    DrawHistos(cer_vals("angle_meridional_Pla", "angle_meridional_Qua"), _ck,
               0, 1, "Cosine of Capture Angle (Meridional)", 1e0, 1e6, "Counts",
               f"angle_meridional_{suffix}", **lpos(_ck),
               **{**args_edep, "mycolors": cer_colors("angle_meridional_Pla", "angle_meridional_Qua")})
    # Arrival time: Plastic and Quartz overlaid
    _ck = cer_keys("time_skew_Pla", "time_skew_Qua")
    DrawHistos(cer_vals("time_skew_Pla", "time_skew_Qua"), _ck,
               *t_range["skew"], "Analytical Arrival Time (Skew) [ns]", 1, 1e4, "Counts",
               f"time_skew_{suffix}", **lpos(_ck),
               **{**args_op, "mycolors": cer_colors("time_skew_Pla", "time_skew_Qua")})
    _ck = cer_keys("time_meridional_Pla", "time_meridional_Qua")
    DrawHistos(cer_vals("time_meridional_Pla", "time_meridional_Qua"), _ck,
               *t_range["meridional"], "Analytical Arrival Time (Meridional) [ns]", 1, 1e4, "Counts",
               f"time_meridional_{suffix}", **lpos(_ck),
               **{**args_op, "mycolors": cer_colors("time_meridional_Pla", "time_meridional_Qua")})

    # --- 2D plots ---
    args2d = {**args_edep, "dology": False, "drawoptions": "colz",
              "dologz": True, "zmax": 1e2, "zmin": 1e-4,
              "doth2": True, "addOverflow": False, "addUnderflow": False}
    
    args2d_op = args2d.copy()
    args2d_op["zmax"] = 2e3
    args2d_op["zmin"] = 1e-2

    for label in labels:
        if label not in histos["truthhit_x_vs_truthhit_y"]:
            continue
        DrawHistos([histos["truthhit_x_vs_truthhit_y"][label]], [],
                   -20, 20, "y [cm]", -30, 30, "x [cm]",
                   f"truthhit_x_vs_y_{label}_{suffix}", **args2d)
        DrawHistos([histos["truthhit_x_vs_truthhit_z"][label]], [],
                   z_lo, z_hi, "z [cm]", -30, 30, "x [cm]",
                   f"truthhit_x_vs_z_{label}_{suffix}", **args2d)
        DrawHistos([histos["truthhit_y_vs_truthhit_z"][label]], [],
                   z_lo, z_hi, "z [cm]", -30, 30, "y [cm]",
                   f"truthhit_y_vs_z_{label}_{suffix}", **args2d)
        DrawHistos([histos["truthhit_r_vs_truthhit_z"][label]], [],
                   z_lo, z_hi, "z [cm]", 0, 40, "r [cm]",
                   f"truthhit_r_vs_z_{label}_{suffix}", **args2d)
        for ftype in ("Pla", "Qua"):
            hk = f"time_skew_vs_truthhit_z_{ftype}"
            if label in histos[hk]:
                DrawHistos([histos[hk][label]], [],
                           z_lo, z_hi, "z [cm]", t_min, t_max, "Time [ns] (Skew)",
                           f"time_skew_vs_z_{ftype}_{label}_{suffix}", **args2d_op)
            hk = f"time_meridional_vs_truthhit_z_{ftype}"
            if label in histos[hk]:
                DrawHistos([histos[hk][label]], [],
                           z_lo, z_hi, "z [cm]", *t_range["meridional"], "Time [ns] (Meridional)",
                           f"time_meridional_vs_z_{ftype}_{label}_{suffix}", **args2d_op)
        
    # ── Profile overlays: mean arrival time vs z, Plastic and Quartz overlaid ──
    args_prof = {**args1d, "dology": False, "donormalize": False}
    markers = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34]

    for mode, col_m, col_q in [
        ("skew",       "time_skew_vs_z_prof_Pla",       "time_skew_vs_z_prof_Qua"),
        ("meridional", "time_meridional_vs_z_prof_Pla", "time_meridional_vs_z_prof_Qua"),
    ]:
        pv = cer_vals(col_m, col_q)
        pk = cer_keys(col_m, col_q)
        pc = cer_colors(col_m, col_q)
        pm = cer_markers(col_m, col_q)
        if pv:
            DrawHistos(
                pv, pk,
                z_lo, z_hi, "Production z [cm]",
                *t_range[mode], f"Mean arrival time [ns] ({mode})",
                f"time_{mode}_vs_z_prof_{suffix}",
                drawoptions=["E1"] * len(pv),
                markerstyles=pm,
                **lpos(pk),
                **{**args_prof, "mycolors": pc})

    # ── ProfileX from 2D time vs z — all labels × Plastic/Quartz on one canvas ──
    # Each profile is fit with pol1; results shown in a TPaveText via extraToDraw.
    fit_group = fit_group or []
    for mode in ("skew", "meridional"):
        hk_pla = f"time_{mode}_vs_truthhit_z_Pla"
        hk_qua = f"time_{mode}_vs_truthhit_z_Qua"
        t_profs, t_keys, t_cols, t_markers = [], [], [], []
        fit_entries = []  # (legend_label, color, p0, p0err, p1, p1err)
        # Fit ranges matched by label prefix: ele* -> ele range, pi* -> pion range
        FIT_RANGES = {"ele": (-100, -50), "pi": (-100, 10)}

        _rotated = "rotate" in suffix

        def _fit_range(lbl):
            if _rotated:
                return None   # rotated geometry: skip per-label fits; use fit_group instead
            for prefix, rng in FIT_RANGES.items():
                if lbl.startswith(prefix):
                    return rng
            return (z_lo, z_hi)

        for i, label in enumerate(labels):
            fit_range = _fit_range(label)
            for hk, fib, col, mkr in [
                (hk_pla, "Plastic", PLA_COLORS[i % len(PLA_COLORS)], PLA_MARKERS[i % len(PLA_MARKERS)]),
                (hk_qua, "Quartz",  QUA_COLORS[i % len(QUA_COLORS)], QUA_MARKERS[i % len(QUA_MARKERS)]),
            ]:
                if label not in histos[hk]:
                    continue
                prof = histos[hk][label].ProfileX(f"{hk}_prof_{label}_{suffix}")
                ROOT.SetOwnership(prof, False)
                if fit_range is not None:
                    prof.Fit("pol1", "QR", "", *fit_range)
                f1 = prof.GetFunction("pol1") if fit_range is not None else None
                if f1:
                    f1.SetLineColor(col)
                    f1.SetLineWidth(2)
                    f1.SetLineStyle(2)
                    fit_entries.append((f"{label} ({fib})", col,
                                        f1.GetParameter(0), f1.GetParameter(1)))
                t_profs.append(prof)
                t_keys.append(f"{label} ({fib})")
                t_cols.append(col)
                t_markers.append(mkr)

        if not t_profs:
            continue

        pave = make_fit_pave(fit_entries)

        # Plot 1: profiles with pol1 fit lines and parameter box
        DrawHistos(t_profs, t_keys,
                   z_lo, z_hi, "z [cm]", *t_range[mode],
                   f"Mean arrival time [ns] ({mode})",
                   f"time_{mode}_vs_z_prof2d_{suffix}",
                   drawoptions=["E1"] * len(t_profs) + ["same"],
                   markerstyles=t_markers,
                   **lpos(t_keys),
                   extraToDraw=[pave],
                   **{**args_prof, "mycolors": t_cols})

        # Plot 2: profiles only, no fit lines
        t_profs_nofit = []
        for prof in t_profs:
            p = prof.Clone(prof.GetName() + "_nofit")
            ROOT.SetOwnership(p, False)
            p.GetListOfFunctions().Clear()
            t_profs_nofit.append(p)
        DrawHistos(t_profs_nofit, t_keys,
                   z_lo, z_hi, "z [cm]", *t_range[mode],
                   f"Mean arrival time [ns] ({mode})",
                   f"time_{mode}_vs_z_prof2d_nofit_{suffix}",
                   drawoptions=["E1"] * len(t_profs_nofit) + ["same"],
                   markerstyles=t_markers,
                   **lpos(t_keys),
                   **{**args_prof, "mycolors": t_cols})

        # Plot 3: combined fit — only if --fit-group labels were given
        grp_profs, grp_keys, grp_cols, grp_mkrs = [], [], [], []
        grp_entries = []
        grp_labels = [l for l in fit_group if l in labels]
        if grp_labels:
            # Group fit always covers the full detector z range
            grp_fit_lo, grp_fit_hi = z_lo, z_hi
            combined = {}  # fib -> (TProfile clone, color)
            for i, label in enumerate(grp_labels):
                for hk, fib, col, mkr in [
                    (hk_pla, "Plastic", PLA_COLORS[i % len(PLA_COLORS)], PLA_MARKERS[i % len(PLA_MARKERS)]),
                    (hk_qua, "Quartz",  QUA_COLORS[i % len(QUA_COLORS)], QUA_MARKERS[i % len(QUA_MARKERS)]),
                ]:
                    if label not in histos[hk]:
                        continue
                    prof = histos[hk][label].ProfileX(f"{hk}_grp_{label}_{suffix}")
                    ROOT.SetOwnership(prof, False)
                    if fib not in combined:
                        cp = prof.Clone(f"grp_comb_{fib}_{mode}_{suffix}")
                        ROOT.SetOwnership(cp, False)
                        combined[fib] = (cp, col, mkr)
                    else:
                        combined[fib][0].Add(prof)
            for fib, (cp, col, mkr) in combined.items():
                cp.Fit("pol1", "QR", "", grp_fit_lo, grp_fit_hi)
                f1 = cp.GetFunction("pol1")
                if f1:
                    f1.SetLineColor(col)
                    f1.SetLineWidth(2)
                    f1.SetLineStyle(2)
                    grp_entries.append((f"Combined ({fib})", col,
                                        f1.GetParameter(0), f1.GetParameter(1)))
                grp_profs.append(cp)
                grp_keys.append(f"Combined ({fib})")
                grp_cols.append(col)
                grp_mkrs.append(mkr)
        if grp_profs:
            grp_pave = make_fit_pave(grp_entries)
            DrawHistos(grp_profs, grp_keys,
                       z_lo, z_hi, "z [cm]", *t_range[mode],
                       f"Mean arrival time [ns] ({mode}) — combined fit",
                       f"time_{mode}_vs_z_prof2d_grpfit_{suffix}",
                       drawoptions=["E1"] * len(grp_profs) + ["same"],
                       markerstyles=grp_mkrs,
                       **lpos(grp_keys),
                       extraToDraw=[grp_pave],
                       **{**args_prof, "mycolors": grp_cols})


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def run(particle_files, suffix, beam_energy_gev,
        from_root=None, calib_json=None, outdir=None, max_files=100, fit_group=None):
    """
    Parameters
    ----------
    particle_files : dict {label: file_list_path}
                     Labels can be particle names ("ele", "pion") or arbitrary
                     position/scan tags ("ele_center", "ele_x5y0", ...).
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
        all_hnames = hnames
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
    draw_histos(histos, draw, suffix, outdir, labels, calib=calib,
                fit_group=fit_group or [])
    print(f"\nDone. Plots in: {outdir}/")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Make energy deposition and Cherenkov z-profile plots.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    # ── legacy single-particle shortcuts (backward compatible) ──────────────
    p.add_argument("--ele",  help="Text file listing electron ROOT files")
    p.add_argument("--pion", help="Text file listing pion ROOT files")
    p.add_argument("--neu",  help="Text file listing neutron ROOT files (optional)")
    p.add_argument("--mu",   help="Text file listing muon ROOT files (optional)")

    # ── generic multi-input: one or more --input LABEL FILE pairs ───────────
    p.add_argument(
        "--input", nargs=2, metavar=("LABEL", "FILE"),
        action="append", default=[],
        help=(
            "Add an input with an arbitrary label (repeatable). "
            "Use this to compare the same particle at different positions: "
            "--input ele_center inputs/ele_center.txt "
            "--input ele_x5y0 inputs/ele_x5.txt"
        ),
    )

    # ── standard options ─────────────────────────────────────────────────────
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
                   help="Max ROOT files to read per label (default: 100)")
    p.add_argument(
        "--fit-group", nargs="+", metavar="LABEL", default=[],
        help=(
            "Labels to combine and fit together in a single pol1 fit "
            "(one fit per fiber type). Saves an additional PDF alongside the "
            "per-label profile plots. Repeatable: --fit-group ele_x0 ele_x5 ele_x10"
        ),
    )
    return p.parse_args()


if __name__ == "__main__":
    import ROOT
    ROOT.gROOT.SetBatch(True)

    args = parse_args()

    # ── Build particle_files dict ────────────────────────────────────────────
    # Legacy shortcuts appear first (preserves legend order), then --input entries.
    particle_files = OrderedDict()

    for label, path in [("ele", args.ele), ("pion", args.pion),
                        ("neu", args.neu), ("mu", args.mu)]:
        if path:
            particle_files[label] = path

    for label, path in args.input:
        if label in particle_files:
            print(f"WARNING: label '{label}' already registered; overwriting with {path}")
        particle_files[label] = path

    if not particle_files and not args.from_root:
        print("ERROR: provide at least one of --ele / --pion / --neu / --mu, "
              "or --input LABEL FILE, or use --from-root.")
        sys.exit(1)

    # If replotting from ROOT we still need labels — infer from what was given,
    # or fall back to ele+pion as the historical default.
    if args.from_root and not particle_files:
        particle_files = OrderedDict([("ele", ""), ("pion", "")])

    run(particle_files,
        suffix          = args.suffix or f"{int(args.energy)}GeV",
        beam_energy_gev = args.energy,
        from_root       = args.from_root,
        calib_json      = args.calib,
        outdir          = args.outdir,
        max_files       = args.max_files,
        fit_group       = args.fit_group)

    print("Done")