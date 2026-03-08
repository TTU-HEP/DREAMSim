"""
utils.py — shared utilities for DREAMSim plotters
"""

import os
import sys
import json
from collections import OrderedDict
import ROOT


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def build_chain(file_list_path, tree_name="tree", max_files=100):
    """
    Read a text file of ROOT file paths (one per line, # = comment) and
    return a filled TChain.  Raises FileNotFoundError if the list is missing.
    Prints a warning and stops if max_files is reached.
    """
    chain = ROOT.TChain(tree_name)
    if not os.path.exists(file_list_path):
        raise FileNotFoundError(f"Input list not found: {file_list_path}")

    n = 0
    with open(file_list_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            n += 1
            if n > max_files:
                print(f"WARNING [{file_list_path}]: reached {max_files}-file cap, "
                      f"remaining entries ignored.")
                break
            chain.Add(line)
            print(f"  + {line}")
    print(f"  → {n} file(s) added to chain from {file_list_path}")
    return chain


def build_rdfs(particle_files, tree_name="tree", max_files=100):
    """
    Parameters
    ----------
    particle_files : dict  {label: file_list_path}
        e.g. {"ele": "inputs/electrons_40GeV.txt", "pion": "inputs/pions_40GeV.txt"}

    Returns
    -------
    rdfs   : OrderedDict {label: RDataFrame}
    chains : OrderedDict {label: TChain}   (kept alive to prevent GC)
    """
    chains = OrderedDict()
    rdfs   = OrderedDict()
    for label, path in particle_files.items():
        print(f"\nReading {label}: {path}")
        chains[label] = build_chain(path, tree_name, max_files)
        rdfs[label]   = ROOT.RDataFrame(chains[label])
    return rdfs, chains


def ensure_dir(path):
    """Create directory (and parents) if it does not exist."""
    os.makedirs(path, exist_ok=True)
    return path


# ---------------------------------------------------------------------------
# Histogram save / reload
# ---------------------------------------------------------------------------

def save_histos_to_root(histos, outpath):
    """
    Write all histograms in a nested dict  {name: {label: TH1/TH2 proxy}}
    to a ROOT file.  Histograms are stored under the key  "name__label".
    """
    ensure_dir(os.path.dirname(outpath) or ".")
    f = ROOT.TFile(outpath, "RECREATE")
    for hname, hdict in histos.items():
        for label, hptr in hdict.items():
            h = hptr.GetValue() if hasattr(hptr, "GetValue") else hptr
            # Preserve TProfile error option across Clone() — Clone() does
            # not reliably copy SetErrorOption("s"), so read and re-apply.
            err_opt = h.GetErrorOption() if isinstance(h, ROOT.TProfile) else None
            h = h.Clone(f"{hname}__{label}")
            if err_opt is not None:
                h.SetErrorOption(err_opt)
            h.SetDirectory(f)
            h.Write()
    f.Close()
    print(f"Saved histograms → {outpath}")


def load_histos_from_root(inpath, hnames, labels):
    """
    Reload histograms previously saved by save_histos_to_root.
    Returns a nested OrderedDict {name: {label: TH1/TH2}}.
    Missing keys are silently skipped.
    """
    f = ROOT.TFile.Open(inpath)
    if not f or f.IsZombie():
        raise IOError(f"Cannot open ROOT file: {inpath}")

    histos = OrderedDict()
    for hname in hnames:
        histos[hname] = OrderedDict()
        for label in labels:
            key = f"{hname}__{label}"
            h = f.Get(key)
            if h:
                # Clone immediately while the file is still open — this gives
                # Python full ownership of a typed copy (critical for TProfile,
                # which is not reliably kept alive by SetDirectory(0) alone).
                h = h.Clone()
                h.SetDirectory(0)
                if isinstance(h, ROOT.TProfile):
                    h.SetErrorOption("")
                histos[hname][label] = h
            else:
                print(f"WARNING: key '{key}' not found in {inpath}")
    f.Close()
    return histos

# ---------------------------------------------------------------------------
# Calibration sidecar
# ---------------------------------------------------------------------------

def load_calibration(json_path):
    """
    Load calibration constants from a JSON file produced alongside the
    simulation.  Returns a dict; returns {} if the file does not exist.

    Expected format (all values optional):
    {
        "calibSen": 1.766,
        "calibCen": 2.764,
        "calibCph": 4659.0,
        "calibSph": 1.766
    }
    """
    if not os.path.exists(json_path):
        print(f"INFO: no calibration file at {json_path}, using defaults.")
        return {}
    with open(json_path) as f:
        data = json.load(f)
    print(f"Loaded calibration from {json_path}: {data}")
    return data


# ---------------------------------------------------------------------------
# Axis-range helpers
# ---------------------------------------------------------------------------

def energy_axis_ranges(beam_energy_gev, detector_half_z_cm):
    """
    Return two dicts:
      book  — ranges used when booking histograms (determine what data is kept)
      draw  — ranges used when calling DrawHistos  (cosmetic only, can differ)
    """
    E    = beam_energy_gev
    z_lo = -detector_half_z_cm
    z_hi =  detector_half_z_cm
    # Generous time upper bound: n/c * full length * 10 % headroom
    #t_max = detector_half_z_cm * 2 / 100.0 * 1.6 * 1.1  # ns
    t_min = 0
    t_max_zoomed = 10.0
    t_max = 20

    book = {
        "eLeak":          (0,         E * 0.30),
        "eCalo":          (0,         E * 1.05),
        "eTotal":         (E * 0.85,  E * 1.10),
        "eTotalGeant":    (E * 0.80,  E * 1.10),
        "eRod":           (E * 0.45,  E * 1.05),
        "eCentruth":      (0,         E * 0.015),
        "eScintruth":     (0,         E * 0.010),
        "nOPsCer_Cer":    (0,         E * 1.4e3),
        "nOPsCer_Sci":    (0,         E * 1.1e3),
        "z":              (z_lo,      z_hi),
        "time":           (t_min,     t_max),
        "time_zoomed":    (t_min,     t_max_zoomed),
    }
    book["nOPsCer_Cer_perGeV"] = (0, book["nOPsCer_Cer"][1] / max(book["eCentruth"][1], 1e-9))
    book["nOPsCer_Sci_perGeV"] = (0, book["nOPsCer_Sci"][1] / max(book["eScintruth"][1], 1e-9))

    # Draw ranges are the same by default; override selectively as needed.
    draw = dict(book)
    draw["eLeak"]  = (0, E * 0.20)   # tighter for display
    draw["z"]      = (z_lo, z_hi)

    return book, draw


# ---------------------------------------------------------------------------
# Color map
# ---------------------------------------------------------------------------
COLORS = {
    "ele":  ROOT.kRed,
    "pion": ROOT.kGreen + 2,
    "neu":  ROOT.kBlue,
    "mu":   ROOT.kOrange + 1,
    "op":   ROOT.kViolet,
}

_COLD_PALETTE = [
    (0.09, 0.11, 0.49),   # deep navy
    #(0.00, 0.30, 0.60),   # deep blue
    (0.12, 0.47, 0.71),   # steel blue
    #(0.00, 0.75, 0.75),   # cyan
    #(0.00, 0.60, 0.50),   # teal
    (0.18, 0.72, 0.38),   # green
    (0.60, 0.80, 0.20),   # lime
    #(0.95, 0.85, 0.10),   # yellow
    (0.95, 0.60, 0.10),   # amber
    (0.90, 0.35, 0.05),   # orange
    (0.80, 0.10, 0.10),   # red
    (0.50, 0.00, 0.00),   # dark red
]


# Cache: label -> ROOT color index (stable across calls)
_COLOR_CACHE = {}


def get_colors(labels):
    """
    Return a ROOT color index for each label.
    Known labels (ele, pion, ...) use the fixed COLORS dict.
    Unknown labels cycle through _COLD_PALETTE using TColor.GetColor()
    which is the reliable PyROOT API for registering custom RGB colors.
    """
    unknown = sorted(set(l for l in labels if l not in COLORS and l not in _COLOR_CACHE))
    for lbl in unknown:
        r, g, b = _COLD_PALETTE[len(_COLOR_CACHE) % len(_COLD_PALETTE)]
        # GetColor(r,g,b) is the correct PyROOT API — it registers the color
        # and returns a stable index without GC / slot-reuse issues.
        _COLOR_CACHE[lbl] = ROOT.TColor.GetColor(r, g, b)

    return [COLORS[l] if l in COLORS else _COLOR_CACHE[l] for l in labels]