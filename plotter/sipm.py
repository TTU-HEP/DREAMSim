"""
sipm.py — SiPM response convolution helpers

Provides make_drs_histos(), which convolves photon arrival-time histograms
with a measured SiPM single-photon response to produce DRSOutput histograms
ready for plotting.

The SiPM response file is a 1-D numpy array sampled at SIPM_DT_NS (50 ps)
per bin, normalized so that the peak equals 1.0.
"""

import os

import numpy as np
import ROOT

# ---------------------------------------------------------------------------
# SiPM response file
# ---------------------------------------------------------------------------
SIPM_RESPONSE_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "data", "AvePulse3mmSIPM_50ps_2200ts_max682ts_Inter4_R1592.npy")
SIPM_DT_NS = 0.05   # 50 ps per bin


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def make_drs_histos(histos, labels, suffix, sipm_path=None):
    """
    Convolve photon arrival-time histograms with the SiPM response.

    Parameters
    ----------
    histos    : nested dict {hname: {label: TH1}}
                Must contain keys "time_{mode}_{fib}" for mode in
                ("skew", "meridional") and fib in ("Pla", "Qua").
    labels    : ordered list of label strings
    suffix    : str — appended to ROOT histogram names to ensure uniqueness
    sipm_path : str or None — override path to the .npy SiPM response file

    Returns
    -------
    dict {hkey: {label: TH1D}} for hkey in:
        "DRSOutput_skew_Pla",       "DRSOutput_skew_Qua",
        "DRSOutput_meridional_Pla", "DRSOutput_meridional_Qua"

    Returns an empty dict if the SiPM file cannot be found.

    Notes
    -----
    The source histograms must have the same bin width as the SiPM response
    (SIPM_DT_NS = 0.05 ns = 50 ps).  The output histograms span
    [t0, t0 + (n_photon_bins + n_sipm_bins - 1) * SIPM_DT_NS].
    """
    if sipm_path is None:
        sipm_path = SIPM_RESPONSE_FILE

    try:
        sipm_full = np.load(sipm_path)
    except FileNotFoundError:
        print(f"WARNING: SiPM response file not found at {sipm_path}; "
              "skipping DRSOutput plots")
        return {}

    # First-peak-only variant: truncate at the first zero-crossing after the
    # peak, removing the negative undershoot and secondary bump.
    peak_idx  = int(np.argmax(sipm_full))
    zero_cross = next(
        (i for i in range(peak_idx, len(sipm_full)) if sipm_full[i] <= 0),
        len(sipm_full))
    sipm_1st = sipm_full[:zero_cross]

    pulses = [("DRSOutput",        sipm_full),
              ("DRSOutput1stPeak", sipm_1st)]

    drs = {}
    for prefix, sipm in pulses:
        for mode in ("skew", "meridional"):
            for fib, hk_src in [("Pla", f"time_{mode}_Pla"),
                                 ("Qua", f"time_{mode}_Qua")]:
                key = f"{prefix}_{mode}_{fib}"
                drs[key] = {}
                for l in labels:
                    if l not in histos.get(hk_src, {}):
                        continue
                    drs[key][l] = _convolve_th1(
                        histos[hk_src][l], sipm,
                        f"{key}_{l}_{suffix}")
    return drs


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _convolve_th1(h_src, sipm, hname):
    """
    Convolve the bin contents of h_src with sipm and return a new TH1D.

    The SiPM pulse peak is shifted to t=0 before building the output time
    axis, so the convolution only smears the distribution without shifting
    its mean.  The output has (n_src + n_sipm - 1) bins of width SIPM_DT_NS.
    """
    n  = h_src.GetNbinsX()
    t0 = h_src.GetXaxis().GetBinLowEdge(1)

    arr  = np.array([h_src.GetBinContent(i + 1) for i in range(n)])
    conv = np.convolve(arr, sipm, mode='full')

    # Shift the time axis so the SiPM pulse peak aligns with t=0,
    # removing the ~34 ns offset introduced by the pulse delay.
    peak_offset = int(np.argmax(sipm)) * SIPM_DT_NS
    n_out  = len(conv)
    t0_out = t0 - peak_offset
    h = ROOT.TH1D(hname, hname, n_out, t0_out, t0_out + n_out * SIPM_DT_NS)
    ROOT.SetOwnership(h, False)
    for i, v in enumerate(conv):
        h.SetBinContent(i + 1, float(v))
    return h
