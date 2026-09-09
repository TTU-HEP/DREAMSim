#!/usr/bin/env python3
"""Extract the DREAMSim per-copper fiber map from the two engineering plots.

The detector plot draws one cell per 3 (x) by 4 (y) block of coppers, as a colour
box stacked on top of a blue box:

    green  over blue  -> 4 plastic Cherenkov + 3 scintillating fibers   ('P')
    orange over blue  -> 4 quartz  Cherenkov + 3 scintillating fibers   ('Q')
    red    over red   -> solid copper with no fibers                    ('.')
    nothing drawn     -> outside the outline, no copper at all          ('_')

The central region (area 3 in CaloXID) is drawn as an inset in the detector plot
and separately, much more legibly, in the zoomed plot.  There one cell is 3 (x)
by 1 (y) coppers and the colour box sits beside the blue box rather than above
it, so the two plots are read with slightly different rules.

Usage:
    fibermap_from_image.py --full image.png --zoom image_zoomed.png -o fibermap.json
"""

import argparse
import json
import sys
from collections import Counter

import numpy as np

try:
    from PIL import Image
except ImportError:
    sys.exit("this tool needs Pillow:  pip install Pillow")

# Reference colours, sampled from the clean (zoomed) plot.
REF = {
    "blue":   (39, 118, 187),
    "green":  (19, 103, 53),
    "orange": (242, 103, 44),
    "red":    (237, 28, 36),
    "white":  (255, 255, 255),
    "black":  (0, 0, 0),
}
KIND = {"green": "P", "orange": "Q", "red": "."}

# Lattice of the full detector plot, in pixels.  30 columns; each row is a colour
# box at Y_COLOUR[r] with its blue partner ~17 px below.
X_COL = [67, 99, 127, 157, 190, 214, 246, 278, 309, 341, 372, 404, 435, 467, 499,
         532, 562, 594, 626, 658, 689, 721, 752, 784, 816, 840, 872, 902, 931, 962]
Y_COLOUR = [35, 71, 108, 144, 181, 217, 254, 292, 332, 372,
            424, 464, 503, 543, 580, 616, 653, 689, 726, 762]

# Lattice of the zoomed central plot: 4 cell columns, each a (blue, colour) pair
# of sub-columns, and 16 rows (4 values of iy x 4 layers within each).
Z_XPAIR = [(52, 84), (148, 181), (248, 280), (344, 376)]
Z_YROW = [44, 82, 117, 152, 186, 220, 254, 288,
          332, 366, 400, 434, 469, 504, 538, 576]

# The central region, in cell coordinates.  Matches CaloXID::findArea().
IX_MIN, IX_MAX, IY_MIN, IY_MAX = 13, 16, 8, 11

RODS_PER_CELL_X = 3
LAYERS_PER_CELL_Y = 4

# The outermost module on each side is drawn shifted by half a cell in y, so its
# boxes cannot be matched to the 4-layer row grid unambiguously -- there is no
# integer row assignment that keeps the map symmetric, which is what a genuine
# half-cell (2 layer) stagger looks like.  What each column *contains* is not in
# doubt: reading its boxes from the top gives one empty cell, two plastic cells
# and two empty cells, identically on both sides.  That is pinned here, aligned
# so the first cell sits on row 7.  The alignment is uncertain by one row.
EDGE_COLUMNS = {0: (7, ".PP.."), 29: (7, ".PP..")}

# Something is drawn across the bottom centre of the detector plot -- a small
# overlay covering the cells at iy = 0, ix = 13..16 -- so those four cannot be
# read.  Every other row of the map mirrors its partner under iy -> 19 - iy
# exactly, and the mirror of this one (iy = 19) is PPPP, so that is what is used
# here.  Delete this entry to leave them empty instead.
OBSCURED_CELLS = {(0, 13): "P", (0, 14): "P", (0, 15): "P", (0, 16): "P"}


def classify(path, tol):
    """Map every pixel onto the nearest reference colour, or -1 if too far."""
    a = np.array(Image.open(path).convert("RGB")).astype(int)
    names = list(REF)
    ref = np.array([REF[n] for n in names])
    d = ((a[:, :, None, :] - ref[None, None, :, :]) ** 2).sum(-1)
    lab = np.where(np.sqrt(d.min(-1)) < tol, d.argmin(-1), -1)
    return lab, names


def components(mask, min_px):
    """Connected components of a boolean mask, as bounding boxes + centroids."""
    h, w = mask.shape
    seen = np.zeros((h, w), bool)
    out = []
    for y0 in range(h):
        for x0 in np.flatnonzero(mask[y0]):
            if seen[y0, x0]:
                continue
            stack, ys, xs = [(y0, x0)], [], []
            seen[y0, x0] = True
            while stack:
                cy, cx = stack.pop()
                ys.append(cy)
                xs.append(cx)
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ny, nx = cy + dy, cx + dx
                    if 0 <= ny < h and 0 <= nx < w and mask[ny, nx] and not seen[ny, nx]:
                        seen[ny, nx] = True
                        stack.append((ny, nx))
            if len(ys) >= min_px:
                ys, xs = np.array(ys), np.array(xs)
                bw = xs.max() - xs.min() + 1
                bh = ys.max() - ys.min() + 1
                out.append((xs.mean(), ys.mean(), bw, bh, len(ys) / float(bw * bh)))
    return out


def boxes_of(path, colours, tol, min_px, wmin, wmax, min_fill=0.55):
    """Every filled box of the given colours, at the expected size.

    Every box in these plots is drawn with a thin green outline, so a blue or
    orange box contributes a green ring that has the size of a real green box.
    A ring covers only its own perimeter, so requiring the component to fill
    most of its bounding box separates fills from outlines.
    """
    lab, names = classify(path, tol)
    found = []
    for c in colours:
        for cx, cy, w, h, fill in components(lab == names.index(c), min_px):
            if wmin <= w <= wmax and wmin <= h <= wmax and fill >= min_fill:
                found.append((cx, cy, c))
    return found


def read_full(path, verbose):
    """The 30 x 20 cell grid, as a list of rows ordered top row first.

    The outermost module on each side is drawn offset by half a cell in y, so
    every column gets its own offset, taken as the median displacement of its
    boxes from the nominal lattice.  Rows are then assigned in that column's
    own frame.
    """
    #  No fill cut here: in this plot the boxes are outlined in black, and they
    #  are small enough that a cut on filled area is not a safe discriminator.
    found = boxes_of(path, ("green", "orange", "red"), tol=60, min_px=25,
                     wmin=6, wmax=16, min_fill=0.0)
    by_col = {}
    for cx, cy, colour in found:
        c = int(np.argmin([abs(cx - x) for x in X_COL]))
        if abs(cx - X_COL[c]) > 14:
            continue                       # a box that belongs to no column
        by_col.setdefault(c, []).append((cy, colour))

    grid = [["_"] * len(X_COL) for _ in Y_COLOUR]
    stagger = {}
    for c, entries in sorted(by_col.items()):
        raw = [cy - min(Y_COLOUR, key=lambda y: abs(y - cy)) for cy, _ in entries]
        off = float(np.median(raw))
        if abs(off) > 4.0:
            stagger[c] = off
        for cy, colour in entries:
            r = int(np.argmin([abs(cy - off - y) for y in Y_COLOUR]))
            #  Half the row pitch: a staggered column sits up to half a
            #  cell from the nominal lattice and must still be matched.
            if abs(cy - off - Y_COLOUR[r]) > 18:
                continue
            want = KIND[colour]
            if grid[r][c] not in ("_", want) and verbose:
                print(f"  clash at row {r} col {c}: {grid[r][c]} vs {want}")
            grid[r][c] = want
    for c, (row0, cells) in EDGE_COLUMNS.items():
        for r in range(len(Y_COLOUR)):
            grid[r][c] = "_"
        for k, ch in enumerate(cells):
            grid[row0 + k][c] = ch

    for (iy, ix), ch in sorted(OBSCURED_CELLS.items()):
        r = len(Y_COLOUR) - 1 - iy
        print(f"  NOTE: cell ix={ix}, iy={iy} is hidden behind an overlay in the "
              f"plot;\n        taking {ch!r} from its mirror at iy={len(Y_COLOUR) - 1 - iy}"
              f" -- please check.")
        grid[r][ix] = ch

    for c, off in sorted(stagger.items()):
        cells = "".join(grid[r][c] for r in range(len(Y_COLOUR)))
        pinned = " (pinned by EDGE_COLUMNS)" if c in EDGE_COLUMNS else ""
        print(f"  NOTE: column ix={c} is drawn {off:+.0f} px from the nominal\n"
              f"        lattice (roughly half a cell of 4 layers), so which row\n"
              f"        each of its cells belongs to is ambiguous by one.\n"
              f"        Read as {cells!r}{pinned} -- please check.")

    return ["".join(r) for r in grid], stagger


def read_zoom(path, verbose):
    """The 4 x 16 central grid, as a list of rows ordered top row first.

    The lattice of this plot is regular and known, so each box is read by
    sampling its interior rather than by finding connected components.  That
    avoids both traps here: every box carries a thin green outline, and the
    outlines of neighbouring boxes touch, so component finding merges the whole
    plot into one blob; and each box has its channel number painted across it,
    which splits a component in two.  A patch well inside the box sees neither.
    """
    lab, names = classify(path, tol=60)
    ib, ig, io = (names.index(c) for c in ("blue", "green", "orange"))
    rows = []
    for y in Z_YROW:
        row = ""
        for xa, xb in Z_XPAIR:
            fills = []
            for x in (xa, xb):
                patch = lab[int(y) - 9:int(y) + 10, int(x) - 9:int(x) + 10].ravel()
                counts = {c: int((patch == c).sum()) for c in (ib, ig, io)}
                best = max(counts, key=counts.get)
                fills.append(names[best] if counts[best] >= 40 else None)
            colour = [f for f in fills if f is not None and f != "blue"]
            if len(colour) != 1:
                print(f"  warning: central cell at y~{y}, x~{xa} reads {fills}; "
                      f"marking unknown")
                row += "?"
            else:
                row += KIND[colour[0]]
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--full", required=True, help="the whole-detector plot")
    ap.add_argument("--zoom", required=True, help="the central-region plot")
    ap.add_argument("-o", "--out", required=True, help="fiber map to write")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    print(f"reading {args.full} ...")
    grid, stagger = read_full(args.full, args.verbose)
    print(f"reading {args.zoom} ...")
    central = read_zoom(args.zoom, args.verbose)

    # The inset in the detector plot covers the central cells, so whatever was
    # read there is unreliable; blank it and let the zoomed plot supply it.
    for r, _ in enumerate(grid):
        iy = len(Y_COLOUR) - 1 - r
        if IY_MIN <= iy <= IY_MAX:
            row = list(grid[r])
            for ix in range(IX_MIN, IX_MAX + 1):
                row[ix] = "*"
            grid[r] = "".join(row)

    doc = {
        "format": "DREAMSim fiber map v1",
        "source": f"{args.full} + {args.zoom}",
        "nRods": len(X_COL) * RODS_PER_CELL_X,
        "nLayers": len(Y_COLOUR) * LAYERS_PER_CELL_Y,
        "rodsPerCellX": RODS_PER_CELL_X,
        "layersPerCellY": LAYERS_PER_CELL_Y,
        "legend": {
            "Q": "4 quartz Cherenkov + 3 scintillating fibers",
            "P": "4 plastic Cherenkov + 3 scintillating fibers",
            ".": "solid copper with no fibers (drawn red)",
            "_": "no copper at all, air (outside the drawn detector outline)",
            "*": "taken from the central grid below",
        },
        "gridNote": ("one row per iy, listed from iy = nLayers/layersPerCellY - 1 "
                     "(top) down to iy = 0; one character per ix = 0 (left) to "
                     "nRods/rodsPerCellX - 1 (right)"),
        "grid": grid,
        "central": {
            "note": ("the cells of area 3, one row per layer, from the top layer "
                     "(iy = iyMax, iyy = 3) down to (iy = iyMin, iyy = 0); one "
                     "character per ix = ixMin (left) to ixMax (right)"),
            "ixMin": IX_MIN, "ixMax": IX_MAX,
            "iyMin": IY_MIN, "iyMax": IY_MAX,
            "grid": central,
        },
    }

    with open(args.out, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")

    print(f"\nwrote {args.out}")
    print(f"  outer grid {len(grid)} rows x {len(grid[0])} cols  "
          f"{dict(Counter(''.join(grid)))}")
    print(f"  central   {len(central)} rows x {len(central[0])} cols  "
          f"{dict(Counter(''.join(central)))}")
    print()
    print("  " + "".join(str(i // 10) for i in range(len(X_COL))))
    print("  " + "".join(str(i % 10) for i in range(len(X_COL))))
    for r, row in enumerate(grid):
        print(f"  {row}   iy={len(Y_COLOUR) - 1 - r}")
    print()
    for r, row in enumerate(central):
        iy = IY_MAX - r // LAYERS_PER_CELL_Y
        iyy = 3 - r % LAYERS_PER_CELL_Y
        print(f"  central {row}   iy={iy} iyy={iyy}")


if __name__ == "__main__":
    main()
