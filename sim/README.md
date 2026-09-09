## The CaloX simulation

Geometry, the fiber map, run parameters and output branches. For setting up the
environment, submitting jobs and making plots, see the [top-level
README](../README.md).

### Building

From inside the container:

```
cd sim
mkdir build && cd build
cmake ..
make -j 8
```

`cmake` copies the `.mac` files and `data/fibermap.json` next to the executable,
so the defaults resolve when you run from the build directory.

### Geometry

The calorimeter is 90 copper rods across (x) by 80 layers (y), each rod 0.4 x
0.4 cm in cross section. The copper is 2 m long and the fibers are 2.5 m: they
start at the front face of the copper and run 50 cm past its back, where the
light is guided out to the SiPMs. So the copper spans z = -100..+100 cm and the
fibers z = -100..+150 cm.

The copper length is not free. `CaloXID` measures z from `z0 = -1000 mm` over a
calorimeter of `zback = 2000 mm`, and the leakage split in `CaloXSteppingAction`
is taken at `z = 100 cm`, so the readout already assumes exactly this. The
calorimeter envelope is air, symmetric about the copper and long enough to reach
the far end of the fibers, which lets it be placed and rotated exactly as the
copper itself. The `Calorimeter` and `Layer` volumes are air throughout — all
the absorber is in the rods.

The nesting is

```
World -> Calorimeter -> Layer -> Rod -> Hole -> fiber cladding -> fiber core
```

which is what the copy numbers in `CaloXSteppingAction` count on: at a step in a
fiber core, `GetCopyNumber(1)` is the fiber, `(2)` the hole, `(3)` the rod and
`(4)` the layer.

Rods are placed explicitly rather than replicated. A `G4PVReplica` shares one
logical volume across all its copies, so with replicas every copper is forced to
be identical and the map below could not be expressed at all. Layers that hold
the same sequence of rod types do share a logical volume, so the 90 x 80
calorimeter is built from about a dozen distinct layer patterns rather than 7200
separate placements.

There are four kinds of rod:

| Rod | Contents |
|-----|----------|
| `RodQuartz`  | a hole with 4 quartz Cherenkov + 3 scintillating fibers |
| `RodPlastic` | a hole with 4 plastic Cherenkov + 3 scintillating fibers |
| `RodEmpty`   | solid copper, no hole |
| `AirGap`     | air — outside the detector outline, no copper at all |
| `FiberTail`  | air — past the back of the copper, holding the last 50 cm of fiber |

The first three are placed as physical volume `Rod`, which is how the stepping
action recognises absorber. `AirGap` and `FiberTail` deliberately are not, so
they fall through to `caloType 0` and contribute to `eCalotruth` but not to
`eRodtruth` — the same treatment the air inside the fiber holes gets.

A `FiberTail` sits at the same nesting depth as a rod and carries the same copy
number, so a step in a tail fiber reports the same rod and layer as one in the
copper. Because of that the fiber is built as two segments, the 2 m inside the
copper and the 50 cm tail, which keeps every volume a simple box or tube and
lets the copper carve its own hole. The segments meet face to face at z = +100
cm with the same material on both sides and no optical surface between them, so
a photon crosses the joint without noticing it.

#### Fibers within a copper

Seven slots: one on the axis and six at 30, 90, 150, 210, 270 and 330 degrees,
on a circle of radius `2 * r_clad + 10 um`. Four carry Cherenkov fibers and
three carry scintillating fibers:

```
            90 S
     150 C        30 C
            0 C            <- on the axis
     210 S        330 S
            270 C
```

Copy numbers are 0..3 for the Cherenkov fibers and 1..3 for the scintillating
ones. Only the Cherenkov flavour varies between coppers; the slot layout never
does. Each fiber is a 0.40 mm radius cladding around a 0.39 mm radius core.

### The fiber map

`data/fibermap.json` says what each copper contains. Its path comes from the
`fiberMapFile` mac parameter. It is required — if it cannot be read the run
stops, rather than silently building a detector nobody asked for. `noRods` and
`noLayers` are taken from the map so the two can never disagree.

The map is a grid of one character per cell:

| Character | Copper contains |
|-----------|-----------------|
| `Q` | 4 quartz Cherenkov + 3 scintillating fibers |
| `P` | 4 plastic Cherenkov + 3 scintillating fibers |
| `.` | nothing — solid copper |
| `_` | nothing, and no copper either — air |

A cell is 3 rods by 4 layers, which is exactly `CaloXID`'s `(ix, iy)`, so the
map is indexed in the same coordinates as the readout. In the central region —
`ix` 13..16, `iy` 8..11, the `area == 3` cells with the 3 mm SiPMs — a cell is 3
rods by 1 layer instead, addressed by `(ix, iy, iyy)` with `iyy = layer % 4`.
That finer block lives in the `central` object; the outer grid marks those cells
`*` to show they are taken from there.

Rows run from the highest `iy` (top) downwards and each row runs from `ix = 0`
(left) rightwards, which makes the file read like a picture of the detector.
The map as actually built is printed at the start of every run, so each output
file records the geometry it used.

#### Regenerating it

`tools/fibermap_from_image.py` extracts the map from the two engineering plots:

```
python3 tools/fibermap_from_image.py \
    --full image.png --zoom image_zoomed.png -o sim/data/fibermap.json
```

It needs `numpy` and `Pillow`. Two features of the plots matter if you touch
this. Every box is drawn with a thin green outline, which looks exactly like a
green fill to a component finder and connects neighbouring boxes into a single
blob; and every box in the central plot has its channel number painted across
it, which splits one box into several components. The whole-detector plot is
therefore read by finding components with a cut on how much of its bounding box
each one fills, and the central plot by sampling box interiors on its known
lattice.

Two parts of the map cannot be read from the plots and are pinned in the tool,
with the reasoning beside them. **Both are worth checking against the real
engineering drawings.**

- `EDGE_COLUMNS` — the outermost column on each side is drawn shifted by half a
  cell in y, and no integer row assignment keeps the map symmetric, which is
  what a genuine 2-layer stagger looks like. What those columns contain is not
  in doubt (one empty cell, two plastic, two empty, identically on both sides);
  which rows they land on is, by one. If the real modules are staggered, a 4
  layer cell cannot express that and the map needs finer granularity there.
- `OBSCURED_CELLS` — the four cells at `iy = 0`, `ix = 13..16` sit behind an
  overlay in the plot. They are taken from their mirror at `iy = 19`. Every
  other row of the map is symmetric under both reflections.

The resulting detector has 1878 quartz coppers, 3210 plastic, 696 bare copper
and 1416 absent.

### Run parameters

Set in the `.mac` file as `#$$$ key value`, and overridable on the command line
as `-key value`. Beyond the beam and bookkeeping parameters, these control the
detector and the output size:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `fiberMapFile` | `data/fibermap.json` | the map described above; required |
| `fiberTailLength` | 50.0 | cm the fibers run past the 2 m copper; 0 gives 2 m fibers |
| `saveOpticalPhotons` | `false` | write the per-photon `OP_*` branches |
| `saveTruthHits` | — | write the per-step `truthhit_*` branches |
| `opSampleRod`, `opSampleLayer` | 45, 40 | the one rod whose optical photons are tracked past their first step |
| `sipmType` | — | 1 = J 6 mm 6.0 V, 2 = J 6 mm 2.5 V |
| `caloRotationX`, `caloRotationY` | — | calorimeter tilt in degrees |

`saveOpticalPhotons` is off by default because one record is written for every
optical photon, which is order 10^5 per 100 GeV shower and dominates the file
size — roughly 13x larger files and 20% slower. The cheap per-event counters
`nOPsCer`, `nOPsCer_Pla`, `nOPsCer_Qua` and `nOPsCer_Sci` are always written, so
the Cherenkov photon yields are available in every file.

### Output

One `TTree` named `tree`. The branches most easily misread:

- `sum3dCC`, `sum3dQQ` — **detected photoelectrons**, in plastic and quartz
  Cherenkov fibers, summed over cells with at least one photoelectron and over
  the `area >= 2` fiducial region. Detection is sampled per photon as a
  Bernoulli trial against the SiPM PDE, so these fluctuate as
  `Binomial(N_trapped, pde)` and carry the correct photostatistics. Accumulating
  the PDE as a weight instead — which is what the code used to do — gives the
  expectation value of the count with a variance smaller by roughly `<pde>`, and
  so far too little Cherenkov resolution smearing.
- `sum3dSS` — scintillation, Birks-suppressed energy in GeV, over the same
  fiducial region and cell logic. Use this rather than a raw sum of
  `truthhit_edep` when comparing against `sum3dCC + sum3dQQ`, or the two sides
  see different acceptances.
- `truthhit_ncer` — Cherenkov photons produced in the step, before any capture
  or PDE.
- `truthhit_ncertrap` — those inside the capture cone, before the PDE.
- `truthhit_ncercap` — sampled photoelectrons, an integer count.
- `eScintruth`, `ePlatruth`, `eQuatruth`, `eRodtruth` — energy deposited in each
  fiber type and in the copper, in GeV, with no acceptance cuts.
- `eInvisible` — nuclear binding energy and the like, from `findInvisible`,
  which reconstructs it per step. It subtracts the rest mass of protons and
  neutrons created out of the vacuum but not that of pions, so an event with a
  photonuclear pion carries about 135 MeV of energy that was never there, and
  `eCalotruth + eWorldtruth + eLeaktruth + eInvisible` overshoots the beam
  energy by that much. It fires in a couple of events in ten at 100 GeV.

The light-trapping model is a fixed cone, `theta < 0.336` about the global z
axis at the point of production, with no propagation, attenuation or cladding
refraction. That is the numerical aperture in air rather than inside the core,
it counts only the forward cone, and it ignores the calorimeter tilt, so the
absolute photoelectron yield should not be taken at face value. It mostly shifts
the mean rather than the fluctuations. A fuller analytic treatment — proper
total internal reflection with skew rays, Beer-Lambert attenuation and arrival
times — is computed per photon and written to the `OP_*` branches when
`saveOpticalPhotons` is on, but it does not feed `sum3dCC` / `sum3dQQ`.
