### HG-DREAM G4 simulation: DREAMSim

#### Environment:

For machines with cvmfs mounted, can directly source the environment
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-dbg/setup.sh
```
Change the path according to the OS.

**On HPCC**, everything (ROOT and GEANT4) is compiled inside the singularity environment. Log into the interactive node ([more information](https://www.depts.ttu.edu/hpcc/userguides/Job_User_Guide.pdf)) with e.g.

```
interactive -p nocona
```
From there run the singularity container with the following command:
```
singularity run --cleanenv --bind /lustre:/lustre /lustre/work/yofeng/SimulationEnv/alma9forgeant4_v3.sif
```
The corresponding docker image can be found [here](https://hub.docker.com/repository/docker/yongbinfeng/alma9geant/general), with the build file [here](https://github.com/TTU-HEP/SimulationEnv).

**If running with Docker** on your local machine, firstly pull the docker image
```
docker pull yongbinfeng/alma9geant:latest
```
Then run the container with the following command:
```
docker run -it --rm -h dreamsim -v /path/to/DREAMSIM:/DREAMSIM yongbinfeng/alma9geant:latest
```

**Note** if you have conda installed, exit the conda environment before running the singularity container, otherwise it might cause conflicts with different ROOT versions etc.


#### Compile:

Inside the singularity environment, build program in "build" area,
```
cd /path/to/DREAMSim/directory
cd sim
mkdir build
cd build
cmake ..
make -j 4
```

Structure of software:

- `sim/CaloXSim.cc`: main program
- `sim/src/CaloXDetectorConstruction.cc`: detector geometry and fiber materials
- `sim/src/CaloXSteppingAction.cc`: per-step hit and photon processing
- `sim/src/CaloXEventAction.cc`: per-event data accumulation
- `sim/src/CaloXTree.cc`: ROOT ntuple output and analysis

#### Fiber types

Three fiber types are simulated:

| Fiber     | Role    | Material            | Elements | Density (g/cm³) | n     | Att. length |
|-----------|---------|---------------------|----------|-----------------|-------|-------------|
| S-fiber   | core    | Polystyrene         | C8H8     | 1.05            | 1.622 | 3.0 m       |
|           | clad    | PMMA_Clad           | C5H8O2   | 1.19            | 1.504 | 5.0 m       |
| C-Plastic | core    | PMMA                | C5H8O2   | 1.19            | 1.504 | 5.0 m       |
|           | clad    | Fluorinated_Polymer | C2F2     | 1.43            | 1.42  | 10.0 m      |
| C-Quartz  | core    | Fused_Silica        | SiO2     | 1.19            | 1.468 | 10.0 m      |
|           | clad    | Hard_Polymer        | C2F2     | 1.43            | 1.42  | 10.0 m      |

S-fibers measure scintillation light; C-Plastic and C-Quartz fibers measure Cherenkov light with different core materials.

#### Run the code

```
./CaloXSim -b paramBatch03_single.mac  \
    -jobName testjob -runNumber 001 -runSeq 003  \
    -numberOfEvents 10 -eventsInNtuple 100    \
    -gun_particle e+ -gun_energy_min 100.0 -gun_energy_max 100.0 \
    -sipmType 1
```

All run parameters are defined in `paramBatch03_single.mac` and may be overridden on the command line or via `runBatch03_single_param.sh`.

The output is a ROOT file containing an ntuple (TTree) with per-event hit information for all fiber types.

Verbosity is quiet by default. To enable verbose G4 output, uncomment the relevant lines in `paramBatch03_single.mac`:
```
#/process/list
#/physics_list/list
#/run/verbose 2
#/event/verbose 2
#/tracking/verbose 1
```

#### Job submission on HPCC

The script `jobs/jobSubmission.py` handles that. Run
```
cd jobs
python jobSubmission.py
# It will produce the shell scripts to submit and run the jobs on HPCC.
bash submit_all.sh
```
to submit the jobs, which can then be monitored with `squeue -u $USER`. More information on the HPCC batch system can be found [here](https://www.depts.ttu.edu/hpcc/userguides/Job_User_Guide.pdf).

#### Analysis

The script `plotter/makePlotsRDF.py` handles this. It runs on multithreads with ROOT's RDataFrame. Therefore, from the login node, log to an interactive node with multiple cores, e.g.
```
interactive -p nocona -c 4
```
Then run the script
```
cd plotter
python makePlotsRDF.py
```
under the singularity environment.

To run `makePlotsOptics.py` for the optics study, compile the helper functions in `plotter/macros` first with
```
cd plotter/macros
root -b -q -e ".L functions.cc+"
```
Then run the script
```
python makePlotsOptics.py
```
