submissiontext = """#!/bin/bash
#SBATCH -J JOBNAME
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o LOGDIR/%x.%j.out
#SBATCH -e LOGDIR/%x.%j.err
#SBATCH -p nocona
"""
singularity_cmd = "singularity run --cleanenv --bind /lustre:/lustre /lustre/work/yofeng/SimulationEnv/alma9forgeant4_v3.sif"

# Gun position/momentum flags are appended at generation time for rotated mode.
# Standard (along Z) mode: no gun_* overrides — mac file defaults apply.
# Rotated (along X) mode:  gun fires in +X, position scanned in Z.
run_cmd = ("BUILDDIR/CaloXSim -b BUILDDIR/paramBatch03_single.mac"
           "      -jobName JOBNAME -runNumber RUNNUM -runSeq RUNSEQ"
           "      -numberOfEvents NEVENTS  -eventsInNtuple NEVENTS"
           "      -gun_particle PARTICLE"
           "      -gun_energy_min ENERGY_MIN -gun_energy_max ENERGY_MAX"
           "      -sipmType 1"
           "      GUN_OVERRIDES")

import os
current_dir = os.getcwd()
print("Current directory: ", current_dir)
output_dir = "/lustre/work/yofeng/SimulationOutputs/20260311_v3_norm/"
output_dir = "/lustre/research/hep/yofeng/SimulationOutputs/20260311_v3"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
sim_build_dir = f"{current_dir}/../sim/build"
log_dir       = f"{current_dir}/log"
njobs             = 100
nevents_per_job   = 50
runnumber         = 123
jobname_prefix    = "dreamsim"

# ---------------------------------------------------------------------------
# Rotated-beam scan settings
# ---------------------------------------------------------------------------
# When rotated_beam=True the particle gun fires in +X (transverse to fibres).
# gun_x is fixed just outside the detector face (half-width = 18 cm).
# gun_z is scanned across the detector depth (fibre half-length = 125 cm).
# gun_y stays at the same transverse position as the standard setup.
ROTATED_BEAM   = True
GUN_X          = -40.0   # cm  — just outside the X face (detector face at -18 cm)
GUN_Y          = -2.5    # cm  — same transverse position as standard setup
GUN_Z_SCAN     = list(range(-120, 121, 60))   # cm: -120, -100, ..., +120 (49 z slices)
#
# To run a single position instead of a scan, set e.g.:
#   GUN_Z_SCAN = [0.0]
# To run the standard (non-rotated) beam, set:
ROTATED_BEAM = False

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def gun_overrides(rotated, gun_z=None):
    """Return the command-line flag string for gun position/momentum."""
    if not rotated:
        return ""   # use mac-file defaults
    return (f" -gun_x_min {GUN_X} -gun_x_max {GUN_X}"
            f" -gun_y_min {GUN_Y} -gun_y_max {GUN_Y}"
            f" -gun_z_min {gun_z} -gun_z_max {gun_z}"
            f" -pMomentum_x_min 1.0 -pMomentum_x_max 1.0"
            f" -pMomentum_y_min 0.0 -pMomentum_y_max 0.0"
            f" -pMomentum_z_min 0.0 -pMomentum_z_max 0.0")


def z_tag(gun_z):
    """Human-readable z label safe for filenames, e.g. -120 → 'm120', 20 → 'p020'."""
    sign = "m" if gun_z < 0 else "p"
    return f"{sign}{abs(int(gun_z)):03d}"


def ensure_dirs(*dirs):
    for d in dirs:
        if not os.path.exists(d):
            print(f"Creating directory: {d}")
            os.makedirs(d)


# ---------------------------------------------------------------------------
# Job generation
# ---------------------------------------------------------------------------

def generate_submission_script(particle, energy_min, energy_max,
                                rotated=ROTATED_BEAM, z_positions=None):
    """
    Generate per-job SLURM scripts and a master submit_all_*.sh.

    Parameters
    ----------
    particle    : str   e.g. "pi+"
    energy_min  : float GeV
    energy_max  : float GeV
    rotated     : bool  True → shoot in +X with z_positions scan
    z_positions : list  gun_z values to scan (only used when rotated=True)
    """
    if rotated and z_positions is None:
        z_positions = GUN_Z_SCAN
    if not rotated:
        z_positions = [None]   # single pass, no z override

    fnames = []
    for gun_z in z_positions:
        ztag  = f"_z{z_tag(gun_z)}" if gun_z is not None else ""
        overrides = gun_overrides(rotated, gun_z)

        for i in range(njobs):
            jobname = (f"{jobname_prefix}_{particle}_{i}"
                       f"_{energy_min}_{energy_max}{ztag}")

            cmd = (run_cmd
                   .replace("BUILDDIR",    sim_build_dir)
                   .replace("JOBNAME",     jobname)
                   .replace("RUNNUM",      str(runnumber))
                   .replace("RUNSEQ",      str(i))
                   .replace("NEVENTS",     str(nevents_per_job))
                   .replace("PARTICLE",    particle)
                   .replace("ENERGY_MIN",  str(energy_min))
                   .replace("ENERGY_MAX",  str(energy_max))
                   .replace("GUN_OVERRIDES", overrides))

            full_cmd = (singularity_cmd
                        + " bash -c \""
                        + "cd " + output_dir + " && "
                        + cmd
                        + " \"")

            header = (submissiontext
                      .replace("JOBNAME", jobname)
                      .replace("LOGDIR",  log_dir))

            fname = f"{log_dir}/submit_{jobname}.sh"
            with open(fname, "w") as f:
                f.write(header)
                f.write("\n\n")
                f.write(full_cmd)
            fnames.append(fname)

    rot_tag = "_rotated" if rotated else ""
    submit_sh = (f"{current_dir}/submit_all"
                 f"_{particle}_{energy_min}_{energy_max}{rot_tag}.sh")
    with open(submit_sh, "w") as f:
        f.write("#!/bin/bash\n")
        for fname in fnames:
            f.write(f"sbatch {fname}\n")
    os.system(f"chmod +x {submit_sh}")

    n_z  = len(z_positions)
    print(f"  {particle} {energy_min}-{energy_max} GeV"
          f"{'  rotated' if rotated else ''}"
          f"  {n_z} z-position(s) × {njobs} jobs = {n_z * njobs} jobs total")
    print(f"  → {submit_sh}")
    print(f"  To submit: bash {submit_sh}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    ensure_dirs(log_dir, output_dir)

    samples = {
        #"pi+": [(40, 40.5), (80, 80.5)],
        #"e+":  [(40.0, 40.5), (80.0, 80.5)],
        "e+": [(40.0, 40.5)],
        "pi+": [(40, 40.5)],
    }

    for particle, energies in samples.items():
        for energy_min, energy_max in energies:
            generate_submission_script(particle, energy_min, energy_max)
