submissiontext = """#!/bin/bash
#SBATCH -J JOBNAME
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o LOGDIR/%x.%j.out
#SBATCH -e LOGDIR/%x.%j.err
#SBATCH -p nocona
"""

singularity_cmd = "singularity run --cleanenv --bind /lustre:/lustre /lustre/work/yofeng/SimulationEnv/alma9forgeant4_v1_sbox/"
#run_cmd = "BUILDDIR/exampleB4b -b BUILDDIR/paramBatch03_single.mac      -jobName JOBNAME -runNumber RUNNUM -runSeq RUNSEQ      -numberOfEvents NEVENTS  -eventsInNtupe NEVENTS -gun_particle PARTICLE -gun_energy_min ENERGY_MIN -gun_energy_max ENERGY_MAX     -sipmType 1"
#run_cmd = "BUILDDIR/exampleB4b -b BUILDDIR/paramBatch03_single_photon.mac      -jobName JOBNAME -runNumber RUNNUM -runSeq RUNSEQ      -numberOfEvents NEVENTS  -eventsInNtupe NEVENTS -gun_particle PARTICLE -sipmType 1"
run_cmd = "BUILDDIR/exampleB4b -b BUILDDIR/paramBatch03_single.mac      -jobName JOBNAME -runNumber RUNNUM -runSeq RUNSEQ      -numberOfEvents NEVENTS  -eventsInNtupe NEVENTS -gun_particle PARTICLE -gun_energy_min ENERGY_MIN -gun_energy_max ENERGY_MAX -sipmType 1"

import os

current_dir = os.getcwd()
print("Current directory: ", current_dir)

output_dir = "/lustre/work/yofeng/SimulationOutputs/20250530/"

#
# change from here
#
sim_build_dir = f"{current_dir}/../sim/build"
log_dir = f"{current_dir}/log"
njobs = 100
nevents_per_job = 30
runnumber = 100
particle = "e+"
energy_min = 100
energy_max = 101
jobname_prefix = "dreamsim"


if not os.path.exists(log_dir):
    print(f"Creating log directory: {log_dir}")
    os.makedirs(log_dir)
    
if not os.path.exists(output_dir):
    print(f"Creating output directory: {output_dir}")
    os.makedirs(output_dir)
    

def generate_submission_script(particle, energy_min, energy_max):
    fnames = []
    for i in range(njobs):
        jobname = f"{jobname_prefix}_{particle}_{i}_{energy_min}_{energy_max}"
        run_cmd_tmp = run_cmd.replace("BUILDDIR", sim_build_dir).\
            replace("JOBNAME", jobname).\
            replace("RUNNUM", str(runnumber)).\
            replace("RUNSEQ", str(i)).\
            replace("NEVENTS", str(nevents_per_job)).\
            replace("PARTICLE", particle).\
            replace("ENERGY_MIN", str(energy_min)).\
            replace("ENERGY_MAX", str(energy_max))

        run_cmd_tmp = singularity_cmd + " bash -c \"" + "cd " + output_dir + " && " + run_cmd_tmp + " \""

        submissiontext_tmp = submissiontext.replace("JOBNAME", jobname).replace("LOGDIR", log_dir)

        # write the submission script
        fname = f"{log_dir}/submit_{jobname}.sh"
        with open(fname, "w") as f:
            f.write(submissiontext_tmp)
            f.write("\n\n")
            f.write(run_cmd_tmp)

        fnames.append(fname)

    submit_sh = f"{current_dir}/submit_all_{particle}_{energy_min}_{energy_max}.sh"
    with open(submit_sh, "w") as f:
        f.write("#!/bin/bash\n")
        for fname in fnames:
            f.write(f"sbatch {fname}\n")

    os.system(f"chmod +x {submit_sh}")

    print(f"Submission script written to {submit_sh}")
    print("To submit jobs, run:")
    print(f"bash {submit_sh}")
    
    
if __name__ == "__main__":
    samples = {
        "pi+": [(10, 10.1), (20, 20.1), (30, 30.1), (50, 50.1), 
                (100, 100.1)],
        "e+": [(10.0, 10.1), (20.0, 20.1), (30.0, 30.1), (50.0, 50.1),
                (100.0, 100.1)],
    }
    
    for particle, energies in samples.items():
        for energy_min, energy_max in energies:
            print(f"Generating jobs for {particle} with energy range {energy_min} - {energy_max}")
            generate_submission_script(particle, energy_min, energy_max)
        
