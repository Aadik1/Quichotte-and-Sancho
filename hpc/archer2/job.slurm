#!/bin/bash --login
# =============================================================================
# hpc/archer2/job.slurm  --  ONE run of QUICHOTTE on ARCHER2
#
#   cd /work/<project>/<project>/$USER/my_run      # must be under /work
#   cp /work/<project>/<project>/$USER/QUICHOTTE/examples/helix_4x2_order2/input.dat .  # then edit it
#   sbatch /work/<project>/<project>/$USER/QUICHOTTE/hpc/archer2/job.slurm
#
# =============================================================================
#SBATCH --job-name=quichotte
#SBATCH --account=e05-nano-lev        # <-- your budget code
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --output=%x.%j.out

EXE=${EXE:-/work/e05/e05/$USER/QUICHOTTE/bin/quichotte.x}   # <-- path to your build

module restore
module load PrgEnv-gnu
module load cray-fftw

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M

cd "$SLURM_SUBMIT_DIR"
[[ -f input.dat || -f inputSOC.dat ]] || { echo "no input.dat in $PWD"; exit 1; }
srun --distribution=block:block --hint=nomultithread "$EXE"
rc=$?
[[ $rc == 3 ]] && echo "Finished, but some voltages did not converge (column 5 = 0)."
exit $rc
