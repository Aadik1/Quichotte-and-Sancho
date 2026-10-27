#!/bin/bash --login
# =============================================================================
# hpc/archer2/job.slurm  --  ONE run of SANCHO on ARCHER2
#
#   cd /work/<project>/<project>/$USER/my_run      # must be under /work
#   cp /work/<project>/<project>/$USER/QUICHOTTE/SANCHO/examples/chain_N20_order2/input.dat .  # then edit it
#   sbatch /work/<project>/<project>/$USER/QUICHOTTE/SANCHO/hpc/archer2/job.slurm
#
# For many runs (parameter scans) use hpc/archer2/sweep/ instead: it packs
# several runs onto each 128-core node and tracks their status.
# =============================================================================
#SBATCH --job-name=sancho
#SBATCH --account=e05-nano-lev        # <-- your budget code
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --output=%x.%j.out

EXE=${EXE:-/work/e05/e05/$USER/QUICHOTTE/SANCHO/bin/sancho.x}   # <-- path to your build

module restore
module load PrgEnv-gnu

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M

cd "$SLURM_SUBMIT_DIR"
[[ -f input.dat ]] || { echo "no input.dat in $PWD"; exit 1; }
srun --distribution=block:block --hint=nomultithread "$EXE"
rc=$?
[[ $rc == 3 ]] && echo "Finished, but some voltages did not converge (column 3 = 0)."
exit $rc
