#!/bin/bash -l
# =============================================================================
# hpc/create/job.slurm  --  one run of SANCHO on KCL CREATE
#
#   mkdir -p /scratch/users/$USER/my_run && cd $_
#   cp ~/QUICHOTTE/SANCHO/examples/chain_N20_order2/input.dat .     # then edit it
#   sbatch ~/QUICHOTTE/SANCHO/hpc/create/job.slurm
#
# The job runs in the directory you submit from. Set EXE below once.
# =============================================================================
#SBATCH --job-name=sancho
#SBATCH --account=kcl
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --output=%x.%j.out

EXE=${EXE:-$HOME/QUICHOTTE/SANCHO/bin/sancho.x}   # <-- path to your build

module purge
module load gcc/13.2.0
module load openblas/0.3.24-gcc-13.2.0-multithreading-openmp

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=1     # BLAS calls are inside OpenMP loops: keep them serial
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M

cd "$SLURM_SUBMIT_DIR"
[[ -f input.dat ]] || { echo "no input.dat in $PWD"; exit 1; }
"$EXE"
rc=$?
[[ $rc == 3 ]] && echo "Finished, but some voltages did not converge (column 3 = 0)."
exit $rc
