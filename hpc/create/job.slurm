#!/bin/bash -l
# =============================================================================
# hpc/create/job.slurm  --  one run of QUICHOTTE on KCL CREATE
#
#   1. build once on a login node:
#        source hpc/create/modules.sh && make PLATFORM=create
#   2. make a run folder with an input.dat in it, and submit FROM that folder:
#        mkdir -p /scratch/users/$USER/my_run && cd $_
#        cp <QUICHOTTE>/examples/helix_4x2_order2/input.dat .     # then edit it
#        sbatch <QUICHOTTE>/hpc/create/job.slurm
#
#   The job runs in the folder you submit from; output goes there too.
#   Binary: EXE below, or override per submission without editing this file:
#        sbatch --export=ALL,EXE=/path/to/quichotte.x <QUICHOTTE>/hpc/create/job.slurm
#   The modules come from hpc/create/modules.sh next to that binary's repo,
#   so building and running always use the same environment.
#
#   Cores: QUICHOTTE is OpenMP only (one node). Set --cpus-per-task to the
#   physical cores of a node:  sinfo -p cpu -o "%c" | sort -u
# =============================================================================
#SBATCH --job-name=quichotte
#SBATCH --account=kcl
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --output=%x.%j.out

EXE=${EXE:-$HOME/QUICHOTTE/bin/quichotte.x}     # <-- path to your build (set once)

#...repo root = the folder above bin/ ; its modules.sh sets up gcc, openblas, fftw
QROOT=$(dirname "$(dirname "$EXE")")
[[ -x "$EXE" ]] || { echo "ERROR: no executable at EXE=$EXE"; exit 1; }
[[ -f "$QROOT/hpc/create/modules.sh" ]] || { echo "ERROR: no $QROOT/hpc/create/modules.sh"; exit 1; }
source "$QROOT/hpc/create/modules.sh"

#...threads: one per physical core, pinned; BLAS serial (it is called inside OpenMP loops)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M

cd "$SLURM_SUBMIT_DIR"
[[ -f input.dat || -f inputSOC.dat ]] || { echo "ERROR: no input.dat in $PWD"; exit 1; }

#...fail early (not after the queue wait) if a shared library cannot be found
if ldd "$EXE" | grep -q "not found"; then
   echo "ERROR: missing shared libraries:"; ldd "$EXE" | grep "not found"; exit 1
fi

echo "Job $SLURM_JOB_ID on $(hostname), $OMP_NUM_THREADS threads, $(date)"
echo "EXE = $EXE"
echo "Run folder = $PWD"

"$EXE"
rc=$?
[[ $rc == 3 ]] && echo "Finished, but some voltages did not converge (column 5 = 0)."
echo "Exit code $rc, $(date)"
exit $rc
