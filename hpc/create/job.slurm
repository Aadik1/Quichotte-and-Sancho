#!/bin/bash -l
# =============================================================================
# hpc/create/job.slurm  --  one run of QUICHOTTE on KCL CREATE, safe on every node
#
#   mkdir -p /scratch/users/$USER/my_run && cd $_
#   cp <QUICHOTTE>/examples/helix_4x2_order2/input.dat .     # then edit it
#   sbatch <QUICHOTTE>/hpc/create/job.slurm
#
# The cpu partition mixes AMD zen2, zen3 and zen5 nodes, and the spack modules are
# built per CPU type. A binary built on the (Intel) login node links Intel-tuned
# libraries and dies on the compute nodes with "Illegal instruction". So this job
# builds QUICHOTTE ON THE NODE IT RUNS ON, with that node's modules:
#   build_<cpu>/ and bin/quichotte_<cpu>.x, e.g. build_znver2, bin/quichotte_znver2.x
# The first job on each CPU type compiles (about a minute); later jobs on the same
# type reuse the binary (make only rebuilds what changed). flock stops two jobs
# that start at once from building into the same folder.
#
# Set QROOT below once: the QUICHOTTE folder (the one with the Makefile).
# =============================================================================
#SBATCH --job-name=quichotte
#SBATCH --account=kcl
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --output=%x.%j.out
#SBATCH --export=NONE

QROOT=/scratch/users/k19020583/QUICHOTTE     # <-- your QUICHOTTE folder

#...start clean (--export=NONE: no login-node module paths), load this node's modules
source "$QROOT/hpc/create/modules.sh"

#...CPU type of this node, as gcc names it (znver2, znver3, znver5, ...)
CPU=$(gcc -march=native -Q --help=target 2>/dev/null | awk '/^ +-march=/{print $2; exit}')
CPU=${CPU:-unknown}
EXE=$QROOT/bin/quichotte_$CPU.x
echo "Node $(hostname), CPU $CPU, job $SLURM_JOB_ID, $(date)"
echo "MODULEPATH: $MODULEPATH"

#...build for this CPU type (no-op if already up to date)
( cd "$QROOT" && flock "build_$CPU.lock" make PLATFORM=create BUILD="build_$CPU" EXE="bin/quichotte_$CPU.x" ) \
   > "build_$CPU.log" 2>&1 || { echo "ERROR: build failed, see build_$CPU.log"; exit 1; }

#...threads: one per physical core, pinned; BLAS serial (it is called inside OpenMP loops)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M

cd "$SLURM_SUBMIT_DIR"
[[ -f input.dat || -f inputSOC.dat ]] || { echo "ERROR: no input.dat in $PWD"; exit 1; }
if ldd "$EXE" | grep -q "not found"; then
   echo "ERROR: missing shared libraries:"; ldd "$EXE" | grep "not found"; exit 1
fi

echo "EXE = $EXE"
"$EXE"
rc=$?
[[ $rc == 3 ]] && echo "Finished, but some voltages did not converge (column 5 = 0)."
echo "Exit code $rc, $(date)"
exit $rc
