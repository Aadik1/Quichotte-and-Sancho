# Parameter sweeps on ARCHER2

Many runs of SANCHO, packed onto whole 128-core nodes as one Slurm
job array. The runs are listed in `manifest.txt`; everything else is
generated from it.

- **Runs are packed onto nodes.** ARCHER2 charges whole nodes, so each array
  task takes one node and runs `RUNS_PER_NODE` runs side by side
  (default 8 runs × 16 threads).
- **Everything must live under `/work`.** Compute nodes cannot see `/home`.
- **The executable is the repo build**, `SANCHO/bin/sancho.x` (set in `sweep.conf`).

| file | what it is |
|---|---|
| `sweep.conf` | budget code, resources, runs per node |
| `defaults.dat` | baseline `input.dat` (every key, including the Broyden settings) |
| `manifest.txt` | the runs: columns override keys of `defaults.dat` |
| `generate.sh` | makes `runs/<name>/input.dat` for every run, writes `runs.list` |
| `submit.sh` | submits runs as one Slurm job array, one node per array task |
| `job_array.sh` | what each node runs (called by `submit.sh`; don't sbatch it) |
| `status.sh` | one line per run: state, job, walltime, voltages done |
| `collect.sh` | all I–V curves into one table, `results_IV.dat` |

## One-time setup

```bash
cd /work/<project>/<project>/$USER
git clone https://github.com/Aadik1/QUICHOTTE.git && cd QUICHOTTE/SANCHO
source hpc/archer2/modules.sh
make PLATFORM=archer2 && make test
cd hpc/archer2/sweep        # or copy this folder somewhere else under /work and fix BIN
```

In `sweep.conf`, set `ACCOUNT` to your budget code. Check `TIME`, and
`RUNS_PER_NODE` × `CPUS` (at most 128).

## Workflow

```bash
nano manifest.txt        # say what to run
./generate.sh            # create run directories
DRY=1 ./submit.sh        # see what would be submitted
./submit.sh              # submit everything not yet submitted
./status.sh              # check progress
./submit.sh failed       # resubmit FAILED, TIMEOUT and vanished runs
./collect.sh             # build results_IV.dat
```

You can pilot a single run by name, e.g. `./submit.sh N20_Ut8.0_o3`. It gets
a node to itself, so it costs a full node for its walltime.

## Run states

| state | meaning |
|---|---|
| NEW | never submitted |
| QUEUED / RUNNING | in a live Slurm job |
| DONE | finished, every voltage converged |
| UNCONV | finished, but some voltages did not converge (exit code 3) |
| FAILED | crashed or other non-zero exit (exit code 1 = bad input.dat, see `run.log`) |
| TIMEOUT | hit the time limit (Slurm warns the script 5 min before) |
| STALE | marked queued/running but the job is gone |

`UNCONV` runs are **not** resubmitted by `./submit.sh failed`, because they
would give the same answer again. Loosen the mixing instead: add e.g.
`alpha_mix` or `maxscf` columns to their manifest line, run
`FORCE=1 ./generate.sh`, then submit them by name.

Each run's convergence history per voltage is in `runs/<name>/scf_V_<V>.dat`.

## The results table

`results_IV.dat` has one row per (run, voltage):

```
# N  U_int  t_hop  U/t  order  delta  V  I  conv  state  run
```

`conv` is 1 if the SCF converged at that voltage and 0 if not. The table
includes DONE and UNCONV runs, so filter on `conv`, e.g.

```bash
awk '$1==12 && $4==2 && $5==3 && $9==1 {print $7, $8}' results_IV.dat
```
