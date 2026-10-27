# SANCHO

**S**pinless **A**nalogue to **N**onequilibrium **C**orrelated **H**opping

Current–voltage characteristics of a tight-binding chain with a Hubbard
interaction between two leads, from self-consistent non-equilibrium Green's
functions. The interaction is treated to 0th, 1st, 2nd or 3rd order in U,
with every order built from factorised particle–hole / particle–particle
bubbles. Spin is not resolved (paramagnetic chain, factor 2 in the current).

SANCHO lives inside the [QUICHOTTE](../) repository but is a **complete,
independent code**: it has its own Makefile, tests, examples and manual, and
builds without anything outside this folder. To give it its own repository
later: `git subtree split --prefix SANCHO -b sancho-only`.

📖 **User manual:** [`docs/manual/sancho_manual.pdf`](docs/manual/sancho_manual.pdf)

---

## Quick start

```bash
cd QUICHOTTE/SANCHO
make                 # build bin/sancho.x   (gfortran + BLAS/LAPACK)
make test            # 6 regression tests, ~10 s; all should PASS

mkdir -p ~/runs/chain && cd ~/runs/chain
cp ~/QUICHOTTE/SANCHO/examples/chain_N20_order2/input.dat .     # then edit it
OMP_NUM_THREADS=4 ~/QUICHOTTE/SANCHO/bin/sancho.x
```

| platform | build | install / load first |
|---|---|---|
| macOS | `make` | `brew install gcc openblas` |
| Linux | `make` | `apt install gfortran libopenblas-dev` |
| KCL CREATE | `make PLATFORM=create` | `source hpc/create/modules.sh` |
| ARCHER2 | `make PLATFORM=archer2` | `source hpc/archer2/modules.sh` |

## Layout

```
src/            SANCHO source (one module per file)
common/         shared with QUICHOTTE -- must stay IDENTICAL to ../common/
make/           compiler and library settings per platform
examples/       annotated input files
tests/          regression tests (Fortran runner)
tools/          dIdV (differential conductance), stm_junction (experimental geometry)
hpc/            job scripts: create/, young/, archer2/ + archer2/sweep/ (parameter scans)
docs/manual/    user manual (LyX source + PDF)
```

| file | module | contents |
|---|---|---|
| `src/sancho.f90` | *(program)* | set-up, grid, voltage loop, output |
| `src/parameters.f90` | `Parameters` | run parameters, sizes, unit numbers |
| `src/hamiltonian.f90` | `Hamiltonian` | chain H (+ experimental xyz H) |
| `src/green_functions.f90` | `GreenFunctions` | GF storage, non-interacting GFs, G<(t=0) |
| `src/self_energies.f90` | `SelfEnergies` | bubbles χ, φ; 2nd- and 3rd-order Σ |
| `src/scf.f90` | `SCF` | SCF driver, residual, Dyson/Keldysh (incl. Hartree term) |
| `src/observables.f90` | `Observables` | current, populations, spectral function |
| `src/input_output.f90` | `InputOutput` | reads `input.dat`, per-site U, printout |

## Input (`input.dat`)

| key | default | meaning |
|---|---|---|
| `V`, `Vf`, `delv` | required | first voltage, last voltage, step |
| `T`, `mu` | required | temperature [K], chemical potential [eV] |
| `Nsites` | required | number of sites (old name `Natoms`) |
| `E_CC`, `t_hop`, `U_int` | required | on-site energy, hopping, Hubbard U [eV] |
| `hubbard_file` | – | one U per site (overrides `U_int`) |
| `order` | required | 0, 1, 2 or 3 |
| `GammaL`, `GammaR` | 1.0, 1.0 | lead couplings to the end sites [eV] |
| `dw`, `up`, `delta` | required | energy grid margins and spacing [eV] |
| `alpha_mix`, `nhist`, `maxscf`, `tolf` | 0.3, 8, 200, 1e-6 | Broyden mixing and convergence |
| `continuation` | T | start each voltage from the previous converged one |

Missing required keys stop the run with a list; unknown keys give a warning;
`restart` and `ver` from old files are accepted and ignored.

## Output

| file | contents |
|---|---|
| `Volt_Current_<order>.dat` | V, J, converged (1/0) |
| `populations.dat` | V, site, n per spin, n total |
| `spectral_function.dat` | V, ω, A_ii(ω) per site; one block per voltage |
| `scf_V_<V>.dat` | SCF history at that voltage |
| `info_Hamiltonian.dat`, `runtime_datasheet.dat` | H and spectrum; code version, wall time |

**Exit codes:** 0 all converged · 1 bad input · 3 some voltages did not converge.

## Cluster runs

* One run: `hpc/create/job.slurm`, `hpc/young/job.sge`, `hpc/archer2/job.slurm`,
  submitted from the run directory.
* Parameter scans on ARCHER2: `hpc/archer2/sweep/` (see its README). It packs
  several runs per node, tracks each run's state and collects every I–V curve.

## Changing the code

`make test` before every commit. After a deliberate physics change:
`make update-references` and commit the references with the change.
The files in `common/` are shared with QUICHOTTE. Change them in both places;
`make check-common` (from the repository root, also run by CI) catches differences.
