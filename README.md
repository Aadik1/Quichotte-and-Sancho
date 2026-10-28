# QUICHOTTE

**QU**antum **I**nteracting **CH**ains: **O**ut-of-equilibrium **T**ranspor**T**

Spin-resolved current–voltage characteristics of a helical molecule with
spin–orbit coupling and a Hubbard interaction, between two leads, from
self-consistent non-equilibrium Green's functions (NEGF). The interaction is
treated to 0th, 1st (Hartree–Fock) or 2nd order in U, the 2nd order by FFT.
Written for studying chirality-induced spin selectivity (CISS).

Its companion **[SANCHO](SANCHO/)** (**S**pinless **A**nalogue to
**N**onequilibrium **C**orrelated **H**opping) is a separate, self-contained
code in `SANCHO/`: a spinless chain to 3rd order in U, using factorised bubbles.

📖 **User manual:** [`docs/manual/quichotte_manual.pdf`](docs/manual/quichotte_manual.pdf)
(model, every equation the code solves, numerics, input/output reference, workflow).

---

## Quick start

```bash
make                 # build bin/quichotte.x   (gfortran + BLAS/LAPACK + FFTW3)
make test            # 9 regression tests (7 physics, 2 input errors); all should PASS

mkdir -p ~/runs/first && cd ~/runs/first
cp ~/QUICHOTTE/examples/helix_4x2_order2/input.dat .     # then edit it
OMP_NUM_THREADS=4 ~/QUICHOTTE/bin/quichotte.x
```

Run from a directory **outside the repo**, one directory per run.

| platform | build | install / load first |
|---|---|---|
| macOS | `make` | `brew install gcc openblas fftw` |
| Linux | `make` | `apt install gfortran libopenblas-dev libfftw3-dev` |
| KCL CREATE | `make PLATFORM=create` | `source hpc/create/modules.sh` |
| ARCHER2 | `make PLATFORM=archer2` | `source hpc/archer2/modules.sh` |

`make help` lists every target: `tools`, `debug` (bounds/NaN checking),
`update-references`, `check-common`, `sancho`, `test-all`, `clean`.

## Repository layout

```
src/            QUICHOTTE source (one module per file, see below)
common/         code shared with SANCHO -- IDENTICAL to SANCHO/common/ (make check-common)
make/           compiler and library settings per platform
examples/       annotated input files to copy: helix_4x2_order2 (helix), helix_from_xyz
                (positions from a file), custom_model (your own Hamiltonian)
tests/          regression tests (cases.txt, cases/<name>/, Fortran runner)
tools/          helix_xyz (geometry -> helix.xyz), dIdV (differential conductance)
hpc/            job scripts and module lists: create/, young/, archer2/
docs/manual/    user manual (LyX source + PDF)
SANCHO/         the SANCHO code, same layout, builds and tests on its own
```

| file | module | contents |
|---|---|---|
| `src/quichotte.f90` | *(program)* | set-up, energy grid, voltage loop, output files |
| `src/parameters.f90` | `Parameters` | run parameters, sizes, Pauli matrices, unit numbers |
| `src/hamiltonian.f90` | `Hamiltonian` | builds H: helix, from `.xyz`, or your own model (`build_custom_hamiltonian`) |
| `src/green_functions.f90` | `GreenFunctions` | storage of the GFs, non-interacting GFs, G<(t=0) |
| `src/self_energies.f90` | `SelfEnergies` | first-order Σ, second-order Σ by FFT |
| `src/scf.f90` | `SCF` | SCF driver, residual, Dyson/Keldysh at one energy |
| `src/observables.f90` | `Observables` | current, spin texture, spectral function |
| `src/input_output.f90` | `InputOutput` | reads `input.dat`, per-site U, Hamiltonian printout |
| `common/constants.f90` | `Constants` | `dp`, ħ, k_B, π, i |
| `common/linear_algebra.f90` | `LinearAlgebra` | LAPACK wrappers, packing of G for Broyden |
| `common/input_parser.f90` | `InputParser` | `key = value` reader with checks |
| `common/broyden_mixing.f90` | `BroydenMixing` | modified Broyden solver (Johnson 1988) |
| `common/scf_checks.f90` | `ScfChecks` | occupation and spectral-weight sanity checks |

## Input (`input.dat`)

`key = value` per line, any order, `#` starts a comment. Missing required keys
are listed and stop the run; unknown keys (typos) give a warning. The old file
name `inputSOC.dat` and old key names (`Hubbard`, `pulay`, `GammaR_both`) still work.

| key | default | meaning |
|---|---|---|
| `V`, `Vf`, `delv` | required | first voltage, last voltage, step |
| `T`, `mu` | required | temperature [K], chemical potential [eV] |
| `hamiltonian` | helix | model: `helix`, `xyz` (positions from `xyz_file`) or `custom` (your own, see below) |
| `xyz_file` | molecule.xyz | site positions, for `hamiltonian = xyz` |
| `N_sites` | from model | number of sites; required for `custom` |
| `N_ions`, `N_turns` | required for helix | sites per turn, turns |
| `hel_length`, `hel_radius`, `hand` | required for helix | helix geometry; `hand` = ±1 |
| `E_CC`, `t_hop` | required | on-site energy, hopping [eV] |
| `lamb` | required for helix, xyz | SOC strength [eV] |
| `U_int` | required | Hubbard U [eV] |
| `hubbard_file` | – | one U per site (overrides `U_int`) |
| `order` | required | 0, 1 or 2 |
| `GammaL_up`, `GammaL_down`, `GammaR` | required | lead couplings [eV] |
| `dw`, `up`, `delta` | required | energy grid margins and spacing [eV] |
| `alpha_mix`, `nhist`, `maxscf`, `tolf` | 0.5, 16, 300, 1e-6 | Broyden mixing and convergence |
| `continuation` | T | start each voltage from the previous converged one |

## Output

| file | contents |
|---|---|
| `Volt_Current_<order>.dat` | V, J_up, J_down, J_total, converged (1/0) |
| `populations.dat` | V, site, ⟨Sx⟩, ⟨Sy⟩, ⟨Sz⟩, n_down, n_up, n_total |
| `spectral_function.dat` | V, ω, A_aa(ω) per spin-orbital; one block per voltage |
| `scf_V_<V>.dat` | SCF history at that voltage |
| `info_Hamiltonian.dat` | H, eigenvalues, eigenvectors |
| `runtime_datasheet.dat` | code version (`git describe`), wall time |

**Exit codes:** 0 all voltages converged · 1 bad input · 3 some voltages did not
converge (they are flagged, the sweep carries on).

## Cluster runs

```bash
cd /scratch/users/$USER/my_run && cp ~/QUICHOTTE/examples/helix_4x2_order2/input.dat .
sbatch ~/QUICHOTTE/hpc/create/job.slurm        # or hpc/young/job.sge, hpc/archer2/job.slurm
```

## Your own Hamiltonian

1. Edit the body of `build_custom_hamiltonian` at the end of `src/hamiltonian.f90`
   (as shipped it is a plain chain, so it runs straight away).
2. In `input.dat`: `hamiltonian = custom` and `N_sites = <number of sites>`
   (ready-made: `examples/custom_model/input.dat`).
   The helix keys are then optional.
3. `make`, check `info_Hamiltonian.dat`, run `make test`.

To keep several models side by side, give each its own routine in
`src/hamiltonian.f90` and add one line to `build_hamiltonian`, e.g.
`case ('ssh'); call build_ssh_hamiltonian()`; then use `hamiltonian = ssh`.

Rules: spin-orbitals ordered `2i-1` = (site i, up), `2i` = (site i, down); fill
the upper triangle and call `make_hermitian()`; energies in eV; no Hubbard terms
in H; leads attach to site 1 and site `N_sites`. The manual's section *Using your
own Hamiltonian* explains each rule. A different molecular shape with the same
terms as the helix needs no code at all: `hamiltonian = xyz`, `xyz_file = ...`.

## Changing the code

Run `make test` before every commit (CI runs it on every push, together with
SANCHO's tests and `check-common`). After a *deliberate* physics change, check
the numbers, run `make update-references`, and commit the new references with
the change. The manual's "Code structure" section explains how to add an input
key, a geometry or a self-energy. Everything (code, tools, test runner) is Fortran;
the only shell code is the cluster job scripts and SANCHO's ARCHER2 sweep scripts.
