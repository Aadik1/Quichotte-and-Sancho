# QUICHOTTE

**QU**antum **I**nteracting **CH**ains: **O**ut-of-equilibrium **T**ranspor**T**

Spin-resolved current–voltage characteristics of a helical molecule with
spin–orbit coupling and a Hubbard interaction, between two leads, from
self-consistent non-equilibrium Green's functions (NEGF). The interaction is
treated to 0th, 1st (Hartree–Fock) or 2nd order in U, the 2nd order by FFT.
Written for studying chirality-induced spin selectivity (CISS).

The Hamiltonian can be a model helix, or a **multi-orbital Hamiltonian from DFTB+**
(s, p_y, p_z, p_x per atom, with spin–orbit coupling); every spatial orbital is
then a Hubbard site. Only the interacting and contact orbitals are stored, so
memory grows with the interacting region, not with the molecule.

Its companion **[SANCHO](SANCHO/)** (**S**pinless **A**nalogue to
**N**onequilibrium **C**orrelated **H**opping) is a separate, self-contained
code in `SANCHO/`: a spinless chain to 3rd order in U, using factorised bubbles.

📖 **User manual:** [`docs/manual/quichotte_manual.pdf`](docs/manual/quichotte_manual.pdf)
(model, every equation the code solves, numerics, input/output reference, workflow).

---

## Quick start

```bash
make                 # build bin/quichotte.x   (gfortran + BLAS/LAPACK + FFTW3)
make test            # 14 regression tests (10 physics, 4 input errors); all should PASS

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
`update-references`, `check-common`, `sancho`, `test-all`, `clean`, `archive`.

## Repository layout

```
src/            QUICHOTTE source (one module per file, see below)
common/         code shared with SANCHO -- IDENTICAL to SANCHO/common/ (make check-common)
make/           compiler and library settings per platform
examples/       annotated input files to copy: helix_4x2_order2 (helix), helix_from_xyz
                (positions from a file), helicene_dftb_hs (DFTB+ H and S, recommended),
                helicene_dftb (soc_dftb block file), custom_model (your own Hamiltonian)
tests/          regression tests (cases.txt, cases/<name>/, Fortran runner)
tools/          helix_xyz (geometry -> helix.xyz), dIdV (differential conductance)
hpc/            job scripts and module lists: create/, young/, archer2/
docs/manual/    user manual (LyX source + PDF)
SANCHO/         the SANCHO code, same layout, builds and tests on its own
archive.sh      makes dated source archives of both codes (make archive)
archive/        where they go (contents not committed)
```

| file | module | contents |
|---|---|---|
| `src/quichotte.f90` | *(program)* | set-up, energy grid, voltage loop, output files |
| `src/parameters.f90` | `Parameters` | run parameters, sizes, Pauli matrices, unit numbers |
| `src/hamiltonian.f90` | `Hamiltonian` | builds H: helix, from `.xyz`, from DFTB+ (`hamsqr1.dat` + `oversqr.dat`, or a block file), or your own model (`build_custom_hamiltonian`) |
| `src/green_functions.f90` | `GreenFunctions` | active subspace, storage of the GFs, G<(t=0) |
| `src/self_energies.f90` | `SelfEnergies` | first-order Σ, second-order Σ by FFT |
| `src/scf.f90` | `SCF` | SCF driver, residuals, final sweep, Dyson/Keldysh at one energy |
| `src/observables.f90` | `Observables` | current, spin texture (per orbital and per atom, with ⟨L⟩), spectral function |
| `src/input_output.f90` | `InputOutput` | reads `input.dat`, atom/orbital map, lead couplings, per-site U, Hamiltonian printout |
| `common/constants.f90` | `Constants` | `dp`, ħ, k_B, π, i |
| `common/linear_algebra.f90` | `LinearAlgebra` | LAPACK wrappers, packing of G for Broyden |
| `common/input_parser.f90` | `InputParser` | `key = value` reader with checks |
| `common/broyden_mixing.f90` | `BroydenMixing` | modified Broyden solver (Johnson 1988), optional single-precision history |
| `common/scf_checks.f90` | `ScfChecks` | occupation and spectral-weight sanity checks |

## Input (`input.dat`)

`key = value` per line, any order, `#` starts a comment. Missing required keys
are listed and stop the run; unknown keys (typos) give a warning. The old file
name `inputSOC.dat` and old key names (`Hubbard`, `pulay`, `GammaR_both`) still work.

| key | default | meaning |
|---|---|---|
| `V`, `Vf`, `delv` | required | first voltage, last voltage, step |
| `T`, `mu` | required | temperature [K], chemical potential [eV] |
| `hamiltonian` | helix | model: `helix`, `xyz` (positions from `xyz_file`), `dftb_hs` (DFTB+ H and S + SOC), `dftb` (block file) or `custom` (your own); see [DFTB+ Hamiltonians](#dftb-hamiltonians) |
| `xyz_file` | molecule.xyz | site positions (`xyz`); element of each atom (`dftb_hs`) |
| `ham_file` | hamsqr1.dat (`dftb_hs`), hamiltonian_blocks_fortran.txt (`dftb`) | DFTB+ Hamiltonian file |
| `overlap_file` | oversqr.dat | DFTB+ overlap, for `dftb_hs` |
| `max_angular_momentum` | required for `dftb_hs` | shells per element, as in DFTB+: `C:p H:s` |
| `soc_xi` | none (no SOC) | p-shell spin–orbit constant per element [eV], e.g. `C:0.00849` |
| `soc_form` | dual | `dual` (½(SΛ+ΛS), as DFTB+) or `onsite` (Λ) |
| `write_hamiltonian` | F | write H in full precision to `hamiltonian_matrix.dat` |
| `N_sites` | from model | number of sites; required for `custom` |
| `N_ions`, `N_turns` | required for helix | sites per turn, turns |
| `hel_length`, `hel_radius`, `hand` | required for helix | helix geometry; `hand` = ±1 |
| `E_CC`, `t_hop` | required (not for dftb, dftb_hs) | on-site energy, hopping [eV] |
| `lamb` | required for helix, xyz | SOC strength [eV] |
| `U_int` | required | Hubbard U [eV] |
| `hubbard_file` | – | one U per site, i.e. per spatial orbital for DFTB (overrides `U_int`); U = 0 sites are not stored |
| `order` | required | 0, 1 or 2 |
| `GammaL_up`, `GammaL_down`, `GammaR` | required | lead couplings [eV], on every orbital of the contact atoms |
| `lead_left_atoms`, `lead_right_atoms` | first, last atom | contact atoms, e.g. `1 2 3` or `24:26` |
| `dw`, `up`, `delta` | required | energy grid margins and spacing [eV] |
| `alpha_mix`, `nhist`, `maxscf`, `tolf` | 0.5, 16, 300, 1e-6 | Broyden mixing and convergence |
| `tolf_order1` | 1e-8 | convergence of order 1 (mixes G<(t=0), so it can afford a tight tolerance) |
| `broyden_single_precision` | T | keep the Broyden history in single precision (halves its memory; iterate stays double) |
| `continuation` | T | start each voltage from the previous converged one |

## Output

| file | contents |
|---|---|
| `Volt_Current_<order>.dat` | V, J_up, J_down, J_total, converged (1/0) |
| `populations.dat` | V, site, ⟨Sx⟩, ⟨Sy⟩, ⟨Sz⟩, n_down, n_up, n_total (site = orbital for DFTB) |
| `populations_atoms.dat` | V, atom, ⟨Sx⟩, ⟨Sy⟩, ⟨Sz⟩, n_down, n_up, n_total, ⟨Lx⟩, ⟨Ly⟩, ⟨Lz⟩ (only multi-orbital runs; spin and L in units of ħ) |
| `spectral_function.dat` | V, ω, A_aa(ω) per spin-orbital; one block per voltage |
| `scf_V_<V>.dat` | SCF history at that voltage |
| `info_Hamiltonian.dat` | H, eigenvalues, eigenvectors |
| `hamiltonian_matrix.dat` | `i j Re Im` of H in full precision (`write_hamiltonian = T`) |
| `runtime_datasheet.dat` | code version (`git describe`), wall time |

All observables come from one final sweep over the energy grid with the
converged self-energies.

**Exit codes:** 0 all voltages converged · 1 bad input · 3 some voltages did not
converge (they are flagged, the sweep carries on).

## Cluster runs

```bash
cd /scratch/users/$USER/my_run && cp ~/QUICHOTTE/examples/helix_4x2_order2/input.dat .
sbatch ~/QUICHOTTE/hpc/create/job.slurm        # or hpc/young/job.sge, hpc/archer2/job.slurm
```

## Memory: the active subspace

Only the **active** states are stored: the interacting orbitals (U ≠ 0) and the
contact orbitals. The full G^r = [ω − H + iΓ/2 − Σ]⁻¹ exists only one energy at
a time, per thread. Order 2 stores G^r, G^< and Σ on the active block
(∝ n_A² n_ω); order 1 mixes only the 2×2 blocks of G^<(t=0) on the interacting
sites (a vector of 8 numbers per site); order 0 stores nothing but the per-energy
observables. The run prints its estimate (`Memory: stored arrays ~ … GB`). A
`hubbard_file` with U on part of the molecule therefore cuts memory and time.

Measured (24-site helix, δ = 0.01 eV): order 2 with U everywhere 0.85 → 0.48 GB,
U on 8 of 24 sites 0.85 → 0.09 GB and 3.5× faster; order 1 0.77 → 0.009 GB;
order 0 0.14 → 0.008 GB. Order-2 SCF histories are bit-identical to the old code
(with `broyden_single_precision = F`).

## DFTB+ Hamiltonians

A **site** is one spatial orbital (C: s, p_y, p_z, p_x in DFTB+ order; H: s), so
`populations.dat`, `spectral_function.dat` and `hubbard_file` are per orbital and
the Hubbard U is intra-orbital. `populations_atoms.dat` sums them per atom and adds
the orbital moment ⟨L⟩ of the p shell. The leads couple to every orbital of
`lead_left_atoms` / `lead_right_atoms`. On reading, the code prints how strongly H
breaks time reversal: in each 2×2 spin block h₀·1 + **h**·σ, Im h₀ and Re **h** must be
0 for a non-magnetic molecule (Im **h** is the spin–orbit coupling).

### `hamiltonian = dftb_hs` (recommended)

1. Run DFTB+ **without** spin–orbit coupling and with `Scc = No` (U is added by
   QUICHOTTE; an SCC Hamiltonian would count the interaction twice), with
   `Options { WriteHS = Yes }`. DFTB+ writes `hamsqr1.dat` and `oversqr.dat`
   (Hartree, Γ point) and stops.
2. Copy them and the geometry (`xyz_file`, for the element of each atom) to the run
   directory; set `max_angular_momentum` as in DFTB+ and `soc_xi` per element in eV.
3. QUICHOTTE builds, in eV,

   H = (S^-½ ⊗ 1) [ H₀ ⊗ 1 + ½(SΛ + ΛS) ] (S^-½ ⊗ 1),   Λ = Σ_atoms ξ **L**·**S** on the p shell,

   the DFTB+ "dual" SOC form (`soc_form = onsite` uses Λ instead). H is exactly
   time-reversal symmetric. The smallest overlap eigenvalue is printed (stop if ≤ 1e-8).

**ξ for carbon:** 0.312 mHa = **8.49 meV** (Jha & Heine). Note that 0.000312 entered
in eV is 27× too small.

Validated on hexahelicene against an independent construction (5e-14 eV), and the
eigenvalues agree with the soc_dftb Hamiltonian to 2e-6 eV.
Example: `examples/helicene_dftb_hs/input.dat`.

### `hamiltonian = dftb` (soc_dftb block file)

Reads the orthogonal spin–orbit Hamiltonian written by the soc_dftb pipeline
(`export_fortran_blocks` in `dftb_soc_hamiltonian.py`): blocks `BLOCK n | Atom_I: a
Atom_J: b | Dist: d | Dim: (2 n_a)x(2 n_b)` followed by one `re im` line per element,
rows over (orbital, spin) of each atom with spin fastest, energies in eV. The number of
atoms and the orbitals per atom come from the block headers; every atom needs its
on-site block. That H is rebuilt from rounded eigen-data and breaks time reversal by
~1e-5 eV (printed on reading), a sizeable fraction of the SOC; prefer `dftb_hs`.
Example: `examples/helicene_dftb/input.dat`.

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
in H; leads attach to atoms `lead_left_atoms` / `lead_right_atoms` (default: site 1
and site `N_sites`). The manual's section *Using your
own Hamiltonian* explains each rule. A different molecular shape with the same
terms as the helix needs no code at all: `hamiltonian = xyz`, `xyz_file = ...`.

## Changing the code

Run `make test` before every commit (CI runs it on every push, together with
SANCHO's tests and `check-common`). After a *deliberate* physics change, check
the numbers, run `make update-references`, and commit the new references with
the change. The manual's "Code structure" section explains how to add an input
key, a geometry or a self-energy. Everything (code, tools, test runner) is Fortran;
the only shell code is the cluster job scripts, SANCHO's ARCHER2 sweep scripts
and `archive.sh`.

## Archiving

```bash
make archive              # or ./archive.sh [all|quichotte|sancho]
```

writes `archive/QUICHOTTE_<YYYY-MM-DD>.tar.gz` (everything except `SANCHO/`) and
`archive/SANCHO_<YYYY-MM-DD>.tar.gz` (the `SANCHO/` folder). A second archive on
the same day gets `_2`, `_3`, ... instead of overwriting. Each holds the working
tree as it is, uncommitted edits included, minus everything in `.gitignore`
(build products, run output), and unpacks into a folder of the same name with an
`ARCHIVE_INFO.txt` recording the date, `git describe`, the commit, and whether
there were uncommitted changes. The archives stay on your machine: `archive/` is
in `.gitignore`.
