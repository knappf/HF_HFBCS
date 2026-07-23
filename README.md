# HF_HFBCS

A Fortran code that performs spherical Hartree-Fock (HF) and Hartree-Fock-BCS
(HF-BCS) mean-field calculations for atomic nuclei, starting from two-body
(2N) and normal-ordered three-body (3N) interaction matrix elements expanded
in a harmonic-oscillator single-particle basis.

For each closed-subshell nucleus (Z or N equal to 2, 8, 20, 28, 50, 82, or
126) the corresponding species is treated with plain Hartree-Fock; otherwise
BCS pairing is switched on automatically.

## Repository layout

```
src/    Fortran source code and makefile
run/    Working directory: input file, interaction files, executable, output
```

### Source files (`src/`)

| File                          | Purpose                                                              |
|-------------------------------|-----------------------------------------------------------------------|
| `HF_BCS.f90`                  | Main program: drives input reading, basis setup, and the HF-BCS iteration |
| `type_defs.f90`                | Derived types for single-particle and 3-body basis states           |
| `declarations.f90`             | Shared module variables (quantum numbers, densities, matrices, energies) |
| `read_input.f90`               | Reads `input_HFBCS.dat` and sets BCS flags per species               |
| `basis.f90`                    | Builds the harmonic-oscillator single-particle and 3-body bases, sets initial occupations |
| `coupling_coefficients.f90`    | Angular-momentum recoupling coefficients (uses `wigxjpf`)             |
| `Tcm_2b.f90`                   | Center-of-mass kinetic energy two-body correction                    |
| `v2b_bin.f90`                  | Reads two-body matrix elements from a binary interaction file        |
| `v3b_no2b_bin.f90`             | Reads normal-ordered three-body matrix elements from a binary stream file |
| `hamiltonian.f90`              | Assembles the one- and two-body Hamiltonian matrix elements           |
| `HFBCS_solver.f90`             | Self-consistent HF(-BCS) iteration loop and output of levels/summary  |

## Building

Requirements:
- Intel Fortran compiler (`ifx`)
- Intel MKL (`libmkl_intel_lp64`, `libmkl_intel_thread`, `libmkl_core`, `iomp5`)
- [`wigxjpf`](https://fy.chalmers.se/subatom/wigxjpf/) library for Wigner symbols

Edit `src/makefile` to point `LDLIBS`/`WIGX` at your local `wigxjpf`
installation, then build:

```sh
cd src
make
```

This produces the executable `run/HF_BCS`.

## Running

The executable is run from the `run/` directory, where it expects to find:

- `input_HFBCS.dat` — calculation parameters (see below)
- `2belem.bin` — binary two-body interaction matrix elements
- `3belem_NO2.stream.bin` — binary normal-ordered three-body matrix elements
- `NN.bin`, `NNN.bin` — additional interaction input

```sh
cd run
./HF_BCS
```

### Input file: `input_HFBCS.dat`

Formatted text, one entry (or comma-separated group) per line, `#` starts a
trailing comment:

```
18,8      # mass and proton number of the core calculated by HF
3,6       # Nmax1_2N, Nmax12_2N   : model space of the 2N interaction file
3,6,9     # Nmax1_3N, Nmax12_3N, Nmax123_3N : model space of the 3N interaction file
3,6,9     # Nmax1, Nmax12, Nmax123 : model space for the HF(BCS) calculation
16.0      # hbar*omega [MeV]
1.d-7     # convergence threshold (epsilon)
```

Both the proton number `Z` and neutron number `N = A - Z` must be even.

### Output

- `HF_levels.dat` — Hartree-Fock single-particle levels
- `HFBCS_levels.dat` — HF-BCS quasiparticle levels (occupations, energies)
- `HFBCS_summary.dat` — summary of the converged calculation (total/pairing energy, etc.)
- `log` — run log
