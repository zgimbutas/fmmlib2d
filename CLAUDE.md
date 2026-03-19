# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FMMLIB2D is a Fortran library implementing the Fast Multipole Method (FMM) for 2D particle interactions governed by Laplace and Helmholtz equations. Authors: Leslie Greengard and Zydrunas Gimbutas. License: BSD 3-clause.

## Build Commands

All builds use GNU Make. Default Fortran compiler is gfortran with `-O2 -std=legacy -finteger-4-integer-8`.

```bash
make lib                # Build static library (src/libfmm2d.a)
make test               # Compile and run Helmholtz FMM test (single-threaded)
make test-openmp         # Compile and run OpenMP Helmholtz test
make clean              # Remove object files
make distclean          # Remove all generated files
```

### Running individual FMM variant tests

From `examples/`:
```bash
make -f hfmm2dpart.make    # Helmholtz FMM
make -f lfmm2dpart.make    # Laplace FMM (complex)
make -f rfmm2dpart.make    # Laplace FMM (real)
make -f zfmm2dpart.make    # Cauchy FMM
make -f cfmm2dpart.make    # Generalized Cauchy FMM
```

Each produces an interactive executable `int2` that prompts for number of sources, runs FMM vs direct O(N^2) evaluation, and reports errors. Example makefiles link against the prebuilt library `../src/libfmm2d.a`.

### OpenMP builds

Append `-openmp` to any example target (e.g., `make -f hfmm2dpart.make HOST=linux-gfortran-openmp`). Requires:
```bash
export OMP_STACKSIZE=4096M
ulimit -s unlimited
```

### MATLAB/Octave MEX

```bash
make mwrap              # Build MWrap wrapper generator (needs bison, flex)
make mex-octave         # Build Octave MEX files
make mex-matlab         # Build MATLAB MEX files
make test-mex-octave    # Run Octave MEX tests
make test-mex-matlab    # Run MATLAB MEX tests
```

## Architecture

### nd-vectorized API (new)

All core FMM routines support an `nd` parameter for processing multiple charge/dipole densities in a single call. The nd-vectorized routines have `_nd` suffix (e.g., `cfmm2dpart_nd`). Legacy wrappers (`cfmm2dpart`, etc.) call the `_nd` routines with `nd=1`, preserving backward compatibility.

### Tree-reuse API (new)

Each FMM variant supports `_tree_build` / `_tree_eval` / `_tree_free` for separating geometry-dependent tree construction from density-dependent evaluation. Tree handles are managed by the `fmm2d_tree_module` in `fmm_tree_manager.f`.

### Core FMM variants (all in `src/`)

| File | Kernel | Notes |
|------|--------|-------|
| `hfmm2dpart.f` | Helmholtz H₀⁽¹⁾(k\|r\|) | Far-field + exponential M2L variant; complex wave number k |
| `lfmm2dpart.f` | Laplace log(r) | Complex charges/dipoles |
| `rfmm2dpart.f` | Laplace log(r) | Real charges/dipoles (~2x faster than lfmm2d) |
| `zfmm2dpart.f` | Cauchy 1/(z-w) | Complex-plane sums |
| `cfmm2dpart.f` | Generalized Cauchy | log(z-w) + dipole/(z-w) |
| `fmm2d_legacy.f` | (all) | Legacy nd=1 wrappers for backward compatibility |

### Alternative variants (in `src/`, not compiled by default)

| File | Notes |
|------|-------|
| `cfmm2dpart_exp.f` | Cauchy with exponential M2L (drop-in for cfmm2dpart.f) |
| `lfmm2dpart_opt.f` | Optimized Laplace using rfmm2d backend (drop-in for lfmm2dpart.f) |

### Supporting routines (`src/`)

- `helmrouts2d.f` / `helmrouts2d_exp.f` / `helmrouts2d_far.f` — Helmholtz kernel evaluation (base, exponential, far-field)
- `laprouts2d.f` / `laprouts2d_exp.f` — Laplace kernel evaluation
- `hfmm2drouts.f` / `hfmm2drouts_exp.f` / `hfmm2drouts_far.f` / `hfmm2drouts_far_exp.f` — Helmholtz FMM translation operators
- `lfmm2drouts.f` / `lfmm2drouts_exp.f` — Laplace FMM translation operators
- `d2tstrcr.f` / `d2tstrcr_omp.f` — adaptive quadtree construction
- `h2dterms.f` / `l2dterms.f` — expansion term count based on precision level
- `hank103.f` — Hankel function H₀⁽¹⁾, H₁⁽¹⁾ evaluation
- `cdjseval2d.f` — Chebyshev-Legendre expansion evaluation
- `dfft.f` / `next235.f` — FFT utilities for exponential variants
- `hwts2e.f` / `lwtsexp2.f` — quadrature weight precomputation (hwts2e.f includes `hwts_*.txt` data files)
- `fmm_tree_manager.f` — F90 module for tree handle management

### Precision levels (iprec parameter)

Values -2 through 5 control accuracy from ~0.5 tolerance to ~0.5e-15 (full double precision). Higher iprec = more expansion terms = slower but more accurate.

### MATLAB/Octave interface (`matlab/`)

- `.m` wrapper files call Fortran via MEX (e.g., `hfmm2dpart.m`)
- `fmm2d.mw` — MWrap interface definition generating `fmm2d.c`
- Pre-compiled MEX binaries included for Linux/macOS/Windows 64-bit
- `test_*_direct.m` — validation tests comparing FMM against direct evaluation

## Fortran conventions

- Fixed-form Fortran 77 with some Fortran 90 features
- Implicit typing: `implicit real*8 (a-h,o-z)`
- GCC 10+ requires `-std=legacy` (already set in makefiles) due to legacy array pointer conventions
- `-finteger-4-integer-8` flag promotes default integers to 8 bytes
- Compiler selection: change `HOST` variable in `src/makefile` (linux-gfortran, macosx, linux-gfortran-openmp, linux-intel, etc.)
