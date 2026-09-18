# fmmlib2d

Fortran FMM library for 2D Laplace, Helmholtz and Cauchy particle interactions; README.md has the overview and `make` with no target lists the build targets.

## Building

- Compiler and flags are selected by the `HOST` variable at the top of `src/makefile` (and, separately, at the top of each `examples/*.make`): `linux-gfortran` (default), `linux-gfortran-openmp`, `macosx`, `linux-intel`, ... Override it on the command line, e.g. from `examples/`: `make -f hfmm2dpart.make HOST=linux-gfortran-openmp`. The `<variant>-openmp` targets in `examples/makefile` (and top-level `make test-openmp`) do exactly that.
- The gfortran hosts compile with `-std=legacy -finteger-4-integer-8`. GCC 10+ needs `-std=legacy` because the code relies on legacy array pointer conventions; `-finteger-4-integer-8` promotes default integers to 8 bytes. Keep both flags when adding a HOST.
- OpenMP runs need a large stack before launching (the OpenMP default is too small):

  ```bash
  export OMP_STACKSIZE=4096M
  ulimit -s unlimited
  ```

## API contracts

### Precision flag `iprec`

Every FMM entry point takes `iprec`, mapped to an FMM tolerance in `src/*fmm2dpart*.f`:

| iprec | -2 | -1 | 0 | 1 | 2 | 3 | 4 | 5 |
|-------|----|----|---|---|---|---|---|---|
| tol | 0.5 | 0.5e-1 | 0.5e-2 | 0.5e-3 | 0.5e-6 | 0.5e-9 | 0.5e-12 | 0.5e-15 |

Higher `iprec` means more expansion terms: slower but more accurate.

### nd-vectorized routines

The core routines carry an `nd` argument and an `_nd` suffix (e.g. `cfmm2dpart_nd`) so several charge/dipole densities are processed in one call. The legacy names (`cfmm2dpart`, `lfmm2dpart`, `rfmm2dpart`, `zfmm2dpart`, and their `self`/`targ` forms) live in `src/fmm2d_legacy.f` and simply call the `_nd` routine with `nd=1`. This backward-compatible wrapper layer must keep working: when an `_nd` signature changes, update the wrapper in the same commit.

### Tree reuse

FMM variants expose `<variant>_tree_build` / `<variant>_tree_eval` / `<variant>_tree_free`, separating geometry-dependent tree construction from density-dependent evaluation. Tree handles are owned by the `fmm2d_tree_module` in `src/fmm_tree_manager.f`; a handle from `_tree_build` may be passed to `_tree_eval` any number of times and must be released with `_tree_free`.

## Fortran conventions

- Fixed-form Fortran 77 with a few Fortran 90 features (the tree manager module).
- Implicit typing throughout: `implicit real*8 (a-h,o-z)`; declare integer arrays and complex variables explicitly.
- Do not add compiler flags or language features that break `-std=legacy` builds.
