# Changelog

## [1.0.0]

### Build system

- **CMake** replaces hand-written shell scripts. CMake auto-detects BLAS (Intel MKL / Apple Accelerate / OpenBLAS), OpenMP, and MPI via `find_package`. Old build scripts moved to `archive/`.
- **inih** (INI config parser) is now a git submodule (`external/inih/`) instead of four identical copies of `INIReader.h`. inih is compiled as a static library.

### Directory layout

- **`qmat_v2/` renamed to `qmat/`** — the `v2` suffix referred to the cache file format, not the program version.
- **`qmat/` (v1) and `qmat_v3_alpha/` moved to `archive/`** — v1 is deprecated (FESim rejects its output), v3_alpha is an incomplete refactoring.

### Fixed

- BLAS header selection no longer depends on compiler identity macros (`__INTEL_COMPILER`, `__APPLE__`). CMake detects the actual BLAS library and defines `BLAS_MKL` / `BLAS_ACCELERATE` / `BLAS_GENERIC` accordingly. Fixes GCC+MKL and Intel compiler+OpenBLAS combinations.
- `vdCos` and `vdLinearFrac` compatibility layers now use `BLAS_*` macros. On Accelerate, `vdCos` maps to `vvcos`; `vdLinearFrac` (always called with scale=1, shift=0) is replaced by element-wise division (`vdDiv` on MKL, `vvdiv` on Accelerate).
- `MAX` / `MIN` macros replaced with `std::max` / `std::min` from `<algorithm>`.
- C99 compound literals (`(double [3]){...}`) replaced with named C++ arrays for GCC compatibility.
- Missing `#include <cmath>` restored for `exp()` on non-Accelerate platforms.

### Added

- Both FESim and QMat accept an optional config file path as the first command-line argument.

### Changed

- `fesim/run-BTO.sh` rewritten: removed cluster-specific `#SBATCH` directives, added `CALC_DIR` / `MAIN_BIN` / `CACHE_Q` variables at the top for user configuration.
- `qmat/run.sh` removed (unnecessary — just `mpirun`).
