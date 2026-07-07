# FESim

## Project overview

FESim simulates the BaTiO3 ferroelectric phase transition using the effective Hamiltonian method from [Zhang et al., 1995](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.52.6301). It consists of two C++ programs:

- **FESim** (`fesim/main.cpp`) — Monte Carlo simulation of a BaTiO3 supercell. Single-file C++ program (~1450 lines).
- **QMat** (`qmat/main.cpp`) — Pre-computes the dipole interaction matrix as a binary cache file that FESim reads at startup. Uses MPI for parallelism.

## Build & run

**CMake** auto-detects BLAS (MKL / Accelerate / OpenBLAS), OpenMP, and MPI.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
# → build/fesim/main    (FESim)
# → build/qmat/qmat     (QMat)
```

Both programs accept an optional config file path as the first argument (default: `./sim.ini` for FESim, `./dipole.ini` for QMat).

**FESim workflow:**
1. Build: `cmake -B build && cmake --build build`
2. Generate dipole matrix: `mpirun -np 4 build/qmat/qmat qmat/dipole.ini`
3. Run simulation: `build/fesim/main ./sim.ini` (config path optional; defaults to `./sim.ini`). See `fesim/run-BTO.sh` for a full cooling-curve example.

**Key constraint:** QMat lattice constants `ax/ay/az` must be set to 1.0 when generating input for FESim, because FESim applies the BaTiO3 lattice constant (7.456 Bohr) internally. The supercell dimensions `Nx/Ny/Nz` must match between QMat and FESim.

**Plotting:** `fesim/plot.py` reads monitor output with numpy and plots strain vs temperature with matplotlib.

## Architecture

### FESim (`fesim/main.cpp`)

Three degrees of freedom per unit cell in an Nx×Ny×Nz supercell with periodic boundaries:

| Variable | Shape | Meaning |
|----------|-------|---------|
| `u` (ux, uy, uz) | `[Nx][Ny][Nz]` each | Local soft-mode amplitude |
| `v` (vx, vy, vz) | `[Nx][Ny][Nz]` each | Inhomogeneous strain (acoustic displacements) |
| `eta_h[6]` | scalar array | Homogeneous strain (Voigt notation: xx, yy, zz, yz, zx, xy) |

**Interaction matrix `J0`:** Six `[Nx][Ny][Nz]` arrays (J0xx, J0yy, J0zz, J0yz, J0zx, J0xy). Stores only the first row (interaction of cell (0,0,0) with all others), exploiting translational symmetry to save memory — for Nx=20 this reduces from 2.8 GB to a tractable size. `init_J0()` loads the QMat cache file and adds self-energy (`keppa2`) and short-range interaction terms (j1–j7) directly into J0.

**Total energy** (`e_global()`, line 1082) has six components:
- `e_u_global()` — local mode: on-site quartic (α, γ) + harmonic coupling via J0 (O(N²), very slow)
- `e_elas_i_global()` — inhomogeneous strain elastic energy (B11, B12, B44)
- `e_elas_h_global()` — homogeneous strain elastic energy
- `e_int_i_global()` — local mode × inhomogeneous strain coupling (B1xx, B1yy, B4yz)
- `e_int_h_global()` — local mode × homogeneous strain coupling
- `e_p_global()` — pV term (pressure × volume)

**Monte Carlo algorithm** (`sim()`, line 1303): Metropolis algorithm with typewriter sweep (sequential left-to-right, per Vanderbilt 1997 convention). Each sweep updates all u → all v → eta_h (2×max(Nx,Ny,Nz) times). Energy changes are computed locally (O(1) per trial move via `e_u_change`, `e_v_change`, `e_eta_h_change`) — the global O(N²) energy is only computed for monitoring.

**Physics parameters** (`Param` struct, line 119): Hardcoded BaTiO3 values from DFT calculations — lattice constant (ax=7.456 Bohr), Born effective charge (z_star=9.956), dielectric constant (ε=5.24), on-site anharmonic coefficients (α=0.320, γ=-0.473), short-range coupling (j1–j7), elastic constants (B11/B12/B44), and mode-strain coupling (B1xx/B1yy/B4yz).

**I/O:** Configuration via INI files (using inih single-header library `INIReader.h`). State files (u, v, eta_h) are raw binary arrays. The cache_q file has a versioned header with supercell dimensions and lattice constants for compatibility checking.

### QMat (`qmat/main.cpp`)

Computes the dipole-dipole interaction matrix in reciprocal space using Ewald summation. Outputs a binary cache file with a version 2 header. Uses MPI to parallelize over k-points. Exploits symmetry to compute ~1/8 of the points.

### QMat versions

- `qmat/` — Current version; lattice constants configurable via INI.
- `archive/qmat_v1/` — Original version with hardcoded BaTiO3 lattice constants.
- `archive/qmat_v3_alpha/` — Experimental (incomplete) refactoring attempt.

## Configuration (sim.ini)

FESim reads `./sim.ini` by default, or the path given as the first command-line argument. Key sections:
- `[sys]` — Nx, Ny, Nz (supercell dimensions)
- `[init]` — Initial conditions: u/v can be "random" (with offset ux0/uy0/uz0 and spread uxd/uyd/uzd) or a path to a binary file; eta_h can be explicit values or a file; q = path to cache_q file
- `[sim]` — T (temperature K), steps, pressure, step sizes (u_step, v_step, eta_h_step)
- `[monitor]` — Output paths for energy and eta_h, with step intervals
- `[out]` — Final state output paths

## Dependencies

- CBLAS: Intel MKL (Linux/HPC) or Apple Accelerate (macOS)
- OpenMP (FESim only, for parallel energy computation)
- OpenMPI (QMat only)
- inih (INI parser, included as `INIReader.h`)
- Python: numpy, matplotlib (plotting only)

## Platform notes

BLAS header selection uses CMake-provided defines (`BLAS_MKL` / `BLAS_ACCELERATE` / `BLAS_GENERIC`) rather than compiler-identity macros. MKL-specific vector math functions (`vdCos`, `vdDiv`) have compatibility wrappers that dispatch to platform-optimized equivalents where available.

Random number generation uses `rand_r` on Unix (thread-safe with per-thread seed) and `rand` on Windows.
