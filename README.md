# Modernized DCC — 2D Transient Potential / Heat Diffusion BEM Solver

Modernization of the legacy Fortran 77 **`DCC`** Boundary Element Method (BEM) program (from the master's thesis by *Teresinha Costa Effren*; legacy files preserved in [`dcc.f`](dcc.f) and [`dcc_orig.f`](dcc_orig.f)) into **modular Fortran 90**, **dynamic memory allocation**, **OpenMP** parallelization, **LAPACK** linear system solving, **Gmsh** mesh generation, and **Python** pre/post-processing.

## Key Modernization Highlights

1. **Modular Fortran 90 (`implicit none` & Dynamic Allocation)**:
   - [`data_mod.f90`](data_mod.f90): global precision parameters and `allocatable` arrays (removing the legacy hardcoded 200-node / time-step limits).
   - [`solver_mod.f90`](solver_mod.f90): explicit LAPACK interfaces (`DGETRF`, `DGETRS`, `DGESV`) replacing manual Gaussian elimination (`SLNPD`). The constant system matrix `C_MAT` is LU-factored **once** (`FACTOR_SYSTEM`) and reused across all `NINTER` time steps (`SOLVE_FACTORED`).
   - [`assembly_mod.f90`](assembly_mod.f90): boundary geometry and outward-normal precomputation (`PRECOMPUTE_GEOMETRY`), parallel influence matrix assembly (`MATGH`), time-marching loop (`SOLVE_TRANSIENT`), and parallel internal point evaluation (`PTINT`) featuring configurable Gauss-Legendre quadrature (`NG`) and the Johnston & Elliott (2005) **`sinh`** transformation (`SINHTRANS`) for near-singular internal points.
   - [`io_mod.f90`](io_mod.f90): reads explicit-connectivity input (`CONEC(NE, 2)`) and time-dependent boundary values from `dcc_modern.dat`, and exports structured results (`dcc_modern.out`, `ubou.dat`, `pbou.dat`, `uint.dat`).
   - [`main.f90`](main.f90): main driver reporting wall-clock runtime (`omp_get_wtime`) and OpenMP thread count.
2. **OpenMP Parallelization**:
   - Outer loops in `MATGH` (`collapse(2)` over element-collocation pairs `(J, I)`) and `PTINT` (over internal evaluation points `KN = 1..NPI`) run in parallel with thread-local accumulators free of race conditions.
3. **Double Precision & Thermal Diffusivity ($K$) Fix**:
   - Compiled with `-fdefault-real-8`.
   - Fixes the transient boundary integral kernel scaling in `IQ` (thesis Eqs. 2.39–2.40) so general diffusivities $K \neq 1$ behave consistently.
4. **Python + Gmsh Workflow**:
   - [`input_data.py`](input_data.py): selects among 3 benchmark validation cases (including time-varying boundary conditions).
   - [`preproc.py`](preproc.py): meshes `.geo` files via the Gmsh Python API (1st-order elements), automatically enforces consistent outward boundary normal orientation using adjacent domain triangles, plots `geometry.png`, and writes `dcc_modern.dat`.
   - [`posproc.py`](posproc.py): validates BEM results against exact analytical solutions (1D/2D Fourier series and Bessel $J_0$ Duhamel superposition), computes relative $L_2$ errors, generates `result_plot.png` and `temperature_map.png`, and exports transient fields to ParaView (`resultado_dcc.vtk`).

## Requirements

- **Fortran Compiler & Math Libraries** (Windows via MSYS2 UCRT64, or Linux):
  - `gfortran` with OpenMP support (`-fopenmp`)
  - `LAPACK` and `BLAS` (`-llapack -lblas`)
- **Python 3** packages:
  - `gmsh`, `meshio`, `numpy`, `matplotlib`

## Building

On Windows (with MSYS2 UCRT64 at `C:\msys64\ucrt64\bin`):
```bat
build.bat
```
Or via `mingw32-make` / `make`:
```bash
mingw32-make
```

## Usage

Select a default benchmark in [`input_data.py`](input_data.py) (`problem = 1`, `2`, or `3`) or pass `--problem`:

```bash
python preproc.py --problem 1
dcc_modern.exe
python posproc.py --problem 1
```

## Benchmark Problems & Validation Results

| Problem | Geometry (`.geo`) | Discretization | Global $L_2$ Error | Final $L_2$ Error | Runtime (8 threads) |
| :--- | :--- | :--- | :---: | :---: | :---: |
| **1** (Thesis Ex. 2, Table 3.1) | [`rect_10x8.geo`](rect_10x8.geo) `[-5,5]x[-4,4]`, $u=1$ | `NE=68`, `NPI=305`, `dt=1.0 s`, `NINTER=20` | `1.9904%` | `0.0714%` | `0.0530 s` |
| **2** (Thesis Ex. 1) | [`placa_6x6.geo`](placa_6x6.geo) `[0,6]x[0,6]`, 1D mixed BCs | `NE=56`, `NPI=202`, `dt=0.5 s`, `NINTER=40` | `0.9267%` | `0.0745%` | `0.0550 s` |
| **3** (Thesis Ex. 3, Thermal Shock) | [`disco_r1.geo`](disco_r1.geo) Disk $R=1$, $u(t)=50 \to 100$ | `NE=48`, `NPI=193`, `dt=0.025 s`, `NINTER=40` | `1.0296%` | `0.0942%` | `0.0480 s` |

On Problem 2, OpenMP parallelization reduces solver runtime from `0.2420 s` (1 thread) to `0.0570 s` (8 threads, $\approx 4.25\times$ speedup), while diffusivity scaling invariance was verified at `K = 2.5`, `dt = 0.2 s` (`0.9267%` global $L_2$ error).

## Generated Outputs

- `geometry.png`: Gmsh triangular mesh, Dirichlet/Neumann boundary segments, outward normals, and internal nodes.
- `result_plot.png`: temperature time histories at validation monitor points and spatial profiles against analytical solutions (plus Table 3.1 comparison on Problem 1).
- `temperature_map.png`: 2D temperature contour snapshots across four time instants.
- `resultado_dcc.vtk`: VTK unstructured triangle mesh with multi-step temperature fields for ParaView.
