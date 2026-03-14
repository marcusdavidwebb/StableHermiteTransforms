# StableHermiteTransforms
Self-contained MATLAB codes for **stable transforms between Hermite function coefficients** and **function values at Gauss–Hermite nodes** (roots of Hermite polynomials).

The main goal is to avoid numerical overflow/instability that appears when assembling Hermite transform matrices directly for moderate/large transform sizes `N`.

## What this repo provides

Given Hermite coefficients `cfs` (length `N`) and values `vals` at Gauss–Hermite nodes (length `N`), the stable algorithms construct a diagonal scaling vector `d` and an orthogonal matrix `Q` such that

- **[coefficients to values]**
  - `vals = d .* (Q' * cfs)`
- **[values to coefficients]**
  - `cfs = Q * (vals ./ d)`

This lets you apply the transform without explicitly forming a dense, ill-conditioned “Hermite matrix”.

## Repository structure

- **[modules/]**
  - **[initialise_Hermite_transform_Golub_Welsch.m]**
    Stable initializer based on the Golub–Welsch eigen-decomposition of the Hermite Jacobi matrix (recommended default).
  - **[initialise_Hermite_transform_Bunck.m]**
    Alternative stabilized recurrence approach (Bunck-style scaling).
  - **[initialise_Hermite_transform_unstable.m]**
    Baseline direct assembly of full transform matrices via a three-term recurrence (intended for comparison; becomes unstable for larger `N`).
  - **[quadrature/]**
    Gauss–Hermite quadrature routines (nodes/weights) based on recurrence coefficients + Gaussian rule construction.
- **[examples/]**
  - **[gross_pitaevskii.m]**
    Strang splitting for the 1D Gross–Pitaevskii equation using Hermite transforms.
  - **[gross_pitaevskii_reference_sln.m]**
    Computes (and caches) a reference solution in Hermite coefficient space.
  - **[gross_pitaevskii_compare_spatial_accuracy.m]**
    Compares spatial accuracy of stable vs direct/unstable transforms vs a reference solution.
- **[performance_evals_run.m / performance_evals_plot.m]**
  Benchmarks/plots comparing assembly time, accuracy, and condition numbers for the different initializers.
- **[data/, images/]**
  Saved benchmark/reference datasets and exported figures.

## Quick start (MATLAB)

From the repo root:

```matlab
addpath('modules/');
addpath('modules/quadrature/');

N = 128;
[d, Q] = initialise_Hermite_transform_Golub_Welsch(N);

% Example: convert coefficients -> values and back
cfs  = randn(N,1);
vals = d .* (Q' * cfs);
cfs2 = Q * (vals ./ d);

norm(cfs - cfs2)
```

## Examples

Run the Gross–Pitaevskii demo:

```matlab
cd examples
gross_pitaevskii
```

This produces a solution plot PDF in `images/`.

To run the spatial accuracy comparison:

```matlab
cd examples
gross_pitaevskii_compare_spatial_accuracy
```

This script computes/loads a cached reference solution under `data/` and exports an accuracy plot PDF.

## Benchmarks

To reproduce the performance study comparing:

- direct unstable assembly (`initialise_Hermite_transform_unstable`)
- Bunck initializer (`initialise_Hermite_transform_Bunck`)
- Golub–Welsch initializer (`initialise_Hermite_transform_Golub_Welsch`)

run:

```matlab
performance_evals_run
performance_evals_plot
```

Plots are exported to `images/`.

## Notes on Gauss–Hermite nodes (`hermpts`)

Some scripts call `hermpts(N)` (primarily for plotting or node generation). Ensure `hermpts` is available on your MATLAB path (commonly via Chebfun) or replace those calls with node generation from this repo’s quadrature routines (e.g. `quad_gauss_hermite(N)`).

## Citation

If you use this code in academic work, please cite the paper listed in `CITATION.cff`.
