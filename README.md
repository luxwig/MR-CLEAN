# MR KLEAN: K-space Local Low-Rank Denoising for MRI

**MR KLEAN** (Magnetic Resonance k-Space Local Low-Rank Estimation for Attenuating Noise) is a generalized local low-rank (LLR) denoising framework that operates directly on complex-valued k-space data prior to image reconstruction. The method is acquisition- and reconstruction-agnostic, and is designed for high-dimensional MRI data.

## Repository Contents

- **Manuscript pipelines** (`paper/`)
  Denoising pipelines used to generate the manuscript results (ASL, cardiac, phantom).

- **Toolbox interface** (`klean.m`)
  General-purpose MATLAB entry point for k-space denoising.

## Quick start

### Expected dimensions

MR KLEAN assumes the following dimension order:

- `ksp_data` : **[READ, PHS1, PHS2, COIL, TAU]**
- `ksp_noise`: **[READ, PHS1, PHS2, COIL]**

Where:
- `READ` = readout
- `PHS1` = phase-encode 1
- `PHS2` = phase-encode 2
- `COIL` = receiver coils
- `TAU` = dynamics / echoes / repetitions

### Basic usage

```matlab
% Patch size in [READ, PHS1, PHS2]
patch_size = [6, 6, 6];

% Denoise k-space
ksp_denoised = klean(ksp_noise, ksp_data, patch_size);
```

## Function interface

```matlab
ksp_denoised = klean(ksp_noise, ksp_data, patch_size)
ksp_denoised = klean(ksp_noise, ksp_data, patch_size, {casorati_shape, sv_threshold})
```

### Inputs

- `ksp_noise`
  Noise-only k-space used to estimate the prewhitening transform.
  Size: `[n_ro, n_pe1, n_pe2, n_coil]`.

- `ksp_data`
  Input k-space data to denoise.
  Size: `[n_ro, n_pe1, n_pe2, n_coil, n_t]`.

- `patch_size`
  Patch size for local low-rank denoising.
  Size: `[patch_ro, patch_pe1, patch_pe2]`.

### Output

- `ksp_denoised`
  Denoised k-space data (same size as `ksp_data`).
  Size: `[n_ro, n_pe1, n_pe2, n_coil, n_t]`.

## Options

### Precomputed singular-value threshold table (optional)

By default, `klean` computes a singular-value threshold table internally based on the data and patch size. If you want to reuse thresholds across datasets (or avoid recomputing), you can pass a precomputed table:

```matlab
[casorati_shape, sv_threshold] = compute_sv_threshold( ...
    [n_ro, n_pe1, n_pe2, n_coil, n_t], patch_size);

ksp_denoised = klean(ksp_noise, ksp_data, patch_size, {casorati_shape, sv_threshold});
```

## Author

Ludwig S. Zhao (ludwigz@seas.upenn.edu)

## License

Released under the license in [`LICENSE`](LICENSE).
