function [ksp_denoised, info] = klean(ksp_noise, ksp_data, patch_size, varargin)
%KLEAN Prewhiten + local low-rank denoising of high-dimensional k-space data.
%
%   ksp_denoised = KLEAN(ksp_noise, ksp_data, patch_size)
%   ksp_denoised = KLEAN(ksp_noise, ksp_data, patch_size, {casorati_shape, sv_threshold})
%   [ksp_denoised, info] = KLEAN(ksp_noise, ksp_data, patch_size)
%   [ksp_denoised, info] = KLEAN(ksp_noise, ksp_data, patch_size, {casorati_shape, sv_threshold})
%
%   This function performs:
%     1) Noise prewhitening using ksp_noise to estimate a whitening transform.
%     2) Local low-rank (LLR) denoising on ksp_data using patch-wise SVD and
%        singular-value thresholding.
%     3) Inverse whitening to return data to the original coil basis.
%
%   Dimensions:
%     ksp_data  : [READ, PHS1, PHS2, COIL, TAU]
%     ksp_noise : [READ, PHS1, PHS2, COIL] or any array with COIL as 4th dim
%     READ  - readout
%     PHS1  - phase-encode 1
%     PHS2  - phase-encode 2
%     COIL  - receiver coils
%     TAU   - time / dynamics / echoes (n_t)
%
%   Inputs:
%     ksp_noise    - Noise-only k-space used to compute the prewhitening matrix.
%                   Size: [n_ro, n_pe1, n_pe2, n_coil] (or equivalent with 4th dim = coils).
%     ksp_data     - Input k-space data to denoise.
%                   Size: [n_ro, n_pe1, n_pe2, n_coil, n_t].
%     patch_size   - Patch size for local low-rank denoising.
%                   Size: [patch_ro, patch_pe1, patch_pe2].
%
%   Optional input (varargin):
%     {casorati_shape, sv_threshold}
%       casorati_shape - Lookup table of Casorati matrix row counts
%                        used to select a threshold for a given patch shape.
%       sv_threshold   - Lookup table of singular-value thresholds 
%                        corresponding to casorati_shape.
%       If not provided, thresholds are computed internally via compute_sv_threshold().
%
%   Outputs:
%     ksp_denoised  - Denoised k-space data (same size as ksp_data).
%                    Size: [n_ro, n_pe1, n_pe2, n_coil, n_t].
%     info          - Optional patch-grid diagnostics, returned when requested.
%                    info.retain_pct: retained singular-value fraction per patch.
%                    info.s: singular values (before thresholded) per patch.
%                    info.thr: singular-value threshold applied per patch.
%
%   Notes:
%     - The coil dimension of ksp_noise and ksp_data must match.
%     - Patches are processed in a non-overlapping grid with step = patch_size along each axis.
%
%   See also: PREWHITEN, COMPUTE_SV_THRESHOLD

    assert(size(ksp_data,4)==size(ksp_noise,4), ...
        'Mismatch in number of coils between data and noise.');

    [~, ~, ~, n_coil] = size(ksp_noise);
    ksp_noise = reshape(ksp_noise, [], n_coil);
    w = prewhiten(ksp_noise);

    [n_ro, n_pe1, n_pe2, ~, n_t] = size(ksp_data);
    
    if ~isempty(varargin)
        casorati_shape = varargin{1}{1};
        sv_threshold   = varargin{1}{2};
    else
        [casorati_shape, sv_threshold] = compute_sv_threshold( ...
            [n_ro, n_pe1, n_pe2, n_coil, n_t], patch_size);
    end

    get_sv_threshold = @(i) lookupTable(casorati_shape, sv_threshold, i);

    ksp_data = reshape(permute(ksp_data, [1 2 3 5 4]), [], n_coil) * w;
    ksp_data = reshape(ksp_data, [n_ro, n_pe1, n_pe2, n_t, n_coil]);
    ksp_data = permute(ksp_data, [1 2 3 5 4]);

    info = [];
    if (nargout > 1)
        info = struct;
        info.retain_pct = zeros(length(1:patch_size(1):n_ro), ...
                       length(1:patch_size(2):n_pe1), ...
                       length(1:patch_size(3):n_pe2));
        info.s = zeros(length(1:patch_size(1):n_ro), ...
                       length(1:patch_size(2):n_pe1), ...
                       length(1:patch_size(3):n_pe2),n_t);
        info.thr = zeros(length(1:patch_size(1):n_ro), ...
                       length(1:patch_size(2):n_pe1), ...
                       length(1:patch_size(3):n_pe2));
    end
    ksp_denoised = zeros(n_ro, n_pe1, n_pe2, n_coil, n_t);
    k_counter = 1;
    for k = 1:patch_size(3):n_pe2
        k_range  = k:min([k+patch_size(3)-1, n_pe2]);
        i_counter = 1;
        for i = 1:patch_size(1):n_ro
            i_range = i:min([i+patch_size(1)-1, n_ro]);
            j_counter = 1;
            for j = 1:patch_size(2):n_pe1            
                j_range = j:min([j+patch_size(2)-1, n_pe1]);
                casorati = reshape( ...
                    ksp_data(i_range, j_range, k_range, :, :), [], n_t);
                [u,s,v] = svd(casorati, 'econ');
                thr = get_sv_threshold(size(casorati,1));
                s(s<thr) = 0;
                if isstruct(info)
                    info.retain_pct(i_counter,j_counter,k_counter) = ...
                        nnz(diag(s) >= thr)/numel(diag(s));
                    info.s(i_counter,j_counter,k_counter,:) = diag(s);
                    info.thr(i_counter,j_counter,k_counter) = thr;
                end
                casorati = u*s*v';
                ksp_denoised(i_range, j_range, k_range, :, :) = ...
                    reshape(casorati, length(i_range), ...
                                      length(j_range), ...
                                      length(k_range), n_coil, n_t);
                j_counter = j_counter + 1;
            end
            i_counter = i_counter + 1;
        end
        k_counter = k_counter + 1;
    end

ksp_denoised = reshape(permute(ksp_denoised, [1 2 3 5 4]), [], n_coil) * inv(w);
ksp_denoised = reshape(ksp_denoised, [n_ro, n_pe1, n_pe2, n_t, n_coil]);
ksp_denoised = permute(ksp_denoised, [1 2 3 5 4]);
