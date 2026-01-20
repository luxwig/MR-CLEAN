function [casorati_shape, sv_threshold] = compute_sv_threshold(data_size, patch_size, varargin)
%COMPUTE_SV_THRESHOLD Compute Casorati row counts and SV thresholds for LLR denoising.
%
%   [casorati_shape, sv_threshold] = COMPUTE_SV_THRESHOLD(data_size, patch_size)
%   [casorati_shape, sv_threshold] = COMPUTE_SV_THRESHOLD(data_size, patch_size, N)
%
%   This function generates a list of possible Casorati matrix size
%   induced by the patch size and boundary patches, and estimates an
%   empirical singular-value threshold for each row count using Monte-Carlo
%   simulations of complex Gaussian random matrices.
%
%   Inputs:
%     data_size     - Size of k-space data in the form:
%                    [n_ro, n_pe1, n_pe2, n_coil, n_t].
%     patch_size    - Patch size for LLR in the form:
%                    [patch_ro, patch_pe1, patch_pe2].
%     N             - (Optional) Number of Monte-Carlo samples used to
%                    estimate thresholds. Default: 1000.
%
%   Outputs:
%     casorati_shape - Row counts (number of spatial samples * n_coil) for
%                     each possible Casorati matrix shape.
%                     Size: [1, n_shapes].
%     sv_threshold   - Estimated thresholds corresponding to casorati_shape.
%                     Each entry is the mean of the largest singular value
%                     of an i.i.d. complex Gaussian matrix of size
%                     [casorati_shape(i), n_t].
%                     Size: [1, n_shapes].

    if ~isempty(varargin)
        N = varargin{1};
        assert(isscalar(N) && N > 0, 'N must be a positive scalar.');
    else
        N = 1000;
    end
    assert(length(patch_size)==3, 'Invalid patch size')
    assert(length(data_size)==5, 'Invalid data size')
    data_size = reshape(data_size, 1, []);
    patch_size = reshape(patch_size, 1, []);

    m = data_size(1:3) >= patch_size;
    m = [m.*patch_size; rem(data_size(1:3), patch_size)];
    [A, B, C] = ndgrid(m(:,1), m(:,2), m(:,3));   
    factors = A(:) .* B(:) .* C(:);               
    casorati_shape = (data_size(4) * factors).';  
    casorati_shape(casorati_shape==0) = [];

    sv_threshold = zeros(size(casorati_shape));
    threshold_mc = zeros(N,1);
    for i = 1:length(casorati_shape)
        for j = 1:N
            [~,tmp,~] = svd(randn([casorati_shape(i) data_size(5)])*sqrt(0.5) ...
                + randn([casorati_shape(i) data_size(5)])*sqrt(0.5)*1j, 'econ');
            threshold_mc(j) = tmp(1, 1);
        end
        sv_threshold(i) = mean(threshold_mc);
    end
end