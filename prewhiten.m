function W = prewhiten(X)
%PREWHITEN Compute a noise prewhitening matrix from noise samples.
%
%   W = PREWHITEN(X)
%
%   Given noise-only data X, this function computes a linear transform W
%   such that multiplying coil vectors by W approximately whitens the noise,
%   i.e., produces unit noise covariance across coils.
%
%   Input:
%     X - Noise samples arranged as [n_samples, n_coil]. Each row is one
%         sample across coils.
%
%   Output:
%     W - Whitening matrix of size [n_coil, n_coil]

     C = X'*X;
     [V, Sig] = eig(C);
     W = V*inv(sqrt(Sig))*sqrt(size(X,1)-1);
end