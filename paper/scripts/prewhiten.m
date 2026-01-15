function Z = prewhiten(X)
%     C = X'*X/(size(X,1)-1);
%     [U,Sig] = svd(C);
%     Z = U*inv(sqrt(Sig));
      C = X'*X;
      [V, Sig] = eig(C);
      Z = V*inv(sqrt(Sig))*sqrt(size(X,1)-1);
    