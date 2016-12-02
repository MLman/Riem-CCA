function uX = unitvec(X)
% UNITVEC converts column vectors in a matrix X into unit column vectors in uX
Z = sqrt(sum(X.^2));
Z(Z < eps) = 1;
uX = X./repmat(Z,size(X,1),1);
