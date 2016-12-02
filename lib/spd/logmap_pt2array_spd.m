function V = logmap_pt2array_spd(p,X)
%LOGMAP_PT2ARRAY_SPD maps X from GL(n)/O(n) to TpM at xi.



[U D] = eig(p);
g = U*sqrt(D);
%invg = inv(g);
invg = diag(1./sqrt(diag(D)))*U'; % 1.3 X faster

V = zeros(size(X));
%% For each data
for i = 1:size(X,3)
    if norm(p-X(:,:,i)) < 1e-18
        V(:,:,i) = zeros(size(p));
        continue
    end
    y = invg*X(:,:,i)*invg';
    [U S] = eig(y);
    H = g*U;
    V(:,:,i) = H*diag(log(diag(S)))*H';
end