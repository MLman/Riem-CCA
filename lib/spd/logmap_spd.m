function v = logmap_spd(P,X)
%LOGMAP_SPD maps Y on GL(n)/O(n) to a tangent space of manifold M at X.
%
%    The tangent space is a set symmetric matrices(n)
if norm(P-X) < 1e-18
    v = zeros(size(P));
    return
end
[U D] = eig(P);
g = U*sqrt(D);
invg = inv(g);
y = invg*X*invg';
[V S] = eig(y);
H = g*V;
v = H*diag(log(diag(S)))*H';
v = (v+v')/2;
%rtX = sqrtm(X);
%invrtX = inv(rtX);
%v = rtX*logm(invrtX*Y*invrtX)*rtX;
%v = (v+v')/2;