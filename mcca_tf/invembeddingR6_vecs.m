function Vnew = invembeddingR6_vecs(p, V)
%INVEMBEDDINGR6 converts v in R6 to vnew in TpM.
%    
%    p is a base point. 
%    vnew is a tangent vector.
%    v is a vector in R6.
%    step 1
%    v = [Sxx Sxy Sxz Syy Syz Szz]';    

   
    
%% Let's assume that input is correct

nmx = size(V, 2);
Vnew = zeros(3,3, nmx);
%    step 1 (common)
%    v = [Sxx Sxy Sxz Syy Syz Szz]';    
w = [1  sqrt(2) sqrt(2)   1 sqrt(2) 1]';
%sqrtp= sqrtm(p);
[U, D] = eig(p);
sqrtp = U*diag(sqrt(diag(D)))*U'; % faster sqrtm in matlab

V = V./repmat(w,1,nmx);

%    step 2 (room for speed up)
for i=1:nmx
    S = vec2symmx(V(:,i));
    Vnew(:,:,i) = sqrtp*S*sqrtp;
    Vnew(:,:,i) = (Vnew(:,:,i)+Vnew(:,:,i)')/2;
end
