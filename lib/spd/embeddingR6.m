function vnew = embeddingR6(p, v) 
%EMBEDDINGR6 embeds v in TpM onto R6.
%    
%    p is a base point. v is a tangent vector.

%    step 1
    try
        invp = inv(p);
    catch
        invp = pinv(p);
        disp('pinv');
    end
    sqrtinvp= sqrtm(invp);
%    step 2    
    S = sqrtinvp*v*sqrtinvp;
%    step 3
%    v = [Sxx Sxy Sxz Syy Syz Szz]';
    v = symmx2vec(S);
    w = [1  sqrt(2) sqrt(2)   1 sqrt(2) 1]';
    vnew = w.*v;
end