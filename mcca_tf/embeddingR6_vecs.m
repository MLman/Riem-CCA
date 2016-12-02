function Vnew = embeddingR6_vecs(p,V)
%EMBEDDINGR6_VECS embeds V in TpM onto R6.
%    
%    p is a base point. V is tangent vectors.
%    V is 3-by-3-n. n is the number of matrices.

if length(size(V)) == 2
    Vnew = embeddingR6(p, V);
    return;
elseif length(size(V)) == 3
    
    % Common part
    nmx = size(V,3);
    Vnew = zeros(6,nmx);
    S = zeros(size(V));
    w = [1  sqrt(2) sqrt(2)   1 sqrt(2) 1]';
    %    step 1 (common)
    [U D] = eig(p);
    sqrtinvp = U*diag(1./sqrt(diag(D)))*U'; % 1.3 X faster
    
%     try
%         invp = inv(p);
%     catch
%         invp = pinv(p);
%         disp('pinv');
%     end
%     sqrtinvp= sqrtm(invp);
    
    %    step 2    
    % For each data set
    for i=1:nmx
       S(:,:,i) = sqrtinvp*V(:,:,i)*sqrtinvp;
    end
    %    step 3
    %    v = [Sxx Sxy Sxz Syy Syz Szz]';

    Vnew(1,:) = S(1,1,:);
    Vnew(2,:) = S(1,2,:);
    Vnew(3,:) = S(1,3,:);
    Vnew(4,:) = S(2,2,:);
    Vnew(5,:) = S(2,3,:);
    Vnew(6,:) = S(3,3,:);

    Vnew = repmat(w , 1,nmx).*Vnew;
else
    error('V is wrong input');
end
