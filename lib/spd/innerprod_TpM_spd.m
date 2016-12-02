function r = innerprod_TpM_spd(U,V,P)
%INNERPROD_TPM_SPD is a inner product in TpM of GL(n)/O(n)
    try
        invP = inv(P);
    catch
        invP = pinv(P);
        disp('pinv');
    end
    sqrtinvP= sqrtm(invP);
    r = trace(sqrtinvP*U*invP*V*sqrtinvP);
%   r = trace(sqrtinvP*U*sqrtinvP*sqrtinvP*V*sqrtinvP);
end