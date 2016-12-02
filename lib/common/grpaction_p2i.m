function VsI = grpaction_p2i(p, Vs)
%GRPACTION_P2I performs group actions for each V of Vs to translate from
%T_{p}M to T_{I}M. p^-1/2 V p^-T/2 = p^-1/2 V p^-1/2 since p^-1/2 is symmetric.
%
%
%   See Also: QUADFUNC, GRPACTION_I2P

%   $ Hyunwoo J. Kim $  $ 2016/04/18 15:07:00 (CDT) $
 
    %pp = inv(sqrtm(p)); % Can be optimized by evd for big matrices.
    [U,S] = eig(p);
    pp = U*diag(1./sqrt(diag(S)))*U';
    VsI = quadfunc(pp, Vs);
end