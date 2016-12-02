function VsP = grpaction_i2p(p, Vs)
%GRPACTION_I2P performs group actions for each V of Vs to translate from
%T_{I}M to T_{p}M. p^1/2 V p^T/2 = p^1/2 V p^1/2 since p^1/2 is symmetric.
%
%
%   See Also: QUADFUNC, GRPACTION_P2I

%   $ Hyunwoo J. Kim $  $ 2016/04/18 15:07:00 (CDT) $
 
    %pp = inv(sqrtm(p)); % Can be optimized by evd for big matrices.
    [U,S] = eig(p);
    pp = U*diag(sqrt(diag(S)))*U';
    VsP = quadfunc(pp, Vs);
end