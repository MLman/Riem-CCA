function Qs = quadfunc(p,Vs)
%QUADFUNC calculates pVp^T for each V in Vs.
%
%
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2016/04/18 15:02:44 (CDT) $

    Qs = zeros(size(p,1), size(p,1), size(Vs,3));
    parfor i=1:size(Vs,3)
        Qs(:,:,i) =  p*Vs(:,:,i)*p';
    end
end