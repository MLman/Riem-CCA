function [wx, wy, r, t,u, status] = rcca_ga(X, Y,varargin)
%RTCCA is a Riem-CCA algorithm 2 with group action in [ECCV2014].
%
%
%   See Also: KARCHER_MEAN_SPD

%   $ Hyunwoo J. Kim $  $ 2016/04/02 17:23:41 (CDT) $


    if nargin == 3
        niter = varargin{1};
    else
        niter = 100;
    end
    
    mu_x = karcher_mean_spd(X,[],niter);
    mu_y = karcher_mean_spd(Y,[],niter);
    
    % Logmap
    Xwr = logmap_pt2array_spd(mu_x, X);
    Ywr = logmap_pt2array_spd(mu_y, Y);
    
    % Embedding in R^d (vectorize tangent vectors or symmetric matrices)
    % Equivalent to group action to bring tangent vectors to Identity.
    %
    % Refactoring is need. 
    % This MUST be a group action and vectoriation for general dimensions.
    
    XwrI = grpaction_p2i(mu_x, Xwr);
    YwrI = grpaction_p2i(mu_y, Ywr);
    Xv = symmx2vec(XwrI)';
    Yv = symmx2vec(YwrI)';
    
    % Xv = embeddingR6_vecs(mu_x,Xwr)';
    % Yv = embeddingR6_vecs(mu_y,Ywr)';
 
    [wxv, wyv, r, t, u, status] = canoncorr(Xv, Yv); % Tangent space output.
    % wx, wy are vectorized tangent vectors in T_{I}M
    wxI = vec2symmx(wxv(:,1));
    wyI = vec2symmx(wyv(:,1));
    wx = grpaction_i2p(mu_x, wxI);
    wy = grpaction_i2p(mu_y, wyI);
    status.mu_x = mu_x;
    status.mu_y = mu_y;
    status.wxI = wxI;
    status.wyI = wyI;
end
    