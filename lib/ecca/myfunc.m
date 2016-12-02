function [fval, t, u]= myfunc(X, Y, wx, wy, mu_x, mu_y, Xwr, Ywr)
%MYFUNC is the objective function for Riem-CCA.
%    X, Y are SPDs
%    t,u are sets of projection results.
%    mu_x and mu_y are intrinsic means of X and Y on SPD.
%    wx, wy are tangent vectors at T_{mu_x}M and T_{mu_y}M.
%    Xwr = logmap_spd(mu_x, X)
%    Ywr = logmap_spd(mu_y, Y)

    % Function value test

    N = size(X,3);
    t = zeros(N,1);
    u = zeros(N,1);

    %if bitand(mode,1) >0    
    for i = 1:N
        t(i) = proj_w(mu_x, wx, X(:,:,i), Xwr(:,:,i));
        u(i) = proj_w(mu_y, wy, Y(:,:,i), Ywr(:,:,i));
    end
    tt = (t - mean(t));
    uu = (u - mean(u));
    T2 = sqrt(mean(tt.^2));
    U2 = sqrt(mean(uu.^2));
    fval = mean(tt.*uu)/(T2*U2);
    %end
end