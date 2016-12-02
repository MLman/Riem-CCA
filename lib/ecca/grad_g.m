function [gg, St, dtSt] = grad_g(mu, w, x, t)
%GRAD_G returns gradient of projection function g.
%    t in R
%    w in T_{\mu}M
%    mu in M. Mean on M
%    x in M. Original proint

    % S(t) = expmap_spd(mu,t*w)
    % St = expmap_spd(mu,t*w);
        P = mu;
        X = t*w;
        if norm(X) < 1e-10
            St = P;
        else
            [U, D] = eig(P);
            g = U*sqrt(D);
            invg = diag(1./sqrt(diag(D)))*U'; % 1.3 X faster

            Y = invg*X*invg';
            [V, S] = eig(Y);
            gv = g*V;
            St = gv*diag(exp(diag(S)))*gv';
        end
    % Derivative of S(t) w.r.t. t

    % A = invg*w*invg';
    % [Va Sa] = eig(A);
    % sqA = Va*sqrt(Sa);d
    % gsqAV = g*sqA*V;
    % dtSt = gsqAV*diag(exp(diag(S)))*gsqAV';


    %% Need to be improved like the above.
    sqrtmu = sqrtm(mu);
    invsqrtmu = inv(sqrtmu);
    A = invsqrtmu*w*invsqrtmu;
    dtSt = sqrtmu*A*expm(t*A)*sqrtmu;
    gg = 2*trace(logm(x\St)*inv(St)*dtSt);
end