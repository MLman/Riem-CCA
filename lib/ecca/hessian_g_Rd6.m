function [htg, dwdt_g, gtg]= hessian_g_Rd6(mu,w,x,t)
%HESSIAN_G returns gradient of projection function g.
%    t in R
%    w in T_{\mu}M
%    mu in M. Mean on M
%    x in M. Original proint

%   See Also:

%   $ Hyunwoo J. Kim $  $ 2016/04/24 17:59:07 (CDT) $

    h = 1e-10;
    gg0 = grad_g(mu, w, x, t);
    gg1 = grad_g(mu, w, x, t+h);
    htg = (gg1-gg0)/h;
    gtg = gg0;

    % For Rd cases.
    dw = invembeddingR6_vecs(mu,h*[1 0 0 0 0 0]');
    ggv1 = grad_g(mu,w+dw,x,t);
    dw = invembeddingR6_vecs(mu,h*[0 1 0 0 0 0]');
    ggv2 = grad_g(mu,w+dw,x,t);
    dw = invembeddingR6_vecs(mu,h*[0 0 1 0 0 0]');
    ggv3 = grad_g(mu,w+dw,x,t);
    dw = invembeddingR6_vecs(mu,h*[0 0 0 1 0 0]');
    ggv4 = grad_g(mu,w+dw,x,t);
    dw = invembeddingR6_vecs(mu,h*[0 0 0 0 1 0]');
    ggv5 = grad_g(mu,w+dw,x,t);
    dw = invembeddingR6_vecs(mu,h*[0 0 0 0 0 1]');
    ggv6 = grad_g(mu,w+dw,x,t);


    dwdt_g = [(ggv1-gg0)/h, (ggv2-gg0)/h, (ggv3-gg0)/h, ...
              (ggv4-gg0)/h, (ggv5-gg0)/h, (ggv6-gg0)/h]';

    dwdt_g = invembeddingR6_vecs(mu, dwdt_g);
end