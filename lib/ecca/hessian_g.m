function [htg, dwdt_g, gtg]= hessian_g(mu,w,x,t)
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
    htg = (gg1-gg0)/h; % d^2g/dt^2
    gtg = gg0;         % dg/dt

    % For Rd cases.
    % d^2g/dwdt
    dimTM = size(w,1)*(size(w,1)+1)/2;
    u_mu = zeros(size(w,1),size(w,1),dimTM);
    ggv = zeros(dimTM,1);
    
    for i=1:dimTM
        uI = zeros(dimTM,1); % basis at T_{I}M
        uI(i)=1;
        u_mu(:,:,i) = invembeddingR6_vecs(mu,uI);% h*uI embedded at T_{mu}M
        ggv(i) = grad_g(mu,w+h*u_mu(:,:,i),x,t);
    end
    dwdt_g = zeros(size(w,1),size(w,1));
    for i=1:dimTM
        dwdt_g = dwdt_g+(ggv(i) - gg0)/h*u_mu(:,:,i);
    end
    
    % For R6 case.
    %
    % dw1 = invembeddingR6_vecs(mu,h*[1 0 0 0 0 0]');
    % ggv1 = grad_g(mu,w+dw1,x,t);
    % dw2 = invembeddingR6_vecs(mu,h*[0 1 0 0 0 0]');
    % ggv2 = grad_g(mu,w+dw2,x,t);
    % dw3 = invembeddingR6_vecs(mu,h*[0 0 1 0 0 0]');
    % ggv3 = grad_g(mu,w+dw3,x,t);
    % dw4 = invembeddingR6_vecs(mu,h*[0 0 0 1 0 0]');
    % ggv4 = grad_g(mu,w+dw4,x,t);
    % dw5 = invembeddingR6_vecs(mu,h*[0 0 0 0 1 0]');
    % ggv5 = grad_g(mu,w+dw5,x,t);
    % dw6 = invembeddingR6_vecs(mu,h*[0 0 0 0 0 1]');
    % ggv6 = grad_g(mu,w+dw6,x,t);
    % 
    %
    % 
    % % d^2g/dwdt
    % dwdt_g = [(ggv1-gg0)/h, (ggv2-gg0)/h, (ggv3-gg0)/h, ...
    %           (ggv4-gg0)/h, (ggv5-gg0)/h, (ggv6-gg0)/h]';
          
    % Double checked. 
    % This is equivalent to  \sum_i ggv_i dw_i.
    % 
    % dwdt_g = invembeddingR6_vecs(mu, dwdt_g);
end