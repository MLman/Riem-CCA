function [f, fLA, feasibility, pt, pu]=myfunc_LA(X, Y, wx, wy, mu_x, mu_y, Xwr, Ywr, t, u, lx, ly, m)
% projected t
% projected u
    [f, pt, pu]= myfunc(X, Y, wx, wy, mu_x, mu_y, Xwr, Ywr);
    N  = size(X,3);

    tt = (t - mean(t));
    uu = (u - mean(t));
    T2 = sum(tt.^2);
    U2 = sum(uu.^2);
    fval = (tt'*uu)/(sqrt(T2)*sqrt(U2));

    % term 1 of x
    % term1 = sum_i \lambda_{x_i} d/dt_i g(wx,t_i,x_i,mu_x)
    term1 = 0;
    term3 = 0;
    for i =1:N
        gg = grad_g(mu_x, wx, X(:,:,i), t(i));
        term1 = term1 +lx(i)*gg;
        term3 = term3 + gg^2;
    end

    % term2 =  sum_i \lambda_{y_i} d/du_i g(wy,u_i,y_i,mu_y)
    term2 = 0;
    for i =1:N
        gg = grad_g(mu_y, wy, Y(:,:,i), u(i));
        term2 = term2 +ly(i)*gg;
        term3 = term3 + gg^2;
    end
    feasibility = term3;
    term3 = term3*m/2;
    fLA = fval + term1 + term2 - term3;
end
