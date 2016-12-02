function [dwx_LA, dwy_LA, dt_LA, du_LA, dlx_LA, dly_LA, lx_kp1, ly_kp1] = ...
    grad_LA(wx, wy, X, Y, t, u, mu_x, mu_y, lx, ly, m)

N  = size(X,3);
dwx_LA = zeros(size(wx));
dwy_LA = zeros(size(wy));
dt_LA = zeros(size(t));
du_LA = zeros(size(u));
dlx_LA = zeros(size(lx));
dly_LA = zeros(size(ly));
lx_kp1 = zeros(size(lx));
ly_kp1 = zeros(size(ly));

[gtf, guf, fval]= grad_f(t, u);

for i = 1:N
    [htg, dwdt_g, gtg] = hessian_g(mu_x,wx,X(:,:,i),t(i));
    lx_kp1(i) = lx(i)-m*gtg;
    dwx_LA = dwx_LA + dwdt_g*(lx(i)-m*gtg);
    dt_LA(i) = htg*(lx(i)-m*gtg); 
    dlx_LA(i) = gtg;
end
dt_LA = gtf + dt_LA;
for i = 1:N
    [hug, dwdu_g, gug] = hessian_g(mu_y,wy,Y(:,:,i),u(i));
    ly_kp1(i) = ly(i)-m*gug;
    dwy_LA = dwy_LA + dwdu_g*(ly(i)-m*gug);
    du_LA(i) = hug*(ly(i)-m*gug); 
    dly_LA(i) = gug;
end

du_LA = guf + du_LA;