function [gtf guf fval]= grad_f(t, u)
%
a = t - mean(t);
b = u - mean(u);
%gf = zeros(2*length(t),1);

abar = mean(a);
bbar = mean(b);

anorm = norm(a);
bnorm = norm(b);

gtf = 1/(anorm*bnorm)*(-bbar + b -a'*b*a/anorm^2);
guf = 1/(anorm*bnorm)*(-abar + a -a'*b*b/bnorm^2);
%gf = [dfdt;dfdu];


fval = a'*b/(anorm*bnorm);