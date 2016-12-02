function [l, ls] = norm_spd_tf(Ms, Vs)

ndim = size(Ms,1);
ls = zeros(ndim, 1);
for i =1:ndim
   ls(i) = innerprod_TpM_spd(Vs{i},Vs{i},Ms{i});
end
l = sum(ls);
ls = ls.^(1/2);