% CCA for tensor field based on MICCAI10 method
clear
close
nsample = 10;
nvoxelsX = 3;
nvoxelsY = 4;
ndim = 3;

%% Data synthesis
X = randspd_for_tensorfield(nvoxelsX, ndim, nsample);
Y = randspd_for_tensorfield(nvoxelsY, ndim, nsample);

%wx,wy are axes. A collection of symmetric matices
[wx, wy, r1, u, v] = mcca_tf(X,Y);
[wx, wy, r2, u, v] = mscca_tf(X,Y);