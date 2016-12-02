%% CCA Implementation
N = 10;
dimx1 = 20;
dimx2 = 20;
X1 = rand(N, dimx1)*2-0.5;
X2 = rand(N, dimx2)*2-0.5;

%% Find threshold
% v values for interval
% idx of constraints

% max iteration
maxiter = 100;
% Sparsity parameters
c1 = sqrt(size(X1,2));
c2 = sqrt(size(X2,2));
c3 = 1e-20;
tic
[ w1, w2, fval, r, status, iter] = scca_ver2(X1, X2, c1, c2, c3, maxiter);
r % correlation


