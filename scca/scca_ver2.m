function [ w1, w2, fval, fvalc, status, iter, u, v] = scca_ver2(X1, X2, c1, c2, varargin)
%SCCA_VER2 is sparse cca of Daniela Witten and Robert Tibshirani.
%    scca_ver1 is L1 and L2 constraints. Centering X1, and X2 is recommended. 
%
%    [ w1, w2 fval fvalc status iter] = SCCA_VER2(X1, X2, c1, c2)
%    [ w1, w2 fval fvalc status iter] = SCCA_VER2(X1, X2, c1, c2, c3)
%    [ w1, w2 fval fvalc status iter] = SCCA_VER2(X1, X2, c1, c2, c3, maxiter)
%
% 
%   See also SUBSOLVER2

%% Find threshold
% v values for interval
% idx of constraints
[N dimx1 ] = size(X1);
[N dimx2 ] = size(X2);

X1 = centering(X1')';
X2 = centering(X2')';
z = std(X1);
X1 = X1./repmat(z,N,1);
z = std(X2);
X2 = X2./repmat(z,N,1);

% Stop condition
if nargin >= 5
    c3 = varargin{1};
else
    c3 = 1e-10;
end

% max iteration
if nargin >= 6
    maxiter = varargin{2};
else
    maxiter = 100;
end

%% SCCA
w1 = unitvec(ones(dimx1,1))/sqrt(dimx1);
w2 = unitvec(ones(dimx2,1))/sqrt(dimx2);
prev_w1 = w1;
prev_w2 = w2;

status = 'maxiter';
for iter = 1:maxiter
    p1 = X1'*X2*w2;
    w1 = subsolver2(p1,c1);
    %% projection 
    l2norm_w1 = norm(w1);
    if l2norm_w1 == 0
        status = 'w1 is zero vector.';
        break
    else
        w1 = w1/l2norm_w1;
    end
    
    p2 = (w1'*X1'*X2)';
    w2 = subsolver2(p2,c2);
    %% projection 
    l2norm_w2 = norm(w2);
    if l2norm_w2 == 0 
        status = 'w2 is zero vector.';
        break
    else
        w2 = w2/l2norm_w2;
    end
    dw1 = sum((prev_w1-w1).^2);
    dw2 = sum((prev_w2-w2).^2);
    
    if dw1 < c3 && dw2 < c3
        status = 'solved';
        break
    end
    prev_w1 = w1;
    prev_w2 = w2;
end
fval = w1'*X1'*X2*w2;
%fvalc = (w1'*X1'*X2*w2)/(sqrt(w1'*X1'*X1*w1)*sqrt(w2'*X2'*X2*w2));
fvalc = corr(X1*w1,X2*w2);
u = X1*w1;
v = X2*w2;
%denom1 = sqrt(w1'*X1'*X1*w1);
%denom2 = sqrt(w2'*X2'*X2*w2);
