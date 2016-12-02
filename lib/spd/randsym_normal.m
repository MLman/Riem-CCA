function S = randsym_normal(n,varargin)
% Random Symmetric matrices
S = vec2symmx(randn(n*(n+1)/2,1));
