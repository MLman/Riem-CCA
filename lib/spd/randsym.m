function S = randsym(n,varargin)
% Random Symmetric matrices
S = rand(n)-0.5;
S = (S +S')/2;
if nargin == 2
    c = varargin{1};
else
    c = 5;
end
S = c*rand(1)*S;