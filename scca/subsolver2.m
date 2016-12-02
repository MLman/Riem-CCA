function [w th ] = subsolver2(p, c)
%SUBSOLVER2 solves a subproblem in Sparse CCA
%    w = argmin_{w} p'w1 
%         subject to ||w||^{2} <= 1
%                       ||w || <= c
%
%          This subsolver give w which may have positive or negative
%          weights
%
%   See also SUBSOLVER

absp = abs(p);
sgn = sign(p);
[v idx] = sort(absp, 'descend');
diffv = [v(1:end-1)-v(2:end);v(end)];
diffv = diffv.*(1:length(diffv))';
cdiff = cumsum(diffv);
n = find(cdiff - c >=0, 1, 'first'); % n is the  number of survivals
if ~isempty(n)
    th = (sum(v(1:n)) -c)/n;
    % Packing the result
    w = zeros(size(p));
    w(idx(1:n)) = sgn(idx(1:n)).*(v(1:n)-th);
else 
    w = p; % Too loose sparsity constraint. No effect of sparsity constraint.
end

