function T = isspd(mx,varargin)
%ISSPD check MX is a symmetric positive definite matrix.
%    This check the smallest eigen value is bigger than c.
%    Default c is epsilon.

if nargin ==2
    c = varargin{1};
else 
    c = eps;
end
% Check matrices are symmetric positive definite.
T = zeros(size(mx,3),1);
for i=1:size(mx,3)
    T(i) = (sum(eig(mx(:,:,i)) <= 0+c ) ==0);
end