function Ps = randspds(n, nmx, varargin)
%RANDSPDS    
%NMX    number of matrix.

% Random SPD
if nargin == 3
    c = varargin{1};
else
    c = 3;
end

Ps = zeros(n,n,nmx);

for k=1:nmx
    P = c*(rand(n)-0.5);
    P = P*P';
    Ps(:,:,k) = P;
end
