function P = randspd(n, varargin)
% Random SPD
if nargin == 2
    c = varargin{1};
else
    c = 3;
end
P = c*(rand(n)-0.5);
P = P*P';
