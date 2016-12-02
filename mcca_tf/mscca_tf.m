function [wx, wy, r, u, v, precomp] = mscca_tf(X,Y, varargin)
%MSCCA_TF Riemannian CCA on tensor field. 
precomp = [];
if nargin >= 3
    c1 = varargin{1};
else 
    c1 = sqrt(size(X,1)*6);
end

if nargin >= 4
    c2 = varargin{2};
else 
    c2 = sqrt(size(Y,1)*6);
end

if nargin >= 5
    c3 = varargin{3};
else 
    c3 =1e-20;
end
if nargin >= 6
    maxiter = varargin{4};
else 
    maxiter =100;
end

if nargin >=7 
     precomp = varargin{5};
end
if ~isempty(precomp)
    Xe = precomp.Xe;
    Ye = precomp.Ye;
    muX = precomp.muX;
    muY = precomp.muY;
    % Logarithm map
%     Xwrs = precomp.Xwrs;
%     Ywrs = precomp.Ywrs;
%     % Embedding product space of Euclidean space (R6)
%     Xvecs = precomp.Xvecs;
%     Yvecs = precomp.Yvecs;
else
    % Calculate the karcher mean of each tensor field
    disp('Karcher Mean Calculation');
    muX = karcher_mean_spd_tf(X);
    muY = karcher_mean_spd_tf(Y);

    % Logarithm map
    disp('Logmap Calculation');
    Xwrs = logmap_pt2array_spd_tf(muX,X);
    Ywrs = logmap_pt2array_spd_tf(muY,Y);

    % Embedding product space of Euclidean space (R6)
    disp('Embedding');
    Xvecs = embeddingR6_vecs_tf(muX,Xwrs);
    Yvecs = embeddingR6_vecs_tf(muY,Ywrs);

    % Embedding in one Euclidean space
    disp('Xe, Xe');
    Xe = cat(1,Xvecs{:});
    Ye = cat(1,Yvecs{:});
end


try
    assert(isreal(Xe)==1);
catch
    fprintf('Numerical Error.');
    Xe = real(Xe);
end
try
    assert(isreal(Ye)==1);
catch
    fprintf('Numerical Error.');
    Ye = real(Ye);
end
disp('Perform SCCA') 
tic;
[Wxe, Wye, fval, r, U, V] = scca_ver2(Xe',Ye', c1, c2, c3, maxiter);
toc;
% Convert Axes
wxvec = num2cell(reshape(Wxe(:,1),6,[]),1)';
wyvec = num2cell(reshape(Wye(:,1),6,[]),1)';

wx = invembeddingR6_vecs_tf(muX, wxvec);
wy = invembeddingR6_vecs_tf(muY, wyvec);

% Correlation of the projection with the first axis
%r = r(1);

%% Need to be implemented.
u = Xe'*Wxe;
v = Ye'*Wye;
r = corr(u,v);

%% Packaging precomp
if isempty(precomp)
    precomp.Xe = Xe;
    precomp.Ye = Ye;
    precomp.muX = muX;
    precomp.muY = muY;
end
