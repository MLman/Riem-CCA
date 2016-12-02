function [wx, wy, r, u, v, muX, muY] = mcca_tf(X,Y)
% Calculate the karcher mean of each tensor field
muX = karcher_mean_spd_tf(X);
muY = karcher_mean_spd_tf(Y);

% Embedding product space of Euclidean space (R6)
Xwrs = embeddingR6_vecs_tf(muX,X);
Ywrs = embeddingR6_vecs_tf(muY,Y);

% Embedding in one Euclidean space
Xe = cat(1,Xwrs{:});
Ye = cat(1,Ywrs{:});
[Wxe, Wye, r, U, V] = canoncorr(Xe',Ye');
assert(isreal(Wxe)==1);
assert(isreal(Wye)==1);

% Convert Axes
wxvec = num2cell(reshape(Wxe(:,1),6,[]),1)';
wyvec = num2cell(reshape(Wye(:,1),6,[]),1)';
%assert(isreal(wxvec)==1);
%assert(isreal(wyvec)==1);
wx = invembeddingR6_vecs_tf(muX, wxvec);
wy = invembeddingR6_vecs_tf(muY, wyvec);

% Correlation of the projection with the first axis
r = r(1);

% projection result of the first axis
u = U(:,1);
v = V(:,1);