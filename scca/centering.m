function Xc = centering(X)
Xm = mean(X,2);
Xc = X - repmat(Xm,1,size(X,2));
