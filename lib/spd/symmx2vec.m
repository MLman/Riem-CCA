function v = symmx2vec(mx)
%SYMMX2VEC converts matrices MX to vectors V. 
%   n by n matrices to n(n+1)/2 dimensional vectors.
    [ nrow ncol ndata ] = size(mx);
    v = zeros(nrow*(nrow+1)/2,ndata);
    k =1;
    for i=1:ncol
        for j=i:ncol
            v(k,:) = squeeze(mx(i,j,:))';
            k = k + 1;
        end

    end
end