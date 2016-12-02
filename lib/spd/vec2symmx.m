function mx = vec2symmx(v)
%VEC2SYMMX converts vectors V to matrices MX.
%   n(n+1)/2 dimensional vectors to n by n matrices.
    [dimv ndata] = size(v);
    n = (-1 + sqrt(1+8*dimv))/2;
    mx = zeros(n,n,ndata);
    k = 1;
    for i=1:n
        for j=i:n
            mx(i,j,:) = v(k,:);
            if i ~=j
                mx(j,i,:) = v(k,:);
            end
            k = k + 1;
        end
    end
end