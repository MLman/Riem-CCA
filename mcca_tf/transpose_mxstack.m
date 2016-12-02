function mxstackT = transpose_mxstack(mxstack)
mxstackT = zeros(size(mxstack));
n  = size(mxstack,3);
for i = 1:n
    mxstackT(:,:,i) = mxstack(:,:,i)';
end