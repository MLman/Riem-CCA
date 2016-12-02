function mat = mxstack2mat(mxstack)
%MXSSTACK2MAT converts a matrix stack to a matrice of row vectors.
%
%    MXSTACK is a matrice stack. n-by-n, m.
%    MAT is m-by-n^2
[nrows, ncols, nslices ] = size(mxstack);
mat = zeros(nslices, nrows*(nrows-1)/2);
i  = 1;
for irow = 1:nrows
    for icol = irow:ncols
        mat(:,i) = squeeze(mxstack(irow,icol,:));
        i = i + 1;
    end
end
