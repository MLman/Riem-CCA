function S = symmxstack(mxstack)
%SYMMXSTACK symmetrize mxstack
    mxstackT = transpose_mxstack(mxstack);
    S = (mxstackT+mxstack)/2;