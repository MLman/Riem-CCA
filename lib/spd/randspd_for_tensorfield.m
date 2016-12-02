function T = randspd_for_tf(nvoxels, ndim, nsample)
%RANDSPD_FOR_TENSORFIELD generates data for tensor field experiments.
%    nvoxels is a number of voxels or cells.
%    ndim is the dimension of SPD matrices.
%    nsample is the number of samples
T = cell(nvoxels,1);
for ic = 1:nvoxels
    tmp = zeros(ndim,ndim,nsample);
    for i=1:nsample
       tmp(:, :, i) = randspd(ndim);
    end
    T{ic} = tmp;
end