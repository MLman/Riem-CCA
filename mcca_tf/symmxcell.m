function Sc =  symmxcell(mxcells)
n = length(mxcells);
Sc = cell(size(mxcells));
for i=1:n
    Sc{i} = symmxstack(mxcells{i});
end