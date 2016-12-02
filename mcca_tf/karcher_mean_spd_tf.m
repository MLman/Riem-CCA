function ms = karcher_mean_spd_tf(Xs)
%KARCHER_MEAN_SPD_TF finds karcher mean of each space.
%    X consists of multiple spaces. Each space is in a cell.

n = length(Xs);
ms = cell(size(Xs));

parfor i=1:n
    ms{i} = karcher_mean_spd(Xs{i},[],100);
end