function Vs = logmap_pt2array_spd_tf(ps,Xs)
N = length(ps);
Vs = cell(size(Xs));
for i = 1:N
    p = ps{i};
    Vs{i} = logmap_pt2array_spd(p,Xs{i});
end
