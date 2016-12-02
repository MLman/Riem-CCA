function Vs = logmap_pt2array_spd_tf_par(ps,Xs)                               
% Hyunwoo 2014.03.05 11:36 am.                                               
N = length(ps);                   
Vs = cell(size(Xs));                                      
parfor i = 1:N                                                                   
    p = ps{i};                                                               
    Vs{i} = logmap_pt2array_spd(p,Xs{i});                                     
end   