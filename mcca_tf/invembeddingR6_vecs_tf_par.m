function Vs = invembeddingR6_vecs_tf_par(ps, Vvs)
% Hyunwoo Kim
% 2014/03/05/ 14:16 pm
N = length(ps);
Vs = cell(length(Vvs),1);

parfor i = 1:N
    p = ps{i};
    Vs{i} = invembeddingR6_vecs(p,Vvs{i});
end
