function Vs = invembeddingR6_vecs_tf(ps, Vvs)
%
N = length(ps);
Vs = cell(length(Vvs),1);

for i = 1:N
    p = ps{i};
    Vs{i} = invembeddingR6_vecs(p,Vvs{i});
end