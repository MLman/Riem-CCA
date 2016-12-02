function Vvs = embeddingR6_vecs_tf(ps,Vs)
%EMBEDDINGR6_VECS_TF embeds Vs in R6 spaces
N = length(ps);
Vvs = cell(length(Vs),1);

for i = 1:N
    p = ps{i};
    Vvs{i} = embeddingR6_vecs(p,Vs{i});
end