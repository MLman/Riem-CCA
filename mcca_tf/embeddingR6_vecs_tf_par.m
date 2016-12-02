function Vvs = embeddingR6_vecs_tf_par(ps,Vs)
% Hyunwoo Kim.
% 2014/03/05 11:47 am.
%EMBEDDINGR6_VECS_TF embeds Vs in R6 spaces
N = length(ps);
Vvs = cell(length(Vs),1);  

parfor i = 1:N
    p = ps{i};
    Vvs{i} = embeddingR6_vecs(p,Vs{i});
end
