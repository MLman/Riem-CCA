function Vnew = embeddingRd_vecs(p,V)
%EMBEDDINGR6_VECS embeds V in TpM onto R6.
%    
%    p is a base point. V is tangent vectors.
%    V is d-by-d-n. n is the number of matrices.

        XwrI = grpaction_p2i(mu_x, Xwr);
        Vnew = symmx2vec(XwrI)';
end
