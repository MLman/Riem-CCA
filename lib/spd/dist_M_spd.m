function d = dist_M_spd(X,Y)
    V = logmap_spd(X,Y);
    d = sqrt(innerprod_TpM_spd(V,V,X));
    
    %% This version will be faster.
    %rtX = sqrtm(X);
    %tmp = logm(rtX*Y*rtX);
    %d = sqrt(trace(tmp*tmp));
end