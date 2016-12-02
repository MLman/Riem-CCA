function r = norm_TpM_spd(P,V)
    r = sqrt(innerprod_TpM_spd(V,V,P));
end