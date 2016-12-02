function fval = geval(t, w, p, x)
%GEVAL returns projection error.
%    projection error evaluation

    zwr = w*t;
    z = expmap_spd(p,zwr);
    fval = dist_M_spd(z,x);
end
