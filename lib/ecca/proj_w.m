function [t, z, fs, ggs, ts] = proj_w(p, w, x, xwr, varargin)

    if nargin == 5
        niter = varargin{1};
    else
        niter = 100;
    end
    
    c1 = 1; % maximum step size

    % initial guess
    normw = sqrt(innerprod_TpM_spd(w,w,p));
    uw = w/normw;

    if isempty(xwr)
        error('Not Implemented.');
        %xwr = logmap_spd(p,x);
    end
    
    t = innerprod_TpM_spd(uw,xwr,p);
    zwr = uw*t;
    z = expmap_spd(p,zwr);
    step_0 = 0.1;
    t = t/normw;
    maxchk = 50;
    fs = [];
    ggs = [];
    ts = [];
    for iter = 1:niter
        ts = [ts;t];
        fval = dist_M_spd(z,x);
        fs = [fs;fval];
        
        [gg, St, dtSt] = grad_g(p, w, x, t);
        ggs = [ggs; gg];

        step = step_0;
        
        if abs(gg) > c1
            gg = gg*c1/abs(gg);
        end
        for ichk = 1:maxchk
            if geval(t-step*gg, w, p, x) < fval
                t = t-step*gg;
                z = expmap_spd(p, t*w);
                break
            else
                step = step*1/2;
            end
            if step*gg < 1e-18
            end
        end
        
        if ichk == maxchk
            break
        end
    end
    
end


