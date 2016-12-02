function [t, u, r, status] = ealcca(X, Y, varargin)
    
    if nargin >= 3
        maxiter = varargin{1};
    else
        maxiter = 10;
    end
    if nargin >= 4
        stepsz0 = varargin{2};
    else
        stepsz0 = 0.01;
    end
    if nargin >= 5
        c4 = varargin{3};
    else
        c4 = 1.2;
    end
   
    mu_x = karcher_mean_spd(X,[],100);
    mu_y = karcher_mean_spd(Y,[],100);

    %% Preprocessing
    Xwr = logmap_pt2array_spd(mu_x, X);
    Ywr = logmap_pt2array_spd(mu_y, Y);
    Xv = embeddingR6_vecs(mu_x,Xwr)';
    Yv = embeddingR6_vecs(mu_y,Ywr)';
    
    %% Random Guess
%    wx = randsym(3);
%    wy = randsym(3);
%    wx = 0.5*wx/sqrt(innerprod_TpM_spd(wx,wx,mu_x));
%    wy = 0.5*wy/sqrt(innerprod_TpM_spd(wy,wy,mu_y));
    
    %% Guess by MICCAI'10 method
     [wxv, wyv] = canoncorr(Xv, Yv); 
     wx = invembeddingR6_vecs(mu_x,wxv(:,1));
     wy = invembeddingR6_vecs(mu_y,wyv(:,1));
    
    %% Optimize La
    N = size(X,3);

    fs = [];
    [fval, t, u] = myfunc(X, Y, wx, wy, mu_x, mu_y, Xwr, Ywr);
    fs = [fs ; fval];
    fprintf('First corr %f \n',fval);    
    %%
    tau_max = 1e+3;
    tau = 1e-2;
    m = 1;
    % From optimality condition c_i(xk) = -\lambda_i^{*}/m_k 
    %(eq.17.35 in Steve Wright's book).
    lx =zeros(size(t));
    ly =zeros(size(u));
    for i = 1:N
        lx(i) = -m*grad_g(mu_x, wx, X(:,:,i), t(i));
        ly(i) = -m*grad_g(mu_y, wy, Y(:,:,i), u(i));
    end

    maxiLA =20;
    maxlinesearch = 20;

    
    %% Main loop
    for iter = 1:maxiter
        % Solve LA
        [f, fLA, feasibility] = myfunc_LA(X, Y, wx, wy, mu_x, mu_y, Xwr, Ywr, t,...
                    u, lx, ly, m);
        fLA_old = fLA;
        f_old = f;
        % Gradient of LA
        [dwx_LA, dwy_LA, dt_LA, du_LA, dlx_LA, dly_LA, lx_kp1, ly_kp1] = ...
        grad_LA(wx, wy, X, Y, t, u, mu_x, mu_y, lx, ly, m);
        % norm check
        norm_dwx_LA = sqrt(innerprod_TpM_spd(dwx_LA,dwx_LA,mu_x));
        norm_dwy_LA = sqrt(innerprod_TpM_spd(dwy_LA,dwy_LA,mu_y));
        norm_dt_LA = norm(dt_LA);
        norm_du_LA = norm(du_LA);
        % Gradient descent with line search
        maxnorm = max([norm_dwx_LA, norm_dwy_LA, norm_dt_LA, norm_du_LA]);
        fprintf('>> Corr %d[%d]: f %f, fLA %f, feasibility %f, |gLA| %f, tau %f m %f\n',...
                iter , 0, f, fLA, feasibility, maxnorm, tau, m); 
            
        for iLA = 1:maxiLA
            % Gradient of LA
            [dwx_LA, dwy_LA, dt_LA, du_LA, dlx_LA, dly_LA, lx_kp1, ly_kp1] = ...
            grad_LA(wx, wy, X, Y, t, u, mu_x, mu_y, lx, ly, m);
            
            % norm check
            norm_dwx_LA = sqrt(innerprod_TpM_spd(dwx_LA,dwx_LA,mu_x));
            norm_dwy_LA = sqrt(innerprod_TpM_spd(dwy_LA,dwy_LA,mu_y));
            norm_dt_LA = norm(dt_LA);
            norm_du_LA = norm(du_LA);
            
            % Gradient descent with line search
            maxnorm = max([norm_dwx_LA, norm_dwy_LA, norm_dt_LA, norm_du_LA]);
            if maxnorm > 3
                stepsz = stepsz0/maxnorm;
            else 
                stepsz = stepsz0;
            end
            moved = false;
            for iLS = 1:maxlinesearch
                iLS
                wx_tmp = wx + stepsz*dwx_LA;
                wy_tmp = wy + stepsz*dwy_LA;
                t_tmp  = t + stepsz*dt_LA;
                u_tmp  = u + stepsz*du_LA;

                [f, fLA, feasibility] = myfunc_LA(X, Y, wx_tmp, wy_tmp, mu_x, mu_y, Xwr, Ywr, t_tmp,...
                    u_tmp, lx, ly, m);
                fprintf('Corr %d[%d][%d]: f %f, fLA %f, feasibility %f, |gLA| %f, tau %f m %f\n',...
                iter , iLA, iLS, f, fLA, feasibility, maxnorm, tau, m);                
                if  f > f_old || fLA > fLA_old 
                    wx = wx_tmp;
                    wy = wy_tmp;
                    t = t_tmp;
                    u = u_tmp;
                    moved = true;
                    fLA_old = fLA;
                    f_old = f;
                    break;
                end
                stepsz = stepsz*0.5;
            end
            fs = [fs ; f];
            fprintf('>> Corr %d[%d]: f %f, fLA %f, feasibility %f, |gLA| %f, tau %f m %f\n',...
                iter , iLA, f, fLA, feasibility, maxnorm, tau, m); 
            % Stop codition
            if norm_dwx_LA+norm_dwy_LA+norm_dt_LA+norm_du_LA < tau ||~moved 
                break;
            end
        end
        lx = lx_kp1;
        ly = ly_kp1;
        m = c4*m;
        if iLA < maxiLA
            tau = min(m*tau, tau_max);
        end
    end
    r = fs(end);
    status.fs = fs;    
end
    
    
    
    
    
    
    
    
