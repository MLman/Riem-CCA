function [X, Y, group1, group2] = generate_data_for_demo1(N, dimX, dimY, c1,c2)    

    nG1 = N/2;    
    nG2 = N/2;

    Xgroup1_mean = randspd(dimX);
    Xgroup2_mean = randspd(dimX);
    X = zeros(dimX,dimX,N);

    for isample = 1:nG1
        % Random additive noise in the tangent space (Sym(dimX)).
        noise = c1*randsym_normal(dimX);
        norm_noise = norm_TpM_spd(Xgroup1_mean, noise);
        if  norm_noise > c2
            noise  = noise / norm_noise;
        end
        X(:,:,isample) = expmap_spd(Xgroup1_mean, noise);
    end
            
    for isample = 1:nG2
        noise = c1*randsym_normal(dimX);
        norm_noise = norm_TpM_spd(Xgroup2_mean, noise);

        if  norm_noise > c2
            noise  = noise / norm_noise;
        end
        X(:,:,nG1+isample) = expmap_spd(Xgroup2_mean, noise);
    end
    % Generate two groups of samples in the space of Y.
	Ygroup1_mean = randspd(dimY);
    Ygroup2_mean = randspd(dimY);
    Y = zeros(dimY,dimY,N);
    for isample = 1:nG1
        % Random additive noise in the tangent space (Sym(dimX)).
        noise = c1*randsym_normal(dimY);
        norm_noise = norm_TpM_spd(Ygroup1_mean, noise);
        if  norm_noise > c2
            noise  = noise / norm_noise;
        end
        Y(:,:,isample) = expmap_spd(Ygroup1_mean, noise);
    end
    for isample = 1:nG2
        noise = c1*randsym_normal(dimY);
        norm_noise = norm_TpM_spd(Ygroup2_mean, noise);
        if  norm_noise > c2
            noise  = noise / norm_noise;
        end
        Y(:,:,nG1+isample) = expmap_spd(Ygroup2_mean, noise);
    end
    group1 = 1:nG1;
    group2 = (nG1+1):(nG1+nG2);
end