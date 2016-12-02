function exp_p_x = expmap_pt2array_spd(p, X)
%EXPMAP_pt2array_spd is exponential map for GL(n)/O(n).


    

    % Step 1 common step
    [U D] = eig(p);
    g = U*sqrt(D);
%    invg = inv(g);
    invg = diag(1./sqrt(diag(D)))*U'; % 1.3 X faster    
    exp_p_x = zeros(size(X));

    % For each data

    for i =1:size(X,3)
        if norm(X(:,:,i)) < 1e-18
            exp_p_x(:,:,i) = p;
            continue
        end
        Y = invg*X(:,:,i)*invg';
        [V S] = eig(Y);
        gv = g*V;
        exp_p_x(:,:,i) = gv*diag(exp(diag(S)))*gv';
    end 
   
end