function exp_p_x = expmap_spd(P,X)
%EXPMAP_SPD is exponential map for GL(n)/O(n).
%    P is a base matrix. It is a symmetric positive definite matrix.
%    V is a tangent vector. It is a symmetric matrix.
    
    if norm(X) < 1e-18
        exp_p_x = P;
        return
    end
    [U D] = eig(P);
    g = U*sqrt(D);
%    invg = inv(g);
    invg = diag(1./sqrt(diag(D)))*U'; % 1.3 X faster
    
    Y = invg*X*invg';
    [V S] = eig(Y);
    gv = g*V;
    exp_p_x = gv*diag(exp(diag(S)))*gv';
    exp_p_x = (exp_p_x+exp_p_x')/2;
%    rtP = sqrtm(P);
%    invrtP = inv(rtP);
%    exp_p_v = rtP*expm(invrtP*V*invrtP)*rtP;
end