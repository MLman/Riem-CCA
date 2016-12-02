function xbar = karcher_mean_spd(X, W, niter)
% KARCHER_MEAN calculates weighted means on d-sphere.
% W is weight
% X is column vectors
% niter is max iteration
xbar = X(:,:,1);

if isempty(W)
    for iter = 1:niter
        phi = mean(logmap_pt2array_spd(xbar,X),3);
        xbar = expmap_spd(xbar, phi);
        if norm(phi) < 1e-10
            break
        end
    end
else
    W = W/norm(W,1);
    for iter = 1:niter
        tmp = logmap_pt2array_spd(xbar,X);
        wtmp = zeros(size(tmp));
        for i = 1:size(tmp,3)
            wtmp(:,:,i) = W(i)*tmp(:,:,i);
        end
        phi = sum(wtmp,3);
        xbar = expmap_spd(xbar, phi);
        if norm(phi) < 1e-18
            break
        end
    end
end
