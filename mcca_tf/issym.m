function r =  issym(mx)
    if ndims(mx) <= 2
        r = issymsub(mx);
    elseif ndims(mx) == 3
        r = zeros(size(mx,3),1);
        for i =1:size(mx,3)
            r(i) = issymsub(mx(:,:,i));
        end
    else
        error('Not Implemented.');
    end
end

function r = issymsub(mx)
    r = (sum(sum(abs(mx -mx'))) == 0);
end
