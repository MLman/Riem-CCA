function r = issym_cell(X)
n = length(X);
r = zeros(size(X));
for i=1:n
    r(i) = prod(issym(X{i}));
end


