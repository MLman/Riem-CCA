function r = isreal_cell(X)
n = length(X);
r = zeros(size(X));
for i=1:n
    r(i) = isreal(X{i});
end
