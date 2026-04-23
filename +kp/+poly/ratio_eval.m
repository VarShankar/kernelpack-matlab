function r = ratio_eval(a, b, x, N)
%RATIO_EVAL Evaluate ratios of successive orthonormal polynomials.

nx = numel(x);

assert(N > 0);
assert(N < length(a));
assert(N < length(b));

r = zeros(nx, N);
xf = x(:);

p0 = ones(nx, 1) / sqrt(b(1));
p1 = ((xf - a(1)) .* p0) / sqrt(b(2));

r1 = p1 ./ p0;
r(:, 1) = r1;

for q = 2:N
    r2 = (xf - a(q)) - sqrt(b(q)) ./ r1;
    r1 = r2 / sqrt(b(q + 1));
    r(:, q) = r1;
end

if nx == 1
    r = r(:);
end
end
