function X = fibonacciSphere(n)
%FIBONACCISPHERE Quasi-uniform points on the unit sphere.

if n <= 0
    X = zeros(0, 3);
    return;
end

k = (0:n-1).';
phi = (1 + sqrt(5)) / 2;
z = 1 - 2 * (k + 0.5) / n;
r = sqrt(max(1 - z.^2, 0));
theta = 2 * pi * k / phi;
X = [r .* cos(theta), r .* sin(theta), z];

