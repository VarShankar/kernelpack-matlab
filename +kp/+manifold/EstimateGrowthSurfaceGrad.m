function [q1, q2, q3, nw] = EstimateGrowthSurfaceGrad(Gx, Gy, Gz, taux, tauy, tauz, X, nr, h)
%ESTIMATEGROWTHSURFACEGRAD Estimate spurious gradient growth exponents.

if isempty(h)
    nw = 10000;
else
    nw = 2 / h;
end

f = zeros(size(X, 1), 1);
gcx = zeros(size(X, 1), 1);
gcy = zeros(size(X, 1), 1);
gcz = zeros(size(X, 1), 1);

if size(X, 2) == 2
    f(:, 1) = exp(1i * nw * (X(:, 1) + X(:, 2)));
    gx_true = 1i * nw * exp(1i * nw * X(:, 1)) .* exp(1i * nw * X(:, 2));
    gy_true = gx_true;
    nx = nr(:, 1);
    ny = nr(:, 2);
    gcx(:, 1) = (1 - nx.^2) .* gx_true(:, 1) - (nx .* ny) .* gy_true(:, 1);
    gcy(:, 1) = -(nx .* ny) .* gx_true(:, 1) + (1 - ny.^2) .* gy_true(:, 1);
elseif size(X, 2) == 3
    f(:, 1) = exp(1i * nw * X(:, 1)) .* exp(1i * nw * X(:, 2)) .* exp(1i * nw * X(:, 3));
    gx_true = 1i * nw * exp(1i * nw * X(:, 1)) .* exp(1i * nw * X(:, 2)) .* exp(1i * nw * X(:, 3));
    gy_true = gx_true;
    gz_true = gx_true;
    nx = nr(:, 1);
    ny = nr(:, 2);
    nz = nr(:, 3);
    gcx(:, 1) = (1 - nx.^2) .* gx_true(:, 1) - (nx .* ny) .* gy_true(:, 1) - (nx .* nz) .* gz_true(:, 1);
    gcy(:, 1) = -(nx .* ny) .* gx_true(:, 1) + (1 - ny.^2) .* gy_true(:, 1) - (ny .* nz) .* gz_true(:, 1);
    gcz(:, 1) = -(nz .* nx) .* gx_true(:, 1) - (nz .* ny) .* gy_true(:, 1) + (1 - nz.^2) .* gz_true(:, 1);
else
    error('kp:manifold:BadDimension', 'EstimateGrowthSurfaceGrad expects 2D or 3D ambient coordinates.');
end

gx = Gx * f;
gy = Gy * f;
gz = Gz * f;

q1 = (log(norm(gx(:, 1) - gcx(:, 1), 2)) - log(taux) - log(norm(f(:, 1), 2))) ./ log(nw);
q2 = (log(norm(gy(:, 1) - gcy(:, 1), 2)) - log(tauy) - log(norm(f(:, 1), 2))) ./ log(nw);
q3 = (log(norm(gz(:, 1) - gcz(:, 1), 2)) - log(tauz) - log(norm(f(:, 1), 2))) ./ log(nw);
q1 = real(q1);
q2 = real(q2);
q3 = real(q3);
end
