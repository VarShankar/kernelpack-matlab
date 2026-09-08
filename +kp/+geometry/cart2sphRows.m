function uv = cart2sphRows(X)
%CART2SPHROWS Row-wise Cartesian to spherical coordinates.

[az, el, r] = cart2sph(X(:, 1), X(:, 2), X(:, 3));
uv = [az, el, r];

