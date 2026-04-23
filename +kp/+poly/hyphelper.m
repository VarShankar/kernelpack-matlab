function [v, hv] = hyphelper(dim, hyppow, ell)
%HYPHELPER Symbolic helper for monomial and repeated-Laplacian evaluation.
%
% This is a direct MATLAB-oriented port of the helper used with the older
% polynomial support code. It returns function handles produced from symbolic
% expressions.

if dim == 2
    syms x y
    cnt = 0;
    npoly = nchoosek(ell + dim, dim);
    hv = sym(zeros(1, npoly));
    v = hv;
    for i = 0:ell
        for j = 0:ell
            if (i + j) <= ell
                cnt = cnt + 1;
                v(:, cnt) = x.^i .* y.^j;
                hv(:, cnt) = v(:, cnt);
                for it = 1:hyppow
                    hv(:, cnt) = laplacian(hv(cnt), [x, y]);
                end
            end
        end
    end
    v = symfun(v, [x, y]);
    hv = symfun(hv, [x, y]);
elseif dim == 3
    syms x y z
    cnt = 1;
    v = sym([]);
    hv = sym([]);
    for k = 0:ell
        for l = k:-1:0
            for m = (k - l):-1:0
                v(:, cnt) = x.^l .* y.^m .* z.^(k - l - m);
                hv(:, cnt) = v(:, cnt);
                for it = 1:hyppow
                    hv(:, cnt) = laplacian(hv(cnt), [x, y, z]);
                end
                cnt = cnt + 1;
            end
        end
    end
    v = symfun(v, [x, y, z]);
    hv = symfun(hv, [x, y, z]);
else
    error('kp:poly:BadDimension', 'hyphelper currently supports only dim = 2 or 3.');
end

v = matlabFunction(v);
hv = matlabFunction(hv);
end
