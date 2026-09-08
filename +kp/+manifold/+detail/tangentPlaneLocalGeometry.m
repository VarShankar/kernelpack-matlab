function geom = tangentPlaneLocalGeometry(x, nrl)
%TANGENTPLANELOCALGEOMETRY Shared tangent-plane frame and local coordinates.

R = tangentPlaneBasis(nrl(:));
xw = (x - x(1, :)) * R;
rd = hypot(xw(:, 1) - xw(:, 1).', xw(:, 2) - xw(:, 2).');
w2 = max(max(abs(xw - xw(1, :))));
if w2 <= eps
    error('kp:manifold:DegenerateStencil', 'Tangent-plane stencil has zero local width.');
end

geom = struct();
geom.R = R;
geom.xw = xw;
geom.rd = rd;
geom.xc = (xw - xw(1, :)) ./ w2;
geom.rdc = rd ./ w2;
geom.diffxe = xw - xw(1, :);
geom.w2 = w2;
end

function R = tangentPlaneBasis(nrl)
[~, argmax] = max(abs(nrl));
if argmax == 1
    e2 = [0; 1; 0];
    e3 = [0; 0; 1];
elseif argmax == 2
    e2 = [1; 0; 0];
    e3 = [0; 0; 1];
else
    e2 = [1; 0; 0];
    e3 = [0; 1; 0];
end
t1 = e2 - (nrl' * e2) * nrl;
t1 = t1 / norm(t1);
t2 = e3 - (nrl' * e3) * nrl - (t1' * e3) * t1;
t2 = t2 / norm(t2);
R = [t1, t2];
end
