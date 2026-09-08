function X = buildPlanarParametricNodes2D(cart, N, equispaced, custom)
%BUILDPLANARPARAMETRICNODES2D Best-fit-plane parameter sites for a 3D patch.

if nargin < 3
    equispaced = 0;
end
if nargin < 4
    custom = 0;
end

if equispaced == 1
    [Xp, ~, ~] = kp.geometry.projectToBestFitPlane(cart);
    x0 = min(Xp(:, 1));
    xm = max(Xp(:, 1));
    y0 = min(Xp(:, 2));
    ym = max(Xp(:, 2));
    m = max(2, floor(sqrt(max(N, 4))));
    [xx, yy] = meshgrid(linspace(x0, xm, m), linspace(y0, ym, m));
    X = [xx(:), yy(:)];
    return;
end

if custom == 1
    [X, ~, ~] = kp.geometry.projectToBestFitPlane(cart);
    return;
end

X = kp.geometry.buildPlanarParametricEvalNodes2D(cart, N);

