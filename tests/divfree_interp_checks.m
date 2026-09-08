function divfree_interp_checks()
%DIVFREE_INTERP_CHECKS Targeted checks for the divergence-free PHS+poly port.

% Check the kernel block builder sizes in 3D.
[dfrbf, lapDfrbf, hessRbf, fullRbf] = kp.divfree.DFPHS(3, 5); %#ok<ASGLU>
X3 = [0 0 0; 1 0 0; 0 1 0];
Y3 = [0 0 1; 1 1 0];
r3 = kp.geometry.distanceMatrix(X3, Y3);
A3 = kp.divfree.DivFreeGram(dfrbf, X3, Y3, r3);
assert(isequal(size(A3), [9, 6]));

% Build a 2D divergence-free polynomial field from a streamfunction and
% confirm exact reproduction off-grid when the polynomial degree is high enough.
theta = linspace(0, 2*pi, 13).';
theta(end) = [];
X = [0.7*cos(theta), 0.5*sin(theta)] + 0.15 * [cos(2*theta), sin(3*theta)];
x = X(:, 1);
y = X(:, 2);
U = [3*x.*y.^2, -(4*x.^3 + y.^3)];

interp = kp.divfree.DivFreePHSInterpolant.fit(X, U, 3, 5);
Pself = interp.PolyData.eval(X);
[Pbuilt, ~] = kp.divfree.df_poly_basis_from_jacobi(X, 3, @(N) kp.poly.jacobi_recurrence(N, 0, 0), struct());
assert(norm(Pself - Pbuilt, inf) < 1e-10);

Xq = [0.10 0.20; -0.20 0.15; 0.25 -0.10; -0.10 -0.20];
xq = Xq(:, 1);
yq = Xq(:, 2);
Uexact = [3*xq.*yq.^2, -(4*xq.^3 + yq.^3)];
Uq = interp.evaluate(Xq);
assert(norm(Uq - Uexact, inf) < 1e-10);

% Check the cached local interpolator path against the same exact field.
localInterp = kp.divfree.LocalDivFreeInterpolator.fit(X, U, 3, 5, 9);
UqLocal = localInterp.evaluate(Xq);
assert(norm(UqLocal - Uexact, inf) < 1e-10);
assert(numel(localInterp.ActiveCenters) == size(X, 1));

disp('divfree interp checks passed');
end
