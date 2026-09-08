function divfree_interp_example()
%DIVFREE_INTERP_EXAMPLE Small divergence-free interpolation demo.

theta = linspace(0, 2*pi, 33).';
theta(end) = [];
X = [0.8*cos(theta), 0.6*sin(theta)];

% A smooth divergence-free field in 2D via a streamfunction.
psi = @(x, y) x.^4 + y.^3 - 3*sin(2*pi*x).*cos(4*pi*y);
u = @(x, y) 3*y.^2 + 12*pi*sin(2*pi*x).*sin(4*pi*y);
v = @(x, y) -(4*x.^3 - 6*pi*cos(2*pi*x).*cos(4*pi*y));
U = [u(X(:, 1), X(:, 2)), v(X(:, 1), X(:, 2))];

interp = kp.divfree.LocalDivFreeInterpolator.fit(X, U, 4, 5, 15);

[xx, yy] = meshgrid(linspace(-1, 1, 61), linspace(-1, 1, 61));
Xq = [xx(:), yy(:)];
Uq = interp.evaluate(Xq);
Umag = reshape(sqrt(sum(Uq.^2, 2)), size(xx));

figure('Name', 'Div-free PHS+poly interpolation');
subplot(1, 2, 1);
quiver(X(:, 1), X(:, 2), U(:, 1), U(:, 2), 0);
axis equal;
title('Stencil samples');

subplot(1, 2, 2);
imagesc(linspace(-1, 1, 61), linspace(-1, 1, 61), Umag);
set(gca, 'YDir', 'normal');
hold on;
quiver(Xq(:, 1), Xq(:, 2), Uq(:, 1), Uq(:, 2), 2.5, 'k');
plot(X(:, 1), X(:, 2), '.', 'MarkerSize', 10);
axis equal tight;
title('Interpolated field magnitude');
colorbar;
end
