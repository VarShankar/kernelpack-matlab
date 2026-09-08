function meshfree_geometric_multilevel_checks()
%MESHFREE_GEOMETRIC_MULTILEVEL_CHECKS Basic checks for MGM plumbing.

[A, X] = gridProblem(6);
solver = kp.solvers.MeshfreeGeometricMultilevelSolver(A, X, [], ...
    'Nmin', 20, ...
    'UseParallel', false);

assert(solver.StencilSize == 3 && solver.PolynomialDegree == 0, ...
    'MGM defaults should match Grady Wright''s stencilSizeT/polyDegT choices.');
assert(solver.DomainVolume > 0, 'MGM plumbing should estimate a domain volume when none is supplied.');
assert(isscalar(solver.LevelsData), 'MGM plumbing should retain only the fine-level descriptor.');
assert(isequal(size(solver.LevelsData.Lh), size(A)), 'MGM fine-level descriptor should retain the input matrix.');
assert(isequal(size(solver.LevelsData.nodes), size(X)), 'MGM fine-level descriptor should retain the input point cloud.');

try
    solver.solve(ones(size(A, 1), 1));
    error('kp:tests:ExpectedMGMUnavailable', 'MGM solve should report that the algorithm is unavailable.');
catch ME
    assert(strcmp(ME.identifier, 'kp:solvers:MGMAlgorithmUnavailable'), ...
        'MGM solve should fail with the explicit unavailable-algorithm identifier.');
end

try
    solver.precondition(ones(size(A, 1), 1));
    error('kp:tests:ExpectedMGMUnavailable', 'MGM precondition should report that the algorithm is unavailable.');
catch ME
    assert(strcmp(ME.identifier, 'kp:solvers:MGMAlgorithmUnavailable'), ...
        'MGM precondition should fail with the explicit unavailable-algorithm identifier.');
end

nullSolver = kp.solvers.MeshfreeGeometricMultilevelSolver(A, X, 1.0, ...
    'HasConstNullspace', true, ...
    'UseParallel', false);
state = [(1:size(A, 1)).'; 2];
constrainedProduct = nullSolver.applyConstrainedOperator(state);
expectedProduct = [A * state(1:end-1) + 2 * ones(size(A, 1), 1); sum(state(1:end-1))];
assert(norm(constrainedProduct - expectedProduct) < 1e-12, ...
    'MGM constrained operator plumbing should preserve the augmented matrix action.');

disp('meshfree geometric multilevel plumbing checks passed');
end

function [A, X] = gridProblem(n)
[xx, yy] = ndgrid(linspace(0, 1, n), linspace(0, 1, n));
X = [xx(:), yy(:)];
A = speye(n * n);
end
