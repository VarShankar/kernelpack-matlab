classdef PoissonSolver < handle
    %POISSONSOLVER Simple fixed-domain Poisson solver in the KernelPack shape.

    properties
        LapAssembler = "fd"
        BCAssembler = "fd"
        LapStencil = "rbf"
        BCStencil = "rbf"
        Domain kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
        xi (1,1) double = 0
        num_omp_threads (1,1) double = 1
        X double = zeros(0, 0)
        Xb double = zeros(0, 0)
        nr double = zeros(0, 0)
        N (1,1) double = 0
        Nf (1,1) double = 0
        Lap double = zeros(0, 0)
        BC double = zeros(0, 0)
        LapStencilProperties kp.rbffd.StencilProperties = kp.rbffd.StencilProperties()
        BCStencilProperties kp.rbffd.StencilProperties = kp.rbffd.StencilProperties()
        LapOpProperties kp.rbffd.OpProperties = kp.rbffd.OpProperties('decompose', false, 'storeWeights', true, 'recordStencils', false)
        BCOpProperties kp.rbffd.OpProperties = kp.rbffd.OpProperties('decompose', false, 'storeWeights', true, 'recordStencils', false)
        last_solve_used_nullspace_ (1,1) logical = false
    end

    methods
        function obj = PoissonSolver(varargin)
            if nargin == 0
                return;
            end
            parser = inputParser();
            parser.addParameter('LapAssembler', obj.LapAssembler);
            parser.addParameter('BCAssembler', obj.BCAssembler);
            parser.addParameter('LapStencil', obj.LapStencil);
            parser.addParameter('BCStencil', obj.BCStencil);
            parser.parse(varargin{:});
            obj.LapAssembler = parser.Results.LapAssembler;
            obj.BCAssembler = parser.Results.BCAssembler;
            obj.LapStencil = parser.Results.LapStencil;
            obj.BCStencil = parser.Results.BCStencil;
        end

        function init(obj, domain, xi, num_omp_threads)
            if nargin < 4 || isempty(num_omp_threads)
                num_omp_threads = 1;
            end

            obj.Domain = domain;
            obj.xi = xi;
            obj.num_omp_threads = num_omp_threads;

            % The solver only needs one cached Laplacian; boundary rows are
            % rebuilt per solve because the BC coefficients may change.
            obj.Domain.buildStructs();
            obj.X = obj.Domain.getIntBdryNodes();
            obj.Xb = obj.Domain.getBdryNodes();
            obj.nr = obj.Domain.getNrmls();
            obj.N = size(obj.X, 1);
            obj.Nf = obj.Domain.getNumTotalNodes();

            obj.LapStencilProperties = buildStencilProperties(obj.Domain, obj.xi, 2, "interior_boundary");
            obj.BCStencilProperties = buildStencilProperties(obj.Domain, obj.xi, 1, "boundary");

            if obj.num_omp_threads > 1
                obj.LapOpProperties.UseParallel = true;
                obj.BCOpProperties.UseParallel = true;
            end

            lapAssembler = makeAssembler(obj.LapAssembler, obj.LapStencil);
            lapAssembler.AssembleOp(obj.Domain, 'lap', obj.LapStencilProperties, obj.LapOpProperties);
            obj.Lap = lapAssembler.getOp();

            obj.BC = zeros(0, obj.Nf);
            obj.last_solve_used_nullspace_ = false;
        end

        function result = solve(obj, forcing, NeuCoeffFunc, DirCoeffFunc, bc, varargin)
            parser = inputParser();
            parser.addOptional('InitialGuess', zeros(0, 1), @(x) isnumeric(x));
            parser.parse(varargin{:});
            initialGuess = parser.Results.InitialGuess(:);

            neuCoeff = evaluateNodeCallback(NeuCoeffFunc, obj.Xb, 'boundary coefficient');
            dirCoeff = evaluateNodeCallback(DirCoeffFunc, obj.Xb, 'boundary coefficient');

            % Assemble the boundary operator for the requested BC mix, then
            % pair it with the cached Laplacian in one square ghost-node
            % system.
            bcAssembler = makeAssembler(obj.BCAssembler, obj.BCStencil);
            bcAssembler.AssembleOp(obj.Domain, 'bc', obj.BCStencilProperties, obj.BCOpProperties, ...
                'NeuCoeff', neuCoeff, ...
                'DirCoeff', dirCoeff);
            obj.BC = bcAssembler.getOp();

            rhsTarget = evaluateNodeCallback(forcing, obj.X, 'forcing');
            rhsBoundary = evaluateBoundaryValues(bc, neuCoeff, dirCoeff, obj.nr, obj.Xb);
            % Pure Neumann problems need the usual one-dimensional
            % nullspace fix, so augment the system with a constant row and
            % column when the Dirichlet coefficient is identically zero.
            pureNeumann = max(abs(dirCoeff)) <= 1e-13;
            obj.last_solve_used_nullspace_ = pureNeumann;

            A = buildSystemMatrix(obj.Lap, obj.BC, obj.Nf, pureNeumann);
            b = buildSystemRHS(rhsTarget, rhsBoundary, pureNeumann);
            guess = buildInitialGuess(initialGuess, obj.N, obj.Nf, rhsBoundary, pureNeumann);

            if isempty(guess)
                sol = A \ b;
            else
                sol = gmresWithFallback(A, b, guess);
            end

            if pureNeumann
                fullState = sol(1:obj.Nf);
                lagrangeMultiplier = sol(end);
            else
                fullState = sol;
                lagrangeMultiplier = [];
            end

            result = struct();
            result.u = fullState(1:obj.N);
            result.FullState = fullState;
            result.L = obj.Lap;
            result.BC = obj.BC;
            result.SystemMatrix = A;
            result.RHS = b;
            result.TargetRHS = rhsTarget;
            result.BoundaryRHS = rhsBoundary;
            result.UsedNullspaceAugmentation = pureNeumann;
            result.LagrangeMultiplier = lagrangeMultiplier;
        end

        function out = getLaplacian(obj)
            out = sparse(obj.Lap);
        end

        function out = getBCOp(obj)
            out = sparse(obj.BC);
        end

        function out = lastSolveUsedNullspace(obj)
            out = obj.last_solve_used_nullspace_;
        end
    end
end

function sp = buildStencilProperties(domain, xi, theta, pointSet)
% Translate a target convergence order into the local stencil settings
% used by the RBF-FD layer.
dim = domain.getDim();
sp = kp.rbffd.StencilProperties();
sp.dim = dim;
sp.ell = max(xi + theta - 1, 2);
sp.npoly = size(kp.poly.total_degree_indices(dim, sp.ell), 1);
sp.n = 2 * sp.npoly + 1;
sp.spline_degree = sp.ell;
if mod(sp.spline_degree, 2) == 0
    sp.spline_degree = sp.spline_degree - 1;
end
sp.spline_degree = max(sp.spline_degree, 5);
sp.treeMode = "all";
sp.pointSet = pointSet;
end

function assembler = makeAssembler(assemblerSpec, stencilSpec)
% The solver stays "puzzle piece" simple by accepting either standard or
% overlapped assemblers, and either RBF or WLS stencils.
stencilFactory = resolveStencilFactory(stencilSpec);
assemblerName = lower(string(assemblerSpec));
switch assemblerName
    case {"fd", "fddiffop", "standard"}
        assembler = kp.rbffd.FDDiffOp(stencilFactory);
    case {"fdo", "fdodiffop", "overlapped", "overlap"}
        assembler = kp.rbffd.FDODiffOp(stencilFactory);
    otherwise
        error('kp:solvers:BadAssembler', 'Unknown assembler "%s".', assemblerName);
end
end

function stencilFactory = resolveStencilFactory(stencilSpec)
if isa(stencilSpec, 'function_handle')
    stencilFactory = stencilSpec;
    return;
end

stencilName = lower(string(stencilSpec));
switch stencilName
    case {"rbf", "rbffd", "rbf-fd"}
        stencilFactory = @() kp.rbffd.RBFStencil();
    case {"wls", "weightedleastsquares", "weighted_least_squares"}
        stencilFactory = @() kp.rbffd.WeightedLeastSquaresStencil();
    otherwise
        error('kp:solvers:BadStencil', 'Unknown stencil backend "%s".', stencilName);
end
end

function values = evaluateNodeCallback(func, X, label)
if isa(func, 'function_handle')
    values = func(X);
else
    values = func;
end
values = values(:);
if isscalar(values)
    values = repmat(values, size(X, 1), 1);
end
if numel(values) ~= size(X, 1)
    error('kp:solvers:BadCallbackSize', '%s values must match the node count.', label);
end
end

function values = evaluateBoundaryValues(func, neuCoeff, dirCoeff, nr, Xb)
% Accept either the full KernelPack-style boundary callback signature or a
% simpler value-at-boundary callback.
if isa(func, 'function_handle')
    try
        values = func(neuCoeff, dirCoeff, nr, Xb);
    catch
        values = func(Xb);
    end
else
    values = func;
end
values = values(:);
if numel(values) ~= size(Xb, 1)
    error('kp:solvers:BadBoundaryValueSize', 'Boundary values must match the boundary row count.');
end
end

function A = buildSystemMatrix(L, BC, nCols, pureNeumann)
% Rows are always physical PDE rows followed by boundary rows. Pure
% Neumann solves append one scalar compatibility constraint.
A = [-L; BC];
if pureNeumann
    A = [A, ones(size(A, 1), 1); ones(1, nCols), 0];
end
end

function b = buildSystemRHS(rhsTarget, rhsBoundary, pureNeumann)
b = [rhsTarget(:); rhsBoundary(:)];
if pureNeumann
    b = [b; 0];
end
end

function guess = buildInitialGuess(initialGuess, nTargets, nCols, rhsBoundary, pureNeumann)
if isempty(initialGuess)
    guess = [];
    return;
end

if pureNeumann
    if numel(initialGuess) == nTargets
        guess = [initialGuess(:); rhsBoundary(:); 0];
    elseif numel(initialGuess) == nCols
        guess = [initialGuess(:); 0];
    elseif numel(initialGuess) == nCols + 1
        guess = initialGuess(:);
    else
        error('kp:solvers:BadInitialGuess', 'Pure-Neumann Poisson guess must have length %d, %d, or %d.', nTargets, nCols, nCols + 1);
    end
else
    if numel(initialGuess) == nTargets
        guess = [initialGuess(:); rhsBoundary(:)];
    elseif numel(initialGuess) == nCols
        guess = initialGuess(:);
    else
        error('kp:solvers:BadInitialGuess', 'Poisson guess must have length %d or %d.', nTargets, nCols);
    end
end
end

function sol = gmresWithFallback(A, b, guess)
[sol, flag] = gmres(sparse(A), b, [], 1e-10, 200, [], [], guess);
if flag ~= 0 || any(~isfinite(sol))
    sol = A \ b;
end
end
