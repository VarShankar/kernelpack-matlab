classdef VariablePoissonSolver < handle
    %VARIABLEPOISSONSOLVER Variable-coefficient ghost-node Poisson solver.

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
        Xf double = zeros(0, 0)
        nr double = zeros(0, 0)
        N (1,1) double = 0
        Nf (1,1) double = 0
        Lap double = zeros(0, 0)
        Grad cell = {}
        PDE double = zeros(0, 0)
        BC double = zeros(0, 0)
        LapStencilProperties kp.rbffd.StencilProperties = kp.rbffd.StencilProperties()
        BCStencilProperties kp.rbffd.StencilProperties = kp.rbffd.StencilProperties()
        LapOpProperties kp.rbffd.OpProperties = kp.rbffd.OpProperties('decompose', false, 'storeWeights', true, 'recordStencils', false)
        BCOpProperties kp.rbffd.OpProperties = kp.rbffd.OpProperties('decompose', false, 'storeWeights', true, 'recordStencils', false)
        last_solve_used_nullspace_ (1,1) logical = false
    end

    methods
        function obj = VariablePoissonSolver(varargin)
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

            % Cache the geometry and the coefficient-independent RBF-FD
            % operators. The PDE matrix itself is rebuilt at solve time.
            obj.Domain.buildStructs();
            obj.X = obj.Domain.getIntBdryNodes();
            obj.Xb = obj.Domain.getBdryNodes();
            obj.Xf = obj.Domain.getAllNodes();
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

            obj.Grad = cell(1, obj.Domain.getDim());
            for d = 1:obj.Domain.getDim()
                gradAssembler = makeAssembler(obj.LapAssembler, obj.LapStencil);
                gradProps = obj.LapOpProperties;
                gradProps.selectdim = d - 1;
                gradAssembler.AssembleOp(obj.Domain, 'grad', obj.LapStencilProperties, gradProps);
                obj.Grad{d} = gradAssembler.getOp();
            end

            obj.PDE = zeros(0, obj.Nf);
            obj.BC = zeros(0, obj.Nf);
            obj.last_solve_used_nullspace_ = false;
        end

        function result = solve(obj, forcing, coeff, NeuCoeffFunc, DirCoeffFunc, bc, varargin)
            parser = inputParser();
            parser.addOptional('InitialGuess', zeros(0, 1), @(x) isnumeric(x));
            parser.parse(varargin{:});
            initialGuess = parser.Results.InitialGuess(:);

            coeffAll = evaluateNodeCallback(coeff, obj.Xf, 'coefficient');
            if any(coeffAll <= 0)
                error('kp:solvers:BadCoefficient', 'VariablePoissonSolver expects a positive scalar coefficient field.');
            end

            obj.PDE = buildVariablePDEOperator(obj.Lap, obj.Grad, coeffAll, obj.N);

            neuCoeff = evaluateNodeCallback(NeuCoeffFunc, obj.Xb, 'boundary coefficient');
            dirCoeff = evaluateNodeCallback(DirCoeffFunc, obj.Xb, 'boundary coefficient');

            bcAssembler = makeAssembler(obj.BCAssembler, obj.BCStencil);
            bcAssembler.AssembleOp(obj.Domain, 'bc', obj.BCStencilProperties, obj.BCOpProperties, ...
                'NeuCoeff', neuCoeff, ...
                'DirCoeff', dirCoeff);
            obj.BC = bcAssembler.getOp();

            rhsTarget = evaluateNodeCallback(forcing, obj.X, 'forcing');
            rhsBoundary = evaluateBoundaryValues(bc, neuCoeff, dirCoeff, obj.nr, obj.Xb);
            pureNeumann = max(abs(dirCoeff)) <= 1e-13;
            obj.last_solve_used_nullspace_ = pureNeumann;

            A = buildSystemMatrix(obj.PDE, obj.BC, obj.Nf, pureNeumann);
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
            result.Coefficient = coeffAll;
            result.L = obj.Lap;
            result.Grad = {obj.Grad{:}};
            result.PDE = obj.PDE;
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

        function out = getGradientOps(obj)
            out = cellfun(@sparse, obj.Grad, 'UniformOutput', false);
        end

        function out = getLastPdeOperator(obj)
            out = sparse(obj.PDE);
        end

        function out = getBCOp(obj)
            out = sparse(obj.BC);
        end

        function out = lastSolveUsedNullspace(obj)
            out = obj.last_solve_used_nullspace_;
        end
    end
end

function pde = buildVariablePDEOperator(Lap, GradOps, coeffAll, nRows)
% Build the physical PDE rows for -div(a grad u) = f in expanded form.
coeffLocal = coeffAll(1:nRows);
pde = -spdiags(coeffLocal, 0, nRows, nRows) * Lap;
for d = 1:numel(GradOps)
    gradCoeff = GradOps{d} * coeffAll;
    pde = pde - spdiags(gradCoeff, 0, nRows, nRows) * GradOps{d};
end
end

function sp = buildStencilProperties(domain, xi, theta, pointSet)
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

function A = buildSystemMatrix(PDE, BC, nCols, pureNeumann)
A = [PDE; BC];
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
        error('kp:solvers:BadInitialGuess', 'Pure-Neumann variable Poisson guess must have length %d, %d, or %d.', nTargets, nCols, nCols + 1);
    end
else
    if numel(initialGuess) == nTargets
        guess = [initialGuess(:); rhsBoundary(:)];
    elseif numel(initialGuess) == nCols
        guess = initialGuess(:);
    else
        error('kp:solvers:BadInitialGuess', 'Variable Poisson guess must have length %d or %d.', nTargets, nCols);
    end
end
end

function sol = gmresWithFallback(A, b, guess)
[sol, flag] = gmres(sparse(A), b, [], 1e-10, 200, [], [], guess);
if flag ~= 0 || any(~isfinite(sol))
    sol = A \ b;
end
end
