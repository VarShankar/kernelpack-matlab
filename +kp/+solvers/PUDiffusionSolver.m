classdef PUDiffusionSolver < handle
    %PUDIFFUSIONSOLVER Simple readable MATLAB PU diffusion stepper.

    properties (SetAccess = private)
        Domain kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
        xi (1,1) double = 0
        dt (1,1) double = NaN
        nu (1,1) double = NaN
        num_omp_threads (1,1) double = 1
        X double = zeros(0, 0)
        Xb double = zeros(0, 0)
        nr double = zeros(0, 0)
        N (1,1) double = 0
        Nf (1,1) double = 0
        patch_data struct = struct('centers', zeros(0, 0), 'radius', 0, 'spacing', 0, 'min_patch_nodes', 0, 'node_ids', {{}});
        Lap double = zeros(0, 0)
        BC double = zeros(0, 0)
        cnm2 double = zeros(0, 1)
        cnm1 double = zeros(0, 1)
        cn double = zeros(0, 1)
        completed_steps_ (1,1) double = 0
        fixed_bc_operator_ready_ (1,1) logical = false
        fixed_bc_coefficients_ready_ (1,1) logical = false
        cached_neu_coeff_ double = zeros(0, 1)
        cached_dir_coeff_ double = zeros(0, 1)
    end

    methods
        function init(obj, domain, xi, dlt, d_coeff, num_omp_threads)
            if nargin < 6 || isempty(num_omp_threads)
                num_omp_threads = 1;
            end
            obj.Domain = domain;
            obj.xi = xi;
            obj.dt = dlt;
            obj.nu = d_coeff;
            obj.num_omp_threads = num_omp_threads;

            obj.Domain.buildStructs();
            obj.X = obj.Domain.getIntBdryNodes();
            obj.Xb = obj.Domain.getBdryNodes();
            obj.nr = obj.Domain.getNrmls();
            obj.N = size(obj.X, 1);
            obj.Nf = obj.Domain.getNumTotalNodes();
            if obj.Nf ~= obj.N + obj.Domain.getNumBdryNodes()
                error('kp:solvers:BadGhostLayout', ...
                    'PUDiffusionSolver expects one ghost node per boundary node.');
            end

            % The PU diffusion solver builds one reusable patch covering
            % and then localizes both the Laplacian and BC operators on
            % that same covering.
            obj.patch_data = kp.solvers.puPatchGeometry(obj.Domain, obj.xi);
            obj.Lap = kp.solvers.puLocalizedOperator(obj.Domain, obj.patch_data, obj.xi, obj.X, 'lap');
            obj.BC = zeros(0, obj.Nf);
            obj.cnm2 = zeros(0, 1);
            obj.cnm1 = zeros(0, 1);
            obj.cn = zeros(0, 1);
            obj.completed_steps_ = 0;
            obj.fixed_bc_operator_ready_ = false;
            obj.fixed_bc_coefficients_ready_ = false;
            obj.cached_neu_coeff_ = zeros(0, 1);
            obj.cached_dir_coeff_ = zeros(0, 1);
        end

        function setStepSize(obj, dlt)
            obj.dt = dlt;
            obj.fixed_bc_operator_ready_ = false;
            obj.fixed_bc_coefficients_ready_ = false;
        end

        function setInitialState(obj, c0)
            obj.cnm2 = validatePhysicalState(c0, obj.N);
            obj.cnm1 = zeros(0, 1);
            obj.cn = zeros(0, 1);
            obj.completed_steps_ = 0;
            obj.fixed_bc_operator_ready_ = false;
            obj.fixed_bc_coefficients_ready_ = false;
        end

        function setStateHistory(obj, varargin)
            % Accept one, two, or three physical states so callers can
            % start directly at BDF1, BDF2, or BDF3.
            switch numel(varargin)
                case 1
                    obj.setInitialState(varargin{1});
                case 2
                    obj.cnm2 = validatePhysicalState(varargin{1}, obj.N);
                    obj.cnm1 = validatePhysicalState(varargin{2}, obj.N);
                    obj.cn = zeros(0, 1);
                    obj.completed_steps_ = 1;
                case 3
                    obj.cnm2 = validatePhysicalState(varargin{1}, obj.N);
                    obj.cnm1 = validatePhysicalState(varargin{2}, obj.N);
                    obj.cn = validatePhysicalState(varargin{3}, obj.N);
                    obj.completed_steps_ = 2;
                otherwise
                    error('kp:solvers:BadStateHistory', 'setStateHistory expects one, two, or three physical states.');
            end
            obj.fixed_bc_operator_ready_ = false;
            obj.fixed_bc_coefficients_ready_ = false;
        end

        function out = currentPhysicalState(obj)
            if obj.completed_steps_ <= 0
                out = obj.cnm2;
            elseif obj.completed_steps_ == 1
                out = obj.cnm1;
            else
                out = obj.cn;
            end
        end

        function out = returnsDistributedState(~)
            out = false;
        end

        function out = getOutputRange(obj)
            out = [1, size(obj.X, 1)];
        end

        function out = getOutputNodes(obj)
            out = obj.X;
        end

        function out = bdf1Step(obj, t, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            if isempty(obj.cnm2)
                error('kp:solvers:MissingInitialState', 'PUDiffusionSolver.bdf1Step requires setInitialState() first.');
            end
            previous = obj.currentPhysicalState();
            % BDF1 is just backward Euler written in the same localized
            % operator language as the higher-order steps.
            rhsPhysical = previous + obj.dt * evaluateForcingCallback(forcing, obj.nu, t, obj.X);
            out = obj.takeStep(rhsPhysical, t, NeuCoeffFunc, DirCoeffFunc, bc, -obj.nu * obj.dt);
        end

        function out = bdf2Step(obj, t, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            if obj.completed_steps_ < 1
                error('kp:solvers:MissingStateHistory', 'PUDiffusionSolver.bdf2Step requires one prior step in the state history.');
            end
            % BDF2 keeps the same solve shape, only changing the history
            % coefficients on the right-hand side.
            rhsPhysical = (4/3) * obj.cnm1 - (1/3) * obj.cnm2 + (2/3) * obj.dt * evaluateForcingCallback(forcing, obj.nu, t, obj.X);
            out = obj.takeStep(rhsPhysical, t, NeuCoeffFunc, DirCoeffFunc, bc, -(2/3) * obj.nu * obj.dt);
        end

        function out = bdf3Step(obj, t, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            if obj.completed_steps_ < 2
                error('kp:solvers:MissingStateHistory', 'PUDiffusionSolver.bdf3Step requires two prior steps in the state history.');
            end
            % BDF3 is the highest-order fixed-step option currently
            % mirrored for the PU diffusion path.
            rhsPhysical = (18/11) * obj.cn - (9/11) * obj.cnm1 + (2/11) * obj.cnm2 + (6/11) * obj.dt * evaluateForcingCallback(forcing, obj.nu, t, obj.X);
            out = obj.takeStep(rhsPhysical, t, NeuCoeffFunc, DirCoeffFunc, bc, -(6/11) * obj.nu * obj.dt);
        end
    end

    methods (Access = private)
        function out = takeStep(obj, rhsPhysical, t, NeuCoeffFunc, DirCoeffFunc, bc, lapScale)
            % Every time step solves one implicit elliptic system with the
            % cached PU Laplacian and the current PU boundary operator.
            [neuCoeff, dirCoeff] = obj.getBoundaryCoefficients(t, NeuCoeffFunc, DirCoeffFunc);
            obj.ensureBoundaryOperator(neuCoeff, dirCoeff);
            rhsBoundary = evaluateTransientBoundaryValues(bc, neuCoeff, dirCoeff, obj.nr, t, obj.Xb);
            A = buildImplicitSystem(obj.Lap, obj.BC, obj.N, lapScale);
            sol = A \ [rhsPhysical(:); rhsBoundary(:)];
            nextState = sol(1:obj.N);
            obj.pushCompletedStep(nextState);
            out = obj.currentPhysicalState();
        end

        function [neuCoeff, dirCoeff] = getBoundaryCoefficients(obj, t, NeuCoeffFunc, DirCoeffFunc)
            if isFixedBoundaryCallback(NeuCoeffFunc) && isFixedBoundaryCallback(DirCoeffFunc)
                if ~obj.fixed_bc_coefficients_ready_
                    obj.cached_neu_coeff_ = evaluateBoundaryCoefficient(NeuCoeffFunc, obj.Xb);
                    obj.cached_dir_coeff_ = evaluateBoundaryCoefficient(DirCoeffFunc, obj.Xb);
                    obj.fixed_bc_coefficients_ready_ = true;
                end
                neuCoeff = obj.cached_neu_coeff_;
                dirCoeff = obj.cached_dir_coeff_;
                return;
            end
            neuCoeff = evaluateBoundaryCoefficient(NeuCoeffFunc, obj.Xb, t);
            dirCoeff = evaluateBoundaryCoefficient(DirCoeffFunc, obj.Xb, t);
        end

        function ensureBoundaryOperator(obj, neuCoeff, dirCoeff)
            % Like the FD diffusion solver, keep the boundary operator
            % cached once the BC coefficients stop changing in time.
            if obj.fixed_bc_operator_ready_ && ~isempty(obj.BC)
                return;
            end
            obj.BC = kp.solvers.puLocalizedOperator(obj.Domain, obj.patch_data, obj.xi, obj.Xb, 'bc', ...
                'Normals', obj.nr, ...
                'NeuCoeff', neuCoeff, ...
                'DirCoeff', dirCoeff);
            if obj.fixed_bc_coefficients_ready_
                obj.fixed_bc_operator_ready_ = true;
            end
        end

        function pushCompletedStep(obj, nextState)
            if obj.completed_steps_ <= 0
                obj.cnm1 = nextState;
            elseif obj.completed_steps_ == 1
                obj.cn = nextState;
            else
                obj.cnm2 = obj.cnm1;
                obj.cnm1 = obj.cn;
                obj.cn = nextState;
            end
            obj.completed_steps_ = obj.completed_steps_ + 1;
        end
    end
end

function state = validatePhysicalState(state, n)
state = state(:);
if numel(state) ~= n
    error('kp:solvers:BadPhysicalStateSize', 'Expected a physical state of length %d.', n);
end
end

function tf = isFixedBoundaryCallback(f)
if ~isa(f, 'function_handle')
    tf = false;
    return;
end
try
    tf = nargin(f) == 1;
catch
    tf = false;
end
end

function values = evaluateBoundaryCoefficient(func, X, varargin)
if isa(func, 'function_handle')
    if isempty(varargin)
        values = func(X);
    else
        try
            values = func(varargin{1}, X);
        catch
            values = func(X);
        end
    end
else
    values = func;
end
values = values(:);
if isscalar(values)
    values = repmat(values, size(X, 1), 1);
end
if numel(values) ~= size(X, 1)
    error('kp:solvers:BadCoeffSize', 'Boundary coefficients must match the boundary row count.');
end
end

function values = evaluateForcingCallback(func, nu, t, X)
if isa(func, 'function_handle')
    try
        values = func(nu, t, X);
    catch
        try
            values = func(t, X);
        catch
            values = func(X);
        end
    end
else
    values = func;
end
values = values(:);
if numel(values) ~= size(X, 1)
    error('kp:solvers:BadForcingSize', 'Forcing values must match the physical row count.');
end
end

function values = evaluateTransientBoundaryValues(func, neuCoeff, dirCoeff, nr, t, Xb)
if isa(func, 'function_handle')
    try
        values = func(neuCoeff, dirCoeff, nr, t, Xb);
    catch
        try
            values = func(t, Xb);
        catch
            values = func(Xb);
        end
    end
else
    values = func;
end
values = values(:);
if numel(values) ~= size(Xb, 1)
    error('kp:solvers:BadBoundaryValueSize', 'Boundary values must match the boundary row count.');
end
end

function A = buildImplicitSystem(L, BC, nPhysical, lapScale)
A = lapScale * L;
A(1:nPhysical, 1:nPhysical) = A(1:nPhysical, 1:nPhysical) + speye(nPhysical);
A = [A; BC];
end
