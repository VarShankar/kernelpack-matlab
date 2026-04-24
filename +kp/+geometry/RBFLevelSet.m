classdef RBFLevelSet < handle
    %RBFLEVELSET KernelPack-like implicit boundary representation.

    properties (SetAccess = private)
        n (1,1) double = 0
        ell (1,1) double = 1
        dim (1,1) double = 0
        m_spline_degree (1,1) double = 3
        npoly (1,1) double = 0
        xd double = zeros(0, 0)
        nrd double = zeros(0, 0)
        ls_xd double = zeros(0, 1)
        mean_potential (1,1) double = 0
    end

    properties (Access = private)
        Centers double = zeros(0, 0)
        Values (:,1) double = zeros(0, 1)
        Weights (:,1) double = zeros(0, 1)
        PolyCoeffs (:,1) double = zeros(0, 1)
    end

    methods
        function BuildLevelSetFromCFI(obj, varargin)
            if nargin == 4
                splineDegree = varargin{1};
                x = varargin{2};
                nr = varargin{3};
            elseif nargin == 3
                splineDegree = 3;
                x = varargin{1};
                nr = varargin{2};
            else
                error('RBFLevelSet:BadInput', ...
                    'Use BuildLevelSetFromCFI(x, nr) or BuildLevelSetFromCFI(splineDegree, x, nr).');
            end

            obj.dim = size(x, 2);
            obj.n = size(x, 1);
            obj.ell = 1;
            obj.m_spline_degree = splineDegree;
            obj.npoly = obj.dim + 1;
            obj.xd = x;
            obj.nrd = kp.geometry.normalizeRows(nr);
            obj.ls_xd = zeros(obj.n, 1);

            % Build a signed offset cloud around the boundary and fit a
            % scalar RBF interpolant through zero, inside, and outside
            % level-set values.
            sep = obj.estimateOffsetDistance(x);
            insidePts = x - sep * obj.nrd;
            outsidePts = x + sep * obj.nrd;

            obj.Centers = [x; insidePts; outsidePts];
            obj.Values = [zeros(obj.n, 1); -sep * ones(obj.n, 1); sep * ones(obj.n, 1)];

            nCenters = size(obj.Centers, 1);
            P = [ones(nCenters, 1), obj.Centers];
            R = kp.geometry.distanceMatrix(obj.Centers, obj.Centers);
            K = kp.geometry.phsKernel(R, obj.m_spline_degree);
            reg = 1e-12 * max(1.0, max(abs(K), [], 'all'));
            A = [K + reg * eye(nCenters), P; P.', zeros(obj.dim + 1, obj.dim + 1)];
            rhs = [obj.Values; zeros(obj.dim + 1, 1)];
            coeffs = A \ rhs;

            obj.Weights = coeffs(1:nCenters);
            obj.PolyCoeffs = coeffs(nCenters + 1:end);
            obj.mean_potential = mean(obj.Evaluate(x), 'all');
        end

        function val = Evaluate(obj, xe)
            val = kp.geometry.RBFLevelSet.evaluateModel(obj.getEvaluationModel(), xe);
        end

        function grad = EvaluateGradient(obj, xe)
            % Differentiate the fitted scalar field analytically by
            % differentiating the radial kernel.
            nPts = size(xe, 1);
            grad = zeros(nPts, obj.dim);
            R = kp.geometry.distanceMatrix(xe, obj.Centers);
            dphi = obj.m_spline_degree * R .^ max(obj.m_spline_degree - 2, 0);
            invR = zeros(size(R));
            mask = R > 0;
            invR(mask) = 1 ./ R(mask);
            for d = 1:obj.dim
                delta = xe(:, d) - obj.Centers(:, d).';
                grad(:, d) = (dphi .* delta .* invR) * obj.Weights + obj.PolyCoeffs(d + 1);
            end
        end

        function result = ProjectToSurfaceNewton(obj, initialPoints, options)
            if nargin < 3 || isempty(options)
                options = struct();
            end
            options = obj.withNewtonDefaults(options);

            x = initialPoints;
            result = struct( ...
                'points', x, ...
                'levelSetValues', zeros(size(x, 1), 1), ...
                'iterations', zeros(size(x, 1), 1), ...
                'converged', zeros(size(x, 1), 1), ...
                'stalled', zeros(size(x, 1), 1));

            % Use a simple Newton projection onto the implicit surface,
            % tracking both stalled points and converged points.
            for it = 1:options.maxIterations
                active = (result.converged == 0) & (result.stalled == 0);
                if ~any(active)
                    break;
                end

                xa = x(active, :);
                phi = obj.Evaluate(xa);
                grad = obj.EvaluateGradient(xa);
                g2 = sum(grad.^2, 2);

                localIdx = find(active);
                for j = 1:numel(localIdx)
                    idx = localIdx(j);
                    result.levelSetValues(idx) = phi(j);
                    result.iterations(idx) = it;
                    if abs(phi(j)) <= options.valueTolerance
                        result.converged(idx) = 1;
                        continue;
                    end
                    if g2(j) <= options.gradientTolerance^2
                        result.stalled(idx) = 1;
                        continue;
                    end
                    step = -(phi(j) / g2(j)) * grad(j, :);
                    stepNorm = norm(step, 2);
                    if isfinite(options.maxStepNorm) && stepNorm > options.maxStepNorm
                        step = step * (options.maxStepNorm / stepNorm);
                    end
                    x(idx, :) = x(idx, :) + step;
                    if norm(step, 2) <= options.stepTolerance
                        result.converged(idx) = 1;
                    end
                end
            end

            result.points = x;
            result.levelSetValues = obj.Evaluate(x);
            result.converged(abs(result.levelSetValues) <= options.valueTolerance) = 1;
        end

        function flags = IsPointInSurface(obj, xe, tol)
            if nargin < 3
                tol = 1e-3;
            end
            flags = uint32(obj.Evaluate(xe) >= 0.5 * tol);
        end

        function flags = IsPointOutsideSurface(obj, xe, tol)
            if nargin < 3
                tol = 1e-3;
            end
            flags = uint32(obj.Evaluate(xe) <= -0.5 * tol);
        end

        function model = getEvaluationModel(obj)
            model = struct( ...
                'Centers', obj.Centers, ...
                'Weights', obj.Weights, ...
                'PolyCoeffs', obj.PolyCoeffs, ...
                'm_spline_degree', obj.m_spline_degree, ...
                'mean_potential', obj.mean_potential);
        end
    end

    methods (Access = private)
        function sep = estimateOffsetDistance(~, x)
            % Pick an offset distance from the nearest-neighbor spacing so
            % the inside/outside cloud scales with the data cloud.
            if size(x, 1) < 2
                sep = 1e-2;
                return;
            end
            D = kp.geometry.distanceMatrix(x, x);
            D(D == 0) = inf;
            sep = max(0.5 * min(min(D, [], 2)), 1e-3);
        end

        function options = withNewtonDefaults(~, options)
            defaults = struct( ...
                'valueTolerance', 1e-12, ...
                'stepTolerance', 1e-12, ...
                'gradientTolerance', 1e-14, ...
                'maxStepNorm', inf, ...
                'maxIterations', 20);
            fields = fieldnames(defaults);
            for k = 1:numel(fields)
                name = fields{k};
                if ~isfield(options, name)
                    options.(name) = defaults.(name);
                end
            end
        end
    end

    methods (Static)
        function val = evaluateModel(model, xe)
            R = kp.geometry.distanceMatrix(xe, model.Centers);
            val = kp.geometry.phsKernel(R, model.m_spline_degree) * model.Weights + ...
                [ones(size(xe, 1), 1), xe] * model.PolyCoeffs - model.mean_potential;
        end
    end
end
