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
        RbfCoeffs double = zeros(0, 0)
        PolyCoeffs (:,1) double = zeros(0, 1)
        ResidualWeights (:,1) double = zeros(0, 1)
        ResidualConstant (1,1) double = 0
        MonomialExponents double = zeros(0, 0)
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
            obj.ell = 2;
            obj.MonomialExponents = kp.poly.total_degree_indices(obj.dim, obj.ell);
            obj.npoly = size(obj.MonomialExponents, 1);
            obj.xd = x;
            obj.nrd = kp.geometry.normalizeRows(nr);
            obj.ls_xd = zeros(obj.n, 1);
            obj.Centers = x;

            % KernelPack builds a curl-free interpolant to the surface
            % normals, then extracts the discretely zero-mean potential.
            cfA = kp.geometry.RBFLevelSet.buildCurlFreeGram(x, x, obj.m_spline_degree, true);
            cfP = kp.geometry.RBFLevelSet.buildCurlFreePolynomialMatrix(x, obj.MonomialExponents);
            reg = 1e-12 * max(1.0, max(abs(cfA), [], 'all'));
            A = [cfA + reg * eye(size(cfA, 1)), cfP; cfP.', zeros(size(cfP, 2), size(cfP, 2))];
            rhs = [obj.nrd(:); zeros(size(cfP, 2), 1)];
            coeffs = A \ rhs;

            obj.RbfCoeffs = reshape(coeffs(1:obj.dim * obj.n), obj.n, obj.dim);
            obj.PolyCoeffs = coeffs(obj.dim * obj.n + 1:end);
            [obj.ResidualWeights, obj.ResidualConstant, pot] = ...
                kp.geometry.RBFLevelSet.computeZeroMeanResidual( ...
                    x, obj.RbfCoeffs, obj.PolyCoeffs, obj.MonomialExponents, obj.m_spline_degree);
            obj.ls_xd = pot;
            obj.mean_potential = mean(pot, 'all');

            % Keep the sign convention aligned with KernelPack's geometry
            % usage by making the inward normal offsets positive.
            sep = obj.estimateOffsetDistance(x);
            insideProbe = x - sep * obj.nrd;
            if mean(obj.Evaluate(insideProbe), 'all') < 0
                obj.RbfCoeffs = -obj.RbfCoeffs;
                obj.PolyCoeffs = -obj.PolyCoeffs;
                obj.ResidualWeights = -obj.ResidualWeights;
                obj.ResidualConstant = -obj.ResidualConstant;
                obj.ls_xd = -obj.ls_xd;
                obj.mean_potential = -obj.mean_potential;
            end
        end

        function val = Evaluate(obj, xe)
            val = kp.geometry.RBFLevelSet.evaluateModel(obj.getEvaluationModel(), xe);
        end

        function grad = EvaluateGradient(obj, xe)
            grad = kp.geometry.RBFLevelSet.evaluateModelGradient(obj.getEvaluationModel(), xe);
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

        function result = FindSegmentRootNewton(obj, insidePoint, outsidePoint, initialParameter, options)
            if nargin < 4 || isempty(initialParameter)
                initialParameter = 0.5;
            end
            if nargin < 5 || isempty(options)
                options = struct();
            end
            options = obj.withNewtonDefaults(options);

            dir = outsidePoint - insidePoint;
            t = min(max(initialParameter, 0.0), 1.0);
            x = insidePoint + t * dir;
            phi = obj.Evaluate(x);
            converged = false;
            stalled = false;
            iterations = 0;

            for it = 1:options.maxIterations
                iterations = it;
                if abs(phi) <= options.valueTolerance
                    converged = true;
                    break;
                end
                grad = obj.EvaluateGradient(x);
                dphi = grad * dir.';
                if abs(dphi) <= options.gradientTolerance
                    stalled = true;
                    break;
                end
                dt = -phi / dphi;
                if isfinite(options.maxStepNorm) && abs(dt) > options.maxStepNorm
                    dt = sign(dt) * options.maxStepNorm;
                end
                tNew = min(max(t + dt, 0.0), 1.0);
                x = insidePoint + tNew * dir;
                phi = obj.Evaluate(x);
                if abs(tNew - t) <= options.stepTolerance
                    t = tNew;
                    converged = abs(phi) <= options.valueTolerance;
                    break;
                end
                t = tNew;
            end

            result = struct( ...
                'point', x, ...
                'parameter', t, ...
                'value', phi, ...
                'iterations', iterations, ...
                'converged', converged, ...
                'stalled', stalled);
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
                'RbfCoeffs', obj.RbfCoeffs, ...
                'PolyCoeffs', obj.PolyCoeffs, ...
                'ResidualWeights', obj.ResidualWeights, ...
                'ResidualConstant', obj.ResidualConstant, ...
                'MonomialExponents', obj.MonomialExponents, ...
                'dim', obj.dim, ...
                'm_spline_degree', obj.m_spline_degree, ...
                'mean_potential', obj.mean_potential);
        end
    end

    methods (Access = private)
        function sep = estimateOffsetDistance(~, x)
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
            [val, ~] = kp.geometry.RBFLevelSet.evaluatePotentialPieces(model, xe);
        end

        function grad = evaluateModelGradient(model, xe)
            % The zero-mean potential gradient is the negative curl-free
            % field plus the scalar residual correction.
            field = kp.geometry.RBFLevelSet.evaluateCurlFreeField(model, xe);
            grad = -field;

            diff = kp.geometry.RBFLevelSet.differenceTensor(xe, model.Centers);
            r = sqrt(sum(diff.^2, 3));
            resFactor = kp.geometry.RBFLevelSet.radialDerivativeFactor(r, 1);
            for d = 1:model.dim
                grad(:, d) = grad(:, d) + (resFactor .* diff(:, :, d)) * model.ResidualWeights;
            end
        end

        function [field, flat] = evaluateCurlFreeField(model, xe)
            gram = kp.geometry.RBFLevelSet.buildCurlFreeGram(xe, model.Centers, model.m_spline_degree, false);
            poly = kp.geometry.RBFLevelSet.buildCurlFreePolynomialMatrix(xe, model.MonomialExponents);
            flat = gram * model.RbfCoeffs(:) + poly * model.PolyCoeffs;
            field = reshape(flat, size(xe, 1), model.dim);
        end

        function [value, rawPotential] = evaluatePotentialPieces(model, xe)
            diff = kp.geometry.RBFLevelSet.differenceTensor(xe, model.Centers);
            r = sqrt(sum(diff.^2, 3));
            derivFactor = kp.geometry.RBFLevelSet.radialDerivativeFactor(r, model.m_spline_degree);
            rawPotential = zeros(size(xe, 1), 1);
            for d = 1:model.dim
                rawPotential = rawPotential - (derivFactor .* diff(:, :, d)) * model.RbfCoeffs(:, d);
            end

            monomials = kp.geometry.RBFLevelSet.evaluateMonomials(xe, model.MonomialExponents);
            if ~isempty(monomials)
                rawPotential = rawPotential + monomials(:, 2:end) * model.PolyCoeffs;
            end

            residual = r * model.ResidualWeights + model.ResidualConstant;
            value = -(rawPotential - residual);
        end

        function gram = buildCurlFreeGram(X, Y, splineDegree, xEqualsY)
            nX = size(X, 1);
            nY = size(Y, 1);
            dim = size(X, 2);
            diff = kp.geometry.RBFLevelSet.differenceTensor(X, Y);
            r2 = sum(diff.^2, 3);
            r = sqrt(r2);
            gram = zeros(dim * nX, dim * nY);
            p = splineDegree;

            diagFactor = p * r .^ max(p - 2, 0);
            offFactor = p * max(p - 2, 0) * r .^ max(p - 4, 0);

            if p < 2
                diagFactor = zeros(size(r));
            end
            if p < 4
                offFactor = zeros(size(r));
            end

            zeroMask = r == 0;
            diagFactor(zeroMask) = 0;
            offFactor(zeroMask) = 0;

            for rowDim = 1:dim
                rowSpan = (rowDim - 1) * nX + (1:nX);
                for colDim = 1:dim
                    colSpan = (colDim - 1) * nY + (1:nY);
                    block = offFactor .* diff(:, :, rowDim) .* diff(:, :, colDim);
                    if rowDim == colDim
                        block = block + diagFactor;
                    end
                    gram(rowSpan, colSpan) = -block;
                end
            end

            if xEqualsY
                gram(1:size(gram, 1)+1:end) = 0;
            end
        end

        function poly = buildCurlFreePolynomialMatrix(X, exponents)
            dim = size(X, 2);
            derivs = kp.geometry.RBFLevelSet.evaluateMonomialDerivativesAll(X, exponents);
            nTerms = size(exponents, 1) - 1;
            poly = zeros(dim * size(X, 1), nTerms);
            for d = 1:dim
                deriv = derivs(:, 2:end, d);
                rows = (d - 1) * size(X, 1) + (1:size(X, 1));
                poly(rows, :) = deriv;
            end
        end

        function [resWeights, resConstant, pot] = computeZeroMeanResidual(X, rbfCoeffs, polyCoeffs, exponents, splineDegree)
            model = struct( ...
                'Centers', X, ...
                'RbfCoeffs', rbfCoeffs, ...
                'PolyCoeffs', polyCoeffs, ...
                'ResidualWeights', zeros(size(X, 1), 1), ...
                'ResidualConstant', 0, ...
                'MonomialExponents', exponents, ...
                'dim', size(X, 2), ...
                'm_spline_degree', splineDegree, ...
                'mean_potential', 0);
            [~, pot] = kp.geometry.RBFLevelSet.evaluatePotentialPieces(model, X);
            r = kp.geometry.distanceMatrix(X, X);
            A = [r, ones(size(X, 1), 1); ones(1, size(X, 1)), 0];
            rhs = [pot; 0];
            coeffs = A \ rhs;
            resWeights = coeffs(1:size(X, 1));
            resConstant = coeffs(end);
        end

        function tensor = differenceTensor(X, Y)
            tensor = reshape(X, size(X, 1), 1, size(X, 2)) - reshape(Y, 1, size(Y, 1), size(Y, 2));
        end

        function monomials = evaluateMonomials(X, exponents)
            monomials = ones(size(X, 1), size(exponents, 1));
            for d = 1:size(X, 2)
                if any(exponents(:, d) ~= 0)
                    monomials = monomials .* bsxfun(@power, X(:, d), exponents(:, d).');
                end
            end
        end

        function deriv = evaluateMonomialDerivatives(X, exponents, dimIndex)
            deriv = kp.geometry.RBFLevelSet.evaluateMonomialDerivativesAll(X, exponents);
            deriv = deriv(:, :, dimIndex);
        end

        function deriv = evaluateMonomialDerivativesAll(X, exponents)
            nX = size(X, 1);
            nTerms = size(exponents, 1);
            dim = size(X, 2);
            deriv = zeros(nX, nTerms, dim);
            for d = 1:dim
                alpha = exponents;
                coeff = alpha(:, d).';
                active = coeff > 0;
                if ~any(active)
                    continue;
                end
                alpha(active, d) = alpha(active, d) - 1;
                term = ones(nX, nTerms);
                for j = 1:dim
                    if any(alpha(:, j) ~= 0)
                        term = term .* bsxfun(@power, X(:, j), alpha(:, j).');
                    end
                end
                term(:, ~active) = 0;
                deriv(:, :, d) = term .* coeff;
            end
        end

        function factor = radialDerivativeFactor(r, degree)
            factor = zeros(size(r));
            if degree == 0
                return;
            end
            if degree == 1
                mask = r > 0;
                factor(mask) = 1 ./ r(mask);
                return;
            end
            factor = degree * r .^ (degree - 2);
            factor(r == 0) = 0;
        end
    end
end
