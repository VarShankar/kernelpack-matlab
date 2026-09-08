classdef WeightedLeastSquaresStencil < handle
    %WEIGHTEDLEASTSQUARESSTENCIL Fixed weighted polynomial least-squares stencil.

    properties (SetAccess = private)
        x_stencil double = zeros(0, 0)
        Xc double = zeros(0, 0)
        Xm double = zeros(1, 0)
        width (1,1) double = 1
        s_dim (1,1) double = 0
        n (1,1) double = 0
        fit_ell (1,1) double = 0
        fit_npoly (1,1) double = 0
        node_weights double = zeros(0, 1)
        interp_metric double = zeros(0, 0)
        reconstructor double = zeros(0, 0)
        coeffs double = zeros(0, 0)
        coeffs_already_computed (1,1) logical = false
        basis kp.poly.PolynomialBasis = kp.poly.PolynomialBasis(zeros(1, 1))
        Wt double = zeros(0, 0)
        L double = zeros(0, 0)
        BC double = zeros(0, 0)
    end

    methods
        function InitializeGeometry(obj, X, sp)
            obj.coeffs_already_computed = false;
            obj.s_dim = size(X, 2);
            obj.n = size(X, 1);
            obj.fit_ell = sp.ell;
            obj.fit_npoly = sp.npoly;
            obj.x_stencil = X;
            obj.Xm = X(1, :);
            r2 = sum((X - obj.Xm).^2, 2);
            obj.width = max(sqrt(max(r2, [], 'all')), 1);
            obj.Xc = (X - obj.Xm) / obj.width;
            obj.basis = kp.poly.PolynomialBasis.fromTotalDegree(obj.s_dim, obj.fit_ell, ...
                'Family', 'legendre', 'Center', zeros(1, obj.s_dim), 'Scale', 1);

            % The WLS stencil works purely in the polynomial space, with a
            % diagonal radial weighting that biases the fit toward the
            % center of the stencil.
            P = obj.basis.evaluate(obj.Xc, zeros(1, obj.s_dim), true);
            obj.node_weights = exp(-4 * r2 / (obj.width^2));
            obj.node_weights = min(max(obj.node_weights, 1e-10), 1);
            sqrtw = sqrt(obj.node_weights);
            weighted_P = P .* sqrtw;
            weighted_identity = diag(sqrtw);
            gram = weighted_P.' * weighted_P;
            % Prefer a direct weighted least-squares solve when the Gram
            % matrix is healthy, but fall back to a pseudoinverse when the
            % fit is rank-deficient or poorly scaled.
            if rank(weighted_P) < obj.fit_npoly || rcond(gram) < 1e-12
                obj.reconstructor = pinv(weighted_P) * weighted_identity;
            else
                obj.reconstructor = weighted_P \ weighted_identity;
            end
            if any(~isfinite(obj.reconstructor), 'all')
                obj.reconstructor = pinv(weighted_P) * weighted_identity;
            end
            obj.reconstructor(~isfinite(obj.reconstructor)) = 0;
            obj.interp_metric = gram;
        end

        function W = ComputeWeights(obj, X, varargin)
            if isnumeric(varargin{1}) && size(varargin{1}, 2) == size(X, 2) && numel(varargin) >= 7
                nr = varargin{1};
                NeuCoeff = varargin{2};
                DirCoeff = varargin{3};
                sp = varargin{4};
                op = varargin{5};
                ApplyOp = varargin{6};
                rhs_indices = varargin{7};
                W = computeWeightsBoundary(obj, X, nr, NeuCoeff, DirCoeff, sp, op, ApplyOp, rhs_indices);
            else
                sp = varargin{1};
                op = varargin{2};
                ApplyOp = varargin{3};
                rhs_indices = varargin{4};
                W = computeWeightsInterior(obj, X, sp, op, ApplyOp, rhs_indices);
            end
        end

        function W = ComputeWeightsAtPoints(obj, X, Xe, sp, op, ApplyOp)
            obj.InitializeGeometry(X, sp);
            if isempty(Xe)
                W = zeros(obj.n, 0);
                return;
            end
            Xec = (Xe - obj.Xm) / obj.width;
            Bpoly = obj.applyOperator(ApplyOp, sp, op, [], Xe, X, Xec, obj.Xc);
            W = obj.reconstructor.' * Bpoly;
        end

        function values = EvalStencil(obj, sp, Xe, f, cache_flag)
            if nargin < 5
                cache_flag = true;
            end
            % Cache the polynomial coefficients when the same local field
            % is queried repeatedly.
            if ~obj.coeffs_already_computed
                obj.coeffs = obj.reconstructor * f;
                if cache_flag
                    obj.coeffs_already_computed = true;
                end
            end
            W = obj.EvalWeights(sp, Xe);
            values = W * f;
        end

        function weights = EvalWeights(obj, ~, Xe)
            weights = zeros(size(Xe, 1), obj.n);
            for row = 1:size(Xe, 1)
                r2 = sum((obj.x_stencil - Xe(row, :)).^2, 2);
                [minr, nearest] = min(r2);
                % Exact stencil-node evaluation should reproduce the data
                % value directly instead of going through the fit.
                if sqrt(minr) <= 1e-12 * max(obj.width, 1)
                    weights(row, nearest) = 1;
                    continue;
                end
                Xeval = (Xe(row, :) - obj.Xm) / obj.width;
                Pe = obj.basis.evaluate(Xeval, zeros(1, obj.s_dim), true);
                weights(row, :) = Pe * obj.reconstructor;
            end
        end

        function setCoeffsStateToFalse(obj), obj.coeffs_already_computed = false; end
        function setCoeffsStateToTrue(obj), obj.coeffs_already_computed = true; end

        function flush(obj, op_name)
            if string(op_name) == "all"
                obj.coeffs = zeros(0, 0);
                obj.reconstructor = zeros(0, 0);
                obj.interp_metric = zeros(0, 0);
                obj.node_weights = zeros(0, 1);
                obj.Wt = zeros(0, 0);
                obj.L = zeros(0, 0);
                obj.BC = zeros(0, 0);
                obj.x_stencil = zeros(0, 0);
                obj.Xc = zeros(0, 0);
                obj.Xm = zeros(1, 0);
                obj.coeffs_already_computed = false;
            end
        end

        function out = getInterpMat(obj), out = obj.interp_metric; end
        function out = getLaplacian(obj), out = obj.L; end
        function out = getBCOp(obj), out = obj.BC; end
        function out = getLagrange(obj), out = obj.Wt; end
        function out = getWidth(obj), out = obj.width; end
        function out = getCentroid(obj), out = obj.Xm; end
        function out = getStencilNodes(obj), out = obj.x_stencil; end
        function out = getScaledStencilNodes(obj), out = obj.Xc; end
        function out = getPolyDegree(obj), out = obj.fit_ell; end
        function out = getNPoly(obj), out = obj.fit_npoly; end
        function out = getReconstructor(obj), out = obj.reconstructor; end

        function B = LapOp(obj, ~, ~, ~, ~, ~, X_at_origin_subset, ~)
            % All WLS differential operators are derivatives of the local
            % polynomial basis, then projected back to nodal weights.
            total = zeros(size(X_at_origin_subset, 1), obj.fit_npoly);
            for d = 1:obj.s_dim
                dd = zeros(1, obj.s_dim);
                dd(d) = 2;
                total = total + obj.basis.evaluate(X_at_origin_subset, dd, true) / (obj.width^2);
            end
            B = total.';
        end

        function B = GradOp(obj, ~, op, ~, ~, ~, X_at_origin_subset, ~)
            B = obj.basis.evaluate(X_at_origin_subset, unitMultiIndex(obj.s_dim, op.selectdim + 1), true).' / obj.width;
        end

        function B = BCOp(obj, ~, ~, NeuCoeff, DirCoeff, ~, ~, ~, X_at_origin_subset, ~, nr_subset)
            total = zeros(obj.fit_npoly, size(X_at_origin_subset, 1));
            if NeuCoeff ~= 0
                for d = 1:obj.s_dim
                    grad = obj.basis.evaluate(X_at_origin_subset, unitMultiIndex(obj.s_dim, d), true);
                    total = total + NeuCoeff * (grad.' .* nr_subset(:, d).') / obj.width;
                end
            end
            if DirCoeff ~= 0
                total = total + DirCoeff * obj.basis.evaluate(X_at_origin_subset, zeros(1, obj.s_dim), true).';
            end
            if DirCoeff == 0 && NeuCoeff == 0
                error('kp:rbffd:ZeroBC', 'Both boundary coefficients cannot be zero.');
            end
            B = total;
        end

        function B = InterpOp(obj, ~, ~, ~, ~, ~, X_at_origin_subset, ~)
            B = obj.basis.evaluate(X_at_origin_subset, zeros(1, obj.s_dim), true).';
        end
    end

    methods (Access = private)
        function W = computeWeightsInterior(obj, X, sp, op, ApplyOp, rhs_indices)
            obj.InitializeGeometry(X, sp);
            rhs_inds = rhs_indices(:);
            X_subset = X(rhs_inds, :);
            X_at_origin_subset = obj.Xc(rhs_inds, :);
            Bpoly = obj.applyOperator(ApplyOp, sp, op, [], X_subset, X, X_at_origin_subset, obj.Xc);
            W = obj.reconstructor.' * Bpoly;
            if isscalar(rhs_indices) && rhs_indices(1) == 1
                switch lower(string(ApplyOp))
                    case {"lap", "laplacian"}
                        obj.L = W;
                    case {"interp", "interpolation"}
                        obj.Wt = W;
                end
            end
        end

        function W = computeWeightsBoundary(obj, X, nr, NeuCoeff, DirCoeff, sp, op, ApplyOp, rhs_indices)
            obj.InitializeGeometry(X, sp);
            rhs_inds = rhs_indices(:);
            X_subset = X(rhs_inds, :);
            X_at_origin_subset = obj.Xc(rhs_inds, :);
            nr_subset = nr(rhs_inds, :);
            Bpoly = obj.applyBoundaryOperator(ApplyOp, sp, op, NeuCoeff, DirCoeff, [], X_subset, X, X_at_origin_subset, obj.Xc, nr_subset);
            W = obj.reconstructor.' * Bpoly;
            if isscalar(rhs_indices) && rhs_indices(1) == 1
                obj.BC = W;
            end
        end

        function B = applyOperator(obj, ApplyOp, sp, op, r_rhs, X_subset, X, X_at_origin_subset, X_at_origin)
            if isa(ApplyOp, 'function_handle')
                B = ApplyOp(obj, sp, op, r_rhs, X_subset, X, X_at_origin_subset, X_at_origin);
                return;
            end
            switch lower(string(ApplyOp))
                case {"lap", "laplacian"}
                    B = obj.LapOp(sp, op, r_rhs, X_subset, X, X_at_origin_subset, X_at_origin);
                case {"grad", "gradient"}
                    B = obj.GradOp(sp, op, r_rhs, X_subset, X, X_at_origin_subset, X_at_origin);
                case {"interp", "interpolation"}
                    B = obj.InterpOp(sp, op, r_rhs, X_subset, X, X_at_origin_subset, X_at_origin);
                otherwise
                    error('kp:rbffd:UnknownOperator', 'Unknown operator "%s".', string(ApplyOp));
            end
        end

        function B = applyBoundaryOperator(obj, ApplyOp, sp, op, NeuCoeff, DirCoeff, r_rhs, X_subset, X, X_at_origin_subset, X_at_origin, nr_subset)
            if isa(ApplyOp, 'function_handle')
                B = ApplyOp(obj, sp, op, NeuCoeff, DirCoeff, r_rhs, X_subset, X, X_at_origin_subset, X_at_origin, nr_subset);
                return;
            end
            switch lower(string(ApplyOp))
                case {"bc", "boundary"}
                    B = obj.BCOp(sp, op, NeuCoeff, DirCoeff, r_rhs, X_subset, X, X_at_origin_subset, X_at_origin, nr_subset);
                otherwise
                    error('kp:rbffd:UnknownBoundaryOperator', 'Unknown boundary operator "%s".', string(ApplyOp));
            end
        end
    end
end

function d = unitMultiIndex(dim, selectdim)
    d = zeros(1, dim);
    d(selectdim) = 1;
end
