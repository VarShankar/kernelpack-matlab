classdef RBFStencil < handle
    %RBFSTENCIL Local PHS-plus-polynomial stencil model.

    properties (SetAccess = private)
        A double = zeros(0, 0)
        solve_lhs double = zeros(0, 0)
        coeffs double = zeros(0, 0)
        x_stencil double = zeros(0, 0)
        Xc double = zeros(0, 0)
        Xm double = zeros(1, 0)
        width (1,1) double = 1
        s_dim (1,1) double = 0
        n (1,1) double = 0
        ell (1,1) double = 0
        npoly (1,1) double = 0
        basis kp.poly.PolynomialBasis = kp.poly.PolynomialBasis(zeros(1, 1))
        coeffs_already_computed (1,1) logical = false
        Wt double = zeros(0, 0)
        L double = zeros(0, 0)
        BC double = zeros(0, 0)
        Gx double = zeros(0, 0)
        Gy double = zeros(0, 0)
        Gz double = zeros(0, 0)
    end

    methods
        function InitializeGeometry(obj, X, sp)
            obj.coeffs_already_computed = false;
            obj.s_dim = size(X, 2);
            obj.n = size(X, 1);
            obj.x_stencil = X;
            obj.ell = sp.ell;
            obj.npoly = sp.npoly;

            % Shift and scale the stencil first so the polynomial block is
            % built on a numerically tame local coordinate system.
            r = kp.geometry.distanceMatrix(X, X);
            obj.width = max(max(r, [], 'all'), 1);
            obj.Xm = mean(X, 1);
            obj.Xc = (X - obj.Xm) / obj.width;
            obj.basis = kp.poly.PolynomialBasis.fromTotalDegree(obj.s_dim, obj.ell, ...
                'Family', 'legendre', 'Center', zeros(1, obj.s_dim), 'Scale', 1);

            P = obj.basis.evaluate(obj.Xc, zeros(1, obj.s_dim), true);
            % Assemble the standard augmented PHS+polynomial saddle-point
            % system used to recover interpolation and differentiation
            % weights.
            obj.A = zeros(obj.n + obj.npoly, obj.n + obj.npoly);
            obj.A(1:obj.n, 1:obj.n) = kp.rbffd.RBFStencil.phsRbf(r, sp.spline_degree);
            obj.A(1:obj.n, obj.n+1:end) = P;
            obj.A(obj.n+1:end, 1:obj.n) = P.';
            obj.solve_lhs = obj.A;
        end

        function W = ComputeWeights(obj, X, varargin)
            if nargin < 3
                error('kp:rbffd:BadCall', 'ComputeWeights requires stencil properties and operator information.');
            end

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
            r_rhs = kp.geometry.distanceMatrix(Xe, X);
            Xec = (Xe - obj.Xm) / obj.width;
            B = obj.applyOperator(ApplyOp, sp, op, r_rhs, Xe, X, Xec, obj.Xc);
            Wfull = kp.rbffd.RBFStencil.stableSolve(obj.solve_lhs, B);
            W = Wfull(1:obj.n, :);
        end

        function values = EvalStencil(obj, sp, Xe, f, cache_flag)
            if nargin < 5
                cache_flag = true;
            end
            re = kp.geometry.distanceMatrix(Xe, obj.x_stencil);
            Xec = (Xe - obj.Xm) / obj.width;
            Pe = obj.basis.evaluate(Xec, zeros(1, obj.s_dim), true);
            Ae = [kp.rbffd.RBFStencil.phsRbf(re, sp.spline_degree), Pe];

            f_padded = [f; zeros(obj.npoly, size(f, 2))];
            % Solving once and caching the interpolation coefficients is
            % useful when the same stencil is evaluated at many points.
            if ~obj.coeffs_already_computed
                obj.coeffs = kp.rbffd.RBFStencil.stableSolve(obj.solve_lhs, f_padded);
                if cache_flag
                    obj.coeffs_already_computed = true;
                end
            end
            values = Ae * obj.coeffs;
        end

        function weights = EvalWeights(obj, sp, Xe)
            if isempty(Xe)
                weights = zeros(0, obj.n);
                return;
            end
            re = kp.geometry.distanceMatrix(Xe, obj.x_stencil);
            Xec = (Xe - obj.Xm) / obj.width;
            Pe = obj.basis.evaluate(Xec, zeros(1, obj.s_dim), true);
            Rt = [kp.rbffd.RBFStencil.phsRbf(re, sp.spline_degree), Pe].';
            % Weight evaluation is just the transposed interpolation
            % problem with point-evaluation right-hand sides.
            lagrange = kp.rbffd.RBFStencil.stableSolve(obj.solve_lhs, Rt);
            weights = lagrange(1:obj.n, :).';
        end

        function setCoeffsStateToFalse(obj), obj.coeffs_already_computed = false; end
        function setCoeffsStateToTrue(obj), obj.coeffs_already_computed = true; end

        function flush(obj, op_name)
            switch string(op_name)
                case "lagrange"
                    obj.Wt = zeros(0, 0);
                case "interp"
                    obj.A = zeros(0, 0);
                case "lap"
                    obj.L = zeros(0, 0);
                case "BC"
                    obj.BC = zeros(0, 0);
                case "grad"
                    obj.Gx = zeros(0, 0);
                    obj.Gy = zeros(0, 0);
                    obj.Gz = zeros(0, 0);
                case "all"
                    obj.A = zeros(0, 0);
                    obj.solve_lhs = zeros(0, 0);
                    obj.coeffs = zeros(0, 0);
                    obj.Wt = zeros(0, 0);
                    obj.L = zeros(0, 0);
                    obj.BC = zeros(0, 0);
                    obj.Gx = zeros(0, 0);
                    obj.Gy = zeros(0, 0);
                    obj.Gz = zeros(0, 0);
                    obj.x_stencil = zeros(0, 0);
                    obj.Xc = zeros(0, 0);
                    obj.Xm = zeros(1, 0);
                    obj.coeffs_already_computed = false;
            end
        end

        function out = getInterpMat(obj), out = obj.A(1:obj.n, 1:obj.n); end
        function out = getLaplacian(obj), out = obj.L; end
        function out = getBCOp(obj), out = obj.BC; end
        function out = getLagrange(obj), out = obj.Wt; end
        function out = getPartialX(obj), out = obj.Gx; end
        function out = getPartialY(obj), out = obj.Gy; end
        function out = getPartialZ(obj), out = obj.Gz; end
        function out = getWidth(obj), out = obj.width; end
        function out = getCentroid(obj), out = obj.Xm; end
        function out = getStencilNodes(obj), out = obj.x_stencil; end
        function out = getScaledStencilNodes(obj), out = obj.Xc; end
        function out = getNPoly(obj), out = obj.npoly; end

        function B = LapOp(obj, sp, ~, r_rhs, ~, ~, X_at_origin_subset, ~)
            % Build the augmented right-hand side for Laplacian weights:
            % radial-kernel Laplacian plus polynomial second derivatives.
            Bpoly = zeros(obj.npoly, size(X_at_origin_subset, 1));
            for d = 1:obj.s_dim
                dd = zeros(1, obj.s_dim);
                dd(d) = 2;
                Bpoly = Bpoly + obj.basis.evaluate(X_at_origin_subset, dd, true).' / (obj.width^2);
            end
            B = [kp.rbffd.RBFStencil.phsLap(r_rhs, sp.spline_degree, obj.s_dim).'; Bpoly];
        end

        function B = GradOp(obj, sp, op, r_rhs, X_subset, X, X_at_origin_subset, ~)
            % Gradient weights use the selected coordinate derivative on
            % both the kernel and polynomial blocks.
            dim = op.selectdim + 1;
            diff = X_subset(:, dim) - X(:, dim).';
            Bpoly = obj.basis.evaluate(X_at_origin_subset, unitMultiIndex(obj.s_dim, dim), true).' / obj.width;
            B = [(diff .* kp.rbffd.RBFStencil.phsDrOverR(r_rhs, sp.spline_degree)).'; Bpoly];
        end

        function B = BCOp(obj, sp, op, NeuCoeff, DirCoeff, r_rhs, X_subset, X, X_at_origin_subset, ~, nr_subset)
            % Boundary rows are just the requested Neumann and Dirichlet
            % pieces superposed in one local augmented system.
            total = zeros(obj.n + obj.npoly, size(X_at_origin_subset, 1));
            if NeuCoeff ~= 0
                for d = 1:obj.s_dim
                    diff = X_subset(:, d) - X(:, d).';
                    gradRbf = (diff .* kp.rbffd.RBFStencil.phsDrOverR(r_rhs, sp.spline_degree)).';
                    gradPoly = obj.basis.evaluate(X_at_origin_subset, unitMultiIndex(obj.s_dim, d), true).' / obj.width;
                    total = total + NeuCoeff * [gradRbf; gradPoly] .* nr_subset(:, d).';
                end
            end
            if DirCoeff ~= 0
                Binterp = obj.InterpOp(sp, op, r_rhs, X_subset, X, X_at_origin_subset, []);
                total = total + DirCoeff * Binterp;
            end
            if DirCoeff == 0 && NeuCoeff == 0
                error('kp:rbffd:ZeroBC', 'Both boundary coefficients cannot be zero.');
            end
            B = total;
        end

        function B = InterpOp(obj, sp, ~, r_rhs, ~, ~, X_at_origin_subset, ~)
            Bpoly = obj.basis.evaluate(X_at_origin_subset, zeros(1, obj.s_dim), true).';
            B = [kp.rbffd.RBFStencil.phsRbf(r_rhs, sp.spline_degree).'; Bpoly];
        end
    end

    methods (Access = private)
        function W = computeWeightsInterior(obj, X, sp, op, ApplyOp, rhs_indices)
            obj.InitializeGeometry(X, sp);
            rhs_inds = rhs_indices(:);
            X_subset = X(rhs_inds, :);
            X_at_origin_subset = obj.Xc(rhs_inds, :);
            r = kp.geometry.distanceMatrix(X, X);
            r_rhs = r(rhs_inds, :);
            B = obj.applyOperator(ApplyOp, sp, op, r_rhs, X_subset, X, X_at_origin_subset, obj.Xc);
            % Some callers want the assembled right-hand side directly; the
            % normal path solves and then keeps only the nodal weights.
            if op.nosolve
                W = B;
            else
                Wfull = kp.rbffd.RBFStencil.stableSolve(obj.solve_lhs, B);
                W = Wfull(1:obj.n, :);
            end
            if isscalar(rhs_indices) && rhs_indices(1) == 1
                obj.captureOperator(ApplyOp, W);
            end
        end

        function W = computeWeightsBoundary(obj, X, nr, NeuCoeff, DirCoeff, sp, op, ApplyOp, rhs_indices)
            obj.InitializeGeometry(X, sp);
            rhs_inds = rhs_indices(:);
            X_subset = X(rhs_inds, :);
            X_at_origin_subset = obj.Xc(rhs_inds, :);
            nr_subset = nr(rhs_inds, :);
            r = kp.geometry.distanceMatrix(X, X);
            r_rhs = r(rhs_inds, :);
            B = obj.applyBoundaryOperator(ApplyOp, sp, op, NeuCoeff, DirCoeff, r_rhs, X_subset, X, X_at_origin_subset, obj.Xc, nr_subset);
            % Boundary operators always solve the local augmented system,
            % because the coefficients live in the BC definition itself.
            Wfull = kp.rbffd.RBFStencil.stableSolve(obj.solve_lhs, B);
            W = Wfull(1:obj.n, :);
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

        function captureOperator(obj, ApplyOp, W)
            switch lower(string(ApplyOp))
                case {"lap", "laplacian"}
                    obj.L = W;
                case {"interp", "interpolation"}
                    obj.Wt = W;
                case {"grad", "gradient"}
                    obj.Gx = W;
            end
        end
    end

    methods (Static, Access = private)
        function Phi = phsRbf(r, degree)
            Phi = r .^ degree;
            if mod(degree, 2) == 0
                Phi = Phi .* log(r + 2e-16);
            end
        end

        function D = phsDrOverR(r, degree)
            if mod(degree, 2) == 0
                D = r .^ (degree - 2) .* (degree * log(r + 2e-16) + 1);
            else
                D = degree * r .^ (degree - 2);
            end
            D(~isfinite(D)) = 0;
        end

        function L = phsLap(r, degree, dim)
            if mod(degree, 2) == 0
                L = r .^ (degree - 2) .* ...
                    (dim + 2 * degree + degree^2 * log(r + 2e-16) - 2 * degree * log(r + 2e-16) + dim * degree * log(r + 2e-16) - 2);
            else
                L = degree * (dim + degree - 2) * r .^ (degree - 2);
            end
            L(~isfinite(L)) = 0;
        end

        function X = stableSolve(A, B)
            warnNear = warning('query', 'MATLAB:nearlySingularMatrix');
            warnSing = warning('query', 'MATLAB:singularMatrix');
            cleanupObj = onCleanup(@() restoreWarnings(warnNear, warnSing)); %#ok<NASGU>
            warning('off', 'MATLAB:nearlySingularMatrix');
            warning('off', 'MATLAB:singularMatrix');
            X = A \ B;
            if any(~isfinite(X), 'all')
                X = pinv(A) * B;
            end
            X(~isfinite(X)) = 0;
        end
    end
end

function d = unitMultiIndex(dim, selectdim)
    d = zeros(1, dim);
    d(selectdim) = 1;
end

function restoreWarnings(warnNear, warnSing)
    warning(warnNear.state, 'MATLAB:nearlySingularMatrix');
    warning(warnSing.state, 'MATLAB:singularMatrix');
end
