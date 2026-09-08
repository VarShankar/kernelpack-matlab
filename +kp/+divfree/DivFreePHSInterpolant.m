classdef DivFreePHSInterpolant < handle
    %DIVFREEPHSINTERPOLANT Local divergence-free PHS+poly interpolant.

    properties (SetAccess = private)
        Dim (1,1) double = 0
        PhsDegree (1,1) double = 0
        PolyDegree (1,1) double = 0
        Nodes double = zeros(0, 0)
        Values double = zeros(0, 0)
        Kernel cell = {}
        PolyData struct = struct()
        RBFCoeffs double = zeros(0, 1)
        PolyCoeffs double = zeros(0, 1)
        SaddleMatrix double = zeros(0, 0)
    end

    methods (Static)
        function obj = fit(X, U, polyDegree, phsDegree, varargin)
            obj = kp.divfree.DivFreePHSInterpolant();
            obj.initialize(X, U, polyDegree, phsDegree, varargin{:});
        end
    end

    methods
        function initialize(obj, X, U, polyDegree, phsDegree, varargin)
            validateattributes(X, {'numeric'}, {'2d', 'finite', 'real'});
            validateattributes(U, {'numeric'}, {'2d', 'finite', 'real', 'nrows', size(X, 1), 'ncols', size(X, 2)});
            validateattributes(polyDegree, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
            validateattributes(phsDegree, {'numeric'}, {'scalar', 'integer', 'positive'});

            obj.Dim = size(X, 2);
            obj.PhsDegree = phsDegree;
            obj.PolyDegree = polyDegree;
            obj.Nodes = X;
            obj.Values = U;

            parser = inputParser();
            parser.addParameter('Recurrence', @(N) kp.poly.jacobi_recurrence(N, 0, 0), @(f) isa(f, 'function_handle'));
            parser.addParameter('PolyOptions', struct(), @(s) isstruct(s));
            parser.parse(varargin{:});

            [obj.Kernel, ~, ~, ~] = kp.divfree.DFPHS(obj.Dim, obj.PhsDegree);
            r = sqrt(max(kp.geometry.distanceMatrix(X, X).^2, 0));
            A = kp.divfree.DivFreeGram(obj.Kernel, X, X, r);
            [P, poly] = kp.divfree.df_poly_basis_from_jacobi(X, polyDegree, parser.Results.Recurrence, parser.Results.PolyOptions);

            rhs = stackField(U);
            saddle = [A, P; P.', zeros(size(P, 2), size(P, 2))];
            coeffs = saddle \ [rhs; zeros(size(P, 2), 1)];

            nRbf = size(A, 2);
            obj.RBFCoeffs = coeffs(1:nRbf);
            obj.PolyCoeffs = coeffs(nRbf + 1:end);
            obj.PolyData = poly;
            obj.SaddleMatrix = saddle;
        end

        function Uq = evaluate(obj, Xq)
            validateattributes(Xq, {'numeric'}, {'2d', 'finite', 'real', 'ncols', obj.Dim});
            r = kp.geometry.distanceMatrix(Xq, obj.Nodes);
            Aeval = kp.divfree.DivFreeGram(obj.Kernel, Xq, obj.Nodes, r);
            Pq = obj.PolyData.eval(Xq);
            stacked = Aeval * obj.RBFCoeffs + Pq * obj.PolyCoeffs;
            Uq = unstackField(stacked, obj.Dim);
        end
    end
end

function rhs = stackField(U)
dim = size(U, 2);
rhs = zeros(numel(U), 1);
for d = 1:dim
    rows = (d - 1) * size(U, 1) + 1:d * size(U, 1);
    rhs(rows) = U(:, d);
end
end

function U = unstackField(rhs, dim)
n = numel(rhs) / dim;
U = zeros(n, dim);
for d = 1:dim
    rows = (d - 1) * n + 1:d * n;
    U(:, d) = rhs(rows);
end
end
