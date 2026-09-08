classdef PolynomialBasis < handle
    %POLYNOMIALBASIS Polynomial basis with optional geometric normalization.
    %
    % By default this object represents a Legendre tensor-product basis.
    % The basis can store a center and scalar scale factor so stencil points
    % can be shifted and scaled into a unit disk/sphere consistently before
    % polynomial evaluation.

    properties (SetAccess = private)
        Family (1, :) char = 'legendre'
        Alpha (1, 1) double = 0
        Beta (1, 1) double = 0
        IndexSet double = zeros(0, 0)
        Center double = zeros(1, 0)
        Scale (1, 1) double = 1
        Dimension (1, 1) double = 0
    end

    methods
        function obj = PolynomialBasis(indexSet, varargin)
            if nargin < 1 || isempty(indexSet)
                indexSet = zeros(1, 1);
            end

            % The basis object stores both the polynomial family and the
            % geometric normalization that maps a local stencil into a
            % unit-sized computational patch.
            parser = inputParser();
            parser.addRequired('indexSet', @(x) validateattributes(x, {'numeric'}, {'2d', 'nonnegative', 'integer'}));
            parser.addParameter('Family', 'legendre', @(x) any(strcmpi(x, {'legendre', 'jacobi', 'chebyshev'})));
            parser.addParameter('Alpha', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite'}));
            parser.addParameter('Beta', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite'}));
            parser.addParameter('Center', [], @(x) isempty(x) || (isnumeric(x) && isvector(x) && all(isfinite(x))));
            parser.addParameter('Scale', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
            parser.parse(indexSet, varargin{:});

            obj.IndexSet = parser.Results.indexSet;
            obj.Dimension = size(obj.IndexSet, 2);
            obj.Family = lower(parser.Results.Family);
            obj.Alpha = parser.Results.Alpha;
            obj.Beta = parser.Results.Beta;
            if strcmp(obj.Family, 'legendre')
                obj.Alpha = 0;
                obj.Beta = 0;
            elseif strcmp(obj.Family, 'chebyshev')
                obj.Alpha = 0.5;
                obj.Beta = 0.5;
            end

            if isempty(parser.Results.Center)
                obj.Center = zeros(1, obj.Dimension);
            else
                obj.Center = reshape(parser.Results.Center, 1, []);
            end
            if numel(obj.Center) ~= obj.Dimension
                error('kp:poly:CenterDimensionMismatch', 'Center must have one entry per dimension.');
            end
            obj.Scale = parser.Results.Scale;
        end

        function setNormalization(obj, center, scale)
            center = reshape(center, 1, []);
            if numel(center) ~= obj.Dimension
                error('kp:poly:CenterDimensionMismatch', 'Center must have one entry per dimension.');
            end
            validateattributes(scale, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
            obj.Center = center;
            obj.Scale = scale;
        end

        function fitNormalizationFromPoints(obj, X)
            % Fit a simple center-and-radius normalization directly from a
            % point cloud so local polynomial work happens on O(1) scales.
            validateattributes(X, {'numeric'}, {'2d', 'finite', 'real', 'ncols', obj.Dimension});
            center = mean(X, 1);
            shifted = X - center;
            radii = sqrt(sum(shifted.^2, 2));
            scale = max(radii);
            if ~(scale > 0)
                scale = 1;
            end
            obj.Center = center;
            obj.Scale = scale;
        end

        function Xn = normalizePoints(obj, X)
            validateattributes(X, {'numeric'}, {'2d', 'finite', 'real', 'ncols', obj.Dimension});
            Xn = (X - obj.Center) / obj.Scale;
        end

        function X = denormalizePoints(obj, Xn)
            validateattributes(Xn, {'numeric'}, {'2d', 'finite', 'real', 'ncols', obj.Dimension});
            X = obj.Scale * Xn + obj.Center;
        end

        function p = evaluate(obj, X, d, assumeNormalized)
            if nargin < 3 || isempty(d)
                d = zeros(1, obj.Dimension);
            end
            if nargin < 4
                assumeNormalized = false;
            end
            % Callers can either provide already-normalized coordinates or
            % let the basis object do the shift-and-scale internally.
            if assumeNormalized
                Xwork = X;
            else
                Xwork = obj.normalizePoints(X);
            end

            % The actual basis evaluation is delegated to the shared
            % Jacobi/Chebyshev tensor-product machinery.
            recurrence = @(N) obj.getRecurrence(N);
            p = kp.poly.mpoly_eval(Xwork, obj.IndexSet, recurrence, d);

            if any(d(:) > 0)
                % Derivatives in normalized coordinates need the usual
                % chain-rule correction back to the physical scale.
                totalOrder = sum(d, 2);
                for q = 1:numel(totalOrder)
                    p(:, :, q) = p(:, :, q) / (obj.Scale^totalOrder(q));
                end
            end
        end

        function [a, b] = getRecurrence(obj, N)
            switch obj.Family
                case 'legendre'
                    [a, b] = kp.poly.jacobi_recurrence(N, 0, 0);
                case 'jacobi'
                    [a, b] = kp.poly.jacobi_recurrence(N, obj.Alpha, obj.Beta);
                case 'chebyshev'
                    [a, b] = kp.poly.chebyshev_recurrence(N);
                otherwise
                    error('kp:poly:UnknownFamily', 'Unknown polynomial family.');
            end
        end
    end

    methods (Static)
        function obj = fromTotalDegree(dim, degree, varargin)
            indexSet = kp.poly.total_degree_indices(dim, degree);
            obj = kp.poly.PolynomialBasis(indexSet, varargin{:});
        end

        function obj = fromHyperbolicCross(dim, degree, varargin)
            indexSet = kp.poly.hyperbolic_cross_indices(dim, degree);
            obj = kp.poly.PolynomialBasis(indexSet, varargin{:});
        end
    end
end
