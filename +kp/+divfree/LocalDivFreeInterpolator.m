classdef LocalDivFreeInterpolator < handle
    %LOCALDIVFREEINTERPOLATOR Cached local divergence-free PHS+poly models.

    properties (SetAccess = private)
        Nodes double = zeros(0, 0)
        Values double = zeros(0, 0)
        PolyDegree (1,1) double = 0
        PhsDegree (1,1) double = 0
        StencilSize (1,1) double = 0
        Tree = []
        CenterModels cell = {}
        StencilIndices cell = {}
        ActiveCenters double = zeros(0, 1)
        Recurrence function_handle = @(N) kp.poly.jacobi_recurrence(N, 0, 0)
        PolyOptions struct = struct()
    end

    methods (Static)
        function obj = fit(X, U, polyDegree, phsDegree, stencilSize, varargin)
            obj = kp.divfree.LocalDivFreeInterpolator();
            obj.initialize(X, U, polyDegree, phsDegree, stencilSize, varargin{:});
        end
    end

    methods
        function initialize(obj, X, U, polyDegree, phsDegree, stencilSize, varargin)
            validateattributes(X, {'numeric'}, {'2d', 'finite', 'real'});
            validateattributes(U, {'numeric'}, {'2d', 'finite', 'real', ...
                'nrows', size(X, 1), 'ncols', size(X, 2)});
            validateattributes(polyDegree, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
            validateattributes(phsDegree, {'numeric'}, {'scalar', 'integer', 'positive'});
            validateattributes(stencilSize, {'numeric'}, {'scalar', 'integer', 'positive'});

            parser = inputParser();
            parser.addParameter('Recurrence', @(N) kp.poly.jacobi_recurrence(N, 0, 0), ...
                @(f) isa(f, 'function_handle'));
            parser.addParameter('PolyOptions', struct(), @(s) isstruct(s));
            parser.addParameter('ActiveCenters', [], ...
                @(v) isempty(v) || (isnumeric(v) && isvector(v)));
            parser.parse(varargin{:});

            obj.Nodes = X;
            obj.Values = U;
            obj.PolyDegree = polyDegree;
            obj.PhsDegree = phsDegree;
            obj.StencilSize = min(stencilSize, size(X, 1));
            obj.Recurrence = parser.Results.Recurrence;
            obj.PolyOptions = parser.Results.PolyOptions;
            obj.Tree = KDTreeSearcher(X);
            obj.CenterModels = cell(size(X, 1), 1);
            obj.StencilIndices = cell(size(X, 1), 1);

            if isempty(parser.Results.ActiveCenters)
                activeCenters = (1:size(X, 1)).';
            else
                activeCenters = unique(parser.Results.ActiveCenters(:), 'stable');
            end
            obj.ActiveCenters = activeCenters;

            for k = 1:numel(activeCenters)
                j = activeCenters(k);
                idx = knnsearch(obj.Tree, X(j, :), 'K', obj.StencilSize);
                obj.StencilIndices{j} = idx(:);
                obj.CenterModels{j} = kp.divfree.DivFreePHSInterpolant.fit( ...
                    X(idx, :), U(idx, :), obj.PolyDegree, obj.PhsDegree, ...
                    'Recurrence', obj.Recurrence, ...
                    'PolyOptions', obj.PolyOptions);
            end
        end

        function centerIdx = assignCenters(obj, Xq)
            validateattributes(Xq, {'numeric'}, {'2d', 'finite', 'real', ...
                'ncols', size(obj.Nodes, 2)});
            centerIdx = knnsearch(obj.Tree, Xq, 'K', 1);
        end

        function Uq = evaluate(obj, Xq, varargin)
            validateattributes(Xq, {'numeric'}, {'2d', 'finite', 'real', ...
                'ncols', size(obj.Nodes, 2)});

            parser = inputParser();
            parser.addParameter('CenterIndices', [], ...
                @(v) isempty(v) || (isnumeric(v) && isvector(v) && numel(v) == size(Xq, 1)));
            parser.parse(varargin{:});

            if isempty(parser.Results.CenterIndices)
                centerIdx = obj.assignCenters(Xq);
            else
                centerIdx = parser.Results.CenterIndices(:);
            end

            Uq = zeros(size(Xq, 1), size(obj.Values, 2));
            activeCenters = unique(centerIdx, 'stable');
            for k = 1:numel(activeCenters)
                j = activeCenters(k);
                mask = (centerIdx == j);
                model = obj.CenterModels{j};
                if isempty(model)
                    idx = knnsearch(obj.Tree, obj.Nodes(j, :), 'K', obj.StencilSize);
                    obj.StencilIndices{j} = idx(:);
                    model = kp.divfree.DivFreePHSInterpolant.fit( ...
                        obj.Nodes(idx, :), obj.Values(idx, :), obj.PolyDegree, obj.PhsDegree, ...
                        'Recurrence', obj.Recurrence, ...
                        'PolyOptions', obj.PolyOptions);
                    obj.CenterModels{j} = model;
                end
                Uq(mask, :) = model.evaluate(Xq(mask, :));
            end
        end
    end
end
