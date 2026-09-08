classdef CrossNodeDiffOp < handle
    %CROSSNODEDIFFOP Assemble cross-domain interpolation/derivative operators.

    properties (SetAccess = private)
        locations double = zeros(0, 2)
        values double = zeros(0, 1)
        N1 (1,1) double = 0
        N2 (1,1) double = 0
        approxFactory function_handle = @() kp.rbffd.RBFStencil()
    end

    methods
        function obj = CrossNodeDiffOp(approxFactory)
            if nargin > 0
                obj.approxFactory = approxFactory;
            end
        end

        function AssembleOp(obj, sourceDomain, targetDomain, opName, stProps, opProps, varargin)
            parser = inputParser();
            parser.addParameter('ActiveRows', [], @(x) isnumeric(x));
            parser.parse(varargin{:});
            opts = parser.Results;

            [targetPoints, targetRowIds] = pickCenters(targetDomain, stProps.pointSet);
            activeRows = opts.ActiveRows(:);
            if isempty(activeRows)
                activeRows = (1:size(targetPoints, 1)).';
            end

            sourceGlobals = sourceDomain.getTreeGlobals(stProps.treeMode);
            obj.N1 = max(targetRowIds);
            obj.N2 = max(sourceGlobals);
            nrows = numel(activeRows);
            allIndices = cell(nrows, 1);
            allWeights = cell(nrows, 1);
            rowIds = zeros(nrows, 1);

            for k = 1:nrows
                localRow = activeRows(k);
                centerPoint = targetPoints(localRow, :);
                rowIds(k) = targetRowIds(localRow);
                [indices, ~] = sourceDomain.queryKnn(stProps.treeMode, centerPoint, stProps.n);
                indices = indices(1, :).';
                if numel(indices) < stProps.n
                    error('kp:rbffd:InsufficientStencilNodes', ...
                        ['Requested stencil size n=%d, but only %d nodes are available in tree mode "%s". ' ...
                         'Decrease the target order or use a larger domain / smaller h.'], ...
                        stProps.n, numel(indices), string(stProps.treeMode));
                end
                stencilPoints = sourceDomain.getTreePoints(stProps.treeMode);
                loc_x = stencilPoints(indices, :);
                stencil = obj.approxFactory();
                W = stencil.ComputeWeightsAtPoints(loc_x, centerPoint, stProps, opProps, opName);
                allIndices{k} = indices;
                allWeights{k} = W(:, 1);
            end

            tripletCount = sum(cellfun(@numel, allIndices));
            obj.locations = zeros(tripletCount, 2);
            obj.values = zeros(tripletCount, 1);
            cursor = 1;
            for k = 1:nrows
                cols = allIndices{k};
                W = allWeights{k};
                rowGlobal = rowIds(k);
                for j = 1:numel(cols)
                    obj.locations(cursor, :) = [rowGlobal, sourceGlobals(cols(j))];
                    obj.values(cursor) = W(j, 1);
                    cursor = cursor + 1;
                end
            end
        end

        function A = getOp(obj)
            if isempty(obj.locations)
                A = sparse(obj.N1, obj.N2);
                return;
            end
            A = sparse(obj.locations(:, 1), obj.locations(:, 2), obj.values, obj.N1, obj.N2);
        end
    end
end

function [points, rowIds] = pickCenters(domain, pointSet)
switch kp.rbffd.StencilProperties.normalizePointSet(pointSet)
    case "all"
        points = domain.getAllNodes();
    case "interior_boundary"
        points = domain.getIntBdryNodes();
    case "boundary"
        points = domain.getBdryNodes();
    otherwise
        error('kp:rbffd:BadPointSet', 'Unknown pointSet.');
end
rowIds = (1:size(points, 1)).';
end
