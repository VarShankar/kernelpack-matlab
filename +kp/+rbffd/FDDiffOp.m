classdef FDDiffOp < handle
    %FDDIFFOP One-row-per-center stencil assembler.

    properties (SetAccess = private)
        locations double = zeros(0, 2)
        values double = zeros(0, 1)
        N1 (1,1) double = 0
        N2 (1,1) double = 0
        stencils cell = {}
        recorded_stencil_centers cell = {}
        recorded_stencil_globals double = zeros(0, 1)
        approxFactory function_handle = @() kp.rbffd.RBFStencil()
    end

    methods
        function obj = FDDiffOp(approxFactory)
            if nargin > 0
                obj.approxFactory = approxFactory;
            end
        end

        function AssembleOp(obj, domain, op, stProps, opProps, varargin)
            parser = inputParser();
            parser.addParameter('NeuCoeff', [], @(x) isnumeric(x));
            parser.addParameter('DirCoeff', [], @(x) isnumeric(x));
            parser.addParameter('ActiveRows', [], @(x) isnumeric(x));
            parser.parse(varargin{:});
            opts = parser.Results;

            % Pick the row centers first, then query the chosen tree for
            % each local stencil neighborhood.
            [centerPoints, centerRowIds, centerColGlobals, centerNormals] = pickCenters(domain, stProps.pointSet);
            activeRows = opts.ActiveRows(:);
            if isempty(activeRows)
                activeRows = (1:size(centerPoints, 1)).';
            end

            stencilGlobals = domain.getTreeGlobals(stProps.treeMode);
            obj.N1 = max(centerRowIds);
            obj.N2 = max(stencilGlobals);
            nrows = numel(activeRows);
            allIndices = cell(nrows, 1);
            allWeights = cell(nrows, 1);
            allStencils = cell(nrows, 1);
            centersRecorded = cell(nrows, 1);
            globalsRecorded = zeros(nrows, 1);

            useBoundary = ~isempty(opts.NeuCoeff) || ~isempty(opts.DirCoeff);
            if useBoundary
                assert(numel(opts.NeuCoeff) == size(centerPoints, 1) && numel(opts.DirCoeff) == size(centerPoints, 1), ...
                    'Boundary coefficients must align with the chosen center cloud.');
            end

            % Every requested row is independent, so standard RBF-FD
            % assembly parallelizes naturally over centers.
            useParallel = opProps.UseParallel && license('test', 'Distrib_Computing_Toolbox') && nrows > 8;
            if useParallel
                parfor k = 1:nrows
                    localRow = activeRows(k);
                    [allIndices{k}, allWeights{k}, allStencils{k}, centersRecorded{k}, rowIds(k), globalsRecorded(k)] = ...
                        assembleOne(obj.approxFactory, domain, centerPoints, centerRowIds, centerColGlobals, centerNormals, localRow, stProps, opProps, op, useBoundary, opts.NeuCoeff, opts.DirCoeff);
                end
            else
                rowIds = zeros(nrows, 1);
                for k = 1:nrows
                    localRow = activeRows(k);
                    [allIndices{k}, allWeights{k}, allStencils{k}, centersRecorded{k}, rowIds(k), globalsRecorded(k)] = ...
                        assembleOne(obj.approxFactory, domain, centerPoints, centerRowIds, centerColGlobals, centerNormals, localRow, stProps, opProps, op, useBoundary, opts.NeuCoeff, opts.DirCoeff);
                end
            end

            % Store the result as triplets first; the public operator is
            % materialized as a sparse matrix on demand.
            tripletCount = sum(cellfun(@numel, allIndices));
            obj.locations = zeros(tripletCount, 2);
            obj.values = zeros(tripletCount, 1);
            cursor = 1;
            obj.stencils = {};
            obj.recorded_stencil_centers = {};
            obj.recorded_stencil_globals = globalsRecorded;
            for k = 1:nrows
                cols = allIndices{k};
                W = allWeights{k};
                rowGlobal = rowIds(k);
                for j = 1:numel(cols)
                    obj.locations(cursor, :) = [rowGlobal, stencilGlobals(cols(j))];
                    obj.values(cursor) = W(j, 1);
                    cursor = cursor + 1;
                end
                obj.recorded_stencil_centers{end+1,1} = centersRecorded{k}; %#ok<AGROW>
                if opProps.recordStencils
                    obj.stencils{end+1,1} = struct('Approx', allStencils{k}, 'Indices', stencilGlobals(cols)); %#ok<AGROW>
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

function [indices, W, stencil, centerPoint, centerRowId, centerColGlobal] = assembleOne(approxFactory, domain, centerPoints, centerRowIds, centerColGlobals, centerNormals, localRow, stProps, opProps, op, useBoundary, NeuCoeff, DirCoeff)
    centerPoint = centerPoints(localRow, :);
    centerRowId = centerRowIds(localRow);
    centerColGlobal = centerColGlobals(localRow);
    [indices, ~] = domain.queryKnn(stProps.treeMode, centerPoint, stProps.n);
    indices = indices(1, :).';
    if numel(indices) < stProps.n
        error('kp:rbffd:InsufficientStencilNodes', ...
            ['Requested stencil size n=%d, but only %d nodes are available in tree mode "%s". ' ...
             'Decrease the target order or use a larger domain / smaller h.'], ...
            stProps.n, numel(indices), string(stProps.treeMode));
    end
    stencilPoints = domain.getTreePoints(stProps.treeMode);
    loc_x = stencilPoints(indices, :);
    stencil = approxFactory();
    if useBoundary
        loc_nr = centerNormals(localRow, :);
        W = stencil.ComputeWeights(loc_x, loc_nr, NeuCoeff(localRow), DirCoeff(localRow), stProps, opProps, op, 1);
    else
        W = stencil.ComputeWeights(loc_x, stProps, opProps, op, 1);
    end
end

function [points, rowIds, colGlobals, normals] = pickCenters(domain, pointSet)
    normals = [];
    switch kp.rbffd.StencilProperties.normalizePointSet(pointSet)
        case "all"
            points = domain.getAllNodes();
            rowIds = (1:size(points, 1)).';
            colGlobals = rowIds;
        case "interior_boundary"
            points = domain.getIntBdryNodes();
            rowIds = (1:size(points, 1)).';
            colGlobals = rowIds;
        case "boundary"
            points = domain.getBdryNodes();
            ni = domain.getNumInteriorNodes();
            rowIds = (1:size(points, 1)).';
            colGlobals = (ni + (1:size(points, 1))).';
            normals = domain.getNrmls();
        otherwise
            error('kp:rbffd:BadPointSet', 'Unknown pointSet.');
    end
    if isempty(normals) && kp.rbffd.StencilProperties.normalizePointSet(pointSet) == "boundary"
        normals = domain.getNrmls();
    end
end
