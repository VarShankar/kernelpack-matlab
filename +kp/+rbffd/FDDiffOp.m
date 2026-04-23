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

            [centerPoints, centerGlobals, centerNormals] = pickCenters(domain, stProps.pointSet);
            [stencilPoints, stencilGlobals] = pickStencilPoints(domain, stProps.treeMode);
            activeRows = opts.ActiveRows(:);
            if isempty(activeRows)
                activeRows = (1:size(centerPoints, 1)).';
            end

            obj.N1 = max(centerGlobals);
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

            useParallel = opProps.UseParallel && license('test', 'Distrib_Computing_Toolbox') && nrows > 8;
            if useParallel
                parfor k = 1:nrows
                    localRow = activeRows(k);
                    [allIndices{k}, allWeights{k}, allStencils{k}, centersRecorded{k}, globalsRecorded(k)] = ...
                        assembleOne(obj.approxFactory, centerPoints, centerGlobals, centerNormals, stencilPoints, stencilGlobals, localRow, stProps, opProps, op, useBoundary, opts.NeuCoeff, opts.DirCoeff);
                end
            else
                for k = 1:nrows
                    localRow = activeRows(k);
                    [allIndices{k}, allWeights{k}, allStencils{k}, centersRecorded{k}, globalsRecorded(k)] = ...
                        assembleOne(obj.approxFactory, centerPoints, centerGlobals, centerNormals, stencilPoints, stencilGlobals, localRow, stProps, opProps, op, useBoundary, opts.NeuCoeff, opts.DirCoeff);
                end
            end

            tripletCount = nrows * stProps.n;
            obj.locations = zeros(tripletCount, 2);
            obj.values = zeros(tripletCount, 1);
            cursor = 1;
            obj.stencils = {};
            obj.recorded_stencil_centers = {};
            obj.recorded_stencil_globals = globalsRecorded;
            for k = 1:nrows
                cols = allIndices{k};
                W = allWeights{k};
                rowGlobal = globalsRecorded(k);
                for j = 1:stProps.n
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

function [indices, W, stencil, centerPoint, centerGlobal] = assembleOne(approxFactory, centerPoints, centerGlobals, centerNormals, stencilPoints, ~, localRow, stProps, opProps, op, useBoundary, NeuCoeff, DirCoeff)
    centerPoint = centerPoints(localRow, :);
    centerGlobal = centerGlobals(localRow);
    d = kp.geometry.distanceMatrix(centerPoint, stencilPoints);
    [~, order] = sort(d, 2, 'ascend');
    indices = order(1:stProps.n);
    loc_x = stencilPoints(indices, :);
    stencil = approxFactory();
    if useBoundary
        loc_nr = centerNormals(localRow, :);
        W = stencil.ComputeWeights(loc_x, loc_nr, NeuCoeff(localRow), DirCoeff(localRow), stProps, opProps, op, 1);
    else
        W = stencil.ComputeWeights(loc_x, stProps, opProps, op, 1);
    end
end

function [points, globals, normals] = pickCenters(domain, pointSet)
    normals = [];
    switch kp.rbffd.StencilProperties.normalizePointSet(pointSet)
        case "all"
            points = domain.getAllNodes();
            globals = (1:size(points, 1)).';
        case "interior_boundary"
            points = domain.getIntBdryNodes();
            globals = (1:size(points, 1)).';
        case "boundary"
            points = domain.getBdryNodes();
            ni = domain.getNumInteriorNodes();
            globals = (ni + (1:size(points, 1))).';
            normals = domain.getNrmls();
        otherwise
            error('kp:rbffd:BadPointSet', 'Unknown pointSet.');
    end
    if isempty(normals) && kp.rbffd.StencilProperties.normalizePointSet(pointSet) == "boundary"
        normals = domain.getNrmls();
    end
end

function [points, globals] = pickStencilPoints(domain, treeMode)
    switch kp.rbffd.StencilProperties.normalizeTreeMode(treeMode)
        case "all"
            points = domain.getAllNodes();
            globals = (1:size(points, 1)).';
        case "interior_boundary"
            points = domain.getIntBdryNodes();
            globals = (1:size(points, 1)).';
        case "boundary"
            points = domain.getBdryNodes();
            ni = domain.getNumInteriorNodes();
            globals = (ni + (1:size(points, 1))).';
        otherwise
            error('kp:rbffd:BadTreeMode', 'Unknown treeMode.');
    end
end
