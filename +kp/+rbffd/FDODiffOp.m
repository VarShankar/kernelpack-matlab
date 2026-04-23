classdef FDODiffOp < handle
    %FDODIFFOP Overlapped stencil assembler.

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
        function obj = FDODiffOp(approxFactory)
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
            activeSet = true(size(centerPoints, 1), 1);
            activeSet(:) = false;
            activeSet(activeRows) = true;

            useBoundary = ~isempty(opts.NeuCoeff) || ~isempty(opts.DirCoeff);
            loc_lim = max(1, floor(opProps.OverlapLoad * stProps.n));
            rowToLocal = containers.Map('KeyType', 'double', 'ValueType', 'double');
            for k = 1:size(centerPoints, 1)
                rowToLocal(centerGlobals(k)) = k;
            end

            obj.N1 = max(centerGlobals);
            obj.N2 = max(stencilGlobals);
            tripletLocations = zeros(numel(activeRows) * stProps.n, 2);
            tripletValues = zeros(numel(activeRows) * stProps.n, 1);
            cursor = 1;
            obj.stencils = {};
            obj.recorded_stencil_centers = {};
            obj.recorded_stencil_globals = zeros(0, 1);

            while any(activeSet)
                localCenter = find(activeSet, 1, 'first');
                centerPoint = centerPoints(localCenter, :);
                d = kp.geometry.distanceMatrix(centerPoint, stencilPoints);
                [~, order] = sort(d, 2, 'ascend');
                indices = order(1:stProps.n);
                rhs_indices = 1:min(loc_lim, numel(indices));
                loc_x = stencilPoints(indices, :);
                stencil = obj.approxFactory();
                if useBoundary
                    W = stencil.ComputeWeights(loc_x, centerNormals(localCenter, :), opts.NeuCoeff(localCenter), opts.DirCoeff(localCenter), stProps, opProps, op, rhs_indices);
                else
                    W = stencil.ComputeWeights(loc_x, stProps, opProps, op, rhs_indices);
                end
                A = stencil.getInterpMat();
                lebesgue = sum(abs(W(1:stProps.n, :)), 1);
                native = zeros(1, size(W, 2));
                for j = 1:size(W, 2)
                    native(j) = abs(W(:, j).' * (A * W(:, j)));
                end
                leb0 = lebesgue(1);
                nat0 = native(1);

                acceptedAny = false;
                for j = 1:min(loc_lim, numel(indices))
                    candidateGlobal = stencilGlobals(indices(j));
                    if ~isKey(rowToLocal, candidateGlobal)
                        continue;
                    end
                    candidateLocal = rowToLocal(candidateGlobal);
                    if ~activeSet(candidateLocal)
                        continue;
                    end
                    if lebesgue(j) > leb0 || native(j) > nat0
                        continue;
                    end
                    for q = 1:stProps.n
                        tripletLocations(cursor, :) = [candidateGlobal, stencilGlobals(indices(q))];
                        tripletValues(cursor) = W(q, j);
                        cursor = cursor + 1;
                    end
                    activeSet(candidateLocal) = false;
                    acceptedAny = true;
                end

                if ~acceptedAny
                    error('kp:rbffd:FDODiffOpStall', 'Current overlapped stencil could not accept any active row.');
                end

                obj.recorded_stencil_centers{end+1,1} = centerPoint; %#ok<AGROW>
                obj.recorded_stencil_globals(end+1,1) = centerGlobals(localCenter); %#ok<AGROW>
                if opProps.recordStencils
                    obj.stencils{end+1,1} = struct('Approx', stencil, 'Indices', stencilGlobals(indices)); %#ok<AGROW>
                end
            end

            obj.locations = tripletLocations(1:cursor-1, :);
            obj.values = tripletValues(1:cursor-1);
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

function [points, globals, normals] = pickCenters(domain, pointSet)
    normals = [];
    switch pointSet
        case 0
            points = domain.getAllNodes();
            globals = (1:size(points, 1)).';
        case 1
            points = domain.getIntBdryNodes();
            globals = (1:size(points, 1)).';
        case 2
            points = domain.getBdryNodes();
            ni = domain.getNumInteriorNodes();
            globals = (ni + (1:size(points, 1))).';
            normals = domain.getNrmls();
        otherwise
            error('kp:rbffd:BadPointSet', 'Unknown pointSet.');
    end
end

function [points, globals] = pickStencilPoints(domain, treeMode)
    switch treeMode
        case 0
            points = domain.getAllNodes();
            globals = (1:size(points, 1)).';
        case 1
            points = domain.getIntBdryNodes();
            globals = (1:size(points, 1)).';
        case 2
            points = domain.getBdryNodes();
            ni = domain.getNumInteriorNodes();
            globals = (ni + (1:size(points, 1))).';
        otherwise
            error('kp:rbffd:BadTreeMode', 'Unknown treeMode.');
    end
end
