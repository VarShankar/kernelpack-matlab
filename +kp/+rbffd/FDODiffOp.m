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

            % Start from the same center cloud as standard RBF-FD, but let
            % one accepted stencil row serve several nearby rows.
            [centerPoints, centerRowIds, centerColGlobals, centerNormals] = pickCenters(domain, stProps.pointSet);
            stencilPoints = domain.getTreePoints(stProps.treeMode);
            stencilGlobals = domain.getTreeGlobals(stProps.treeMode);
            activeRows = opts.ActiveRows(:);
            if isempty(activeRows)
                activeRows = (1:size(centerPoints, 1)).';
            end
            activeSet = true(size(centerPoints, 1), 1);
            activeSet(:) = false;
            activeSet(activeRows) = true;

            useBoundary = ~isempty(opts.NeuCoeff) || ~isempty(opts.DirCoeff);
            % KernelPack uses a fixed half-stencil acceptance neighborhood
            % for the overlapped path.
            loc_lim = max(1, floor(0.5 * stProps.n));
            rowToLocal = containers.Map('KeyType', 'double', 'ValueType', 'double');
            for k = 1:size(centerPoints, 1)
                rowToLocal(centerColGlobals(k)) = k;
            end

            obj.N1 = max(centerRowIds);
            obj.N2 = max(stencilGlobals);
            tripletLocations = zeros(numel(activeRows) * stProps.n, 2);
            tripletValues = zeros(numel(activeRows) * stProps.n, 1);
            cursor = 1;
            obj.stencils = {};
            obj.recorded_stencil_centers = {};
            obj.recorded_stencil_globals = zeros(0, 1);

            % Greedily accept overlapped rows until every requested row has
            % been covered by some accepted stencil.
            while any(activeSet)
                localCenter = find(activeSet, 1, 'first');
                centerPoint = centerPoints(localCenter, :);
                [indices, ~] = domain.queryKnn(stProps.treeMode, centerPoint, stProps.n);
                indices = indices(1, :).';
                if numel(indices) < stProps.n
                    error('kp:rbffd:InsufficientStencilNodes', ...
                        ['Requested stencil size n=%d, but only %d nodes are available in tree mode "%s". ' ...
                         'Decrease the target order or use a larger domain / smaller h.'], ...
                        stProps.n, numel(indices), string(stProps.treeMode));
                end
                rhs_indices = 1:min(loc_lim, numel(indices));
                loc_x = stencilPoints(indices, :);
                stencil = obj.approxFactory();
                if useBoundary
                    loc_nrmls = zeros(size(loc_x));
                    for q = 1:numel(indices)
                        candidateGlobal = stencilGlobals(indices(q));
                        if isKey(rowToLocal, candidateGlobal)
                            loc_nrmls(q, :) = centerNormals(rowToLocal(candidateGlobal), :);
                        end
                    end
                    W = stencil.ComputeWeights(loc_x, loc_nrmls, opts.NeuCoeff(localCenter), opts.DirCoeff(localCenter), stProps, opProps, op, rhs_indices);
                else
                    W = stencil.ComputeWeights(loc_x, stProps, opProps, op, rhs_indices);
                end
                A = stencil.getInterpMat();
                lebesgue = sum(abs(W(1:stProps.n, :)), 1);
                native = zeros(1, size(W, 2));
                useEnergyMetric = ~isempty(A) && size(A, 1) == size(W, 1) && size(A, 2) == size(W, 1);
                for j = 1:size(W, 2)
                    if useEnergyMetric
                        native(j) = abs(W(:, j).' * (A * W(:, j)));
                    else
                        native(j) = sum(abs(W(:, j)).^2);
                    end
                end
                % The center row defines the acceptance thresholds for the
                % other candidate rows harvested from the same stencil.
                leb0 = lebesgue(1);
                nat0 = native(1);

                acceptedAny = false;
                acceptedGlobals = zeros(0, 1);
                for j = 1:min(loc_lim, numel(indices))
                    candidateColGlobal = stencilGlobals(indices(j));
                    if ~isKey(rowToLocal, candidateColGlobal)
                        continue;
                    end
                    candidateLocal = rowToLocal(candidateColGlobal);
                    if ~activeSet(candidateLocal)
                        continue;
                    end
                    if lebesgue(j) > leb0 || native(j) > nat0
                        continue;
                    end
                    for q = 1:stProps.n
                        tripletLocations(cursor, :) = [centerRowIds(candidateLocal), stencilGlobals(indices(q))];
                        tripletValues(cursor) = W(q, j);
                        cursor = cursor + 1;
                    end
                    activeSet(candidateLocal) = false;
                    acceptedGlobals(end+1, 1) = candidateColGlobal; %#ok<AGROW>
                    acceptedAny = true;
                end

                if ~acceptedAny
                    error('kp:rbffd:FDODiffOpStall', 'Current overlapped stencil could not accept any active row.');
                end

                obj.recorded_stencil_centers{end+1,1} = centerPoint; %#ok<AGROW>
                obj.recorded_stencil_globals(end+1,1) = centerColGlobals(localCenter); %#ok<AGROW>
                if opProps.recordStencils
                    obj.stencils{end+1,1} = struct( ...
                        'Approx', stencil, ...
                        'Indices', stencilGlobals(indices), ...
                        'AcceptedGlobals', acceptedGlobals); %#ok<AGROW>
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
end
