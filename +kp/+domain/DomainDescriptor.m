classdef DomainDescriptor < handle
    %DOMAINDESCRIPTOR Pared-down KernelPack-style domain storage object.

    properties (SetAccess = private)
        Xi double = zeros(0, 0)
        Xb double = zeros(0, 0)
        Xg double = zeros(0, 0)
        X double = zeros(0, 0)
        Xf double = zeros(0, 0)
        Nrmls double = zeros(0, 0)
        sepRad (1,1) double = NaN
        tallTree struct = struct('Points', [], 'Searcher', [], 'HasSearcher', false)
        intBdryTree struct = struct('Points', [], 'Searcher', [], 'HasSearcher', false)
        bdryTree struct = struct('Points', [], 'Searcher', [], 'HasSearcher', false)
        outerLevelSet = []
        boundaryLevelSets cell = {}
    end

    methods
        function setNodes(obj, int_nodes, bdry_nodes, ghost_nodes)
            % Store the three node classes separately, then rebuild the
            % packed node views that the assemblers query later.
            if nargin < 4 || isempty(ghost_nodes)
                ghost_nodes = zeros(0, size(int_nodes, 2));
            end
            obj.Xi = int_nodes;
            obj.Xb = bdry_nodes;
            obj.Xg = ghost_nodes;
            obj.setTotalNodes();
        end

        function setNormals(obj, nrmls)
            % Boundary normals must stay aligned with Xb because boundary
            % operators index into those two arrays row-for-row.
            if size(nrmls, 1) ~= size(obj.Xb, 1)
                error('kp:domain:BadNormals', 'Boundary normals must match the number of boundary nodes.');
            end
            obj.Nrmls = nrmls;
        end

        function setOuterLevelSet(obj, levelSet)
            obj.outerLevelSet = levelSet;
        end

        function setSepRad(obj, sepRad)
            validateattributes(sepRad, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
            obj.sepRad = sepRad;
        end

        function setBoundaryLevelSets(obj, levelSets)
            obj.boundaryLevelSets = levelSets;
        end

        function buildStructs(obj)
            % Build the three search layers that mirror the C++ split:
            % all nodes, interior+boundary nodes, and boundary nodes.
            obj.buildTallTree();
            obj.buildIntBdryTree();
            obj.buildBdryTree();
        end

        function buildBdryTree(obj)
            obj.bdryTree = obj.buildTreeStruct(obj.Xb);
        end

        function buildTallTree(obj)
            obj.tallTree = obj.buildTreeStruct(obj.Xf);
        end

        function buildIntBdryTree(obj)
            obj.intBdryTree = obj.buildTreeStruct(obj.X);
        end

        function [indices, distances] = queryKnn(obj, treeMode, queryPoints, k)
            treeMode = string(treeMode);
            [tree, points] = obj.getTreeData(treeMode);
            if isempty(points)
                indices = zeros(size(queryPoints, 1), 0);
                distances = zeros(size(queryPoints, 1), 0);
                return;
            end

            k = min(k, size(points, 1));
            if tree.HasSearcher
                % Prefer the KD-tree path when MATLAB has it, but keep a
                % deterministic dense fallback when that toolbox is absent.
                [indices, distances] = knnsearch(tree.Searcher, queryPoints, 'K', k);
                if isvector(indices)
                    indices = reshape(indices, size(queryPoints, 1), []);
                    distances = reshape(distances, size(queryPoints, 1), []);
                end
                return;
            end

            d = kp.geometry.distanceMatrix(queryPoints, points);
            [distances, order] = sort(d, 2, 'ascend');
            indices = order(:, 1:k);
            distances = distances(:, 1:k);
        end

        function [indices, distances] = queryBall(obj, treeMode, queryPoints, radius)
            treeMode = string(treeMode);
            [tree, points] = obj.getTreeData(treeMode);
            if isempty(points)
                indices = cell(size(queryPoints, 1), 1);
                distances = cell(size(queryPoints, 1), 1);
                return;
            end

            if tree.HasSearcher
                [indices, distances] = rangesearch(tree.Searcher, queryPoints, radius);
                return;
            end

            % The fallback path reproduces the same query contract using a
            % dense distance matrix and explicit masking.
            d = kp.geometry.distanceMatrix(queryPoints, points);
            indices = cell(size(queryPoints, 1), 1);
            distances = cell(size(queryPoints, 1), 1);
            for q = 1:size(queryPoints, 1)
                mask = d(q, :) <= radius;
                indices{q} = find(mask);
                distances{q} = d(q, mask);
            end
        end

        function points = getTreePoints(obj, treeMode)
            [~, points] = obj.getTreeData(string(treeMode));
        end

        function globals = getTreeGlobals(obj, treeMode)
            treeMode = string(treeMode);
            switch treeMode
                case "all"
                    globals = (1:size(obj.Xf, 1)).';
                case "interior_boundary"
                    globals = (1:size(obj.X, 1)).';
                case "boundary"
                    ni = obj.getNumInteriorNodes();
                    globals = (ni + (1:size(obj.Xb, 1))).';
                otherwise
                    error('kp:domain:BadTreeMode', 'Unknown tree mode "%s".', treeMode);
            end
        end

        function out = getInteriorNodes(obj), out = obj.Xi; end
        function out = getBdryNodes(obj), out = obj.Xb; end
        function out = getGhostNodes(obj), out = obj.Xg; end
        function out = getIntBdryNodes(obj), out = obj.X; end
        function out = getAllNodes(obj), out = obj.Xf; end
        function out = getNrmls(obj), out = obj.Nrmls; end
        function out = getTallTree(obj), out = obj.tallTree; end
        function out = getIntBdryTree(obj), out = obj.intBdryTree; end
        function out = getBdryTree(obj), out = obj.bdryTree; end
        function out = getOuterLevelSet(obj), out = obj.outerLevelSet; end
        function out = getBoundaryLevelSets(obj), out = obj.boundaryLevelSets; end
        function out = getSepRad(obj), out = obj.sepRad; end
        function out = getDim(obj), out = size(obj.Xf, 2); end
        function out = getNumTotalNodes(obj), out = size(obj.Xf, 1); end
        function out = getNumIntBdryNodes(obj), out = size(obj.X, 1); end
        function out = getNumInteriorNodes(obj), out = size(obj.Xi, 1); end
        function out = getNumBdryNodes(obj), out = size(obj.Xb, 1); end
        function out = getNumGhostNodes(obj), out = size(obj.Xg, 1); end
    end

    methods (Access = private)
        function setTotalNodes(obj)
            % Keep the packed node views synchronized with the separated
            % Xi/Xb/Xg storage expected by the solvers.
            if isempty(obj.Xi)
                dim = size(obj.Xb, 2);
            else
                dim = size(obj.Xi, 2);
            end
            if isempty(obj.Xg)
                obj.Xg = zeros(0, dim);
            end
            obj.X = [obj.Xi; obj.Xb];
            obj.Xf = [obj.X; obj.Xg];
        end

        function [tree, points] = getTreeData(obj, treeMode)
            % Tree mode names are the semantic MATLAB replacement for the
            % old integer tree selectors in KernelPack C++.
            switch treeMode
                case "all"
                    tree = obj.tallTree;
                    points = obj.Xf;
                case "interior_boundary"
                    tree = obj.intBdryTree;
                    points = obj.X;
                case "boundary"
                    tree = obj.bdryTree;
                    points = obj.Xb;
                otherwise
                    error('kp:domain:BadTreeMode', 'Unknown tree mode "%s".', treeMode);
            end
        end

        function tree = buildTreeStruct(~, points)
            tree = struct('Points', points, 'Searcher', [], 'HasSearcher', false);
            if isempty(points)
                return;
            end
            % A missing KD-tree toolbox should slow us down, not break the
            % assembly pipeline.
            if exist('KDTreeSearcher', 'class') == 8 && exist('knnsearch', 'file') == 2
                tree.Searcher = KDTreeSearcher(points);
                tree.HasSearcher = true;
            end
        end
    end
end
