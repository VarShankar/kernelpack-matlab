classdef DomainDescriptor < handle
    %DOMAINDESCRIPTOR Pared-down KernelPack-style domain storage object.

    properties (SetAccess = private)
        Xi double = zeros(0, 0)
        Xb double = zeros(0, 0)
        Xg double = zeros(0, 0)
        X double = zeros(0, 0)
        Xf double = zeros(0, 0)
        Nrmls double = zeros(0, 0)
        tallTree struct = struct('Points', [])
        intBdryTree struct = struct('Points', [])
        bdryTree struct = struct('Points', [])
        outerLevelSet = []
        boundaryLevelSets cell = {}
    end

    methods
        function setNodes(obj, int_nodes, bdry_nodes, ghost_nodes)
            if nargin < 4 || isempty(ghost_nodes)
                ghost_nodes = zeros(0, size(int_nodes, 2));
            end
            obj.Xi = int_nodes;
            obj.Xb = bdry_nodes;
            obj.Xg = ghost_nodes;
            obj.setTotalNodes();
        end

        function setNormals(obj, nrmls)
            if size(nrmls, 1) ~= size(obj.Xb, 1)
                error('kp:domain:BadNormals', 'Boundary normals must match the number of boundary nodes.');
            end
            obj.Nrmls = nrmls;
        end

        function setOuterLevelSet(obj, levelSet)
            obj.outerLevelSet = levelSet;
        end

        function setBoundaryLevelSets(obj, levelSets)
            obj.boundaryLevelSets = levelSets;
        end

        function buildStructs(obj)
            obj.buildTallTree();
            obj.buildIntBdryTree();
            obj.buildBdryTree();
        end

        function buildBdryTree(obj)
            obj.bdryTree = struct('Points', obj.Xb);
        end

        function buildTallTree(obj)
            obj.tallTree = struct('Points', obj.Xf);
        end

        function buildIntBdryTree(obj)
            obj.intBdryTree = struct('Points', obj.X);
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
    end

    methods (Access = private)
        function setTotalNodes(obj)
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
    end
end
