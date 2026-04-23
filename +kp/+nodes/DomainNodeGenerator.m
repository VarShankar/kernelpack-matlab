classdef DomainNodeGenerator < handle
    %DOMAINNODEGENERATOR KernelPack-shaped box Poisson node generator.

    properties (SetAccess = private)
        Xi double = zeros(0, 0)
        Xb double = zeros(0, 0)
        Xg double = zeros(0, 0)
        Nrmls double = zeros(0, 0)
        Xi_orig double = zeros(0, 0)
        Xi_pds_raw double = zeros(0, 0)
        s_dim (1,1) double = 0
        last_info struct = struct()
        descriptor kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
    end

    methods
        function generatePoissonNodes(obj, radius, x_min, x_max, varargin)
            [X, info] = kp.nodes.generatePoissonNodesInBox(radius, x_min, x_max, varargin{:});
            obj.Xi = X;
            obj.Xb = zeros(0, size(X, 2));
            obj.Xg = zeros(0, size(X, 2));
            obj.Nrmls = zeros(0, size(X, 2));
            obj.Xi_orig = X;
            obj.Xi_pds_raw = X;
            obj.s_dim = size(X, 2);
            obj.last_info = info;
            obj.descriptor = kp.domain.DomainDescriptor();
        end

        function generateInteriorNodesFromGeometry(obj, geometry, radius, varargin)
            parser = inputParser();
            parser.KeepUnmatched = true;
            parser.addParameter('DoOuterRefinement', false, @(x) islogical(x) || isnumeric(x));
            parser.addParameter('OuterFractionOfh', 0.5, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite', 'positive', '<=', 1}));
            parser.addParameter('OuterRefinementZoneSizeAsMultipleOfh', 2.0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
            parser.parse(varargin{:});
            localOpts = parser.Results;
            samplerArgs = namedArgsToCell(parser.Unmatched);

            [x_min, x_max] = kp.nodes.boundingBoxExtents(geometry, true);
            obj.generatePoissonNodes(radius, x_min, x_max, samplerArgs{:});

            if ~localOpts.DoOuterRefinement
                obj.clipToGeometry(geometry, 'Keep', 'inside', 'BoundaryClearance', radius);
                obj.last_info.min_active_radius = radius;
                return;
            end

            fineRadius = localOpts.OuterFractionOfh * radius;
            zoneSize = localOpts.OuterRefinementZoneSizeAsMultipleOfh * radius;
            minRadiusUsed = min(radius, fineRadius);

            levelSet = geometry.getLevelSet();
            if isempty(levelSet) || levelSet.n == 0
                if ismethod(geometry, 'buildLevelSetFromGeometricModel')
                    geometry.buildLevelSetFromGeometricModel([]);
                else
                    geometry.buildLevelSet();
                end
                levelSet = geometry.getLevelSet();
            end

            coarsePhi = levelSet.Evaluate(obj.Xi_pds_raw);
            coarseMask = (coarsePhi <= -zoneSize) & (coarsePhi <= -minRadiusUsed);
            coarseNodes = obj.Xi_pds_raw(coarseMask, :);

            [fineRaw, fineInfo] = kp.nodes.generatePoissonNodesInBox(fineRadius, x_min, x_max, samplerArgs{:});
            [fineNodes, fineMask, finePhi] = kp.nodes.clipPointsByGeometry(fineRaw, geometry, ...
                'Keep', 'inside', ...
                'BoundaryClearance', minRadiusUsed, ...
                'MinSignedDistance', -zoneSize, ...
                'MaxSignedDistance', -minRadiusUsed, ...
                'UseParallel', true);

            obj.Xi_pds_raw = [obj.Xi_pds_raw; fineRaw];
            obj.Xi = [coarseNodes; fineNodes];
            obj.Xi_orig = obj.Xi;
            obj.s_dim = size(obj.Xi, 2);
            obj.last_info.outer_refinement = struct( ...
                'enabled', true, ...
                'coarse_radius', radius, ...
                'fine_radius', fineRadius, ...
                'zone_size', zoneSize, ...
                'min_active_radius', minRadiusUsed, ...
                'coarse_points_kept', size(coarseNodes, 1), ...
                'fine_raw_points', size(fineRaw, 1), ...
                'fine_points_kept', size(fineNodes, 1), ...
                'fine_info', fineInfo, ...
                'fine_mask', fineMask, ...
                'fine_phi', finePhi(fineMask));
            obj.last_info.min_active_radius = minRadiusUsed;
        end

        function [X, mask, phi] = clipToGeometry(obj, geometry, varargin)
            [X, mask, phi] = kp.nodes.clipPointsByGeometry(obj.Xi_pds_raw, geometry, varargin{:});
            obj.Xi = X;
            obj.Xi_orig = X;
            obj.s_dim = size(X, 2);
            obj.last_info.clip_mask = mask;
            obj.last_info.clip_count = size(X, 1);
            obj.last_info.clip_levelset_values = phi(mask);
        end

        function descriptor = buildDomainDescriptorFromGeometry(obj, geometry, radius, varargin)
            obj.generateInteriorNodesFromGeometry(geometry, radius, varargin{:});
            [boundaryNodes, boundaryNormals, boundaryLevelSet] = buildBoundaryState(geometry);
            minRadiusUsed = radius;
            if isfield(obj.last_info, 'min_active_radius')
                minRadiusUsed = obj.last_info.min_active_radius;
            end
            ghostNodes = boundaryNodes + 0.5 * minRadiusUsed * boundaryNormals;

            descriptor = kp.domain.DomainDescriptor();
            descriptor.setNodes(obj.Xi, boundaryNodes, ghostNodes);
            descriptor.setNormals(boundaryNormals);
            descriptor.setSepRad(minRadiusUsed);
            descriptor.setOuterLevelSet(boundaryLevelSet);
            descriptor.setBoundaryLevelSets({boundaryLevelSet});
            descriptor.buildStructs();

            obj.Xb = boundaryNodes;
            obj.Xg = ghostNodes;
            obj.Nrmls = boundaryNormals;
            obj.descriptor = descriptor;
            obj.last_info.boundary_node_count = size(boundaryNodes, 1);
            obj.last_info.ghost_node_count = size(ghostNodes, 1);
            obj.last_info.min_active_radius = minRadiusUsed;
        end

        function out = getInteriorNodes(obj)
            out = obj.Xi;
        end

        function out = getBdryNodes(obj)
            out = obj.Xb;
        end

        function out = getGhostNodes(obj)
            out = obj.Xg;
        end

        function out = getNrmls(obj)
            out = obj.Nrmls;
        end

        function out = getRawPoissonInteriorNodes(obj)
            out = obj.Xi_pds_raw;
        end

        function out = getDomainDescriptor(obj)
            out = obj.descriptor;
        end
    end
end

function args = namedArgsToCell(s)
    names = fieldnames(s);
    args = cell(1, 2 * numel(names));
    for k = 1:numel(names)
        args{2 * k - 1} = names{k};
        args{2 * k} = s.(names{k});
    end
end

function [Xb, Nrmls, levelSet] = buildBoundaryState(geometry)
    if ismethod(geometry, 'getUniformBdryNodes')
        Xb = geometry.getUniformBdryNodes();
        Nrmls = geometry.getUniformBdryNrmls();
    elseif ismethod(geometry, 'getUniformSampleSites')
        Xb = geometry.getUniformSampleSites();
        Nrmls = geometry.getUniformNrmls();
    elseif ismethod(geometry, 'getBdryNodes')
        Xb = geometry.getBdryNodes();
        Nrmls = geometry.getBdryNrmls();
    else
        Xb = geometry.getSampleSites();
        Nrmls = geometry.getNrmls();
    end

    levelSet = geometry.getLevelSet();
    if isempty(levelSet) || ~isa(levelSet, 'kp.geometry.RBFLevelSet') || levelSet.n == 0
        if ismethod(geometry, 'buildLevelSetFromGeometricModel')
            geometry.buildLevelSetFromGeometricModel([]);
        else
            geometry.buildLevelSet();
        end
        levelSet = geometry.getLevelSet();
    end
end
