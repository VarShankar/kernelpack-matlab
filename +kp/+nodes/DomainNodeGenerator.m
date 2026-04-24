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
            % Box-only generation is the lowest-level entry point; no
            % geometry clipping or boundary bookkeeping happens here.
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
            parser.addParameter('RadiusFunction', [], @(x) isempty(x) || isa(x, 'function_handle'));
            parser.addParameter('DoOuterRefinement', false, @(x) islogical(x) || isnumeric(x));
            parser.addParameter('OuterFractionOfh', 1.0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite', 'positive', '<=', 1}));
            parser.addParameter('OuterRefinementZoneSizeAsMultipleOfh', 2.0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
            parser.parse(varargin{:});
            localOpts = parser.Results;
            samplerArgs = namedArgsToCell(parser.Unmatched);

            [boundaryNodes, ~, boundaryLevelSet] = buildBoundaryState(geometry);
            [x_min, x_max] = kp.nodes.boundingBoxExtents(geometry, true);
            % Choose between fixed-radius and variable-radius sampling, but
            % always thread the geometry-aware refinement options through
            % the same sampler call.
            if isempty(localOpts.RadiusFunction)
                samplerInput = radius;
                samplerArgs = [samplerArgs, {'MinRadius', radius}];
            else
                samplerInput = localOpts.RadiusFunction;
                samplerArgs = [samplerArgs, {'MinRadius', radius}];
            end

            if localOpts.DoOuterRefinement
                zoneSize = localOpts.OuterRefinementZoneSizeAsMultipleOfh * radius;
                samplerArgs = [samplerArgs, ...
                    {'BoundaryPoints', boundaryNodes, ...
                     'BoundaryRefinementFraction', localOpts.OuterFractionOfh, ...
                     'BoundaryDistance', zoneSize}];
            end

            [Xraw, info] = kp.nodes.generatePoissonNodesInBox(samplerInput, x_min, x_max, samplerArgs{:});
            obj.Xi = Xraw;
            obj.Xb = zeros(0, size(Xraw, 2));
            obj.Xg = zeros(0, size(Xraw, 2));
            obj.Nrmls = zeros(0, size(Xraw, 2));
            obj.Xi_orig = Xraw;
            obj.Xi_pds_raw = Xraw;
            obj.s_dim = size(Xraw, 2);
            obj.last_info = info;

            % KernelPack-style clipping keeps a clearance band relative to
            % the smallest active spacing near the boundary.
            clearance = localOpts.OuterFractionOfh * radius;
            obj.clipToGeometry(geometry, 'Keep', 'inside', 'BoundaryClearance', clearance);
            obj.last_info.min_active_radius = clearance;
            obj.last_info.outer_refinement = struct( ...
                'enabled', logical(localOpts.DoOuterRefinement), ...
                'refinement_fraction', localOpts.OuterFractionOfh, ...
                'zone_size', localOpts.OuterRefinementZoneSizeAsMultipleOfh * radius, ...
                'boundary_distance', localOpts.OuterRefinementZoneSizeAsMultipleOfh * radius);
            obj.last_info.boundary_node_count = size(boundaryNodes, 1);
            obj.last_info.boundary_level_set_built = ~isempty(boundaryLevelSet);
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
            % This is the high-level path: generate interior nodes, pull
            % boundary nodes from the geometry object, then append ghost
            % nodes in the outward normal direction.
            obj.generateInteriorNodesFromGeometry(geometry, radius, varargin{:});
            [boundaryNodes, boundaryNormals, boundaryLevelSet] = buildBoundaryState(geometry);
            ghostNodes = boundaryNodes + 0.5 * radius * boundaryNormals;

            descriptor = kp.domain.DomainDescriptor();
            descriptor.setNodes(obj.Xi, boundaryNodes, ghostNodes);
            descriptor.setNormals(boundaryNormals);
            descriptor.setSepRad(radius);
            descriptor.setOuterLevelSet(boundaryLevelSet);
            descriptor.setBoundaryLevelSets({boundaryLevelSet});
            descriptor.buildStructs();

            obj.Xb = boundaryNodes;
            obj.Xg = ghostNodes;
            obj.Nrmls = boundaryNormals;
            obj.descriptor = descriptor;
            obj.last_info.boundary_node_count = size(boundaryNodes, 1);
            obj.last_info.ghost_node_count = size(ghostNodes, 1);
            obj.last_info.min_active_radius = obj.last_info.min_active_radius;
            obj.last_info.sep_radius = radius;
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
% Accept either EmbeddedSurface-style or piecewise geometry-style
% interfaces and normalize them to one boundary cloud contract.
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
