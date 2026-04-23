classdef DomainNodeGenerator < handle
    %DOMAINNODEGENERATOR KernelPack-shaped box Poisson node generator.

    properties (SetAccess = private)
        Xi double = zeros(0, 0)
        Xi_orig double = zeros(0, 0)
        Xi_pds_raw double = zeros(0, 0)
        s_dim (1,1) double = 0
        last_info struct = struct()
    end

    methods
        function generatePoissonNodes(obj, radius, x_min, x_max, varargin)
            [X, info] = kp.nodes.generatePoissonNodesInBox(radius, x_min, x_max, varargin{:});
            obj.Xi = X;
            obj.Xi_orig = X;
            obj.Xi_pds_raw = X;
            obj.s_dim = size(X, 2);
            obj.last_info = info;
        end

        function generateInteriorNodesFromGeometry(obj, geometry, radius, varargin)
            [x_min, x_max] = kp.nodes.boundingBoxExtents(geometry, true);
            obj.generatePoissonNodes(radius, x_min, x_max, varargin{:});
            obj.clipToGeometry(geometry, 'Keep', 'inside');
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

        function out = getInteriorNodes(obj)
            out = obj.Xi;
        end

        function out = getRawPoissonInteriorNodes(obj)
            out = obj.Xi_pds_raw;
        end
    end
end
