classdef StencilProperties
    %STENCILPROPERTIES KernelPack-style stencil metadata.

    properties
        n (1,1) double = 0
        dim (1,1) double = 0
        ell (1,1) double = 0
        spline_degree (1,1) double = 3
        npoly (1,1) double = 0
        width (1,1) double = 1
        treeMode (1,1) double = 0
        pointSet (1,1) double = 1
    end

    methods
        function obj = StencilProperties(varargin)
            if nargin == 0
                return;
            end
            parser = inputParser();
            parser.addParameter('n', obj.n);
            parser.addParameter('dim', obj.dim);
            parser.addParameter('ell', obj.ell);
            parser.addParameter('spline_degree', obj.spline_degree);
            parser.addParameter('npoly', obj.npoly);
            parser.addParameter('width', obj.width);
            parser.addParameter('treeMode', obj.treeMode);
            parser.addParameter('pointSet', obj.pointSet);
            parser.parse(varargin{:});

            fields = fieldnames(parser.Results);
            for k = 1:numel(fields)
                obj.(fields{k}) = parser.Results.(fields{k});
            end
            if obj.npoly == 0 && obj.dim > 0
                obj.npoly = size(kp.poly.total_degree_indices(obj.dim, obj.ell), 1);
            end
        end
    end
end
