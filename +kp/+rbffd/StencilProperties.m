classdef StencilProperties
    %STENCILPROPERTIES KernelPack-style stencil metadata.

    properties
        n (1,1) double = 0
        dim (1,1) double = 0
        ell (1,1) double = 0
        spline_degree (1,1) double = 3
        npoly (1,1) double = 0
        width (1,1) double = 1
        treeMode = "all"
        pointSet = "interior_boundary"
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
            parser.addParameter('treeMode', obj.treeMode, @isModeSpec);
            parser.addParameter('pointSet', obj.pointSet, @isModeSpec);
            parser.parse(varargin{:});

            fields = fieldnames(parser.Results);
            for k = 1:numel(fields)
                obj.(fields{k}) = parser.Results.(fields{k});
            end
            obj.treeMode = kp.rbffd.StencilProperties.normalizeTreeMode(obj.treeMode);
            obj.pointSet = kp.rbffd.StencilProperties.normalizePointSet(obj.pointSet);
            if obj.npoly == 0 && obj.dim > 0
                obj.npoly = size(kp.poly.total_degree_indices(obj.dim, obj.ell), 1);
            end
        end
    end

    methods (Static)
        function mode = normalizeTreeMode(mode)
            if isnumeric(mode)
                switch mode
                    case 0
                        mode = "all";
                    case 1
                        mode = "interior_boundary";
                    case 2
                        mode = "boundary";
                    otherwise
                        error('kp:rbffd:BadTreeMode', 'Unknown numeric treeMode.');
                end
                return;
            end

            mode = lower(string(mode));
            switch mode
                case {"all", "all_nodes", "full"}
                    mode = "all";
                case {"interior_boundary", "interior+boundary", "int_bdry", "intboundary", "owned"}
                    mode = "interior_boundary";
                case {"boundary", "bdry", "boundary_only"}
                    mode = "boundary";
                otherwise
                    error('kp:rbffd:BadTreeMode', 'Unknown treeMode "%s".', mode);
            end
        end

        function mode = normalizePointSet(mode)
            if isnumeric(mode)
                switch mode
                    case 0
                        mode = "all";
                    case 1
                        mode = "interior_boundary";
                    case 2
                        mode = "boundary";
                    otherwise
                        error('kp:rbffd:BadPointSet', 'Unknown numeric pointSet.');
                end
                return;
            end

            mode = lower(string(mode));
            switch mode
                case {"all", "all_nodes", "full"}
                    mode = "all";
                case {"interior_boundary", "interior+boundary", "int_bdry", "intboundary", "interior and boundary"}
                    mode = "interior_boundary";
                case {"boundary", "bdry", "boundary_only"}
                    mode = "boundary";
                otherwise
                    error('kp:rbffd:BadPointSet', 'Unknown pointSet "%s".', mode);
            end
        end
    end
end

function tf = isModeSpec(x)
    tf = isnumeric(x) || ischar(x) || isStringScalar(x);
end
