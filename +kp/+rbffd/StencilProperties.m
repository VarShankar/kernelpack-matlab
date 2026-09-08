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
        function obj = fromAccuracy(varargin)
            parser = inputParser();
            parser.addParameter('Operator', 'lap', @(x) ischar(x) || isStringScalar(x));
            parser.addParameter('ConvergenceOrder', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
            parser.addParameter('DiffOpOrder', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0));
            parser.addParameter('Dimension', [], @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
            parser.addParameter('Approximation', 'rbf', @(x) ischar(x) || isStringScalar(x));
            parser.addParameter('StencilFactor', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 1));
            parser.addParameter('SplineDegree', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
            parser.addParameter('treeMode', "all", @isModeSpec);
            parser.addParameter('pointSet', "interior_boundary", @isModeSpec);
            parser.parse(varargin{:});
            opts = parser.Results;

            q = opts.DiffOpOrder;
            if isempty(q)
                q = kp.rbffd.StencilProperties.defaultDiffOrder(opts.Operator);
            end

            p = opts.ConvergenceOrder;
            if isempty(p)
                p = 2;
            end

            ell = p + q - 1;
            ell = max(ell, 0);
            dim = opts.Dimension;
            npoly = size(kp.poly.total_degree_indices(dim, ell), 1);

            approx = lower(string(opts.Approximation));
            if isempty(opts.StencilFactor)
                switch approx
                    case {"rbf", "rbf-fd", "rbffd"}
                        stencilFactor = 2.0;
                    case {"wls", "weighted_least_squares", "weightedleastsquares"}
                        stencilFactor = 1.5;
                    otherwise
                        error('kp:rbffd:BadApproximation', 'Unknown approximation "%s".', approx);
                end
            else
                stencilFactor = opts.StencilFactor;
            end

            if isempty(opts.SplineDegree)
                splineDegree = max(2 * q + 1, 3);
                if mod(splineDegree, 2) == 0
                    splineDegree = splineDegree + 1;
                end
            else
                splineDegree = opts.SplineDegree;
            end

            n = max(npoly + 1, ceil(stencilFactor * npoly));

            obj = kp.rbffd.StencilProperties( ...
                'n', n, ...
                'dim', dim, ...
                'ell', ell, ...
                'spline_degree', splineDegree, ...
                'npoly', npoly, ...
                'treeMode', opts.treeMode, ...
                'pointSet', opts.pointSet);
        end

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

        function q = defaultDiffOrder(op)
            op = lower(string(op));
            switch op
                case {"interp", "interpolation", "identity"}
                    q = 0;
                case {"grad", "gradient", "dx", "dy", "dz"}
                    q = 1;
                case {"lap", "laplacian", "bc", "boundary"}
                    q = 2;
                otherwise
                    error('kp:rbffd:BadOperator', 'Unknown operator "%s". Supply DiffOpOrder explicitly.', op);
            end
        end
    end
end

function tf = isModeSpec(x)
    tf = isnumeric(x) || ischar(x) || isStringScalar(x);
end
