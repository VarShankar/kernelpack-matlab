classdef MeanCurvatureFlowSolver < handle
    %MEANCURVATUREFLOWSOLVER Tangent-plane RBF-FD mean-curvature flow.
    %   Evolves a point-cloud embedding by X_t = Delta_Gamma X. The
    %   Laplace-Beltrami operator is rebuilt on the current surface at each
    %   step. Normals may be supplied analytically through NormalProvider or
    %   estimated from the point cloud by local PCA.

    properties
        X double = zeros(0, 3)
        Normals double = zeros(0, 3)
        xi (1,1) double = 4
        dt (1,1) double = 1.0e-3
        Mode (1,1) string = "semiimplicit"
        ProjectNormal (1,1) logical = true
        NormalProvider = []
        NormalStencilSize (1,1) double = 20
        CompletedSteps (1,1) double = 0
        CurrentTime (1,1) double = 0
    end

    methods
        function obj = MeanCurvatureFlowSolver(varargin)
            if nargin == 0
                return;
            end

            parser = inputParser();
            parser.addParameter('Mode', obj.Mode);
            parser.addParameter('ProjectNormal', obj.ProjectNormal);
            parser.addParameter('NormalProvider', obj.NormalProvider);
            parser.addParameter('NormalStencilSize', obj.NormalStencilSize);
            parser.parse(varargin{:});

            obj.Mode = string(parser.Results.Mode);
            obj.ProjectNormal = parser.Results.ProjectNormal;
            obj.NormalProvider = parser.Results.NormalProvider;
            obj.NormalStencilSize = parser.Results.NormalStencilSize;
        end

        function init(obj, X0, normals0, xi, dt)
            validateNodes(X0);
            if nargin < 3 || isempty(normals0)
                normals0 = obj.estimateNormals(X0, zeros(0, 3));
            end
            if nargin >= 4 && ~isempty(xi)
                obj.xi = xi;
            end
            if nargin >= 5 && ~isempty(dt)
                obj.dt = dt;
            end

            obj.X = X0;
            obj.Normals = kp.geometry.normalizeRows(normals0);
            obj.CompletedSteps = 0;
            obj.CurrentTime = 0;
        end

        function setStepSize(obj, dt)
            obj.dt = dt;
        end

        function X = nodes(obj)
            X = obj.X;
        end

        function nr = normals(obj)
            nr = obj.Normals;
        end

        function V = meanCurvatureVelocity(obj)
            [L, ~] = obj.assembleLaplacian(obj.X, obj.Normals);
            V = L * obj.X;
            if obj.ProjectNormal
                V = sum(V .* obj.Normals, 2) .* obj.Normals;
            end
        end

        function Xnext = step(obj)
            obj.refreshNormals();
            [L, I] = obj.assembleLaplacian(obj.X, obj.Normals);

            switch lower(obj.Mode)
                case "explicit"
                    V = L * obj.X;
                    if obj.ProjectNormal
                        V = sum(V .* obj.Normals, 2) .* obj.Normals;
                    end
                    Xnext = obj.X + obj.dt * V;

                case "semiimplicit"
                    Xnext = (I - obj.dt * L) \ obj.X;
                    if obj.ProjectNormal
                        rawStep = Xnext - obj.X;
                        normalStep = sum(rawStep .* obj.Normals, 2) .* obj.Normals;
                        Xnext = obj.X + normalStep;
                    end

                otherwise
                    error('kp:manifold:BadMCFMode', ...
                        'Unknown mean-curvature-flow mode "%s".', obj.Mode);
            end

            obj.X = Xnext;
            obj.CompletedSteps = obj.CompletedSteps + 1;
            obj.CurrentTime = obj.CurrentTime + obj.dt;
            obj.refreshNormals();
        end
    end

    methods (Access = private)
        function refreshNormals(obj)
            if ~isempty(obj.NormalProvider)
                obj.Normals = kp.geometry.normalizeRows(obj.NormalProvider( ...
                    obj.X, obj.CurrentTime, obj.CompletedSteps));
            else
                obj.Normals = obj.estimateNormals(obj.X, obj.Normals);
            end
        end

        function [L, I] = assembleLaplacian(obj, X, normals)
            op = kp.manifold.rbffdop(2, obj.xi, 2, 0);
            tree = KDTreeSearcher(X);
            [L, ~, ~, ~] = kp.manifold.FormSurfaceDiffOpsTP( ...
                X, op.rbf, op.drbfor, op.d2rbf, normals, tree, ...
                op.stencilSize, op.ell);
            I = speye(size(X, 1), size(X, 1));
        end

        function normals = estimateNormals(obj, X, previousNormals)
            validateNodes(X);
            n = size(X, 1);
            k = min(max(obj.NormalStencilSize, 4), n);
            tree = KDTreeSearcher(X);
            idx = knnsearch(tree, X, 'K', k);
            normals = zeros(n, 3);

            for i = 1:n
                cloud = X(idx(i, :), :);
                cloud = cloud - mean(cloud, 1);
                [~, ~, V] = svd(cloud, 0);
                normals(i, :) = V(:, end).';
            end

            if isempty(previousNormals)
                center = mean(X, 1);
                radial = X - center;
                flip = sum(normals .* radial, 2) < 0;
            else
                flip = sum(normals .* previousNormals, 2) < 0;
            end
            normals(flip, :) = -normals(flip, :);
            normals = kp.geometry.normalizeRows(normals);
        end
    end
end

function validateNodes(X)
if size(X, 2) ~= 3
    error('kp:manifold:BadMCFNodes', ...
        'MeanCurvatureFlowSolver expects an N-by-3 node array.');
end
end
