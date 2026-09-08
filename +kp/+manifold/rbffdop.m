classdef rbffdop
    %RBFFDOP Tangent-plane RBF-FD operator settings for manifold PDEs.

    properties
        ell (1,1) double = 0
        polyM (1,1) double = 0
        stencilSize (1,1) double = 0
        rbfexp (1,1) double = 0
        rbf function_handle = @(ep, r) r
        drbfor function_handle = @(ep, r) r
        d2rbf function_handle = @(ep, r) r
        s_dim (1,1) double = 0
        theta (1,1) double = 0
        xi (1,1) double = 0
        hyppow (1,1) double = 0
    end

    methods
        function obj = rbffdop(s_dim, xi, theta, flag)
            %RBFFDOP Match the source tangent-plane operator-selection rules.
            if nargin == 0
                return;
            end

            if flag == 0
                obj.xi = xi;
                obj.theta = theta;
                obj.s_dim = s_dim;
                obj.ell = obj.xi + obj.theta - 1;
                if mod(obj.ell, 2) == 0
                    obj.rbfexp = obj.ell - 1;
                else
                    obj.rbfexp = obj.ell;
                end
                obj.rbfexp = max(obj.rbfexp, 5);
                obj.rbfexp = min(obj.rbfexp, 11);
                obj.rbf = @(~, r) (r + eps).^obj.rbfexp;
                obj.drbfor = @(~, r) obj.rbfexp * (r + eps).^(obj.rbfexp - 2);
                obj.d2rbf = @(~, r) obj.rbfexp * (obj.rbfexp - 1) * (r + eps).^(obj.rbfexp - 2);
                obj.polyM = nchoosek(obj.ell + obj.s_dim, obj.s_dim);
                obj.stencilSize = 2 * obj.polyM + 1;
            else
                obj.s_dim = s_dim;
                obj.xi = xi;
                obj.ell = obj.xi + 1;
                obj.polyM = nchoosek(obj.ell + obj.s_dim, obj.s_dim);
                obj.stencilSize = 2 * obj.polyM + 1;
                obj.hyppow = floor(1.5 * log(obj.stencilSize));
                obj.rbfexp = 2 * obj.hyppow + 1;
                obj.ell = obj.hyppow;
                obj.polyM = nchoosek(obj.ell + obj.s_dim, obj.s_dim);
                obj.stencilSize = 2 * obj.polyM + 1;
            end
        end
    end
end
