classdef EmbeddedSurface < handle
    %EMBEDDEDSURFACE KernelPack-like geometric description of one boundary object.

    properties
        ubox struct = struct('p', [], 'V', [], 'D', [])
        cbox struct = struct('p', [], 'V', [], 'D', [])
        ebox struct = struct('p', [], 'V', [], 'D', [])
        data_sites double = zeros(0, 0)
        data_site_nrmls double = zeros(0, 0)
        uniform_sample_sites double = zeros(0, 0)
        sample_sites double = zeros(0, 0)
        sample_sites_s double = zeros(0, 0)
        thick_copy double = zeros(0, 0)
        uniform_Nrmls double = zeros(0, 0)
        tangents1 double = zeros(0, 0)
        tangents2 double = zeros(0, 0)
        tangents1_s double = zeros(0, 0)
        tangents2_s double = zeros(0, 0)
        Nrmls double = zeros(0, 0)
        Nrmls_s double = zeros(0, 0)
        geom_model struct = struct()
        Nd (1,1) double = 0
        N (1,1) double = 0
        surf_dim (1,1) double = 0
        sep_rad (1,1) double = NaN
        uniformTree struct = struct('Points', [])
        tallTree struct = struct('Points', [])
        thick_tree struct = struct('Points', [])
        levelSet kp.geometry.RBFLevelSet = kp.geometry.RBFLevelSet()
        tangent_step (1,1) double = 1e-5
    end

    methods
        function setDataSites(obj, dataSites)
            obj.data_sites = dataSites;
            obj.Nd = size(dataSites, 1);
            obj.surf_dim = size(dataSites, 2) - 1;
        end

        function setSampleSites(obj, sampleSites)
            obj.sample_sites = sampleSites;
            obj.N = size(sampleSites, 1);
        end

        function setUniformSampleSites(obj, sampleSites)
            obj.uniform_sample_sites = sampleSites;
        end

        function buildClosedGeometricModelPS(obj, dim, rad, Nb, Ne, eval_mode, method, supersample_fac, mode, chart_mode)
            if nargin < 5 || isempty(Ne), Ne = []; end
            if nargin < 6 || isempty(eval_mode), eval_mode = 'eval + eval_first_ders'; end
            if nargin < 7 || isempty(method), method = 1; end
            if nargin < 8 || isempty(supersample_fac), supersample_fac = 2; end
            if nargin < 9 || isempty(mode), mode = 2; end
            if nargin < 10 || isempty(chart_mode), chart_mode = 1; end
            unused_opts = {eval_mode, mode, chart_mode}; %#ok<NASGU>
            obj.sep_rad = rad;
            obj.Nd = min(Nb, size(obj.data_sites, 1));
            obj.surf_dim = dim - 1;

            % KernelPack builds a smooth parametric model from the raw
            % data sites first, then derives sample sites, tangents, and
            % normals from that model.
            dataSites = obj.data_sites(1:obj.Nd, :);
            degree = 7;
            Ntarget = obj.estimateEvaluationCount(dim, rad, true, Ne);

            if dim == 2
                % Closed curves use a periodic SBF fit on a chord-length
                % parameter. This is the path that drives the smooth 2D
                % boundary sampling used elsewhere in the package.
                t = kp.geometry.chordLengthParam(dataSites, true);
                theta = 2 * pi * t;
                [r, ~] = kp.geometry.periodicChordDistance(theta, theta);
                K = kp.geometry.phsKernel(r, degree);
                reg = 1e-12 * max(1.0, max(abs(K), [], 'all'));
                weights = (K + reg * eye(size(K))) \ dataSites;

                obj.geom_model = struct( ...
                    'type', 'closed-curve-sbf', ...
                    'degree', degree, ...
                    'theta', theta, ...
                    'weights', weights);
                [~, dn] = obj.evalClosedCurveFrame(t);
                obj.data_site_nrmls = dn;

                if method == 1
                    % Oversample the curve, then resample by arc length to
                    % get a cleaner spacing-controlled boundary cloud.
                    Ns = max(round(supersample_fac * 1.5 * Ntarget), Ntarget);
                    ts = linspace(0, 1, Ns + 1).';
                    ts = ts(1:end-1);
                    ptss = obj.evalClosedCurve(ts);
                    [tanS, nrmlS] = obj.evalClosedCurveFrame(ts);
                    obj.sample_sites_s = ptss;
                    obj.tangents1_s = tanS;
                    obj.Nrmls_s = nrmlS;
                    curveLength = obj.closedPolylineLength(ptss);
                    targetSpacing = max(sqrt(2.0) * rad, eps);
                    targetCount = max(2, round(curveLength / targetSpacing));
                    keepInds = kp.geometry.resampleClosedCurveByArcLength(ptss, targetCount);
                    obj.sample_sites = ptss(keepInds, :);
                    obj.tangents1 = tanS(keepInds, :);
                    obj.Nrmls = nrmlS(keepInds, :);
                    obj.uniform_sample_sites = obj.sample_sites;
                    obj.uniform_Nrmls = obj.Nrmls;
                else
                    tEval = linspace(0, 1, Ntarget + 1).';
                    tEval = tEval(1:end-1);
                    obj.sample_sites = obj.evalClosedCurve(tEval);
                    [obj.tangents1, obj.Nrmls] = obj.evalClosedCurveFrame(tEval);
                    obj.uniform_sample_sites = obj.sample_sites;
                    obj.uniform_Nrmls = obj.Nrmls;
                end
            elseif dim == 3
                % Smooth closed surfaces are fit coordinate-wise on the
                % sphere, then sampled from a Fibonacci sphere parameter
                % cloud before optional thinning.
                center = mean(dataSites, 1);
                local = dataSites - center;
                uvw = kp.geometry.cart2sphRows(local);
                uv = uvw(:, 1:2);
                unitCenters = local ./ max(uvw(:, 3), eps);
                [r, ~] = kp.geometry.sphereChordDistance(unitCenters, unitCenters);
                K = kp.geometry.phsKernel(r, degree);
                reg = 1e-12 * max(1.0, max(abs(K), [], 'all'));
                weights = (K + reg * eye(size(K))) \ dataSites;
                obj.geom_model = struct( ...
                    'type', 'closed-surface-sbf', ...
                    'degree', degree, ...
                    'center', center, ...
                    'unitCenters', unitCenters, ...
                    'weights', weights);

                if method == 1
                    Ns = max(round(supersample_fac * Ntarget), Ntarget);
                else
                    Ns = Ntarget;
                end
                xyzParam = kp.geometry.fibonacciSphere(Ns);
                uvEval = kp.geometry.cart2sphRows(xyzParam);
                uvEval = uvEval(:, 1:2);
                ptss = obj.evalClosedSurface(uvEval);
                [tan1S, tan2S, nrmlS] = obj.evalClosedSurfaceFrame(uvEval);
                if method == 1
                    % The 3D path still uses sample thinning after a
                    % supersampled surface evaluation.
                    keep = kp.geometry.weightedSampleEliminationMIS(ptss, rad);
                    if nnz(keep) < max(8, floor(Ntarget / 2))
                        keep = false(size(ptss, 1), 1);
                        keep(round(linspace(1, size(ptss, 1), Ntarget))) = true;
                    end
                    obj.sample_sites_s = ptss;
                    obj.sample_sites = ptss(keep, :);
                    obj.tangents1_s = tan1S;
                    obj.tangents2_s = tan2S;
                    obj.Nrmls_s = nrmlS;
                    obj.tangents1 = tan1S(keep, :);
                    obj.tangents2 = tan2S(keep, :);
                    obj.Nrmls = nrmlS(keep, :);
                else
                    obj.sample_sites = ptss;
                    obj.tangents1 = tan1S;
                    obj.tangents2 = tan2S;
                    obj.Nrmls = nrmlS;
                end
                obj.uniform_sample_sites = obj.sample_sites;
                obj.uniform_Nrmls = obj.Nrmls;
                [~, ~, dn] = obj.evalClosedSurfaceFrame(uv);
                obj.data_site_nrmls = dn;
                obj.orientClosedSurfaceNormalsOutward();
            else
                error('EmbeddedSurface:BadDimension', 'Unsupported dimension.');
            end

            obj.N = size(obj.sample_sites, 1);
            obj.buildTree();
            obj.buildUniformTree();
        end

        function buildGeometricModelPS(obj, dim, rad, Nb, Ne, eval_mode, method, supersample_fac, mode, chart_mode)
            if nargin < 5 || isempty(Ne), Ne = []; end
            if nargin < 6 || isempty(eval_mode), eval_mode = 'eval + eval_first_ders'; end
            if nargin < 7 || isempty(method), method = 1; end
            if nargin < 8 || isempty(supersample_fac), supersample_fac = 2; end
            if nargin < 9 || isempty(mode), mode = 2; end
            if nargin < 10 || isempty(chart_mode), chart_mode = 1; end
            unused_opts = {eval_mode, mode, chart_mode}; %#ok<NASGU>
            obj.sep_rad = rad;
            obj.Nd = min(Nb, size(obj.data_sites, 1));
            obj.surf_dim = dim - 1;

            % Open segments and patches are treated as local chart models
            % rather than periodic closed objects.
            dataSites = obj.data_sites(1:obj.Nd, :);
            degree = 7;
            Ntarget = obj.estimateEvaluationCount(dim, rad, false, Ne);

            if dim == 2
                % Open 2D segments are fit on a 1D parameter line with a
                % polynomial augmentation for stable derivative recovery.
                u = kp.geometry.chordLengthParam(dataSites, false);
                K = kp.geometry.phsKernel(kp.geometry.distanceMatrix(u, u), degree);
                P = [ones(numel(u), 1), u];
                reg = 1e-12 * max(1.0, max(abs(K), [], 'all'));
                A = [K + reg * eye(numel(u)), P; P.', zeros(2, 2)];
                coeffs = A \ [dataSites; zeros(2, 2)];

                obj.geom_model = struct( ...
                    'type', 'open-curve-rbf', ...
                    'degree', degree, ...
                    'u', u, ...
                    'rbfWeights', coeffs(1:numel(u), :), ...
                    'polyCoeffs', coeffs(numel(u) + 1:end, :));
                [~, dn] = obj.evalOpenCurveFrame(u);
                obj.data_site_nrmls = dn;

                if method == 1
                    % The open-curve path still uses geometric thinning on
                    % the supersampled cloud instead of arc-length
                    % redistribution.
                    Ns = max(round(supersample_fac * 1.5 * Ntarget), Ntarget);
                    us = linspace(0, 1, Ns).';
                    ptss = obj.evalOpenCurve(us);
                    [tanS, nrmlS] = obj.evalOpenCurveFrame(us);
                    keep = kp.geometry.weightedSampleEliminationMIS(ptss, rad);
                    keep = obj.ensureNonemptyKeep(keep, size(ptss, 1));
                    obj.sample_sites_s = ptss;
                    obj.sample_sites = ptss(keep, :);
                    obj.tangents1_s = tanS;
                    obj.Nrmls_s = nrmlS;
                    obj.tangents1 = tanS(keep, :);
                    obj.Nrmls = nrmlS(keep, :);
                    obj.uniform_sample_sites = obj.sample_sites;
                    obj.uniform_Nrmls = obj.Nrmls;
                else
                    uEval = linspace(0, 1, Ntarget).';
                    obj.sample_sites = obj.evalOpenCurve(uEval);
                    [obj.tangents1, obj.Nrmls] = obj.evalOpenCurveFrame(uEval);
                    obj.uniform_sample_sites = obj.sample_sites;
                    obj.uniform_Nrmls = obj.Nrmls;
                end
            elseif dim == 3
                % Open 3D patches are parameterized on a best-fit plane and
                % then fit with an RBF patch model in that local chart.
                uv = kp.geometry.buildPlanarParametricNodes2D(dataSites, obj.Nd, 0, 1);
                K = kp.geometry.phsKernel(kp.geometry.distanceMatrix(uv, uv), degree);
                P = [ones(size(uv, 1), 1), uv];
                reg = 1e-12 * max(1.0, max(abs(K), [], 'all'));
                A = [K + reg * eye(size(uv, 1)), P; P.', zeros(3, 3)];
                coeffs = A \ [dataSites; zeros(3, 3)];
                [~, origin, basis] = kp.geometry.projectToBestFitPlane(dataSites);
                obj.geom_model = struct( ...
                    'type', 'surface-patch-rbf', ...
                    'degree', degree, ...
                    'uv', uv, ...
                    'rbfWeights', coeffs(1:size(uv, 1), :), ...
                    'polyCoeffs', coeffs(size(uv, 1) + 1:end, :), ...
                    'origin', origin, ...
                    'basis', basis);

                if method == 1
                    uvEval = kp.geometry.buildPlanarParametricEvalNodes2D(dataSites, max(round(supersample_fac * Ntarget), Ntarget));
                    if size(uvEval, 1) < Ntarget
                        uvEval = uv;
                    end
                else
                    uvEval = kp.geometry.buildPlanarParametricEvalNodes2D(dataSites, Ntarget);
                    if isempty(uvEval)
                        uvEval = uv;
                    end
                end
                ptss = obj.evalOpenSurface(uvEval);
                [tan1S, tan2S, nrmlS] = obj.evalOpenSurfaceFrame(uvEval);
                if method == 1
                    keep = kp.geometry.weightedSampleEliminationMIS(ptss, rad);
                    keep = obj.ensureNonemptyKeep(keep, size(ptss, 1));
                    obj.sample_sites_s = ptss;
                    obj.tangents1_s = tan1S;
                    obj.tangents2_s = tan2S;
                    obj.Nrmls_s = nrmlS;
                    obj.sample_sites = ptss(keep, :);
                    obj.tangents1 = tan1S(keep, :);
                    obj.tangents2 = tan2S(keep, :);
                    obj.Nrmls = nrmlS(keep, :);
                else
                    obj.sample_sites = ptss;
                    obj.tangents1 = tan1S;
                    obj.tangents2 = tan2S;
                    obj.Nrmls = nrmlS;
                end
                obj.uniform_sample_sites = obj.sample_sites;
                obj.uniform_Nrmls = obj.Nrmls;
                [~, ~, dn] = obj.evalOpenSurfaceFrame(uv);
                obj.data_site_nrmls = dn;
            else
                error('EmbeddedSurface:BadDimension', 'Unsupported dimension.');
            end

            obj.N = size(obj.sample_sites, 1);
            obj.buildTree();
            obj.buildUniformTree();
        end

        function buildLevelSetFromGeometricModel(obj, lambda)
            % Build the implicit representation from whichever boundary
            % cloud is currently regarded as the uniform geometric model.
            if nargin < 2 || isempty(lambda)
                pts = obj.uniform_sample_sites;
                nrmls = obj.uniform_Nrmls;
            else
                if strcmp(obj.geom_model.type, 'closed-curve-sbf')
                    pts = obj.evalClosedCurve(lambda);
                    [~, nrmls] = obj.evalClosedCurveFrame(lambda);
                elseif strcmp(obj.geom_model.type, 'closed-surface-sbf')
                    pts = obj.evalClosedSurface(lambda);
                    [~, ~, nrmls] = obj.evalClosedSurfaceFrame(lambda);
                else
                    if obj.surf_dim == 1
                        pts = obj.evalOpenCurve(lambda);
                        [~, nrmls] = obj.evalOpenCurveFrame(lambda);
                    else
                        pts = obj.evalOpenSurface(lambda);
                        [~, ~, nrmls] = obj.evalOpenSurfaceFrame(lambda);
                    end
                end
            end
            obj.levelSet = kp.geometry.RBFLevelSet();
            obj.levelSet.BuildLevelSetFromCFI(pts, nrmls);
        end

        function computeBoundingBox(obj)
            obj.cbox = kp.geometry.pcaOrientedBoundingBox(obj.sample_sites);
        end

        function computeUniformBoundingBox(obj)
            obj.ubox = kp.geometry.pcaOrientedBoundingBox(obj.uniform_sample_sites);
        end

        function computeExtendedBoundingBox(obj, h)
            % The extended box is built from a thickened copy in the normal
            % direction, mirroring the C++ "thick copy" idea.
            obj.thick_copy = obj.uniform_sample_sites + h * obj.uniform_Nrmls;
            obj.thick_tree = struct('Points', obj.thick_copy);
            obj.ebox = kp.geometry.pcaOrientedBoundingBox(obj.thick_copy);
        end

        function buildTree(obj)
            obj.tallTree = struct('Points', obj.sample_sites);
        end

        function buildUniformTree(obj)
            obj.uniformTree = struct('Points', obj.uniform_sample_sites);
        end

        function flipNormals(obj)
            obj.Nrmls = -obj.Nrmls;
            obj.uniform_Nrmls = -obj.uniform_Nrmls;
            obj.Nrmls_s = -obj.Nrmls_s;
            obj.data_site_nrmls = -obj.data_site_nrmls;
        end

        function out = getSampleSites(obj), out = obj.sample_sites; end
        function out = getUniformSampleSites(obj), out = obj.uniform_sample_sites; end
        function out = getNrmls(obj), out = obj.Nrmls; end
        function out = getUniformNrmls(obj), out = obj.uniform_Nrmls; end
        function out = getBoundingBox(obj), out = obj.cbox.p; end
        function out = getUniformBoundingBox(obj), out = obj.ubox.p; end
        function out = getExtendedBoundingBox(obj), out = obj.ebox.p; end
        function out = getKDTree(obj), out = obj.tallTree; end
        function out = getUniformKDTree(obj), out = obj.uniformTree; end
        function out = getThickTree(obj), out = obj.thick_tree; end
        function out = getLevelSet(obj), out = obj.levelSet; end
        function out = getN(obj), out = obj.N; end
    end

    methods (Access = private)
        function keep = ensureNonemptyKeep(~, keep, nPts)
            if ~any(keep) && nPts > 0
                keep = false(nPts, 1);
                keep(1) = true;
            end
        end

        function Ntarget = estimateEvaluationCount(obj, dim, rad, isClosed, NeHint)
            if ~isempty(NeHint)
                Ntarget = max(ceil(double(NeHint)), 1);
                return;
            end

            X = obj.data_sites(1:obj.Nd, :);
            if dim == 2
                if isClosed
                    measure = obj.closedPolylineLength(X);
                else
                    measure = obj.openPolylineLength(X);
                end
                Ntarget = max(ceil(measure / max(rad, eps)), 8);
            elseif dim == 3
                mins = min(X, [], 1);
                maxs = max(X, [], 1);
                ext = max(maxs - mins, eps);
                area = 2 * (ext(1) * ext(2) + ext(1) * ext(3) + ext(2) * ext(3));
                if ~isClosed
                    area = 0.5 * area;
                end
                Ntarget = max(ceil(area / max(rad * rad, eps)), 16);
            else
                error('EmbeddedSurface:BadDimension', 'Unsupported dimension.');
            end
        end

        function L = closedPolylineLength(~, X)
            if size(X, 1) < 2
                L = 0;
                return;
            end
            shifted = X([2:end, 1], :);
            L = sum(sqrt(sum((shifted - X) .^ 2, 2)));
        end

        function L = openPolylineLength(~, X)
            if size(X, 1) < 2
                L = 0;
                return;
            end
            dX = diff(X, 1, 1);
            L = sum(sqrt(sum(dX .^ 2, 2)));
        end

        function orientClosedSurfaceNormalsOutward(obj)
            if isempty(obj.sample_sites) || size(obj.sample_sites, 2) ~= 3
                return;
            end

            center = mean(obj.sample_sites, 1);
            orientation = mean(sum((obj.sample_sites - center) .* obj.Nrmls, 2), 'omitnan');
            if orientation < 0
                obj.Nrmls = -obj.Nrmls;
                obj.uniform_Nrmls = -obj.uniform_Nrmls;
                obj.Nrmls_s = -obj.Nrmls_s;
                obj.data_site_nrmls = -obj.data_site_nrmls;
            end
        end
    end

    methods (Access = private)
        function x = evalClosedCurve(obj, t)
            % Evaluate the periodic closed-curve fit at arbitrary
            % parameter locations in [0,1).
            tw = kp.geometry.wrapPeriodicParameter(t(:));
            theta = 2 * pi * tw;
            [r, ~] = kp.geometry.periodicChordDistance(theta, obj.geom_model.theta);
            K = kp.geometry.phsKernel(r, obj.geom_model.degree);
            x = K * obj.geom_model.weights;
        end

        function [xt, n] = evalClosedCurveFrame(obj, t)
            % Recover the tangent and outward normal from the first
            % derivative of the periodic curve model.
            tw = kp.geometry.wrapPeriodicParameter(t(:));
            theta = 2 * pi * tw;
            [r, delta] = kp.geometry.periodicChordDistance(theta, obj.geom_model.theta);
            degree = obj.geom_model.degree;
            if degree == 1
                dphi = sin(delta) ./ max(r, eps);
                dphi(r == 0) = 0;
            else
                dphi = degree * sin(delta) .* (r .^ (degree - 2));
                dphi(r == 0) = 0;
            end
            xt = (2 * pi) * (dphi * obj.geom_model.weights);
            n = kp.geometry.normalizeRows([xt(:, 2), -xt(:, 1)]);
        end

        function x = evalOpenCurve(obj, u)
            % Evaluate the open-curve RBF chart at arbitrary parameter
            % locations on the unit interval.
            u = min(max(u(:), 0), 1);
            K = kp.geometry.phsKernel(kp.geometry.distanceMatrix(u, obj.geom_model.u), obj.geom_model.degree);
            x = K * obj.geom_model.rbfWeights + [ones(size(u, 1), 1), u] * obj.geom_model.polyCoeffs;
        end

        function [xt, n] = evalOpenCurveFrame(obj, u)
            % Differentiate the open curve model directly and rotate the
            % tangent to get a consistent outward normal.
            u = min(max(u(:), 0), 1);
            du = u - obj.geom_model.u.';
            r = abs(du);
            degree = obj.geom_model.degree;
            if degree == 1
                dphi = sign(du);
            else
                dphi = degree * sign(du) .* (r .^ (degree - 1));
            end
            dphi(r == 0 & degree > 1) = 0;
            xt = dphi * obj.geom_model.rbfWeights + repmat(obj.geom_model.polyCoeffs(2, :), size(u, 1), 1);
            n = kp.geometry.normalizeRows([xt(:, 2), -xt(:, 1)]);
        end

        function x = evalClosedSurface(obj, uv)
            % Evaluate the smooth closed surface from the spherical chart.
            unitQuery = [cos(uv(:, 2)) .* cos(uv(:, 1)), ...
                         cos(uv(:, 2)) .* sin(uv(:, 1)), ...
                         sin(uv(:, 2))];
            [r, ~] = kp.geometry.sphereChordDistance(unitQuery, obj.geom_model.unitCenters);
            K = kp.geometry.phsKernel(r, obj.geom_model.degree);
            x = K * obj.geom_model.weights;
        end

        function [tu, tv, n] = evalClosedSurfaceFrame(obj, uv)
            % Differentiate the spherical-basis surface fit analytically
            % rather than using finite differences in parameter space.
            unitQuery = [cos(uv(:, 2)) .* cos(uv(:, 1)), ...
                         cos(uv(:, 2)) .* sin(uv(:, 1)), ...
                         sin(uv(:, 2))];
            du = [-cos(uv(:, 2)) .* sin(uv(:, 1)), ...
                   cos(uv(:, 2)) .* cos(uv(:, 1)), ...
                   zeros(size(uv, 1), 1)];
            dv = [-sin(uv(:, 2)) .* cos(uv(:, 1)), ...
                  -sin(uv(:, 2)) .* sin(uv(:, 1)), ...
                   cos(uv(:, 2))];
            diff = kp.geometry.RBFLevelSet.differenceTensor(unitQuery, obj.geom_model.unitCenters);
            r = sqrt(sum(diff.^2, 3));
            phiFactor = kp.geometry.RBFLevelSet.radialDerivativeFactor(r, obj.geom_model.degree);
            dotDu = diff(:, :, 1) .* du(:, 1) + diff(:, :, 2) .* du(:, 2) + diff(:, :, 3) .* du(:, 3);
            dotDv = diff(:, :, 1) .* dv(:, 1) + diff(:, :, 2) .* dv(:, 2) + diff(:, :, 3) .* dv(:, 3);
            tu = (phiFactor .* dotDu) * obj.geom_model.weights;
            tv = (phiFactor .* dotDv) * obj.geom_model.weights;
            n = kp.geometry.normalizeRows(cross(tu, tv, 2));
        end

        function x = evalOpenSurface(obj, uv)
            % Evaluate the local open surface patch in its planar chart.
            K = kp.geometry.phsKernel(kp.geometry.distanceMatrix(uv, obj.geom_model.uv), obj.geom_model.degree);
            x = K * obj.geom_model.rbfWeights + [ones(size(uv, 1), 1), uv] * obj.geom_model.polyCoeffs;
        end

        function [tu, tv, n] = evalOpenSurfaceFrame(obj, uv)
            % Differentiate the patch chart with respect to its two local
            % parameters, then take the cross product for normals.
            diffU = uv(:, 1) - obj.geom_model.uv(:, 1).';
            diffV = uv(:, 2) - obj.geom_model.uv(:, 2).';
            r = sqrt(diffU.^2 + diffV.^2);
            degree = obj.geom_model.degree;
            dphi = degree * r .^ max(degree - 2, 0);
            if degree == 1
                dphi = ones(size(r));
            end
            dphi(r == 0 & degree > 1) = 0;
            tu = (dphi .* diffU) * obj.geom_model.rbfWeights + repmat(obj.geom_model.polyCoeffs(2, :), size(uv, 1), 1);
            tv = (dphi .* diffV) * obj.geom_model.rbfWeights + repmat(obj.geom_model.polyCoeffs(3, :), size(uv, 1), 1);
            n = kp.geometry.normalizeRows(cross(tu, tv, 2));
        end
    end
end
