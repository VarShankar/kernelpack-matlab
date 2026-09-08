classdef PiecewiseSmoothEmbeddedSurface < handle
    %PIECEWISESMOOTHEMBEDDEDSURFACE KernelPack-like piecewise boundary object.

    properties (SetAccess = private)
        segments cell = {}
        Xb double = zeros(0, 0)
        Xb_uniform double = zeros(0, 0)
        Nrmls double = zeros(0, 0)
        Nrmls_uniform double = zeros(0, 0)
        Nrmls_raw double = zeros(0, 0)
        Nrmls_uniform_raw double = zeros(0, 0)
        Nrmls_smoothed double = zeros(0, 0)
        Nrmls_uniform_smoothed double = zeros(0, 0)
        thick_copy double = zeros(0, 0)
        N (1,1) double = 0
        segment_map (:,1) double = zeros(0, 1)
        corner_flags (:,1) double = zeros(0, 1)
        cbox struct = struct('p', [], 'V', [], 'D', [])
        ubox struct = struct('p', [], 'V', [], 'D', [])
        ebox struct = struct('p', [], 'V', [], 'D', [])
        bdryTree struct = struct('Points', [])
        uniform_bdryTree struct = struct('Points', [])
        thick_tree struct = struct('Points', [])
        levelSet kp.geometry.RBFLevelSet = kp.geometry.RBFLevelSet()
    end

    methods
        function generatePiecewiseSmoothSurfaceBySegment(obj, bdry_segments, flip_normal, radius, method, supersample_fac, mode, smooth_normals, smooth_neighborhood)
            if nargin < 5 || isempty(method), method = 1; end
            if nargin < 6 || isempty(supersample_fac), supersample_fac = 2; end
            if nargin < 7 || isempty(mode), mode = 2; end
            if nargin < 8 || isempty(smooth_normals), smooth_normals = false; end
            if nargin < 9 || isempty(smooth_neighborhood), smooth_neighborhood = 0; end

            nSeg = numel(bdry_segments);
            obj.segments = cell(size(bdry_segments));
            obj.segment_map = zeros(0, 1);
            obj.corner_flags = zeros(0, 1);

            % Each segment or patch gets its own EmbeddedSurface model
            % first; the global piecewise object is assembled afterward.
            dim = size(bdry_segments{1}, 2);
            if ~isscalar(radius)
                radius = min(radius(:));
            end

            localSegments = cell(size(bdry_segments));
            if obj.shouldUseParfor(nSeg)
                parfor k = 1:nSeg
                    localSegments{k} = kp.geometry.PiecewiseSmoothEmbeddedSurface.buildSegmentModel( ...
                        bdry_segments{k}, flip_normal(k), dim, radius, method, supersample_fac, mode);
                end
            else
                for k = 1:nSeg
                    localSegments{k} = kp.geometry.PiecewiseSmoothEmbeddedSurface.buildSegmentModel( ...
                        bdry_segments{k}, flip_normal(k), dim, radius, method, supersample_fac, mode);
                end
            end
            obj.segments = localSegments;

            % Concatenate the per-segment sample clouds into the global
            % piecewise boundary representation.
            counts = zeros(numel(obj.segments), 1);
            uniformCounts = zeros(numel(obj.segments), 1);
            for k = 1:numel(obj.segments)
                counts(k) = size(obj.segments{k}.getSampleSites(), 1);
                uniformCounts(k) = size(obj.segments{k}.getUniformSampleSites(), 1);
            end
            totalCount = sum(counts);
            totalUniformCount = sum(uniformCounts);
            Xb_s = zeros(totalCount, dim);
            Nrmls_s = zeros(totalCount, dim);
            Xb_su = zeros(totalUniformCount, dim);
            Nrmls_su = zeros(totalUniformCount, dim);
            obj.segment_map = zeros(totalCount, 1);
            obj.corner_flags = zeros(totalCount, 1);
            uniformSegmentMap = zeros(totalUniformCount, 1);
            uniformCornerFlags = zeros(totalUniformCount, 1);
            row0 = 1;
            rowu0 = 1;
            for k = 1:numel(obj.segments)
                seg = obj.segments{k};
                pts = seg.getSampleSites();
                nrmls = seg.getNrmls();
                m = size(pts, 1);
                rows = row0:(row0 + m - 1);
                Xb_s(rows, :) = pts;
                Nrmls_s(rows, :) = nrmls;
                obj.segment_map(rows) = k;
                flags = zeros(m, 1);
                flags(1) = 1;
                flags(end) = 1;
                obj.corner_flags(rows) = flags;
                row0 = row0 + m;

                ptsu = seg.getUniformSampleSites();
                nrmlsu = seg.getUniformNrmls();
                mu = size(ptsu, 1);
                rowsu = rowu0:(rowu0 + mu - 1);
                Xb_su(rowsu, :) = ptsu;
                Nrmls_su(rowsu, :) = nrmlsu;
                uniformSegmentMap(rowsu) = k;
                flagsu = zeros(mu, 1);
                flagsu(1) = 1;
                flagsu(end) = 1;
                uniformCornerFlags(rowsu) = flagsu;
                rowu0 = rowu0 + mu;
            end

            % First remove corner duplicates in the same spirit as the C++
            % code, then enforce a true global minimum spacing so distinct
            % corner-adjacent samples do not survive too close together.
            [obj.Xb, obj.Nrmls, obj.segment_map, obj.corner_flags] = ...
                obj.deduplicateBoundary(Xb_s, Nrmls_s, obj.segment_map, obj.corner_flags, 0.2 * radius);
            [obj.Xb, obj.Nrmls, obj.segment_map, obj.corner_flags] = ...
                obj.enforceMinimumSpacing(obj.Xb, obj.Nrmls, obj.segment_map, obj.corner_flags, radius);
            [obj.Xb_uniform, obj.Nrmls_uniform, uniformSegmentMap, uniformCornerFlags] = ...
                obj.deduplicateBoundary(Xb_su, Nrmls_su, uniformSegmentMap, uniformCornerFlags, 0.2 * radius);
            [obj.Xb_uniform, obj.Nrmls_uniform, ~, ~] = ...
                obj.enforceMinimumSpacing(obj.Xb_uniform, obj.Nrmls_uniform, uniformSegmentMap, uniformCornerFlags, radius);

            obj.segment_map = obj.segment_map(:);
            obj.corner_flags = obj.corner_flags(:);
            obj.Nrmls_raw = obj.Nrmls;
            obj.Nrmls_uniform_raw = obj.Nrmls_uniform;
            obj.Nrmls_smoothed = obj.Nrmls_raw;
            obj.Nrmls_uniform_smoothed = obj.Nrmls_uniform_raw;
            if smooth_normals && smooth_neighborhood > 0
                % Optional normal smoothing is kept as a postprocess so the
                % raw normals remain available for inspection.
                obj.Nrmls_smoothed = obj.smoothNormals(obj.Xb, obj.Nrmls_raw, smooth_neighborhood);
                obj.Nrmls_uniform_smoothed = obj.smoothNormals(obj.Xb_uniform, obj.Nrmls_uniform_raw, smooth_neighborhood);
                obj.Nrmls = obj.Nrmls_smoothed;
                obj.Nrmls_uniform = obj.Nrmls_uniform_smoothed;
            else
                obj.Nrmls = obj.Nrmls_raw;
                obj.Nrmls_uniform = obj.Nrmls_uniform_raw;
            end

            duplicateCornerMask = false(size(obj.corner_flags));
            duplicateCornerMask(obj.corner_flags ~= 0) = true;
            obj.corner_flags = double(duplicateCornerMask);
            obj.N = size(obj.Xb, 1);

            obj.buildTree();
            obj.buildUniformTree();
        end

        function buildLevelSet(obj)
            % The piecewise level set is built from the assembled uniform
            % boundary cloud rather than from the per-segment data sites.
            obj.levelSet = kp.geometry.RBFLevelSet();
            obj.levelSet.BuildLevelSetFromCFI(obj.Xb_uniform, obj.Nrmls_uniform);
        end

        function buildTree(obj)
            obj.bdryTree = struct('Points', obj.Xb);
        end

        function buildUniformTree(obj)
            obj.uniform_bdryTree = struct('Points', obj.Xb_uniform);
        end

        function buildThickTree(obj)
            obj.thick_tree = struct('Points', obj.thick_copy);
        end

        function computeBoundingBox(obj)
            obj.cbox = kp.geometry.pcaOrientedBoundingBox(obj.Xb);
        end

        function computeUniformBoundingBox(obj)
            obj.ubox = kp.geometry.pcaOrientedBoundingBox(obj.Xb_uniform);
        end

        function computeExtendedBoundingBox(obj, h)
            obj.thick_copy = obj.Xb_uniform + h * obj.Nrmls_uniform;
            obj.buildThickTree();
            obj.ebox = kp.geometry.pcaOrientedBoundingBox(obj.thick_copy);
        end

        function out = getSampleSites(obj), out = obj.Xb; end
        function out = getUniformSampleSites(obj), out = obj.Xb_uniform; end
        function out = getNrmls(obj), out = obj.Nrmls; end
        function out = getUniformNrmls(obj), out = obj.Nrmls_uniform; end
        function out = getCornerFlags(obj), out = obj.corner_flags; end
        function out = getKDTree(obj), out = obj.bdryTree; end
        function out = getUniformKDTree(obj), out = obj.uniform_bdryTree; end
        function out = getThickTree(obj), out = obj.thick_tree; end
        function out = getLevelSet(obj), out = obj.levelSet; end
        function out = getBdryNodes(obj), out = obj.Xb; end
        function out = getUniformBdryNodes(obj), out = obj.Xb_uniform; end
        function out = getBdryNrmls(obj), out = obj.Nrmls; end
        function out = getUniformBdryNrmls(obj), out = obj.Nrmls_uniform; end
        function out = getBoundingBox(obj), out = obj.cbox.p; end
        function out = getUniformBoundingBox(obj), out = obj.ubox.p; end
        function out = getExtendedBoundingBox(obj), out = obj.ebox.p; end
    end

    methods (Access = private)
        function tf = shouldUseParfor(~, nSeg)
            tf = false;
            if nSeg < 2
                return;
            end
            if ~license('test', 'Distrib_Computing_Toolbox')
                return;
            end
            tf = true;
        end

        function [Xout, Nout, segOut, cornerOut] = deduplicateBoundary(~, X, N, segMap, cornerFlags, tol)
            n = size(X, 1);
            if n == 0
                Xout = X;
                Nout = N;
                segOut = segMap;
                cornerOut = cornerFlags;
                return;
            end

            searcher = buildSearcher(X);
            visited = false(n, 1);
            Xkeep = zeros(n, size(X, 2));
            Nkeep = zeros(n, size(N, 2));
            segKeep = zeros(n, 1);
            cornerKeep = zeros(n, 1);
            outCount = 0;

            for i = 1:n
                if visited(i)
                    continue;
                end
                % Merge a duplicate cluster into one representative point,
                % while aligning normals before averaging them.
                cluster = rangeQueryIndices(searcher, X, X(i, :), tol);
                visited(cluster) = true;
                outCount = outCount + 1;
                Xkeep(outCount, :) = mean(X(cluster, :), 1);

                nrmls = N(cluster, :);
                nref = nrmls(1, :);
                for j = 2:size(nrmls, 1)
                    if dot(nrmls(j, :), nref) < 0
                        nrmls(j, :) = -nrmls(j, :);
                    end
                end
                Navg = mean(nrmls, 1);
                Nkeep(outCount, :) = Navg ./ max(norm(Navg, 2), eps);
                segKeep(outCount) = segMap(cluster(1));

                isCorner = any(cornerFlags(cluster) ~= 0);
                if numel(unique(segMap(cluster))) > 1
                    isCorner = true;
                end
                if size(nrmls, 1) > 1
                    angles = nrmls * nrmls.';
                    angles = angles(triu(true(size(angles)), 1));
                    if any(angles < cosd(35))
                        isCorner = true;
                    end
                end
                cornerKeep(outCount) = double(isCorner);
            end

            Xout = Xkeep(1:outCount, :);
            Nout = kp.geometry.normalizeRows(Nkeep(1:outCount, :));
            segOut = segKeep(1:outCount);
            cornerOut = cornerKeep(1:outCount);
        end

        function [Xout, Nout, segOut, cornerOut] = enforceMinimumSpacing(~, X, N, segMap, cornerFlags, radius)
            n = size(X, 1);
            if n == 0
                Xout = X;
                Nout = N;
                segOut = segMap;
                cornerOut = cornerFlags;
                return;
            end

            % After duplicate cleanup, do a conservative spacing pass on
            % the assembled cloud. Corner-tagged points get priority.
            searcher = buildSearcher(X);
            keep = true(n, 1);
            cornerMask = cornerFlags(:) ~= 0;
            order = [(find(cornerMask)).' (find(~cornerMask)).'];
            spacingTol = radius * (1 - 1e-12);

            for idx = order
                if ~keep(idx)
                    continue;
                end
                nbrs = rangeQueryIndices(searcher, X, X(idx, :), spacingTol);
                nbrs = nbrs(nbrs ~= idx);
                if isempty(nbrs)
                    continue;
                end
                for jj = 1:numel(nbrs)
                    j = nbrs(jj);
                    if ~keep(j)
                        continue;
                    end
                    if cornerMask(j) && ~cornerMask(idx)
                        keep(idx) = false;
                        break;
                    end
                    keep(j) = false;
                end
            end

            Xout = X(keep, :);
            Nout = N(keep, :);
            segOut = segMap(keep);
            cornerOut = cornerFlags(keep);
        end

        function NrmlsOut = smoothNormals(~, X, NrmlsIn, neighborhood)
            % Smooth by averaging a small nearest-neighbor normal cloud,
            % flipping signs first so opposite orientations do not cancel.
            n = size(X, 1);
            NrmlsOut = NrmlsIn;
            if n == 0 || neighborhood <= 1
                return;
            end
            searcher = buildSearcher(X);
            for i = 1:n
                take = knnQueryIndices(searcher, X, X(i, :), min(neighborhood, n));
                nri = NrmlsIn(take, :);
                nref = nri(1, :);
                for j = 2:size(nri, 1)
                    if dot(nri(j, :), nref) < 0
                        nri(j, :) = -nri(j, :);
                    end
                end
                navg = mean(nri, 1);
                NrmlsOut(i, :) = navg ./ max(norm(navg, 2), eps);
            end
        end
    end

    methods (Static, Access = private)
        function boundary = buildSegmentModel(bdryPts, flipNormal, dim, radius, method, supersample_fac, mode)
            boundary = kp.geometry.EmbeddedSurface();
            boundary.setDataSites(bdryPts);
            boundary.setSampleSites(bdryPts);
            boundary.computeBoundingBox();
            box = boundary.getBoundingBox();
            if dim == 2
                w = abs(max(box(:, 1)) - min(box(:, 1)));
                h = abs(max(box(:, 2)) - min(box(:, 2)));
                bdry_size = 2 * (w + h);
                seg_N = max(8, round(bdry_size / radius));
                boundary.buildGeometricModelPS(dim, radius, size(bdryPts, 1), seg_N, ...
                    'eval + eval_first_ders', method, supersample_fac, mode, 1);
            elseif dim == 3
                w = abs(max(box(:, 1)) - min(box(:, 1)));
                h = abs(max(box(:, 2)) - min(box(:, 2)));
                l = abs(max(box(:, 3)) - min(box(:, 3)));
                bdry_size = 2 * (w * h + h * l + l * w);
                seg_N = max(16, round(bdry_size / max(radius^2, eps)));
                boundary.buildGeometricModelPS(dim, radius, size(bdryPts, 1), seg_N, ...
                    'eval + eval_first_ders', method, supersample_fac, mode, 1);
            else
                error('PiecewiseSmoothEmbeddedSurface:NotImplemented', ...
                    'Unsupported segment dimension.');
            end
            if flipNormal
                boundary.flipNormals();
            end
        end
    end
end

function searcher = buildSearcher(X)
searcher = [];
if isempty(X)
    return;
end
if exist('KDTreeSearcher', 'class') == 8 && exist('rangesearch', 'file') == 2
    searcher = KDTreeSearcher(X);
end
end

function idx = rangeQueryIndices(searcher, X, xq, radius)
if ~isempty(searcher)
    idxCell = rangesearch(searcher, xq, radius);
    idx = idxCell{1}(:);
else
    d = sqrt(sum((X - xq).^2, 2));
    idx = find(d <= radius);
end
end

function idx = knnQueryIndices(searcher, X, xq, k)
if ~isempty(searcher)
    idx = knnsearch(searcher, xq, 'K', k);
    idx = idx(:);
else
    d = sqrt(sum((X - xq).^2, 2));
    [~, order] = sort(d, 'ascend');
    idx = order(1:k);
end
end
