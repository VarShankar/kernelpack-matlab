classdef IBAMRSurfaceTrajectory < handle
    %IBAMRSURFACETRAJECTORY Sphere-homeomorphic IBAMR material trajectory.
    %   The native IB vertices retain fixed material coordinates on S^2.
    %   Cubic Hermite interpolation uses both the exported positions and IB
    %   velocities. A cached spherical-basis interpolation factorization then
    %   evaluates the evolving surface, tangents, normals, and quadrature at
    %   either the native material sites or rearranged marker sites.

    properties (SetAccess = private)
        Folder string
        Times double
        Steps double
        MaterialSites double
        Faces double
        Positions double
        Velocities double
        Diagnostics table
        Degree (1,1) double
        HarmonicDegree (1,1) double
        Xi (1,1) double
        ControlIds double
        ControlInfo struct
    end

    properties (Access = private)
        InterpolationFactor
        ControlSites double
    end

    methods
        function obj = IBAMRSurfaceTrajectory(folder, varargin)
            parser = inputParser();
            parser.addParameter('Degree', 7, @(x) isnumeric(x) && isscalar(x));
            parser.addParameter('Xi', 6, @(x) isnumeric(x) && isscalar(x));
            parser.addParameter('ControlPointCount', NaN, ...
                @(x) isnumeric(x) && isscalar(x));
            parser.parse(varargin{:});

            obj.Folder = string(folder);
            obj.Degree = parser.Results.Degree;
            obj.HarmonicDegree = max(5, round((obj.Degree + 1) / 2));
            obj.Xi = parser.Results.Xi;
            obj.loadFiles();
            obj.buildGeometryFactor(parser.Results.ControlPointCount);
        end

        function material = initialMaterial(obj)
            material = struct('U', obj.MaterialSites);
        end

        function geom = geometry(obj, material, t)
            U = normalizeRows(material.U);
            [nativeX, nativeV] = obj.interpolateNative(t);
            coefficients = obj.solveCoefficients(nativeX(obj.ControlIds, :));
            velocityCoefficients = obj.solveCoefficients(nativeV(obj.ControlIds, :));
            [X, tangentAzimuth, tangentElevation] = obj.evaluateModel(U, coefficients);
            velocity = obj.evaluateValues(U, velocityCoefficients);

            normals = normalizeRows(cross(tangentAzimuth, tangentElevation, 2));
            radial = X - mean(X, 1);
            if median(sum(normals .* radial, 2)) < 0
                normals = -normals;
            end

            sphericalWeights = sphericalLumpedWeights(U);
            uvw = kp.geometry.cart2sphRows(U);
            sphereJacobian = max(cos(uvw(:, 2)), 1.0e-10);
            surfaceJacobian = vecnorm(cross(tangentAzimuth, tangentElevation, 2), 2, 2);
            weights = sphericalWeights .* surfaceJacobian ./ sphereJacobian;

            geom = struct();
            geom.X = X;
            geom.velocity = velocity;
            geom.normals = normals;
            geom.weights = weights;
            geom.area = sum(weights);
            geom.h = sqrt(geom.area / size(X, 1));
            geom.material = struct('U', U, 'X', X);
            % The spherical labels preserve membrane topology when the
            % biconcave surface brings distinct physical sheets close.
            geom.neighborCoordinates = U;
        end

        function material = sampleMaterial(obj, N, t)
            surfaceMap = @(U) obj.positionAt(U, t);
            areaDensity = @(U) obj.areaRatioAt(U, t);
            [U, info] = kp.manifold.sampleSphereSurfaceFPS(N, surfaceMap, ...
                'CandidateFactor', 8, ...
                'RawCandidateFactor', 3, ...
                'AreaDensityFunction', areaDensity, ...
                'UseAreaWeightedCandidates', true);
            material = struct('U', U, 'samplingInfo', info);
        end

        function X = positionAt(obj, U, t)
            nativeX = obj.interpolateNative(t);
            coefficients = obj.solveCoefficients(nativeX(obj.ControlIds, :));
            X = obj.evaluateValues(normalizeRows(U), coefficients);
        end

        function ratio = areaRatioAt(obj, U, t)
            nativeX = obj.interpolateNative(t);
            coefficients = obj.solveCoefficients(nativeX(obj.ControlIds, :));
            U = normalizeRows(U);
            [~, tangentAzimuth, tangentElevation] = obj.evaluateModel(U, coefficients);
            uvw = kp.geometry.cart2sphRows(U);
            ratio = vecnorm(cross(tangentAzimuth, tangentElevation, 2), 2, 2) ./ ...
                max(cos(uvw(:, 2)), 1.0e-10);
        end

        function info = validateVelocityInterpolation(obj)
            if numel(obj.Times) < 3
                info = struct('relativeRMS', NaN, 'maximumRelative', NaN);
                return;
            end
            errors = zeros(numel(obj.Times) - 2, 1);
            relative = zeros(size(errors));
            for k = 2:numel(obj.Times) - 1
                dt = obj.Times(k + 1) - obj.Times(k - 1);
                finiteDifference = (obj.Positions(:, :, k + 1) - ...
                    obj.Positions(:, :, k - 1)) ./ dt;
                residual = finiteDifference - obj.Velocities(:, :, k);
                errors(k - 1) = norm(residual, 'fro');
                relative(k - 1) = errors(k - 1) / ...
                    max(norm(obj.Velocities(:, :, k), 'fro'), eps);
            end
            info = struct('relativeRMS', sqrt(mean(relative .^ 2)), ...
                'maximumRelative', max(relative));
        end

        function targetValues = interpolateMaterialField(obj, sourceMaterial, ...
                sourceValues, targetMaterial)
            %INTERPOLATEMATERIALFIELD Local tangent-plane PHS interpolation.
            sourceSites = normalizeRows(sourceMaterial.U);
            targetSites = normalizeRows(targetMaterial.U);
            if size(sourceValues, 1) ~= size(sourceSites, 1)
                error('kp:manifold:BadIBAMRMaterialField', ...
                    'Source values must have one row per source material site.');
            end

            ell = max(obj.Xi - 1, 1);
            splineDegree = ell;
            if mod(splineDegree, 2) == 0
                splineDegree = splineDegree - 1;
            end
            splineDegree = max(splineDegree, 3);
            polynomialCount = size(kp.poly.total_degree_indices(2, ell), 1);
            stencilSize = min(size(sourceSites, 1), 2 * polynomialCount + 1);
            if stencilSize <= polynomialCount
                error('kp:manifold:InsufficientIBAMRMaterialSites', ...
                    'Local material interpolation needs more than %d source sites.', ...
                    polynomialCount);
            end
            properties = kp.rbffd.StencilProperties( ...
                'n', stencilSize, 'dim', 2, 'ell', ell, ...
                'spline_degree', splineDegree, 'npoly', polynomialCount);
            tree = KDTreeSearcher(sourceSites);
            neighbors = knnsearch(tree, targetSites, 'K', stencilSize);
            targetValues = zeros(size(targetSites, 1), size(sourceValues, 2));
            parfor targetId = 1:size(targetSites, 1)
                frame = tangentBasis(targetSites(targetId, :));
                ids = neighbors(targetId, :);
                localSites = (sourceSites(ids, :) - targetSites(targetId, :)) * frame; %#ok<PFBNS>
                stencil = kp.rbffd.RBFStencil();
                stencil.InitializeGeometry(localSites, properties);
                targetValues(targetId, :) = stencil.EvalStencil( ...
                    properties, zeros(1, 2), sourceValues(ids, :), false); %#ok<PFBNS>
            end
        end
    end

    methods (Access = private)
        function loadFiles(obj)
            required = ["material.csv", "faces.csv", "diagnostics.csv"];
            for filename = required
                path = fullfile(obj.Folder, filename);
                if ~isfile(path)
                    error('kp:manifold:MissingIBAMRTrajectoryFile', ...
                        'Missing IBAMR trajectory file: %s', path);
                end
            end

            materialTable = readtable(fullfile(obj.Folder, 'material.csv'));
            [ids, order] = sort(materialTable.node_id);
            if ~isequal(ids, (0:height(materialTable) - 1).')
                error('kp:manifold:BadIBAMRNodeIds', ...
                    'IBAMR material node IDs must be contiguous and zero based.');
            end
            obj.MaterialSites = normalizeRows([materialTable.u_x(order), ...
                materialTable.u_y(order), materialTable.u_z(order)]);

            faceTable = readtable(fullfile(obj.Folder, 'faces.csv'));
            obj.Faces = [faceTable.node_1, faceTable.node_2, faceTable.node_3] + 1;
            obj.Diagnostics = readtable(fullfile(obj.Folder, 'diagnostics.csv'));
            obj.Steps = obj.Diagnostics.step(:);
            obj.Times = obj.Diagnostics.time(:);

            N = size(obj.MaterialSites, 1);
            numFrames = numel(obj.Times);
            obj.Positions = zeros(N, 3, numFrames);
            obj.Velocities = zeros(N, 3, numFrames);
            for k = 1:numFrames
                filename = sprintf('frame_%06d.csv', obj.Steps(k));
                frame = readtable(fullfile(obj.Folder, filename));
                [frameIds, frameOrder] = sort(frame.node_id);
                if height(frame) ~= N || ~isequal(frameIds, (0:N - 1).')
                    error('kp:manifold:BadIBAMRFrame', ...
                        'Frame %s does not contain every material node once.', filename);
                end
                obj.Positions(:, :, k) = [frame.x(frameOrder), frame.y(frameOrder), frame.z(frameOrder)];
                obj.Velocities(:, :, k) = [frame.u(frameOrder), frame.v(frameOrder), frame.w(frameOrder)];
            end
        end

        function buildGeometryFactor(obj, requestedCount)
            initial = struct();
            initial.X = obj.Positions(:, :, 1);
            initial.material = struct('U', obj.MaterialSites);
            initial.area = obj.Diagnostics.area(1);
            initial.h = sqrt(initial.area / size(initial.X, 1));
            if isnan(requestedCount)
                [count, calibration] = kp.manifold.calibratedSBFControlPointCount( ...
                    initial, obj.Xi, 'Degree', obj.Degree, ...
                    'ControlPointScale', 1 / 3);
            elseif isfinite(requestedCount)
                count = min(size(initial.X, 1), max(1, round(requestedCount)));
                calibration = struct();
            else
                count = size(initial.X, 1);
                calibration = struct();
            end
            [obj.ControlIds, selection] = kp.manifold.selectSBFControlPoints(initial, 'Count', count);
            obj.ControlInfo = selection;
            fields = fieldnames(calibration);
            for k = 1:numel(fields)
                obj.ControlInfo.(fields{k}) = calibration.(fields{k});
            end
            obj.ControlSites = obj.MaterialSites(obj.ControlIds, :);
            [r, ~] = kp.geometry.sphereChordDistance(obj.ControlSites, obj.ControlSites);
            kernel = kp.geometry.phsKernel(r, obj.Degree);
            harmonics = realSphericalHarmonics(obj.ControlSites, obj.HarmonicDegree);
            regularization = 1.0e-13 * max(1.0, max(abs(kernel), [], 'all'));
            augmented = [kernel + regularization * eye(size(kernel)), harmonics; ...
                harmonics.', zeros(size(harmonics, 2))];
            obj.InterpolationFactor = decomposition( ...
                augmented, 'lu');
        end

        function [X, velocity] = interpolateNative(obj, t)
            t = min(max(t, obj.Times(1)), obj.Times(end));
            exact = find(abs(obj.Times - t) <= 10 * eps(max(1, abs(t))), 1);
            if ~isempty(exact)
                X = obj.Positions(:, :, exact);
                velocity = obj.Velocities(:, :, exact);
                return;
            end
            right = find(obj.Times > t, 1);
            left = right - 1;
            dt = obj.Times(right) - obj.Times(left);
            s = (t - obj.Times(left)) / dt;
            h00 = 2 * s^3 - 3 * s^2 + 1;
            h10 = s^3 - 2 * s^2 + s;
            h01 = -2 * s^3 + 3 * s^2;
            h11 = s^3 - s^2;
            X = h00 * obj.Positions(:, :, left) + ...
                h10 * dt * obj.Velocities(:, :, left) + ...
                h01 * obj.Positions(:, :, right) + ...
                h11 * dt * obj.Velocities(:, :, right);

            dh00 = (6 * s^2 - 6 * s) / dt;
            dh10 = 3 * s^2 - 4 * s + 1;
            dh01 = (-6 * s^2 + 6 * s) / dt;
            dh11 = 3 * s^2 - 2 * s;
            velocity = dh00 * obj.Positions(:, :, left) + ...
                dh10 * obj.Velocities(:, :, left) + ...
                dh01 * obj.Positions(:, :, right) + ...
                dh11 * obj.Velocities(:, :, right);
        end

        function kernel = evaluateKernel(obj, U)
            [r, ~] = kp.geometry.sphereChordDistance(U, obj.ControlSites);
            kernel = kp.geometry.phsKernel(r, obj.Degree);
        end

        function coefficients = solveCoefficients(obj, values)
            harmonicCount = (obj.HarmonicDegree + 1) ^ 2;
            coefficients = obj.InterpolationFactor \ ...
                [values; zeros(harmonicCount, size(values, 2))];
        end

        function values = evaluateValues(obj, U, coefficients)
            controlCount = numel(obj.ControlIds);
            kernel = obj.evaluateKernel(U);
            harmonics = realSphericalHarmonics(U, obj.HarmonicDegree);
            values = kernel * coefficients(1:controlCount, :) + ...
                harmonics * coefficients(controlCount + 1:end, :);
        end

        function [X, tangentAzimuth, tangentElevation] = evaluateModel(obj, U, coefficients)
            U = normalizeRows(U);
            [r, ~] = kp.geometry.sphereChordDistance(U, obj.ControlSites);
            kernel = kp.geometry.phsKernel(r, obj.Degree);
            [harmonics, harmonicAzimuth, harmonicElevation] = ...
                realSphericalHarmonics(U, obj.HarmonicDegree);
            controlCount = numel(obj.ControlIds);
            radialCoefficients = coefficients(1:controlCount, :);
            harmonicCoefficients = coefficients(controlCount + 1:end, :);
            X = kernel * radialCoefficients + harmonics * harmonicCoefficients;

            uvw = kp.geometry.cart2sphRows(U);
            azimuth = uvw(:, 1);
            elevation = uvw(:, 2);
            dAzimuth = [-cos(elevation) .* sin(azimuth), ...
                cos(elevation) .* cos(azimuth), zeros(size(azimuth))];
            dElevation = [-sin(elevation) .* cos(azimuth), ...
                -sin(elevation) .* sin(azimuth), cos(elevation)];
            difference = kp.geometry.RBFLevelSet.differenceTensor(U, obj.ControlSites);
            factor = kp.geometry.RBFLevelSet.radialDerivativeFactor(r, obj.Degree);
            dotAzimuth = difference(:, :, 1) .* dAzimuth(:, 1) + ...
                difference(:, :, 2) .* dAzimuth(:, 2) + ...
                difference(:, :, 3) .* dAzimuth(:, 3);
            dotElevation = difference(:, :, 1) .* dElevation(:, 1) + ...
                difference(:, :, 2) .* dElevation(:, 2) + ...
                difference(:, :, 3) .* dElevation(:, 3);
            tangentAzimuth = (factor .* dotAzimuth) * radialCoefficients + ...
                harmonicAzimuth * harmonicCoefficients;
            tangentElevation = (factor .* dotElevation) * radialCoefficients + ...
                harmonicElevation * harmonicCoefficients;
        end
    end
end


function [Y, dAzimuth, dElevation] = realSphericalHarmonics(U, maxDegree)
uvw = kp.geometry.cart2sphRows(normalizeRows(U));
azimuth = uvw(:, 1);
elevation = uvw(:, 2);
z = sin(elevation);
cosElevation = cos(elevation);
numTerms = (maxDegree + 1) ^ 2;
Y = zeros(size(U, 1), numTerms);
dAzimuth = zeros(size(Y));
dElevation = zeros(size(Y));
column = 0;
for ell = 0:maxDegree
    associated = legendre(ell, z.', 'unnorm');
    if ell > 0
        previous = legendre(ell - 1, z.', 'unnorm');
    else
        previous = zeros(0, numel(z));
    end
    for order = 0:ell
        normalization = sqrt((2 * ell + 1) / (4 * pi) * ...
            exp(gammaln(ell - order + 1) - gammaln(ell + order + 1)));
        P = associated(order + 1, :).';
        if ell == 0
            derivativeP = zeros(size(P));
        else
            if order <= ell - 1
                previousP = previous(order + 1, :).';
            else
                previousP = zeros(size(P));
            end
            denominator = z .^ 2 - 1;
            safeDenominator = signNonzero(denominator) .* ...
                max(abs(denominator), 1.0e-12);
            derivativeP = (ell * z .* P - (ell + order) .* previousP) ./ ...
                safeDenominator;
        end
        derivativeElevation = derivativeP .* cosElevation;
        if order == 0
            column = column + 1;
            Y(:, column) = normalization * P;
            dElevation(:, column) = normalization * derivativeElevation;
        else
            scale = sqrt(2) * normalization;
            cosine = cos(order * azimuth);
            sine = sin(order * azimuth);
            column = column + 1;
            Y(:, column) = scale * P .* cosine;
            dAzimuth(:, column) = -scale * order * P .* sine;
            dElevation(:, column) = scale * derivativeElevation .* cosine;
            column = column + 1;
            Y(:, column) = scale * P .* sine;
            dAzimuth(:, column) = scale * order * P .* cosine;
            dElevation(:, column) = scale * derivativeElevation .* sine;
        end
    end
end
end

function s = signNonzero(x)
s = sign(x);
s(s == 0) = 1;
end

function X = normalizeRows(X)
X = X ./ max(vecnorm(X, 2, 2), eps);
end

function frame = tangentBasis(normal)
normal = normal(:) / max(norm(normal), eps);
if abs(normal(3)) < 0.9
    seed = [0; 0; 1];
else
    seed = [1; 0; 0];
end
tangent1 = cross(normal, seed);
tangent1 = tangent1 / max(norm(tangent1), eps);
tangent2 = cross(normal, tangent1);
frame = [tangent1, tangent2];
end

function weights = sphericalLumpedWeights(U)
faces = convhull(U(:, 1), U(:, 2), U(:, 3));
weights = zeros(size(U, 1), 1);
for k = 1:size(faces, 1)
    a = U(faces(k, 1), :);
    b = U(faces(k, 2), :);
    c = U(faces(k, 3), :);
    numerator = abs(dot(a, cross(b, c)));
    denominator = 1 + dot(a, b) + dot(b, c) + dot(c, a);
    area = 2 * atan2(numerator, max(denominator, eps));
    weights(faces(k, :)) = weights(faces(k, :)) + area / 3;
end
end
