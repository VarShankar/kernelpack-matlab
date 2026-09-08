function showcase = moving_surface_adr_tp_showcase_bumpy_sphere(varargin)
%MOVING_SURFACE_ADR_TP_SHOWCASE_BUMPY_SPHERE Source-free hard moving surface.
%   Demonstrates the method on a rotating, breathing, bumpy sphere with no
%   manufactured solution.  By default the run uses global parametric SBF
%   normals and reports source-free mass drift.

parser = inputParser();
parser.addParameter('N', 1024, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('Xi', 4, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('Mu', 0.02, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('FinalTimes', [0.0, 0.20, 0.40], @(x) isnumeric(x) && isvector(x));
parser.addParameter('DtScale', 0.05, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('NormalMode', "globalSbf");
parser.addParameter('NormalNeighborCount', 32, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('GlobalSBFNormalDegree', 7, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('GlobalSBFControlPointCount', Inf, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('UseRearrangement', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('QualityThreshold', 1.75, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('PredictiveLookaheadSteps', 3, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('MinStepsBetweenRearrangements', 8, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('RearrangementTransferMode', "localTp", @(x) isstring(x) || ischar(x));
parser.addParameter('MassCorrectionMode', "balance", @(x) isstring(x) || ischar(x));
parser.addParameter('ImageDir', fullfile(pwd, 'docs', 'figures'));
parser.addParameter('OutputPrefix', fullfile(pwd, 'moving_surface_adr_tp_showcase_bumpy_sphere'));
parser.addParameter('RedrawOnly', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('WriteFigure', true, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});
opts = parser.Results;

if ~exist(opts.ImageDir, 'dir')
    mkdir(opts.ImageDir);
end

if opts.RedrawOnly
    savedPath = string(opts.OutputPrefix) + ".mat";
    if ~isfile(savedPath)
        error('kp:examples:MissingSavedShowcase', ...
            'Cannot redraw because saved showcase does not exist: %s', savedPath);
    end
    S = load(savedPath, 'showcase');
    showcase = S.showcase;
    showcase.figurePath = writeSnapshotFigure(showcase.snapshots, opts.ImageDir);
    save(savedPath, 'showcase');
    return;
end

problem = bumpySphereProblem();
times = opts.FinalTimes(:).';
snapshots = repmat(struct('time', 0, 'X', [], 'C', [], 'massDrift', NaN, ...
    'balanceResidual', NaN), 1, numel(times));

geom0 = problem.geometry(opts.N, 0);
c0 = problem.initial(geom0.material);
m0 = sum(surfaceWeights(geom0) .* c0);

for it = 1:numel(times)
    T = times(it);
    geom = problem.geometry(opts.N, T);
    if T == 0
        C = c0;
        massDrift = 0;
        balanceResidual = 0;
    else
        result = kp.manifold.runLagrangianMovingADRConvergence(problem, opts.Xi, opts.N, opts.DtScale, ...
            'Mu', opts.Mu, ...
            'FinalTime', T, ...
            'DiffMatUpdateMethod', 'defect', ...
            'DefectTolerance', 1.0e-6, ...
            'MaxDefectIterations', 4, ...
            'SpectrumCheck', false, ...
            'HyperviscosityUpdateMode', 'adaptiveGeometry', ...
            'NormalMode', opts.NormalMode, ...
            'NormalNeighborCount', opts.NormalNeighborCount, ...
            'GlobalSBFNormalDegree', opts.GlobalSBFNormalDegree, ...
            'GlobalSBFControlPointCount', opts.GlobalSBFControlPointCount, ...
            'UseRearrangement', opts.UseRearrangement, ...
            'QualityThreshold', opts.QualityThreshold, ...
            'PredictiveLookaheadSteps', opts.PredictiveLookaheadSteps, ...
            'MinStepsBetweenRearrangements', opts.MinStepsBetweenRearrangements, ...
            'RearrangementTransferMode', opts.RearrangementTransferMode, ...
            'MassCorrectionMode', opts.MassCorrectionMode, ...
            'ReturnSolution', true);
        C = result.solution{1};
        massDrift = abs(result.finalMass - m0) / max(abs(m0), 1.0e-14);
        balanceResidual = result.balanceRelResidual;
    end
    snapshots(it).time = T;
    snapshots(it).X = geom.X;
    snapshots(it).C = C;
    snapshots(it).massDrift = massDrift;
    snapshots(it).balanceResidual = balanceResidual;
end

showcase = struct();
showcase.problem = problem;
showcase.N = opts.N;
showcase.xi = opts.Xi;
showcase.mu = opts.Mu;
showcase.snapshots = snapshots;
showcase.table = snapshotTable(snapshots);
showcase.figurePath = "";
if opts.WriteFigure
    showcase.figurePath = writeSnapshotFigure(snapshots, opts.ImageDir);
end

save(string(opts.OutputPrefix) + ".mat", 'showcase');
writetable(showcase.table, string(opts.OutputPrefix) + ".csv");
end

function problem = bumpySphereProblem()
problem.label = "Source-free ADR on a rotating bumpy sphere";
problem.title = "Source-free ADR on a rotating bumpy sphere";
problem.finalTime = 0.40;
problem.geometry = @bumpySphereGeometry;
problem.geometryFromMaterial = @bumpySphereGeometryFromMaterial;
problem.sampleMaterial = @(N, t) bumpySphereMaterial(N);
problem.backtraceMaterial = @(material, tNow, tPast) material;
problem.initial = @initialCondition;
problem.forcing = @(t, material, mu) zeros(size(material.U, 1), 1);
problem.balanceSource = @(t, material, mu) zeros(size(material.U, 1), 1);
end

function geom = bumpySphereGeometry(N, t)
geom = bumpySphereGeometryFromMaterial(bumpySphereMaterial(N), t);
end

function material = bumpySphereMaterial(N)
U = kp.geometry.fibonacciSphere(N);
U = U ./ vecnorm(U, 2, 2);
theta = acos(min(max(U(:, 3), -1), 1));
phi = atan2(U(:, 2), U(:, 1));
material = struct('U', U, 'theta', theta, 'phi', phi);
end

function geom = bumpySphereGeometryFromMaterial(material, t)
U = material.U;
theta = material.theta;
phi = material.phi;
N = size(U, 1);
Q = rotationZ(0.8 * t);
Y = U * Q.';
thetaY = acos(min(max(Y(:, 3), -1), 1));
phiY = atan2(Y(:, 2), Y(:, 1));

breath = 1.0 + 0.10 * sin(1.2 * t);
a = 0.16 * sin(1.7 * t + 0.25);
b = 0.08 * cos(1.1 * t);
s = sin(thetaY);
ct = cos(thetaY);
bump = a .* (s.^4) .* cos(5 * phiY) ...
    + b .* (s.^3) .* sin(3 * phiY);
r = breath + bump;
geom.X = r .* Y;
drTheta = a .* 4 .* s.^3 .* ct .* cos(5 * phiY) ...
    + b .* 3 .* s.^2 .* ct .* sin(3 * phiY);
drPhi = -5 .* a .* s.^4 .* sin(5 * phiY) ...
    + 3 .* b .* s.^3 .* cos(3 * phiY);
grad2 = drTheta.^2 + (drPhi.^2 ./ max(s.^2, 1.0e-14));
areaDensity = r .* sqrt(max(r.^2 + grad2, 0));
geom.weights = (4 * pi / N) .* areaDensity;
geom.area = sum(geom.weights);
geom.h = sqrt(geom.area / N);
geom.material = struct('U', U, 'theta', theta, 'phi', phi, 'X', geom.X);
% No analytic normals are supplied here; globalSbf uses these material sites.
end

function c = initialCondition(material)
U = material.U;
a = [0.45, -0.25, 0.8573214099741123];
a = a ./ norm(a);
b = [-0.55, 0.60, 0.5809475019311126];
b = b ./ norm(b);
c = exp(-18 * (1 - U * a.')) + 0.65 * exp(-28 * (1 - U * b.'));
end

function Q = rotationZ(theta)
c = cos(theta);
s = sin(theta);
Q = [c, -s, 0; s, c, 0; 0, 0, 1];
end

function w = surfaceWeights(geom)
w = geom.weights(:);
end

function T = snapshotTable(snapshots)
n = numel(snapshots);
time = zeros(n, 1);
massDrift = zeros(n, 1);
balanceResidual = zeros(n, 1);
for i = 1:n
    time(i) = snapshots(i).time;
    massDrift(i) = snapshots(i).massDrift;
    balanceResidual(i) = snapshots(i).balanceResidual;
end
T = table(time, massDrift, balanceResidual);
end

function imagePath = writeSnapshotFigure(snapshots, imageDir)
fig = figure('Color', 'w', 'Position', [100, 100, 1380, 460]);
tl = tiledlayout(fig, 1, numel(snapshots), 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(tl, 'Source-free transport on a moving bumpy sphere', ...
    'FontSize', 16, 'FontWeight', 'bold');
for i = 1:numel(snapshots)
    nexttile(tl);
    X = snapshots(i).X;
    C = snapshots(i).C;
    K = convhull(X(:, 1), X(:, 2), X(:, 3));
    trisurf(K, X(:, 1), X(:, 2), X(:, 3), C, ...
        'EdgeColor', 'none', 'FaceColor', 'interp');
    axis equal off;
    view(35, 20);
    camlight headlight;
    lighting gouraud;
    title(sprintf('t = %.2f', snapshots(i).time));
end
colormap(turbo);
cb = colorbar;
cb.Layout.Tile = 'east';
imagePath = fullfile(imageDir, 'moving_surface_adr_tp_showcase_bumpy_sphere_snapshots.png');
kp.plot.exportPaperFigure(fig, imagePath);
close(fig);
end
