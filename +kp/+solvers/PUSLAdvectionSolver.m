classdef PUSLAdvectionSolver < handle
    %PUSLADVECTIONSOLVER MATLAB PU semi-Lagrangian transport solver.

    properties (SetAccess = private)
        Domain kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
        xi (1,1) double = 0
        dt (1,1) double = 0
        output_nodes double = zeros(0, 0)
        output_range (1,2) double = [1 0]
        patch_centers double = zeros(0, 0)
        patch_radius (1,1) double = 0
        patch_spacing (1,1) double = 0
        min_patch_nodes (1,1) double = 0
        patch_node_ids cell = {}
        domain_measure (1,1) double = 0
        boundary_condition struct = struct('mode', "unspecified", 'normal_velocity_tolerance', 1.0e-10, 'periodic_patches', [], 'inflow_value', [])
    end

    methods
        function clearAdvectionBoundaryCondition(obj)
            obj.boundary_condition = struct('mode', "unspecified", 'normal_velocity_tolerance', 1.0e-10, 'periodic_patches', [], 'inflow_value', []);
        end

        function setTangentialFlowBoundary(obj, normal_velocity_tolerance)
            if nargin < 2 || isempty(normal_velocity_tolerance)
                normal_velocity_tolerance = 1.0e-10;
            end
            obj.boundary_condition = struct('mode', "tangential", ...
                'normal_velocity_tolerance', normal_velocity_tolerance, ...
                'periodic_patches', [], ...
                'inflow_value', []);
        end

        function setPeriodicBoundary(obj, periodic_patches, normal_velocity_tolerance)
            if nargin < 3 || isempty(normal_velocity_tolerance)
                normal_velocity_tolerance = 1.0e-10;
            end
            obj.boundary_condition = struct('mode', "periodic", ...
                'normal_velocity_tolerance', normal_velocity_tolerance, ...
                'periodic_patches', periodic_patches, ...
                'inflow_value', []);
        end

        function setInflowDirichletBoundary(obj, inflow_value, normal_velocity_tolerance)
            if nargin < 3 || isempty(normal_velocity_tolerance)
                normal_velocity_tolerance = 1.0e-10;
            end
            obj.boundary_condition = struct('mode', "inflow_dirichlet", ...
                'normal_velocity_tolerance', normal_velocity_tolerance, ...
                'periodic_patches', [], ...
                'inflow_value', inflow_value);
        end

        function out = getAdvectionBoundaryCondition(obj)
            out = obj.boundary_condition;
        end

        function out = returnsDistributedState(~)
            out = false;
        end

        function out = getOutputRange(obj)
            out = obj.output_range;
        end

        function init(obj, domain, xi, dlt, patch_spacing_factor, patch_radius_factor)
            if nargin < 5 || isempty(patch_spacing_factor)
                patch_spacing_factor = 0.0;
            end
            if nargin < 6 || isempty(patch_radius_factor)
                patch_radius_factor = 0.0;
            end

            obj.Domain = domain;
            obj.Domain.buildStructs();
            obj.xi = xi;
            obj.dt = dlt;
            obj.output_nodes = obj.Domain.getIntBdryNodes();
            obj.output_range = [1, size(obj.output_nodes, 1)];

            h = obj.Domain.getSepRad();
            % The PU-SL solver is driven by a patch covering of the
            % physical node cloud rather than by one global interpolant.
            obj.patch_spacing = choosePatchSpacing(h, patch_spacing_factor);
            obj.patch_radius = choosePatchRadius(h, patch_radius_factor);
            obj.min_patch_nodes = chooseMinimumPatchNodes(size(obj.output_nodes, 2), xi);
            obj.patch_centers = choosePatchCenters(obj.output_nodes, obj.patch_spacing);
            obj.patch_node_ids = buildPatchNodeIds(obj.Domain, obj.patch_centers, obj.patch_radius, obj.min_patch_nodes);
            obj.domain_measure = estimateDomainMeasure(obj.Domain);
        end

        function out = getOutputNodes(obj)
            out = obj.output_nodes;
        end

        function coeffs = projectInitial(obj, rho0)
            % The coefficient vector is just the nodal field sampled at the
            % output nodes in this MATLAB mirror.
            X = obj.output_nodes;
            coeffs = zeros(size(X, 1), 1);
            for i = 1:size(X, 1)
                coeffs(i, 1) = rho0(X(i, :));
            end
        end

        function coeffs = projectConstant(obj, value, n_cols)
            if nargin < 3 || isempty(n_cols)
                n_cols = 1;
            end
            coeffs = value * ones(size(obj.output_nodes, 1), n_cols);
        end

        function coeffs = projectSamples(~, nodal_samples)
            coeffs = nodal_samples;
        end

        function values = evaluateAtNodes(~, coeffs)
            values = coeffs;
        end

        function values = evaluateAtPoints(obj, coeffs, X)
            % Evaluation at arbitrary points is done by localized PU
            % reconstruction over the patch covering.
            values = localizedEvaluate(obj, coeffs, X);
        end

        function coeffs_next = backwardSLStep(obj, tn, coeffs_old, velocity, rk)
            % Backward SL: trace departure points, evaluate there, then
            % write the transported values back at the Eulerian nodes.
            X = obj.output_nodes;
            Xdep = tracePointsBackward(X, tn, obj.dt, velocity, rk);
            values = localizedEvaluate(obj, coeffs_old, Xdep);
            values = applyBoundaryCondition(obj, values, Xdep, tn + obj.dt, velocity);
            coeffs_next = values;
        end

        function coeffs_next = forwardSLStep(obj, tn, coeffs_old, velocity, rk)
            % Forward SL is kept for parity with KernelPack, though the
            % current verification work has focused on the backward path.
            X = obj.output_nodes;
            Xarr = tracePointsForward(X, tn, obj.dt, velocity, rk);
            coeffs_arr = coeffs_old;
            values = localizedEvaluate(obj, coeffs_arr, Xarr);
            coeffs_next = localizedEvaluate(obj, values, X);
        end

        function setStepSize(obj, dlt)
            obj.dt = dlt;
        end

        function resetSolveStats(~)
        end

        function out = getSolveStats(~)
            out = struct();
        end

        function out = getNumPatches(obj)
            out = size(obj.patch_centers, 1);
        end

        function out = getNumDofs(obj)
            out = size(obj.output_nodes, 1);
        end

        function out = getPatchRadius(obj)
            out = obj.patch_radius;
        end

        function out = getPatchSpacing(obj)
            out = obj.patch_spacing;
        end

        function out = totalMass(obj, coeffs, col)
            if nargin < 3 || isempty(col)
                col = 1;
            end
            values = coeffs(:, col);
            out = obj.domain_measure * mean(values);
        end

        function out = getDomainMeasure(obj)
            out = obj.domain_measure;
        end
    end
end

function spacing = choosePatchSpacing(h, patch_spacing_factor)
if patch_spacing_factor > 0
    spacing = patch_spacing_factor * h;
else
    spacing = 2.0 * h;
end
end

function radius = choosePatchRadius(h, patch_radius_factor)
if patch_radius_factor > 0
    radius = patch_radius_factor * h;
else
    radius = 3.0 * h;
end
end

function centers = choosePatchCenters(X, spacing)
if isempty(X)
    centers = zeros(0, size(X, 2));
    return;
end

remaining = true(size(X, 1), 1);
centers = zeros(0, size(X, 2));
for i = 1:size(X, 1)
    if ~remaining(i)
        continue;
    end
    xi = X(i, :);
    centers = [centers; xi]; %#ok<AGROW>
    d = sqrt(sum((X - xi).^2, 2));
    remaining(d <= spacing) = false;
end
end

function patch_node_ids = buildPatchNodeIds(domain, centers, radius, min_patch_nodes)
% Make each PU patch robust by forcing a minimum node count, falling back
% to KNN if the geometric ball query is too sparse.
patch_node_ids = cell(size(centers, 1), 1);
for p = 1:size(centers, 1)
    ids = domain.queryBall("interior_boundary", centers(p, :), radius);
    patch_node_ids{p} = ids{1}(:);
    if numel(patch_node_ids{p}) < min_patch_nodes
        [idsKnn, ~] = domain.queryKnn("interior_boundary", centers(p, :), min(min_patch_nodes, domain.getNumIntBdryNodes()));
        patch_node_ids{p} = idsKnn(1, :).';
    end
end
end

function min_nodes = chooseMinimumPatchNodes(dim, xi)
ell = max(xi + 1, 2);
npoly = size(kp.poly.total_degree_indices(dim, ell), 1);
min_nodes = 2 * npoly + 1;
end

function values = localizedEvaluate(obj, coeffs, Xq)
if isempty(Xq)
    values = zeros(0, size(coeffs, 2));
    return;
end

Xnodes = obj.output_nodes;
dim = size(Xnodes, 2);
nq = size(Xq, 1);
nc = size(coeffs, 2);
values = zeros(nq, nc);
weight_sum = zeros(nq, 1);

sp = kp.rbffd.StencilProperties();
sp.dim = dim;
sp.ell = max(obj.xi + 1, 2);
sp.npoly = size(kp.poly.total_degree_indices(dim, sp.ell), 1);
sp.spline_degree = max(5, 2 * floor((sp.ell + 1) / 2) - 1);

% Evaluate every active patch contribution and blend them with compact
% Wendland weights.
for p = 1:size(obj.patch_centers, 1)
    center = obj.patch_centers(p, :);
    diff = Xq - center;
    dist = sqrt(sum(diff.^2, 2));
    mask = dist < obj.patch_radius;
    if ~any(mask)
        continue;
    end
    ids = obj.patch_node_ids{p};
    Xloc = Xnodes(ids, :);
    floc = coeffs(ids, :);
    stencil = kp.rbffd.RBFStencil();
    stencil.InitializeGeometry(Xloc, sp);
    qloc = Xq(mask, :);
    vloc = stencil.EvalStencil(sp, qloc, floc, false);
    wloc = wendlandC6(dist(mask) ./ obj.patch_radius);
    values(mask, :) = values(mask, :) + wloc .* vloc;
    weight_sum(mask) = weight_sum(mask) + wloc;
end

% If no patch contributes, fall back to the nearest nodal value so the
% routine stays total over the query set.
missing = weight_sum <= 1.0e-14;
if any(missing)
    [idx, ~] = obj.Domain.queryKnn("interior_boundary", Xq(missing, :), 1);
    values(missing, :) = coeffs(idx(:, 1), :);
    weight_sum(missing) = 1.0;
end

values = values ./ weight_sum;
end

function w = wendlandC6(r)
w = zeros(size(r));
mask = r < 1;
t = 1 - r(mask);
rm = r(mask);
w(mask) = t.^8 .* (32 * rm.^3 + 25 * rm.^2 + 8 * rm + 1);
end

function Xdep = tracePointsBackward(X, tn, dt, velocity, rk)
Xdep = tracePointsSigned(X, tn, -dt, velocity, rk);
end

function Xarr = tracePointsForward(X, tn, dt, velocity, rk)
Xarr = tracePointsSigned(X, tn, dt, velocity, rk);
end

function Xnext = tracePointsSigned(X, tn, dt, velocity, rk)
if isa(rk, 'function_handle')
    Xnext = rk(tn, X, dt, velocity);
else
    Xnext = rk4step(tn, X, dt, velocity);
end
end

function Xnext = rk4step(t, X, dt, velocity)
X1 = X;
K1 = velocity(t, X1);
X2 = X + 0.5 * dt * K1;
K2 = velocity(t + 0.5 * dt, X2);
X3 = X + 0.5 * dt * K2;
K3 = velocity(t + 0.5 * dt, X3);
X4 = X + dt * K3;
K4 = velocity(t + dt, X4);
Xnext = X + (dt / 6) * (K1 + 2 * K2 + 2 * K3 + K4);
end

function values = applyBoundaryCondition(obj, values, Xdep, tnext, ~)
bc = obj.boundary_condition;
if bc.mode == "unspecified" || bc.mode == "tangential"
    return;
end
if bc.mode == "inflow_dirichlet"
    phi = obj.Domain.getOuterLevelSet().Evaluate(Xdep);
    outside = phi > 0;
    if any(outside)
        inflow = bc.inflow_value(tnext, Xdep(outside, :));
        values(outside, :) = inflow;
    end
    return;
end
if bc.mode == "periodic"
    error('kp:solvers:PeriodicNotYetSupported', 'Periodic PU-SL transport is not yet supported in MATLAB.');
end
end

function measure = estimateDomainMeasure(domain)
phi = domain.getOuterLevelSet();
if isempty(phi)
    measure = size(domain.getInteriorNodes(), 1) * domain.getSepRad()^domain.getDim();
    return;
end
Xall = domain.getIntBdryNodes();
xmin = min(Xall, [], 1) - domain.getSepRad();
xmax = max(Xall, [], 1) + domain.getSepRad();
dim = size(Xall, 2);
ngrid = 32;
switch dim
    case 2
        x = linspace(xmin(1), xmax(1), ngrid);
        y = linspace(xmin(2), xmax(2), ngrid);
        [Xg, Yg] = ndgrid(x, y);
        pts = [Xg(:), Yg(:)];
    case 3
        x = linspace(xmin(1), xmax(1), 18);
        y = linspace(xmin(2), xmax(2), 18);
        z = linspace(xmin(3), xmax(3), 18);
        [Xg, Yg, Zg] = ndgrid(x, y, z);
        pts = [Xg(:), Yg(:), Zg(:)];
    otherwise
        measure = size(domain.getInteriorNodes(), 1) * domain.getSepRad()^dim;
        return;
end
vals = phi.Evaluate(pts);
frac = mean(vals <= 0);
measure = frac * prod(xmax - xmin);
end
