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
        patch_radii double = zeros(0, 1)
        patch_spacing (1,1) double = 0
        min_patch_nodes (1,1) double = 0
        patch_node_ids cell = {}
        patch_center_tree struct = struct('Points', [], 'Searcher', [], 'HasSearcher', false)
        patch_stencil_props kp.rbffd.StencilProperties = kp.rbffd.StencilProperties()
        patch_stencils cell = {}
        patch_mass_rows cell = {}
        coeff_offsets double = zeros(0, 1)
        total_coeffs (1,1) double = 0
        constant_mode double = zeros(0, 1)
        fixed_collocation double = zeros(0, 0)
        fixed_collocation_mass_augmented double = zeros(0, 0)
        nodal_mass_weights double = zeros(0, 1)
        local_constant_column double = zeros(0, 1)
        forward_cardinal_coeffs cell = {}
        domain_measure (1,1) double = 0
        inflow_forward_fallback_reported (1,1) logical = false
        boundary_condition struct = struct('mode', "unspecified", 'normal_velocity_tolerance', 1.0e-10, 'periodic_patches', [], 'inflow_value', [])
        solve_stats struct = defaultSolveStats()
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

            if obj.boundary_condition.mode == "periodic"
                error('kp:solvers:PeriodicNotYetSupported', ...
                    ['PUSL periodic advection is not yet supported robustly. ' ...
                     'The localized PU basis still needs true periodic patch support, not just wrapped traces.']);
            end

            Xi = obj.Domain.getInteriorNodes();
            Xb = obj.Domain.getBdryNodes();
            h = obj.Domain.getSepRad();
            tuning = autoTunePatchGeometry(Xi, Xb, obj.Domain.getOuterLevelSet(), obj.xi, h, patch_radius_factor, patch_spacing_factor);
            obj.patch_spacing = tuning.spacing;
            obj.patch_radius = tuning.radius;
            obj.min_patch_nodes = tuning.target_points;
            [obj.patch_centers, obj.patch_radii, obj.patch_node_ids] = buildPatchData(Xi, Xb, tuning);
            obj.patch_center_tree = buildCenterTree(obj.patch_centers);
            obj.patch_stencil_props = buildPatchStencilProperties(size(obj.output_nodes, 2), obj.xi);
            obj.patch_stencils = buildPatchStencils(obj.output_nodes, obj.patch_node_ids, obj.patch_stencil_props);
            obj.coeff_offsets = buildCoeffOffsets(obj.patch_node_ids, obj.patch_stencil_props.npoly);
            obj.total_coeffs = obj.coeff_offsets(end) - 1;

            [obj.constant_mode, obj.patch_mass_rows, obj.domain_measure] = precomputePatchMassData(obj);
            obj.forward_cardinal_coeffs = cell(size(obj.patch_node_ids));
            for p = 1:numel(obj.patch_node_ids)
                nloc = numel(obj.patch_node_ids{p});
                rhs = [eye(nloc); zeros(obj.patch_stencil_props.npoly, nloc)];
                obj.forward_cardinal_coeffs{p} = solveWithStencilCPQR(obj.patch_stencils{p}, rhs);
            end

            obj.fixed_collocation = nodalEvaluationWeights(obj, obj.output_nodes);
            [obj.fixed_collocation_mass_augmented, obj.nodal_mass_weights, obj.local_constant_column] = ...
                buildForwardMassAugmentedData(obj);
            obj.solve_stats = defaultSolveStats();
        end

        function out = getOutputNodes(obj)
            out = obj.output_nodes;
        end

        function coeffs = projectInitial(obj, rho0)
            coeffs = obj.projectSamples(sampleScalarFunctionAtOutputNodes(obj, rho0));
        end

        function coeffs = projectConstant(obj, value, n_cols)
            if nargin < 3 || isempty(n_cols)
                n_cols = 1;
            end
            coeffs = value * repmat(obj.constant_mode, 1, n_cols);
        end

        function coeffs = projectSamples(obj, nodal_samples)
            if size(nodal_samples, 1) ~= size(obj.output_nodes, 1)
                error('kp:solvers:BadSampleCount', ...
                    'Localized PU projectSamples expects one sample per output node.');
            end
            coeffs = zeros(obj.total_coeffs, size(nodal_samples, 2));
            for p = 1:numel(obj.patch_node_ids)
                ids = obj.patch_node_ids{p};
                coeffs(localCoeffRange(obj, p), :) = solveLocalCoeffs(obj, p, obj.output_nodes(ids, :), nodal_samples(ids, :));
            end
        end

        function values = evaluateAtNodes(obj, coeffs)
            values = obj.evaluateAtPoints(coeffs, obj.output_nodes);
        end

        function values = evaluateAtPoints(obj, coeffs, X)
            if size(coeffs, 1) ~= obj.total_coeffs
                error('kp:solvers:BadCoeffCount', ...
                    'Localized PU evaluateAtPoints received coefficient matrix with wrong row count.');
            end
            values = localizedEvaluateCoefficients(obj, coeffs, X);
        end

        function coeffs_next = backwardSLStep(obj, tn, coeffs_old, velocity, rk)
            validateBoundaryConditionConfiguration(obj);
            bc = obj.boundary_condition;
            preserve_mass = bc.mode ~= "inflow_dirichlet";
            if bc.mode == "tangential"
                validateTangentialBoundaryFlow(obj, tn, velocity);
                validateTangentialBoundaryFlow(obj, tn + obj.dt, velocity);
            end

            X = obj.output_nodes;
            Xdep = tracePointsBackward(X, tn + obj.dt, obj.dt, velocity, rk);
            values = applyBackwardBoundaryCondition(obj, Xdep, X, coeffs_old, tn + obj.dt);
            coeffs_next = obj.projectSamples(values);
            if preserve_mass
                coeffs_next = enforceMassCorrection(obj, coeffs_next, totalMasses(obj, coeffs_old));
            end
        end

        function coeffs_next = forwardSLStep(obj, tn, coeffs_old, velocity, rk)
            validateBoundaryConditionConfiguration(obj);
            bc = obj.boundary_condition;
            if bc.mode == "inflow_dirichlet"
                if ~obj.inflow_forward_fallback_reported
                    warning('kp:solvers:ForwardInflowFallback', ...
                        ['PUSLAdvectionSolver.forwardSLStep requested with inflow-Dirichlet advection BCs. ' ...
                         'Falling back to the backward update because forward SL does not yet inject ' ...
                         'the full boundary-in-time inflow data robustly.']);
                    obj.inflow_forward_fallback_reported = true;
                end
                coeffs_next = obj.backwardSLStep(tn, coeffs_old, velocity, rk);
                return;
            end

            preserve_mass = bc.mode ~= "inflow_dirichlet";
            if bc.mode == "tangential"
                validateTangentialBoundaryFlow(obj, tn, velocity);
                validateTangentialBoundaryFlow(obj, tn + obj.dt, velocity);
            end

            if preserve_mass
                target_mass = totalMasses(obj, coeffs_old);
            else
                target_mass = zeros(1, size(coeffs_old, 2));
            end

            nodal_old = obj.evaluateAtNodes(coeffs_old);
            rhs_local = nodal_old;
            nodal_current = nodal_old;
            traced_points = tracePointsForward(obj.output_nodes, tn, obj.dt, velocity, rk);
            if bc.mode == "periodic"
                error('kp:solvers:PeriodicNotYetSupported', 'Periodic PU-SL transport is not yet supported in MATLAB.');
            end
            if ~isempty(obj.Domain.getOuterLevelSet())
                boundary_start = obj.Domain.getNumInteriorNodes() + 1;
                Xb = obj.Domain.getBdryNodes();
                nr = obj.Domain.getNrmls();
                for local_row = 1:size(traced_points, 1)
                    is_boundary_row = local_row >= boundary_start;
                    x = obj.output_nodes(local_row, :);
                    if bc.mode == "inflow_dirichlet" && is_boundary_row
                        boundary_idx = local_row - boundary_start + 1;
                        bn = boundaryNormalVelocity(Xb(boundary_idx, :), nr(boundary_idx, :), tn + obj.dt, velocity);
                        if bn < -bc.normal_velocity_tolerance
                            traced_points(local_row, :) = x;
                            rhs_local(local_row, :) = bc.inflow_value(tn + obj.dt, x);
                            continue;
                        end
                    end

                    if isInsideLevelSetPoint(obj, traced_points(local_row, :))
                        continue;
                    end
                    if bc.mode == "tangential"
                        error('kp:solvers:TangentialForwardExit', ...
                            'PUSL forward SL exited the domain under a tangential-flow boundary condition.');
                    end
                    if bc.mode == "unspecified"
                        error('kp:solvers:ForwardOutsideDomain', ...
                            'PUSL forward SL traced outside the domain without an advection boundary condition.');
                    end
                    if bc.mode == "inflow_dirichlet"
                        hit = boundaryHitOnSegment(obj, x, traced_points(local_row, :));
                        traced_points(local_row, :) = hit.point;
                    end
                end
            end
            traced_weights = nodalEvaluationWeights(obj, traced_points);

            max_defect_sweeps = 4;
            for iter = 1:max_defect_sweeps
                forward_values = traced_weights * nodal_current;
                residual = rhs_local - forward_values;
                if isempty(residual)
                    resid_norm = 0.0;
                else
                    resid_norm = norm(residual, inf);
                end
                mass_resid_norm = 0.0;
                mass_resid = [];
                if preserve_mass
                    current_mass = obj.nodal_mass_weights.' * nodal_current;
                    mass_resid = target_mass - current_mass;
                    if ~isempty(mass_resid)
                        mass_resid_norm = norm(mass_resid, inf);
                    end
                end
                if resid_norm <= 1.0e-10 && (~preserve_mass || mass_resid_norm <= 1.0e-12)
                    break;
                end

                if preserve_mass
                    rhs_aug = [residual; mass_resid];
                    delta_aug = obj.fixed_collocation_mass_augmented \ rhs_aug;
                    delta = delta_aug(1:size(obj.output_nodes, 1), :);
                else
                    delta = obj.fixed_collocation \ residual;
                end
                nodal_current = nodal_current + delta;
                obj.solve_stats.defect_correction_solves = obj.solve_stats.defect_correction_solves + 1;
                if iter == 1
                    obj.solve_stats.moved_solves_one_defect = obj.solve_stats.moved_solves_one_defect + 1;
                elseif iter == 2
                    obj.solve_stats.moved_solves_two_defect = obj.solve_stats.moved_solves_two_defect + 1;
                end
            end

            coeffs_next = obj.projectSamples(nodal_current);
            if preserve_mass
                coeffs_next = enforceMassCorrection(obj, coeffs_next, target_mass);
            end
        end

        function setStepSize(obj, dlt)
            obj.dt = dlt;
        end

        function resetSolveStats(obj)
            obj.solve_stats = defaultSolveStats();
        end

        function out = getSolveStats(obj)
            out = obj.solve_stats;
        end

        function out = getNumPatches(obj)
            out = size(obj.patch_centers, 1);
        end

        function out = getNumDofs(obj)
            out = obj.total_coeffs;
        end

        function out = getPatchRadius(obj)
            out = obj.patch_radius;
        end

        function out = getPatchRadii(obj)
            out = obj.patch_radii;
        end

        function out = getPatchSpacing(obj)
            out = obj.patch_spacing;
        end

        function out = getLocalizedPatchData(obj)
            out = struct( ...
                'centers', obj.patch_centers, ...
                'radii', obj.patch_radii, ...
                'spacing', obj.patch_spacing, ...
                'node_ids', {obj.patch_node_ids}, ...
                'stencil_properties', obj.patch_stencil_props, ...
                'coeff_offsets', obj.coeff_offsets);
        end

        function out = totalMass(obj, coeffs, col)
            if nargin < 3 || isempty(col)
                col = 1;
            end
            out = 0.0;
            for p = 1:numel(obj.patch_mass_rows)
                if isempty(obj.patch_mass_rows{p})
                    continue;
                end
                out = out + obj.patch_mass_rows{p}.' * coeffs(localCoeffRange(obj, p), col);
            end
        end

        function out = getDomainMeasure(obj)
            out = obj.domain_measure;
        end
    end
end

function values = sampleScalarFunctionAtOutputNodes(obj, rho0)
X = obj.output_nodes;
values = zeros(size(X, 1), 1);
for i = 1:size(X, 1)
    values(i, 1) = rho0(X(i, :));
end
end

function tuning = autoTunePatchGeometry(Xi, Xb, level_set, xi, h, patch_radius_factor, patch_spacing_factor)
X = [Xi; Xb];
dim = size(X, 2);
target_points = localizedTargetPoints(xi, dim);
if size(X, 1) >= 2
    target_points = min(target_points, max(1, floor(size(X, 1) / 2)));
end
eligible = interiorSafeCenters(Xi, Xb, h);
tuning = farthestPointPatchGeometry(X, eligible, Xb, target_points, h, ...
    level_set, patch_radius_factor, patch_spacing_factor);
end

function [centers, patch_radii, patch_node_ids] = buildPatchData(Xi, Xb, tuning)
X = [Xi; Xb];
center_ids = tuning.center_ids(:);
if ~isempty(center_ids)
    [center_ids, order] = sort(center_ids);
    centers = X(center_ids, :);
    patch_radii = tuning.center_radii(order);
else
    centers = vertcat(tuning.center_positions{:});
    patch_radii = tuning.center_radii(:);
end
if isempty(centers)
    patch_radii = zeros(0, 1);
    patch_node_ids = {};
    return;
end
patch_node_ids = cell(size(centers, 1), 1);
for p = 1:size(centers, 1)
    d = sqrt(sum((X - centers(p, :)).^2, 2));
    patch_node_ids{p} = find(d <= patch_radii(p) * (1.0 + 1.0e-12));
end
end

function [avg_nnz, max_nnz] = estimateCollocationSparsity(X, centers, node_ids, radius)
if isempty(X) || isempty(centers)
    avg_nnz = 0;
    max_nnz = 0;
    return;
end
total_nnz = 0;
max_nnz = 0;
for i = 1:size(X, 1)
    marker = false(size(X, 1), 1);
    row_nnz = 0;
    for p = 1:size(centers, 1)
        if norm(X(i, :) - centers(p, :), 2) > radius
            continue;
        end
        ids = node_ids{p};
        new_ids = ids(~marker(ids));
        marker(new_ids) = true;
        row_nnz = row_nnz + numel(new_ids);
    end
    total_nnz = total_nnz + row_nnz;
    max_nnz = max(max_nnz, row_nnz);
end
avg_nnz = total_nnz / size(X, 1);
end

function min_nodes = chooseMinimumPatchNodes(dim, xi)
ell = max(xi + 1, 2);
npoly = size(kp.poly.total_degree_indices(dim, ell), 1);
min_nodes = 2 * npoly + 1;
end

function n = localizedTargetPoints(xi, dim)
ell = max(xi + 1, 2);
npoly = size(kp.poly.total_degree_indices(dim, ell), 1);
n = max(ceil(1.5 * npoly) + 1, npoly + 1);
end

function eligible = interiorSafeCenters(Xi, Xb, h)
if isempty(Xi)
    eligible = zeros(0, 1);
    return;
end
clearance_levels = [1.25, 1.0, 0.75, 0.5];
for fac = clearance_levels
    clearance = fac * h;
    keep = false(size(Xi, 1), 1);
    for i = 1:size(Xi, 1)
        if nearestBoundaryDistance(Xb, Xi(i, :)) >= clearance
            keep(i) = true;
        end
    end
    eligible = find(keep);
    if ~isempty(eligible)
        return;
    end
end
eligible = (1:size(Xi, 1)).';
end

function dist = nearestBoundaryDistance(Xb, x)
if isempty(Xb)
    dist = inf;
    return;
end
d = sqrt(sum((Xb - x).^2, 2));
dist = min(d);
end

function tuning = farthestPointPatchGeometry(X, eligible_centers, Xb, min_points, h, level_set, radius_factor, spacing_factor)
tuning = struct('radius', h, 'spacing', h, 'target_points', 0, 'max_target_points', 0, ...
    'center_ids', [], 'center_positions', {{}}, 'center_radii', [], ...
    'min_patch_nodes', 0, 'max_patch_nodes', 0, 'num_centers', 0, ...
    'avg_patch_nodes', 0.0, 'avg_overlap', 0.0, 'max_overlap', 0.0, 'objective', 0.0);
if isempty(X)
    return;
end

tuning.target_points = min(max(min_points, 1), size(X, 1));
tuning.max_target_points = min(size(X, 1), max(floor(3 * tuning.target_points / 2) + 1, tuning.target_points + 1));
separation_factor = 0.9;

if isempty(eligible_centers)
    eligible = (1:size(X, 1)).';
else
    eligible = eligible_centers(:);
end

candidate_core = zeros(numel(eligible), 1);
candidate_patch = zeros(numel(eligible), 1);
for k = 1:numel(eligible)
    idx = eligible(k);
    [candidate_core(k), candidate_patch(k)] = localizedCoreAndPatchRadiusForCenter( ...
        X, Xb, idx, tuning.target_points, tuning.max_target_points, level_set, h, radius_factor);
end

active = true(numel(eligible), 1);
seed_idx = pointNearestCentroid(X, eligible);
seed_k = find(eligible == seed_idx, 1, 'first');
if isempty(seed_k)
    seed_k = 1;
end

core_radii = add_center(seed_k);
suppress_nearby(seed_k);
core_cover = false(size(X, 1), 1);
update_core_cover(size(tuning.center_positions, 1), core_radii(end));

[~, order] = sort(candidate_core, 'descend');
for pos = 1:numel(order)
    ek = order(pos);
    if ~active(ek)
        continue;
    end
    core_radii(end + 1, 1) = add_center(ek); %#ok<AGROW>
    suppress_nearby(ek);
    update_core_cover(size(tuning.center_positions, 1), core_radii(end));
end

while any(~core_cover)
    next_idx = 0;
    worst_dist = -inf;
    for i = 1:size(X, 1)
        if core_cover(i)
            continue;
        end
        nearest = inf;
        for j = 1:size(tuning.center_positions, 1)
            nearest = min(nearest, norm(X(i, :) - tuning.center_positions{j}, 2) / core_radii(j));
        end
        if nearest > worst_dist
            worst_dist = nearest;
            next_idx = i;
        end
    end
    if next_idx == 0
        break;
    end
    d = sqrt(sum((X(eligible, :) - X(next_idx, :)).^2, 2));
    [~, best_k] = min(d);
    if any(tuning.center_ids == eligible(best_k))
        core_cover(next_idx) = true;
        continue;
    end
    core_radii(end + 1, 1) = add_center(best_k); %#ok<AGROW>
    suppress_nearby(best_k);
    update_core_cover(size(tuning.center_positions, 1), core_radii(end));
end

j = numel(tuning.center_ids);
while j >= 1
    cj = tuning.center_positions{j};
    rj = core_radii(j);
    removable = true;
    for i = 1:size(X, 1)
        if norm(X(i, :) - cj, 2) > rj * (1.0 + 1.0e-12)
            continue;
        end
        covered_by_other = false;
        for k = 1:numel(tuning.center_positions)
            if k == j
                continue;
            end
            if norm(X(i, :) - tuning.center_positions{k}, 2) <= core_radii(k) * (1.0 + 1.0e-12)
                covered_by_other = true;
                break;
            end
        end
        if ~covered_by_other
            removable = false;
            break;
        end
    end
    if removable
        tuning.center_ids(j) = [];
        tuning.center_positions(j) = [];
        tuning.center_radii(j) = [];
        core_radii(j) = [];
    end
    j = j - 1;
end

tuning.center_radii = enforceWeightFloorSupport(X, tuning.center_positions, tuning.center_radii);
tuning.num_centers = numel(tuning.center_ids);
if isempty(tuning.center_radii)
    tuning.radius = h;
else
    tuning.radius = max(tuning.center_radii);
end
if radius_factor > 0
    tuning.center_radii = tuning.radius * ones(size(tuning.center_radii));
end
if spacing_factor > 0
    tuning.spacing = spacing_factor * h;
else
    tuning.spacing = averageNearestCenterDistancePositions(tuning.center_positions);
end

counts = occupancyCountsForCenterPositions(X, tuning.center_positions, tuning.center_radii);
if isempty(counts)
    tuning.min_patch_nodes = 0;
    tuning.max_patch_nodes = 0;
    tuning.avg_patch_nodes = 0.0;
    tuning.avg_overlap = 0.0;
    tuning.max_overlap = 0.0;
    tuning.objective = 0.0;
else
    tuning.min_patch_nodes = min(counts);
    tuning.max_patch_nodes = max(counts);
    tuning.avg_patch_nodes = mean(counts);
    cover = zeros(size(X, 1), 1);
    for i = 1:numel(tuning.center_positions)
        d = sqrt(sum((X - tuning.center_positions{i}).^2, 2));
        ids = find(d <= tuning.center_radii(i) * (1.0 + 1.0e-12));
        cover(ids) = cover(ids) + 1;
    end
    tuning.avg_overlap = mean(cover);
    tuning.max_overlap = max(cover);
    tuning.objective = tuning.avg_patch_nodes + tuning.avg_overlap;
end

    function rcore = add_center(ek)
        idx = eligible(ek);
        tuning.center_ids(end + 1, 1) = idx; %#ok<AGROW>
        tuning.center_positions{end + 1, 1} = X(idx, :); %#ok<AGROW>
        tuning.center_radii(end + 1, 1) = candidate_patch(ek); %#ok<AGROW>
        rcore = candidate_core(ek);
    end

    function suppress_nearby(chosen_k)
        active(chosen_k) = false;
        for ek = 1:numel(eligible)
            if ~active(ek)
                continue;
            end
            sep = separation_factor * min(candidate_core(ek), candidate_core(chosen_k));
            if norm(X(eligible(ek), :) - X(eligible(chosen_k), :), 2) < sep
                active(ek) = false;
            end
        end
    end

    function update_core_cover(center_idx, rc)
        c = tuning.center_positions{center_idx};
        for ii = 1:size(X, 1)
            if norm(X(ii, :) - c, 2) <= rc * (1.0 + 1.0e-12)
                core_cover(ii) = true;
            end
        end
    end
end

function idx = pointNearestCentroid(X, eligible)
centroid = mean(X, 1);
idx = eligible(1);
best = inf;
for k = 1:numel(eligible)
    d = norm(X(eligible(k), :) - centroid, 2);
    if d < best
        best = d;
        idx = eligible(k);
    end
end
end

function radius = kthNeighborRadius(X, center_idx, min_points)
d = sqrt(sum((X - X(center_idx, :)).^2, 2));
[ds, ~] = sort(d, 'ascend');
kth = min(max(min_points, 1), numel(ds));
radius = max(ds(kth), 1.0e-12);
end

function radius = defaultRadiusField(X, Xb, center, level_set)
diameter = cloudDiameter(X);
rmin = 0.08 * diameter;
rmax = 0.20 * diameter;
if ~isempty(level_set)
    boundary_dist = abs(levelSetValue(level_set, center));
elseif ~isempty(Xb)
    boundary_dist = nearestBoundaryDistance(Xb, center);
else
    boundary_dist = inf;
end
if ~isfinite(boundary_dist)
    radius = rmin;
else
    radius = max(rmin, min(rmax, 0.35 * boundary_dist));
end
end

function val = levelSetValue(level_set, x)
val = level_set.Evaluate(x);
end

function diameter = cloudDiameter(X)
if size(X, 1) <= 1
    diameter = 1.0;
    return;
end
xmin = min(X, [], 1);
xmax = max(X, [], 1);
diameter = max(norm(xmax - xmin, 2), 1.0e-8);
end

function [r_core, patch_radius] = localizedCoreAndPatchRadiusForCenter(X, Xb, idx, min_points, max_points, level_set, h, radius_factor)
center = X(idx, :);
r_geom = defaultRadiusField(X, Xb, center, level_set);
r_kmin = kthNeighborRadius(X, idx, min_points);
r_kmax = kthNeighborRadius(X, idx, max_points);
r_core = max(1.0e-12, min(r_kmax, max(r_geom, r_kmin)));
if radius_factor > 0
    patch_radius = radius_factor * h;
else
    patch_radius = 1.25 * r_core;
end
if ~isempty(level_set)
    boundary_radius = boundaryLimitedRadius(level_set, center, patch_radius);
    patch_radius = max(r_core, min(patch_radius, boundary_radius));
elseif ~isempty(Xb)
    bdry_dist = nearestBoundaryDistance(Xb, center);
    if isfinite(bdry_dist)
        patch_radius = min(patch_radius, max(r_core, 1.05 * bdry_dist));
    end
end
patch_radius = max(r_core, patch_radius);
end

function radius = boundaryLimitedRadius(level_set, center, candidate_radius)
radius = candidate_radius;
if ~(candidate_radius > 0)
    return;
end
phi0 = levelSetValue(level_set, center);
if abs(phi0) <= 1.0e-14
    return;
end
inside_positive = phi0 > 0;
dirs = unitDirections(numel(center));
for k = 1:numel(dirs)
    x1 = center + candidate_radius * dirs{k};
    phi1 = levelSetValue(level_set, x1);
    if inside_positive
        crosses = (phi1 < 0);
    else
        crosses = (phi1 > 0);
    end
    if ~crosses
        continue;
    end
    lo = 0.0;
    hi = candidate_radius;
    for iter = 1:30
        mid = 0.5 * (lo + hi);
        xm = center + mid * dirs{k};
        phim = levelSetValue(level_set, xm);
        if inside_positive
            if phim >= 0
                lo = mid;
            else
                hi = mid;
            end
        else
            if phim <= 0
                lo = mid;
            else
                hi = mid;
            end
        end
    end
    radius = min(radius, hi);
end
end

function dirs = unitDirections(dim)
dirs = {};
if dim == 2
    n = 32;
    dirs = cell(n, 1);
    for k = 1:n
        theta = 2 * pi * (k - 1) / n;
        dirs{k} = [cos(theta), sin(theta)];
    end
elseif dim == 3
    n = 48;
    phi = 0.5 * (1.0 + sqrt(5.0));
    dirs = cell(n, 1);
    for k = 1:n
        z = 1.0 - 2.0 * ((k - 1) + 0.5) / n;
        r = sqrt(max(0.0, 1.0 - z * z));
        theta = 2.0 * pi * (k - 1) / phi;
        dirs{k} = [r * cos(theta), r * sin(theta), z];
    end
end
end

function radii = enforceWeightFloorSupport(X, centers, patch_radii, tau)
if nargin < 4 || isempty(tau)
    tau = 1.0e-10;
end
radii = patch_radii(:);
if isempty(X) || isempty(centers) || numel(radii) ~= numel(centers)
    return;
end
r_safe = max(radiusRatioForWeightFloor(tau), 1.0e-8);
for i = 1:size(X, 1)
    safe = false;
    best_k = 0;
    best_required = inf;
    for k = 1:numel(centers)
        dist = norm(X(i, :) - centers{k}, 2);
        if dist <= r_safe * radii(k)
            safe = true;
            break;
        end
        if dist <= radii(k) * (1.0 + 1.0e-12)
            required = dist / r_safe;
            if required < best_required
                best_required = required;
                best_k = k;
            end
        end
    end
    if ~safe && isfinite(best_required)
        radii(best_k) = max(radii(best_k), best_required);
    end
end
end

function ratio = radiusRatioForWeightFloor(tau)
if ~(tau > 0)
    ratio = 1.0;
    return;
end
if wendlandC6(0.0) <= tau
    ratio = 0.0;
    return;
end
if wendlandC6(1.0) >= tau
    ratio = 1.0;
    return;
end
lo = 0.0;
hi = 1.0;
for iter = 1:80
    mid = 0.5 * (lo + hi);
    if wendlandC6(mid) >= tau
        lo = mid;
    else
        hi = mid;
    end
end
ratio = lo;
end

function avg = averageNearestCenterDistancePositions(centers)
if numel(centers) <= 1
    avg = 0.0;
    return;
end
nearest_sum = 0.0;
for i = 1:numel(centers)
    nearest = inf;
    for j = 1:numel(centers)
        if i == j
            continue;
        end
        nearest = min(nearest, norm(centers{i} - centers{j}, 2));
    end
    nearest_sum = nearest_sum + nearest;
end
avg = nearest_sum / numel(centers);
end

function counts = occupancyCountsForCenterPositions(X, centers, radii)
counts = zeros(numel(centers), 1);
for c = 1:numel(centers)
    d = sqrt(sum((X - centers{c}).^2, 2));
    counts(c) = sum(d <= radii(c) * (1.0 + 1.0e-12));
end
end


function tree = buildCenterTree(centers)
tree = struct('Points', centers, 'Searcher', [], 'HasSearcher', false);
if isempty(centers)
    return;
end
if exist('KDTreeSearcher', 'class') == 8 && exist('rangesearch', 'file') == 2
    tree.Searcher = KDTreeSearcher(centers);
    tree.HasSearcher = true;
end
end

function sp = buildPatchStencilProperties(dim, xi)
sp = kp.rbffd.StencilProperties();
sp.dim = dim;
sp.ell = max(xi + 1, 2);
sp.npoly = size(kp.poly.total_degree_indices(dim, sp.ell), 1);
sp.spline_degree = sp.ell;
if mod(sp.spline_degree, 2) == 0
    sp.spline_degree = sp.spline_degree + 1;
end
end

function stencils = buildPatchStencils(Xnodes, node_ids, sp)
stencils = cell(size(node_ids));
for p = 1:numel(node_ids)
    nodes = Xnodes(node_ids{p}, :);
    r = kp.geometry.distanceMatrix(nodes, nodes);
    width = max(max(r, [], 'all'), 1);
    centroid = mean(nodes, 1);
    Xc = (nodes - centroid) / width;
    poly_indices = kp.poly.total_degree_indices(size(nodes, 2), sp.ell);
    P = monomialBasis(Xc, poly_indices);
    solve_lhs = zeros(size(nodes, 1) + sp.npoly, size(nodes, 1) + sp.npoly);
    solve_lhs(1:size(nodes, 1), 1:size(nodes, 1)) = phsRbf(r, sp.spline_degree);
    solve_lhs(1:size(nodes, 1), size(nodes, 1)+1:end) = P;
    solve_lhs(size(nodes, 1)+1:end, 1:size(nodes, 1)) = P.';
    [Q, R, perm] = qr(solve_lhs, 'vector');
    stencils{p} = struct( ...
        'nodes', nodes, ...
        'centroid', centroid, ...
        'width', width, ...
        'poly_indices', poly_indices, ...
        'solve_lhs', solve_lhs, ...
        'Q', Q, ...
        'R', R, ...
        'perm', perm, ...
        'rank', cpqrNumericalRank(R));
end
end

function offsets = buildCoeffOffsets(node_ids, npoly)
offsets = ones(numel(node_ids) + 1, 1);
for p = 1:numel(node_ids)
    offsets(p + 1) = offsets(p) + numel(node_ids{p}) + npoly;
end
end

function range = localCoeffRange(obj, patch_id)
range = obj.coeff_offsets(patch_id):(obj.coeff_offsets(patch_id + 1) - 1);
end

function values = localizedEvaluateCoefficients(obj, coeffs, Xq)
if isempty(Xq)
    values = zeros(0, size(coeffs, 2));
    return;
end

nq = size(Xq, 1);
nc = size(coeffs, 2);
values = zeros(nq, nc);
 patch_ids_per_query = queryPatchIds(obj.patch_center_tree, obj.patch_centers, obj.patch_radii, Xq);

for q = 1:nq
    patch_ids = patch_ids_per_query{q};
    center_dist = sqrt(sum((obj.patch_centers(patch_ids, :) - Xq(q, :)).^2, 2));
    alpha = wendlandC6(center_dist ./ obj.patch_radii(patch_ids));
    alpha_sum = sum(alpha);
    if alpha_sum <= 1.0e-14
        alpha = ones(size(alpha));
        alpha_sum = sum(alpha);
    end
    alpha = alpha / alpha_sum;

    for k = 1:numel(patch_ids)
        p = patch_ids(k);
        basis_row = localBasisRow(obj.patch_stencils{p}, obj.patch_stencil_props, Xq(q, :));
        values(q, :) = values(q, :) + alpha(k) * (basis_row * coeffs(localCoeffRange(obj, p), :));
    end
end
end

function values = localizedEvaluateNodal(obj, nodal_values, Xq)
if isempty(Xq)
    values = zeros(0, size(nodal_values, 2));
    return;
end
W = nodalEvaluationWeights(obj, Xq);
values = W * nodal_values;
end

function W = nodalEvaluationWeights(obj, Xq)
if isempty(Xq)
    W = sparse(0, size(obj.output_nodes, 1));
    return;
end

rows = cell(size(Xq, 1), 1);
cols = cell(size(Xq, 1), 1);
vals = cell(size(Xq, 1), 1);
patch_ids_per_query = queryPatchIds(obj.patch_center_tree, obj.patch_centers, obj.patch_radii, Xq);

for q = 1:size(Xq, 1)
    patch_ids = patch_ids_per_query{q};
    center_dist = sqrt(sum((obj.patch_centers(patch_ids, :) - Xq(q, :)).^2, 2));
    alpha = wendlandC6(center_dist ./ obj.patch_radii(patch_ids));
    alpha_sum = sum(alpha);
    if alpha_sum <= 1.0e-14
        alpha = ones(size(alpha));
        alpha_sum = sum(alpha);
    end
    alpha = alpha / alpha_sum;

    row = zeros(1, size(obj.output_nodes, 1));
    for k = 1:numel(patch_ids)
        p = patch_ids(k);
        ids = obj.patch_node_ids{p};
        basis_row = localBasisRow(obj.patch_stencils{p}, obj.patch_stencil_props, Xq(q, :));
        row(ids) = row(ids) + alpha(k) * (basis_row * obj.forward_cardinal_coeffs{p});
    end

    nz = find(abs(row) > 1.0e-14);
    rows{q} = q * ones(numel(nz), 1);
    cols{q} = nz(:);
    vals{q} = row(nz(:)).';
end

W = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), size(Xq, 1), size(obj.output_nodes, 1));
end

function patch_ids_per_query = queryPatchIds(tree, centers, patch_radii, Xq)
if isempty(centers)
    patch_ids_per_query = cell(size(Xq, 1), 1);
    return;
end
max_radius = max(patch_radii);
if tree.HasSearcher
    patch_ids_per_query = rangesearch(tree.Searcher, Xq, max_radius);
else
    D = kp.geometry.distanceMatrix(Xq, centers);
    patch_ids_per_query = cell(size(Xq, 1), 1);
    for q = 1:size(Xq, 1)
        patch_ids_per_query{q} = find(D(q, :) <= max_radius);
    end
end

for q = 1:numel(patch_ids_per_query)
    if isempty(patch_ids_per_query{q})
        d = sqrt(sum((centers - Xq(q, :)).^2, 2));
        [~, nearest_patch] = min(d);
        patch_ids_per_query{q} = nearest_patch;
    else
        ids = patch_ids_per_query{q}(:);
        d = sqrt(sum((centers(ids, :) - Xq(q, :)).^2, 2));
        keep = d <= patch_radii(ids) * (1.0 + 1.0e-12);
        patch_ids_per_query{q} = ids(keep).';
        if isempty(patch_ids_per_query{q})
            d = sqrt(sum((centers - Xq(q, :)).^2, 2));
            [~, nearest_patch] = min(d);
            patch_ids_per_query{q} = nearest_patch;
        end
    end
end
end

function basis_row = localBasisRow(stencil, sp, xq)
re = kp.geometry.distanceMatrix(xq, stencil.nodes);
Xec = (xq - stencil.centroid) / stencil.width;
Pe = monomialBasis(Xec, stencil.poly_indices);
basis_row = [phsRbf(re, sp.spline_degree), Pe];
end

function coeffs = solveLocalCoeffs(obj, patch_id, sample_sites, sample_values)
stencil = obj.patch_stencils{patch_id};
basis = localBasisMatrix(stencil, obj.patch_stencil_props, sample_sites);
Xnodes = obj.output_nodes(obj.patch_node_ids{patch_id}, :);
Xc = (Xnodes - stencil.centroid) / stencil.width;
P = monomialBasis(Xc, stencil.poly_indices);
lhs = [basis; [P.', zeros(obj.patch_stencil_props.npoly, obj.patch_stencil_props.npoly)]];
rhs = [sample_values; zeros(obj.patch_stencil_props.npoly, size(sample_values, 2))];
fixed_sites = isequal(size(sample_sites), size(Xnodes)) && norm(sample_sites - Xnodes, inf) <= 1.0e-14;

if fixed_sites
    obj.solve_stats.fixed_site_solves = obj.solve_stats.fixed_site_solves + 1;
    coeffs = solveWithStencilCPQR(stencil, rhs);
    return;
end

coeffs = solveWithStencilCPQR(stencil, rhs);
rhs_norm = max(norm(rhs, inf), 1.0);
corrections_used = 0;
prev_resid_norm = inf;
for iter = 1:4
    resid = rhs - lhs * coeffs;
    resid_norm = norm(resid, inf);
    if resid_norm <= 1.0e-10 * rhs_norm
        break;
    end
    if iter > 1 && isfinite(prev_resid_norm) && prev_resid_norm > 0
        obj.solve_stats.residual_ratio_sums(iter - 1) = obj.solve_stats.residual_ratio_sums(iter - 1) + resid_norm / prev_resid_norm;
        obj.solve_stats.residual_ratio_counts(iter - 1) = obj.solve_stats.residual_ratio_counts(iter - 1) + 1;
    end
    if iter > 1 && resid_norm > 0.7 * prev_resid_norm
        coeffs = solveWithCPQROnTheFly(lhs, rhs);
        obj.solve_stats.moved_solves_fallback = obj.solve_stats.moved_solves_fallback + 1;
        return;
    end
    delta = solveWithStencilCPQR(stencil, resid);
    coeffs = coeffs + delta;
    prev_resid_norm = resid_norm;
    corrections_used = corrections_used + 1;
    obj.solve_stats.defect_correction_solves = obj.solve_stats.defect_correction_solves + 1;
end

final_resid = rhs - lhs * coeffs;
if norm(final_resid, inf) <= 1.0e-8 * rhs_norm
    if corrections_used == 0
        obj.solve_stats.moved_solves_zero_defect = obj.solve_stats.moved_solves_zero_defect + 1;
    elseif corrections_used == 1
        obj.solve_stats.moved_solves_one_defect = obj.solve_stats.moved_solves_one_defect + 1;
    else
        obj.solve_stats.moved_solves_two_defect = obj.solve_stats.moved_solves_two_defect + 1;
    end
else
    coeffs = solveWithCPQROnTheFly(lhs, rhs);
    obj.solve_stats.moved_solves_fallback = obj.solve_stats.moved_solves_fallback + 1;
end
end

function rank = cpqrNumericalRank(R)
if isempty(R)
    rank = 0;
    return;
end
diagR = abs(diag(R));
if isempty(diagR) || ~(diagR(1) > 0)
    rank = 0;
    return;
end
tol = max(1.0e-12, 1.0e-10 * diagR(1));
rank = find(diagR > tol, 1, 'last');
if isempty(rank)
    rank = 0;
end
end

function X = solveWithStencilCPQR(stencil, rhs)
X = solveWithCPQR(stencil.Q, stencil.R, stencil.perm, stencil.rank, rhs);
if isempty(X)
    X = stableSolve(stencil.solve_lhs, rhs);
end
end

function X = solveWithCPQROnTheFly(lhs, rhs)
[Q, R, perm] = qr(lhs, 'vector');
rank = cpqrNumericalRank(R);
X = solveWithCPQR(Q, R, perm, rank, rhs);
if isempty(X)
    X = stableSolve(lhs, rhs);
end
end

function X = solveWithCPQR(Q, R, perm, rank, rhs)
ncols = size(R, 2);
X = zeros(ncols, size(rhs, 2));
if rank == 0
    return;
end
y = Q' * rhs;
z = R(1:rank, 1:rank) \ y(1:rank, :);
if any(~isfinite(z), 'all')
    X = [];
    return;
end
X(perm(1:rank), :) = z;
end

function basis = localBasisMatrix(stencil, sp, Xe)
if isempty(Xe)
    basis = zeros(0, size(stencil.nodes, 1) + sp.npoly);
    return;
end
re = kp.geometry.distanceMatrix(Xe, stencil.nodes);
Xec = (Xe - stencil.centroid) / stencil.width;
Pe = monomialBasis(Xec, stencil.poly_indices);
basis = [phsRbf(re, sp.spline_degree), Pe];
end

function P = monomialBasis(X, indices)
if isempty(X)
    P = zeros(0, size(indices, 1));
    return;
end
P = ones(size(X, 1), size(indices, 1));
for j = 1:size(indices, 1)
    for d = 1:size(X, 2)
        e = indices(j, d);
        if e ~= 0
            P(:, j) = P(:, j) .* (X(:, d) .^ e);
        end
    end
end
end

function [constant_mode, patch_mass_rows, domain_measure] = precomputePatchMassData(obj)
constant_mode = zeros(obj.total_coeffs, 1);
patch_mass_rows = cell(numel(obj.patch_node_ids), 1);
domain_measure = 0.0;
quad = patchQuadratureTemplate(size(obj.output_nodes, 2), max(4, obj.patch_stencil_props.ell + 2));

for p = 1:numel(obj.patch_node_ids)
    ids = obj.patch_node_ids{p};
    patch_nodes = obj.output_nodes(ids, :);
    ones_vals = ones(size(patch_nodes, 1), 1);
    constant_mode(localCoeffRange(obj, p)) = solveLocalCoeffs(obj, p, patch_nodes, ones_vals);

    mass_row = zeros(numel(localCoeffRange(obj, p)), 1);
    patch_measure = 0.0;
    center = obj.patch_centers(p, :);
    radius = obj.patch_radii(p);
    for q = 1:numel(quad.weights)
        x = center + radius * quad.points(q, :);
        if ~isInsideDomain(obj, x)
            continue;
        end
        overlaps = queryPatchIds(obj.patch_center_tree, obj.patch_centers, obj.patch_radii, x);
        overlaps = overlaps{1};
        w = normalizedPatchWeight(obj, p, overlaps, x);
        if w <= 0
            continue;
        end
        basis_row = localBasisRow(obj.patch_stencils{p}, obj.patch_stencil_props, x).';
        weight = quad.weights(q) * (radius ^ size(obj.output_nodes, 2)) * w;
        mass_row = mass_row + weight * basis_row;
        patch_measure = patch_measure + weight;
    end
    patch_mass_rows{p} = mass_row;
    domain_measure = domain_measure + patch_measure;
end
end

function [aug, nodal_mass_weights, local_constant_column] = buildForwardMassAugmentedData(obj)
n = size(obj.output_nodes, 1);
nodal_mass_weights = zeros(n, 1);
for p = 1:numel(obj.patch_node_ids)
    if isempty(obj.patch_mass_rows{p})
        continue;
    end
    ids = obj.patch_node_ids{p};
    local_weights = obj.forward_cardinal_coeffs{p}.' * obj.patch_mass_rows{p};
    nodal_mass_weights(ids) = nodal_mass_weights(ids) + local_weights;
end
local_constant_column = sum(obj.fixed_collocation, 2);
aug = [obj.fixed_collocation, local_constant_column; nodal_mass_weights.', 0];
end

function inside = isInsideDomain(obj, x)
if isempty(obj.Domain.getOuterLevelSet())
    inside = true;
    return;
end
inside = obj.Domain.getOuterLevelSet().Evaluate(x) >= -1.0e-12;
end

function inside = isInsideLevelSetPoint(obj, x)
if isempty(obj.Domain.getOuterLevelSet())
    inside = true;
    return;
end
inside = obj.Domain.getOuterLevelSet().Evaluate(x) >= -1.0e-12;
end

function w = normalizedPatchWeight(obj, patch_id, patch_ids, x)
denom = 0.0;
numer = 0.0;
for idx = patch_ids
    r = norm(x - obj.patch_centers(idx, :), 2) / obj.patch_radii(idx);
    psi = wendlandC6(r);
    denom = denom + psi;
    if idx == patch_id
        numer = psi;
    end
end
if denom <= 1.0e-14
    w = 0.0;
else
    w = numer / denom;
end
end

function quad = patchQuadratureTemplate(dim, order)
if dim == 2
    nrad = max(8, order + 2);
    nth = max(24, 2 * order + 8);
    r = ((1:nrad)' - 0.5) / nrad;
    th = 2 * pi * ((1:nth)' - 0.5) / nth;
    pts = zeros(nrad * nth, 2);
    w = zeros(nrad * nth, 1);
    idx = 1;
    for i = 1:nrad
        for j = 1:nth
            pts(idx, :) = [r(i) * cos(th(j)), r(i) * sin(th(j))];
            w(idx) = (2 * pi / nth) * (1 / nrad) * r(i);
            idx = idx + 1;
        end
    end
else
    nrad = max(6, order + 1);
    nth = max(12, 2 * order + 4);
    nph = max(24, 4 * order + 8);
    r = ((1:nrad)' - 0.5) / nrad;
    th = pi * ((1:nth)' - 0.5) / nth;
    ph = 2 * pi * ((1:nph)' - 0.5) / nph;
    pts = zeros(nrad * nth * nph, 3);
    w = zeros(nrad * nth * nph, 1);
    idx = 1;
    for i = 1:nrad
        for j = 1:nth
            for k = 1:nph
                st = sin(th(j));
                pts(idx, :) = [r(i) * st * cos(ph(k)), r(i) * st * sin(ph(k)), r(i) * cos(th(j))];
                w(idx) = (2 * pi / nph) * (pi / nth) * (1 / nrad) * (r(i)^2) * st;
                idx = idx + 1;
            end
        end
    end
end
quad = struct('points', pts, 'weights', w);
end

function Phi = phsRbf(r, degree)
Phi = r .^ degree;
if mod(degree, 2) == 0
    Phi = Phi .* log(r + 2e-16);
end
end

function Xdep = tracePointsBackward(X, t_arrival, dt, velocity, rk)
Xdep = tracePointsSigned(X, t_arrival, -dt, velocity, rk);
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

function validateBoundaryConditionConfiguration(obj)
bc = obj.boundary_condition;
if bc.mode == "periodic" && isempty(bc.periodic_patches)
    error('kp:solvers:PeriodicNeedsPatches', ...
        'PUSL periodic advection boundary condition requires periodic patch rules.');
end
if bc.mode == "periodic"
    error('kp:solvers:PeriodicNotYetSupported', ...
        ['PUSL periodic advection is not yet supported robustly. ' ...
         'The localized PU basis still needs true periodic patch support, not just wrapped traces.']);
end
if bc.mode == "inflow_dirichlet" && ~isa(bc.inflow_value, 'function_handle')
    error('kp:solvers:MissingInflowCallback', ...
        'PUSL inflow-Dirichlet boundary condition requires an inflow callback.');
end
end

function validateTangentialBoundaryFlow(obj, time, velocity)
Xb = obj.Domain.getBdryNodes();
if isempty(Xb)
    return;
end
nr = obj.Domain.getNrmls();
U = velocity(time, Xb);
normal_speed = abs(sum(U .* nr, 2));
if max(normal_speed) > obj.boundary_condition.normal_velocity_tolerance
    error('kp:solvers:TangentialFlowMismatch', ...
        'PUSL tangential-flow boundary condition is incompatible with the supplied velocity field.');
end
end

function traced = adjustForwardTraces(obj, traced, tn, velocity)
bc = obj.boundary_condition;
if bc.mode == "periodic"
    error('kp:solvers:PeriodicNotYetSupported', 'Periodic PU-SL transport is not yet supported in MATLAB.');
end

if isempty(obj.Domain.getOuterLevelSet())
    return;
end
phi = obj.Domain.getOuterLevelSet();
inside = phi.Evaluate(traced) >= -1.0e-12;
if all(inside)
    return;
end

if bc.mode == "tangential"
    error('kp:solvers:TangentialForwardExit', ...
        'PUSL forward SL exited the domain under a tangential-flow boundary condition.');
end
if bc.mode == "unspecified"
    error('kp:solvers:ForwardOutsideDomain', ...
        'PUSL forward SL traced outside the domain without an advection boundary condition.');
end

if bc.mode == "inflow_dirichlet"
    X = obj.output_nodes;
    traced(~inside, :) = X(~inside, :);
    return;
end

boundary_start = obj.Domain.getNumInteriorNodes() + 1;
if boundary_start > size(obj.output_nodes, 1)
    return;
end
Xb = obj.Domain.getBdryNodes();
nr = obj.Domain.getNrmls();
for local_row = boundary_start:size(obj.output_nodes, 1)
    boundary_idx = local_row - boundary_start + 1;
    bn = boundaryNormalVelocity(Xb(boundary_idx, :), nr(boundary_idx, :), tn + obj.dt, velocity);
    if bn < -bc.normal_velocity_tolerance
        traced(local_row, :) = obj.output_nodes(local_row, :);
    end
end
end

function values = applyBackwardBoundaryCondition(obj, Xdep, X, coeffs_old, tnext)
values = obj.evaluateAtPoints(coeffs_old, Xdep);
bc = obj.boundary_condition;
if bc.mode == "unspecified" || bc.mode == "tangential"
    return;
end
if bc.mode == "periodic"
    error('kp:solvers:PeriodicNotYetSupported', 'Periodic PU-SL transport is not yet supported in MATLAB.');
end
if isempty(obj.Domain.getOuterLevelSet())
    return;
end
phi = obj.Domain.getOuterLevelSet().Evaluate(Xdep);
outside = phi < -1.0e-12;
if ~any(outside)
    return;
end
if bc.mode == "inflow_dirichlet"
    outside_ids = find(outside);
    for k = 1:numel(outside_ids)
        i = outside_ids(k);
        hit = boundaryHitOnSegment(obj, X(i, :), Xdep(i, :));
        inflow = bc.inflow_value(tnext - hit.parameter * obj.dt, hit.point);
        values(i, :) = inflow(1, :);
    end
end
end

function hit = boundaryHitOnSegment(obj, inside_point, outside_point)
phi = obj.Domain.getOuterLevelSet();
opts = struct('valueTolerance', 1.0e-12, ...
    'stepTolerance', 1.0e-12, ...
    'gradientTolerance', 1.0e-14, ...
    'maxIterations', 30, ...
    'maxStepNorm', inf);
hit = phi.FindSegmentRootNewton(inside_point, outside_point, 0.5, opts);
if ~hit.converged
    error('kp:solvers:BoundaryTraceFailed', ...
        'PUSL boundary tracing failed to locate a boundary hit on the level set.');
end
end

function bn = boundaryNormalVelocity(point, normal, time, velocity)
X = point;
U = velocity(time, X);
bn = dot(U(1, :), normal);
end

function corrected = enforceMassCorrection(obj, coeffs, target_mass)
corrected = coeffs;
for col = 1:size(coeffs, 2)
    current = obj.totalMass(corrected, col);
    alpha = (target_mass(col) - current) / max(abs(obj.totalMass(obj.constant_mode, 1)), 1.0e-14);
    corrected(:, col) = corrected(:, col) + alpha * obj.constant_mode;
end
end

function masses = totalMasses(obj, coeffs)
masses = zeros(1, size(coeffs, 2));
for col = 1:size(coeffs, 2)
    masses(col) = obj.totalMass(coeffs, col);
end
end

function stats = defaultSolveStats()
stats = struct( ...
    'fixed_site_solves', 0, ...
    'moved_solves_zero_defect', 0, ...
    'moved_solves_one_defect', 0, ...
    'moved_solves_two_defect', 0, ...
    'moved_solves_fallback', 0, ...
    'defect_correction_solves', 0, ...
    'residual_ratio_sums', zeros(1, 5), ...
    'residual_ratio_counts', zeros(1, 5));
end

function w = wendlandC6(r)
w = zeros(size(r));
mask = r < 1;
t = 1 - r(mask);
rm = r(mask);
w(mask) = t.^8 .* (32 * rm.^3 + 25 * rm.^2 + 8 * rm + 1);
end

function X = stableSolve(A, B)
warnNear = warning('query', 'MATLAB:nearlySingularMatrix');
warnSing = warning('query', 'MATLAB:singularMatrix');
cleanupObj = onCleanup(@() restoreWarnings(warnNear, warnSing)); %#ok<NASGU>
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
X = A \ B;
if any(~isfinite(X), 'all')
    X = pinv(A) * B;
end
X(~isfinite(X)) = 0;
end

function restoreWarnings(warnNear, warnSing)
warning(warnNear.state, 'MATLAB:nearlySingularMatrix');
warning(warnSing.state, 'MATLAB:singularMatrix');
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
frac = mean(vals >= 0);
measure = frac * prod(xmax - xmin);
end
