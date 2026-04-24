classdef PUSLPUAdvectionDiffusionReactionSolver < handle
    %PUSLPUADVECTIONDIFFUSIONREACTIONSOLVER PU-SL transport + PU diffusion + explicit reaction.

    properties (SetAccess = private)
        Domain kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
        xi_sl (1,1) double = 0
        xi_pu (1,1) double = 0
        dt (1,1) double = 0
        nu (1,1) double = 0
        direction string = "backward"
        advection kp.solvers.PUSLAdvectionSolver = kp.solvers.PUSLAdvectionSolver()
        diffusion kp.solvers.PUDiffusionSolver = kp.solvers.PUDiffusionSolver()
        coeff_nm2 double = zeros(0, 0)
        coeff_nm1 double = zeros(0, 0)
        coeff_n double = zeros(0, 0)
        completed_steps (1,1) double = 0
        enforce_mass_constraint (1,1) logical = false
        mass_constraint_target (1,1) double = 0
        has_explicit_mass_constraint_target (1,1) logical = false
        constant_mode double = zeros(0, 0)
        constant_mode_mass (1,1) double = 0
        local_begin (1,1) double = 1
        local_end (1,1) double = 0
        physical_nodes double = zeros(0, 0)
    end

    methods
        function init(obj, domain, xi_sl, xi_pu, dt, nu, direction, num_omp_threads)
            if nargin < 7 || isempty(direction)
                direction = "backward";
            end
            if nargin < 8 || isempty(num_omp_threads)
                num_omp_threads = 1;
            end
            obj.Domain = domain;
            obj.xi_sl = xi_sl;
            obj.xi_pu = xi_pu;
            obj.dt = dt;
            obj.nu = nu;
            obj.direction = lower(string(direction));

            obj.advection.init(domain, xi_sl, dt);
            obj.diffusion.init(domain, xi_pu, dt, nu, num_omp_threads);

            output_range = obj.advection.getOutputRange();
            obj.local_begin = output_range(1);
            obj.local_end = output_range(2);
            obj.physical_nodes = obj.advection.getOutputNodes();
            requireMatchingNodes(obj.physical_nodes, obj.diffusion.getOutputNodes(), ...
                'PUSLPUAdvectionDiffusionReactionSolver::init');
            obj.constant_mode = obj.advection.projectConstant(1.0, 1);
            obj.constant_mode_mass = obj.advection.totalMass(obj.constant_mode, 1);

            obj.coeff_nm2 = zeros(0, 0);
            obj.coeff_nm1 = zeros(0, 0);
            obj.coeff_n = zeros(0, 0);
            obj.completed_steps = 0;
        end

        function setStepSize(obj, dt)
            obj.dt = dt;
            obj.advection.setStepSize(dt);
            obj.diffusion.setStepSize(dt);
        end

        function setTransportDirection(obj, direction)
            obj.direction = lower(string(direction));
        end

        function enableHomogeneousNeumannMassConservation(obj, enable)
            obj.enableMassConstraint(enable);
        end

        function enableMassConstraint(obj, enable)
            if nargin < 2
                enable = true;
            end
            obj.enforce_mass_constraint = logical(enable);
            if enable && ~obj.has_explicit_mass_constraint_target
                if obj.completed_steps >= 2 && ~isempty(obj.coeff_n)
                    obj.mass_constraint_target = massOfCoefficients(obj, obj.coeff_n);
                elseif obj.completed_steps >= 1 && ~isempty(obj.coeff_nm1)
                    obj.mass_constraint_target = massOfCoefficients(obj, obj.coeff_nm1);
                elseif ~isempty(obj.coeff_nm2)
                    obj.mass_constraint_target = massOfCoefficients(obj, obj.coeff_nm2);
                end
            end
        end

        function disableMassConstraint(obj)
            obj.enforce_mass_constraint = false;
        end

        function setMassConstraintTarget(obj, target_mass)
            obj.mass_constraint_target = target_mass;
            obj.has_explicit_mass_constraint_target = true;
            obj.enforce_mass_constraint = true;
        end

        function setInitialState(obj, state)
            local_state = extractLocalPhysicalState(obj, state);
            obj.coeff_nm2 = obj.advection.projectSamples(local_state(:));
            obj.coeff_nm1 = zeros(0, 0);
            obj.coeff_n = zeros(0, 0);
            obj.completed_steps = 0;
            initializeMassConstraintTarget(obj, obj.coeff_nm2);
        end

        function setStateHistory(obj, varargin)
            switch numel(varargin)
                case 1
                    obj.setInitialState(varargin{1});
                case 2
                    obj.coeff_nm2 = obj.advection.projectSamples(extractLocalPhysicalState(obj, varargin{1}));
                    obj.coeff_nm1 = obj.advection.projectSamples(extractLocalPhysicalState(obj, varargin{2}));
                    obj.coeff_n = zeros(0, 0);
                    obj.completed_steps = 1;
                    initializeMassConstraintTarget(obj, obj.coeff_nm1);
                case 3
                    obj.coeff_nm2 = obj.advection.projectSamples(extractLocalPhysicalState(obj, varargin{1}));
                    obj.coeff_nm1 = obj.advection.projectSamples(extractLocalPhysicalState(obj, varargin{2}));
                    obj.coeff_n = obj.advection.projectSamples(extractLocalPhysicalState(obj, varargin{3}));
                    obj.completed_steps = 2;
                    initializeMassConstraintTarget(obj, obj.coeff_n);
                otherwise
                    error('kp:solvers:BadStateHistory', 'Expected one, two, or three physical states.');
            end
        end

        function out = getOutputNodes(obj)
            out = obj.physical_nodes;
        end

        function out = getOutputRange(obj)
            out = [obj.local_begin, obj.local_end];
        end

        function out = returnsDistributedState(obj)
            out = obj.advection.returnsDistributedState();
        end

        function out = advectionSolver(obj)
            out = obj.advection;
        end

        function out = diffusionSolver(obj)
            out = obj.diffusion;
        end

        function next_state = bdf1Step(obj, t_next, velocity, rk, reaction, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            assert(~isempty(obj.coeff_nm2), 'PUSLPUAdvectionDiffusionReactionSolver::bdf1Step requires setInitialState first.');
            transported_coeff_nm2 = transportCoefficients(obj, t_next - obj.dt, obj.coeff_nm2, 1, velocity, rk);
            transported_nm2 = evaluateSingleState(obj, transported_coeff_nm2);
            reaction_nm2 = evaluateReactionCallback(reaction, t_next, transported_nm2, obj.physical_nodes);
            reaction_extrapolated = extrapolateReactionBDF1(reaction_nm2);
            obj.diffusion.setStateHistory(transported_nm2);
            next_state = obj.diffusion.bdf1Step(t_next, combineForcing(obj, forcing, reaction_extrapolated), ...
                NeuCoeffFunc, DirCoeffFunc, bc);
            coeff_next = obj.advection.projectSamples(next_state(:));
            enforceCoefficientMassTargetIfNeeded(obj, coeff_next);
            obj.coeff_nm1 = coeff_next;
            obj.completed_steps = 1;
            next_state = currentState(obj);
        end

        function next_state = bdf2Step(obj, t_next, velocity, rk, reaction, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            assert(obj.completed_steps >= 1 && ~isempty(obj.coeff_nm1), ...
                'PUSLPUAdvectionDiffusionReactionSolver::bdf2Step requires one prior step.');
            velocity_fn = velocity;
            transported_coeff_nm2 = transportCoefficients(obj, t_next - 2 * obj.dt, obj.coeff_nm2, 2, velocity_fn, rk);
            transported_coeff_nm1 = transportCoefficients(obj, t_next - obj.dt, obj.coeff_nm1, 1, velocity_fn, rk);
            transported_nm2 = evaluateSingleState(obj, transported_coeff_nm2);
            transported_nm1 = evaluateSingleState(obj, transported_coeff_nm1);
            reaction_nm2 = evaluateReactionCallback(reaction, t_next, transported_nm2, obj.physical_nodes);
            reaction_nm1 = evaluateReactionCallback(reaction, t_next, transported_nm1, obj.physical_nodes);
            reaction_extrapolated = extrapolateReactionBDF2(reaction_nm2, reaction_nm1);
            obj.diffusion.setStateHistory(transported_nm2, transported_nm1);
            next_state = obj.diffusion.bdf2Step(t_next, combineForcing(obj, forcing, reaction_extrapolated), ...
                NeuCoeffFunc, DirCoeffFunc, bc);
            coeff_next = obj.advection.projectSamples(next_state(:));
            enforceCoefficientMassTargetIfNeeded(obj, coeff_next);
            obj.coeff_n = coeff_next;
            obj.completed_steps = 2;
            next_state = currentState(obj);
        end

        function next_state = bdf3Step(obj, t_next, velocity, rk, reaction, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            assert(obj.completed_steps >= 2 && ~isempty(obj.coeff_n), ...
                'PUSLPUAdvectionDiffusionReactionSolver::bdf3Step requires two prior steps.');
            velocity_fn = velocity;
            transported_coeff_nm2 = transportCoefficients(obj, t_next - 3 * obj.dt, obj.coeff_nm2, 3, velocity_fn, rk);
            transported_coeff_nm1 = transportCoefficients(obj, t_next - 2 * obj.dt, obj.coeff_nm1, 2, velocity_fn, rk);
            transported_coeff_n = transportCoefficients(obj, t_next - obj.dt, obj.coeff_n, 1, velocity_fn, rk);
            transported_nm2 = evaluateSingleState(obj, transported_coeff_nm2);
            transported_nm1 = evaluateSingleState(obj, transported_coeff_nm1);
            transported_n = evaluateSingleState(obj, transported_coeff_n);
            reaction_nm2 = evaluateReactionCallback(reaction, t_next, transported_nm2, obj.physical_nodes);
            reaction_nm1 = evaluateReactionCallback(reaction, t_next, transported_nm1, obj.physical_nodes);
            reaction_n = evaluateReactionCallback(reaction, t_next, transported_n, obj.physical_nodes);
            reaction_extrapolated = extrapolateReactionBDF3(reaction_nm2, reaction_nm1, reaction_n);
            obj.diffusion.setStateHistory(transported_nm2, transported_nm1, transported_n);
            next_state = obj.diffusion.bdf3Step(t_next, combineForcing(obj, forcing, reaction_extrapolated), ...
                NeuCoeffFunc, DirCoeffFunc, bc);
            coeff_next = obj.advection.projectSamples(next_state(:));
            enforceCoefficientMassTargetIfNeeded(obj, coeff_next);
            obj.coeff_nm2 = obj.coeff_nm1;
            obj.coeff_nm1 = obj.coeff_n;
            obj.coeff_n = coeff_next;
            obj.completed_steps = max(obj.completed_steps, 3);
            next_state = currentState(obj);
        end

        function out = currentState(obj)
            coeffs = currentCoefficients(obj);
            if isempty(coeffs)
                out = zeros(0, 1);
                return;
            end
            out = evaluateSingleState(obj, coeffs);
        end

        function out = totalMass(obj, local_state)
            coeffs = obj.advection.projectSamples(extractLocalPhysicalState(obj, local_state));
            out = obj.advection.totalMass(coeffs, 1);
        end
    end
end

function out = extractLocalPhysicalState(obj, state)
state = state(:);
local_count = obj.local_end - obj.local_begin + 1;
if numel(state) == local_count
    out = state;
elseif obj.local_end <= numel(state)
    out = state(obj.local_begin:obj.local_end);
else
    error('kp:solvers:BadPhysicalStateSize', 'Received neither a local nor a full physical nodal state.');
end
end

function requireMatchingNodes(lhs, rhs, label)
if ~isequal(size(lhs), size(rhs)) || any(abs(lhs(:) - rhs(:)) > 1.0e-12)
    error('kp:solvers:NodeLayoutMismatch', '%s requires matching nodal coordinates.', label);
end
end

function transported = transportCoefficients(obj, t_start, coeff_state, num_steps, velocity, rk)
if num_steps <= 0
    error('kp:solvers:BadStepCount', 'PUSLPUAdvectionDiffusionReactionSolver::transportCoefficients requires a positive step count.');
end
original_dt = obj.dt;
horizon = num_steps * obj.dt;
obj.advection.setStepSize(horizon);
if obj.direction == "forward"
    coeffs_next = obj.advection.forwardSLStep(t_start, coeff_state, velocity, rk);
else
    coeffs_next = obj.advection.backwardSLStep(t_start, coeff_state, velocity, rk);
end
obj.advection.setStepSize(original_dt);
transported = coeffs_next;
enforceCoefficientMassTargetIfNeeded(obj, transported);
end

function out = evaluateSingleState(obj, coeffs)
out = obj.advection.evaluateAtPoints(coeffs, obj.physical_nodes);
if size(out, 2) ~= 1
    error('kp:solvers:BadStateShape', 'Expected a single-column transported state.');
end
out = out(:, 1);
end

function out = massOfCoefficients(obj, coeffs)
out = obj.advection.totalMass(coeffs, 1);
end

function initializeMassConstraintTarget(obj, reference_coeffs)
if ~obj.enforce_mass_constraint || obj.has_explicit_mass_constraint_target
    return;
end
obj.mass_constraint_target = massOfCoefficients(obj, reference_coeffs);
end

function target = massConstraintTarget(obj)
target = obj.mass_constraint_target;
end

function enforceCoefficientMassTargetIfNeeded(obj, coeffs)
if ~obj.enforce_mass_constraint
    return;
end
current_mass = massOfCoefficients(obj, coeffs);
alpha = (massConstraintTarget(obj) - current_mass) / max(abs(obj.constant_mode_mass), 1.0e-14);
coeffs(:, :) = coeffs + alpha * obj.constant_mode;
end

function coeffs = currentCoefficients(obj)
if obj.completed_steps <= 0
    coeffs = obj.coeff_nm2;
elseif obj.completed_steps == 1
    coeffs = obj.coeff_nm1;
else
    coeffs = obj.coeff_n;
end
end

function combined = combineForcing(obj, forcing, reaction_extrapolated)
combined = @(nuValue, t, X) evaluateCombinedForcing(obj, forcing, reaction_extrapolated, t, X); %#ok<NASGU>
end

function values = evaluateCombinedForcing(obj, forcing, reaction_extrapolated, t, X) %#ok<INUSD>
values = evaluateForcing(obj, forcing, t) + reaction_extrapolated(:);
end

function values = evaluateForcing(obj, forcing, t)
X = obj.physical_nodes;
try
    values = forcing(obj.nu, t, X);
catch
    try
        values = forcing(t, X);
    catch
        values = forcing(X);
    end
end
values = values(:);
end

function values = evaluateReactionCallback(reaction, t, state, X)
try
    values = reaction(t, state(:), X);
catch
    try
        values = reaction(state(:), X);
    catch
        values = reaction(t, X);
    end
end
values = values(:);
if numel(values) ~= size(X, 1)
    error('kp:solvers:BadReactionSize', 'Reaction values must match the physical node count.');
end
end

function values = extrapolateReactionBDF1(r_nm2)
values = r_nm2(:);
end

function values = extrapolateReactionBDF2(r_nm2, r_nm1)
values = 2 * r_nm1(:) - r_nm2(:);
end

function values = extrapolateReactionBDF3(r_nm2, r_nm1, r_n)
values = 3 * r_n(:) - 3 * r_nm1(:) + r_nm2(:);
end
