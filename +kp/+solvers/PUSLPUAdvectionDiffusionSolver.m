classdef PUSLPUAdvectionDiffusionSolver < handle
    %PUSLPUADVECTIONDIFFUSIONSOLVER PU-SL transport coupled with PU diffusion.

    properties (SetAccess = private)
        Domain kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
        xi_sl (1,1) double = 0
        xi_pu (1,1) double = 0
        dt (1,1) double = 0
        nu (1,1) double = 0
        direction string = "backward"
        advection kp.solvers.PUSLAdvectionSolver = kp.solvers.PUSLAdvectionSolver()
        diffusion kp.solvers.PUDiffusionSolver = kp.solvers.PUDiffusionSolver()
        state_nm2 double = zeros(0, 1)
        state_nm1 double = zeros(0, 1)
        state_n double = zeros(0, 1)
        completed_steps (1,1) double = 0
        enforce_mass_constraint (1,1) logical = false
        mass_constraint_target (1,1) double = 0
        has_explicit_mass_constraint_target (1,1) logical = false
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

            % This wrapper mirrors the FD-coupled version, but swaps in PU
            % diffusion for the elliptic part.
            obj.advection.init(domain, xi_sl, dt);
            obj.diffusion.init(domain, xi_pu, dt, nu, num_omp_threads);

            obj.state_nm2 = zeros(0, 1);
            obj.state_nm1 = zeros(0, 1);
            obj.state_n = zeros(0, 1);
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
                if obj.completed_steps >= 2 && ~isempty(obj.state_n)
                    obj.mass_constraint_target = obj.totalMass(obj.state_n);
                elseif obj.completed_steps >= 1 && ~isempty(obj.state_nm1)
                    obj.mass_constraint_target = obj.totalMass(obj.state_nm1);
                elseif ~isempty(obj.state_nm2)
                    obj.mass_constraint_target = obj.totalMass(obj.state_nm2);
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
            obj.state_nm2 = state(:);
            obj.state_nm1 = zeros(0, 1);
            obj.state_n = zeros(0, 1);
            obj.completed_steps = 0;
            initializeMassConstraintTarget(obj, obj.state_nm2);
        end

        function setStateHistory(obj, varargin)
            switch numel(varargin)
                case 1
                    obj.setInitialState(varargin{1});
                case 2
                    obj.state_nm2 = varargin{1}(:);
                    obj.state_nm1 = varargin{2}(:);
                    obj.state_n = zeros(0, 1);
                    obj.completed_steps = 1;
                    initializeMassConstraintTarget(obj, obj.state_nm1);
                case 3
                    obj.state_nm2 = varargin{1}(:);
                    obj.state_nm1 = varargin{2}(:);
                    obj.state_n = varargin{3}(:);
                    obj.completed_steps = 2;
                    initializeMassConstraintTarget(obj, obj.state_n);
                otherwise
                    error('kp:solvers:BadStateHistory', 'Expected one, two, or three physical states.');
            end
        end

        function out = getOutputNodes(obj)
            out = obj.advection.getOutputNodes();
        end

        function out = getOutputRange(obj)
            out = obj.advection.getOutputRange();
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

        function next_state = bdf1Step(obj, t_next, velocity, rk, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            assert(~isempty(obj.state_nm2), 'PUSLPUAdvectionDiffusionSolver::bdf1Step requires setInitialState first.');
            % Transport first, then diffuse the transported state with the
            % corresponding PU diffusion step.
            transported_nm2 = transportState(obj, t_next - obj.dt, obj.state_nm2, 1, velocity, rk);
            obj.diffusion.setStateHistory(transported_nm2);
            next_state = obj.diffusion.bdf1Step(t_next, forcing, NeuCoeffFunc, DirCoeffFunc, bc);
            if obj.enforce_mass_constraint
                forcing_mass = totalMass(obj, evaluateForcing(obj, forcing, t_next));
                target = totalMass(obj, transported_nm2) + obj.dt * forcing_mass;
                next_state = enforceMassTarget(obj, next_state, target);
            end
            obj.state_nm1 = next_state;
            obj.completed_steps = 1;
        end

        function next_state = bdf2Step(obj, t_next, velocity, rk, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            assert(obj.completed_steps >= 1 && ~isempty(obj.state_nm1), ...
                'PUSLPUAdvectionDiffusionSolver::bdf2Step requires one prior step.');
            velocity_fn = velocity;
            transported_nm2 = transportState(obj, t_next - 2 * obj.dt, obj.state_nm2, 2, velocity_fn, rk);
            transported_nm1 = transportState(obj, t_next - obj.dt, obj.state_nm1, 1, velocity_fn, rk);
            obj.diffusion.setStateHistory(transported_nm2, transported_nm1);
            next_state = obj.diffusion.bdf2Step(t_next, forcing, NeuCoeffFunc, DirCoeffFunc, bc);
            if obj.enforce_mass_constraint
                forcing_mass = totalMass(obj, evaluateForcing(obj, forcing, t_next));
                target = (4 * totalMass(obj, transported_nm1) - totalMass(obj, transported_nm2) + 2 * obj.dt * forcing_mass) / 3;
                next_state = enforceMassTarget(obj, next_state, target);
            end
            obj.state_n = next_state;
            obj.completed_steps = 2;
        end

        function next_state = bdf3Step(obj, t_next, velocity, rk, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            assert(obj.completed_steps >= 2 && ~isempty(obj.state_n), ...
                'PUSLPUAdvectionDiffusionSolver::bdf3Step requires two prior steps.');
            velocity_fn = velocity;
            transported_nm2 = transportState(obj, t_next - 3 * obj.dt, obj.state_nm2, 3, velocity_fn, rk);
            transported_nm1 = transportState(obj, t_next - 2 * obj.dt, obj.state_nm1, 2, velocity_fn, rk);
            transported_n = transportState(obj, t_next - obj.dt, obj.state_n, 1, velocity_fn, rk);
            obj.diffusion.setStateHistory(transported_nm2, transported_nm1, transported_n);
            next_state = obj.diffusion.bdf3Step(t_next, forcing, NeuCoeffFunc, DirCoeffFunc, bc);
            if obj.enforce_mass_constraint
                forcing_mass = totalMass(obj, evaluateForcing(obj, forcing, t_next));
                target = (18 * totalMass(obj, transported_n) - 9 * totalMass(obj, transported_nm1) + ...
                    2 * totalMass(obj, transported_nm2) + 6 * obj.dt * forcing_mass) / 11;
                next_state = enforceMassTarget(obj, next_state, target);
            end
            obj.state_nm2 = obj.state_nm1;
            obj.state_nm1 = obj.state_n;
            obj.state_n = next_state;
            obj.completed_steps = max(obj.completed_steps, 3);
        end

        function out = currentState(obj)
            if obj.completed_steps <= 0
                out = obj.state_nm2;
            elseif obj.completed_steps == 1
                out = obj.state_nm1;
            else
                out = obj.state_n;
            end
        end

        function out = totalMass(obj, local_state)
            out = obj.advection.totalMass(local_state(:), 1);
        end
    end
end

function transported = transportState(obj, t_start, state, num_steps, velocity, rk)
% Reuse the PU-SL advection solver over a longer time horizon when older
% BDF states need to be transported to the current step.
original_dt = obj.dt;
horizon = num_steps * obj.dt;
obj.advection.setStepSize(horizon);
if obj.direction == "forward"
    coeffs_next = obj.advection.forwardSLStep(t_start, state(:), velocity, rk);
else
    coeffs_next = obj.advection.backwardSLStep(t_start, state(:), velocity, rk);
end
obj.advection.setStepSize(original_dt);
transported = coeffs_next(:);
if obj.enforce_mass_constraint
    transported = enforceMassTarget(obj, transported, massConstraintTarget(obj));
end
end

function values = evaluateForcing(obj, forcing, t)
X = obj.advection.getOutputNodes();
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

function initializeMassConstraintTarget(obj, reference_state)
if ~obj.enforce_mass_constraint || obj.has_explicit_mass_constraint_target
    return;
end
obj.mass_constraint_target = obj.totalMass(reference_state);
end

function target = massConstraintTarget(obj)
target = obj.mass_constraint_target;
end

function corrected = enforceMassTarget(obj, state, target_mass)
current_mass = obj.totalMass(state);
alpha = (target_mass - current_mass) / max(obj.advection.getDomainMeasure(), 1.0e-14);
corrected = state + alpha;
end
