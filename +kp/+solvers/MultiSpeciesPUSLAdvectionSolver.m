classdef MultiSpeciesPUSLAdvectionSolver < handle
    %MULTISPECIESPUSLADVECTIONSOLVER Thin multi-species wrapper over PUSLAdvectionSolver.

    properties (SetAccess = private)
        num_species (1,1) double = 0
        solver kp.solvers.PUSLAdvectionSolver = kp.solvers.PUSLAdvectionSolver()
    end

    methods
        function clearAdvectionBoundaryCondition(obj)
            obj.solver.clearAdvectionBoundaryCondition();
        end

        function setTangentialFlowBoundary(obj, normal_velocity_tolerance)
            if nargin < 2
                normal_velocity_tolerance = [];
            end
            obj.solver.setTangentialFlowBoundary(normal_velocity_tolerance);
        end

        function setPeriodicBoundary(obj, periodic_patches, normal_velocity_tolerance)
            if nargin < 3
                normal_velocity_tolerance = [];
            end
            obj.solver.setPeriodicBoundary(periodic_patches, normal_velocity_tolerance);
        end

        function setInflowDirichletBoundary(obj, inflow_value, normal_velocity_tolerance)
            if nargin < 3
                normal_velocity_tolerance = [];
            end
            obj.solver.setInflowDirichletBoundary(inflow_value, normal_velocity_tolerance);
        end

        function out = getAdvectionBoundaryCondition(obj)
            out = obj.solver.getAdvectionBoundaryCondition();
        end

        function out = returnsDistributedState(obj)
            out = obj.solver.returnsDistributedState();
        end

        function out = getOutputRange(obj)
            out = obj.solver.getOutputRange();
        end

        function init(obj, domain, xi, dlt, num_species, patch_spacing_factor, patch_radius_factor)
            if nargin < 5 || isempty(num_species)
                error('kp:solvers:BadSpeciesCount', 'MultiSpeciesPUSLAdvectionSolver requires a positive species count.');
            end
            if nargin < 6
                patch_spacing_factor = [];
            end
            if nargin < 7
                patch_radius_factor = [];
            end
            obj.num_species = num_species;
            obj.solver.init(domain, xi, dlt, patch_spacing_factor, patch_radius_factor);
        end

        function out = getNumSpecies(obj)
            out = obj.num_species;
        end

        function out = getOutputNodes(obj)
            out = obj.solver.getOutputNodes();
        end

        function coeffs = projectInitial(obj, rho0)
            ensureInitialized(obj);
            X = obj.solver.getOutputNodes();
            nodal_samples = zeros(size(X, 1), obj.num_species);
            for i = 1:size(X, 1)
                values = rho0(X(i, :));
                if numel(values) ~= obj.num_species
                    error('kp:solvers:BadSpeciesCount', 'Initial callback returned the wrong number of species values.');
                end
                nodal_samples(i, :) = reshape(values, 1, []);
            end
            coeffs = obj.solver.projectSamples(nodal_samples);
        end

        function coeffs = projectConstant(obj, value)
            ensureInitialized(obj);
            coeffs = obj.solver.projectConstant(value, obj.num_species);
        end

        function coeffs = projectConstants(obj, values)
            ensureInitialized(obj);
            values = reshape(values, 1, []);
            if numel(values) ~= obj.num_species
                error('kp:solvers:BadSpeciesCount', 'projectConstants received the wrong number of species values.');
            end
            coeffs = repmat(values, size(obj.solver.getOutputNodes(), 1), 1);
        end

        function coeffs = projectSamples(obj, nodal_samples)
            ensureInitialized(obj);
            validateSpeciesMatrix(obj, nodal_samples, 'projectSamples');
            coeffs = obj.solver.projectSamples(nodal_samples);
        end

        function values = evaluateAtNodes(obj, coeffs)
            ensureInitialized(obj);
            validateSpeciesMatrix(obj, coeffs, 'evaluateAtNodes');
            values = obj.solver.evaluateAtNodes(coeffs);
        end

        function values = evaluateAtPoints(obj, coeffs, X)
            ensureInitialized(obj);
            validateSpeciesMatrix(obj, coeffs, 'evaluateAtPoints');
            values = obj.solver.evaluateAtPoints(coeffs, X);
        end

        function coeffs_next = forwardSLStep(obj, tn, coeffs_old, velocity, rk)
            ensureInitialized(obj);
            validateSpeciesMatrix(obj, coeffs_old, 'forwardSLStep');
            coeffs_next = obj.solver.forwardSLStep(tn, coeffs_old, velocity, rk);
        end

        function coeffs_next = backwardSLStep(obj, tn, coeffs_old, velocity, rk)
            ensureInitialized(obj);
            validateSpeciesMatrix(obj, coeffs_old, 'backwardSLStep');
            coeffs_next = obj.solver.backwardSLStep(tn, coeffs_old, velocity, rk);
        end

        function setStepSize(obj, dlt)
            obj.solver.setStepSize(dlt);
        end

        function resetSolveStats(obj)
            obj.solver.resetSolveStats();
        end

        function out = getSolveStats(obj)
            out = obj.solver.getSolveStats();
        end

        function out = getNumPatches(obj)
            out = obj.solver.getNumPatches();
        end

        function out = getNumDofsPerSpecies(obj)
            out = obj.solver.getNumDofs();
        end

        function out = getTotalDofs(obj)
            ensureInitialized(obj);
            out = obj.num_species * obj.solver.getNumDofs();
        end

        function out = getPatchRadius(obj)
            out = obj.solver.getPatchRadius();
        end

        function out = getPatchSpacing(obj)
            out = obj.solver.getPatchSpacing();
        end

        function out = totalMass(obj, coeffs, species)
            ensureInitialized(obj);
            validateSpeciesMatrix(obj, coeffs, 'totalMass');
            out = obj.solver.totalMass(coeffs, species);
        end

        function out = totalMasses(obj, coeffs)
            ensureInitialized(obj);
            validateSpeciesMatrix(obj, coeffs, 'totalMasses');
            out = zeros(obj.num_species, 1);
            for j = 1:obj.num_species
                out(j) = obj.solver.totalMass(coeffs, j);
            end
        end

        function out = getDomainMeasure(obj)
            out = obj.solver.getDomainMeasure();
        end
    end
end

function ensureInitialized(obj)
if obj.num_species <= 0
    error('kp:solvers:Uninitialized', 'MultiSpeciesPUSLAdvectionSolver must be initialized before use.');
end
end

function validateSpeciesMatrix(obj, values, where)
if size(values, 2) ~= obj.num_species
    error('kp:solvers:BadSpeciesCount', ...
        'MultiSpeciesPUSLAdvectionSolver %s received the wrong number of species columns.', where);
end
end
