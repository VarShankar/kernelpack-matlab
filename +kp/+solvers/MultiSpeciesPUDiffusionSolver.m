classdef MultiSpeciesPUDiffusionSolver < handle
    %MULTISPECIESPUDIFFUSIONSOLVER Thin multi-column wrapper for PUDiffusionSolver.

    properties (SetAccess = private)
        Domain kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
        xi (1,1) double = 0
        dt (1,1) double = 0
        nu (1,1) double = 0
        num_omp_threads (1,1) double = 1
        solvers cell = {}
        num_species (1,1) double = 0
    end

    methods
        function init(obj, domain, xi, dlt, d_coeff, num_omp_threads)
            if nargin < 6 || isempty(num_omp_threads)
                num_omp_threads = 1;
            end
            obj.Domain = domain;
            obj.xi = xi;
            obj.dt = dlt;
            obj.nu = d_coeff;
            obj.num_omp_threads = num_omp_threads;
            obj.solvers = {};
            obj.num_species = 0;
        end

        function setStepSize(obj, dlt)
            obj.dt = dlt;
            for i = 1:numel(obj.solvers)
                obj.solvers{i}.setStepSize(dlt);
            end
        end

        function setInitialState(obj, U0)
            ensureSolvers(obj, size(U0, 2));
            for species = 1:size(U0, 2)
                obj.solvers{species}.setInitialState(U0(:, species));
            end
        end

        function setStateHistory(obj, varargin)
            numSpecies = size(varargin{1}, 2);
            ensureSolvers(obj, numSpecies);
            for arg = 2:numel(varargin)
                if size(varargin{arg}, 2) ~= numSpecies
                    error('kp:solvers:BadSpeciesHistory', 'All state-history matrices must have the same species count.');
                end
            end
            for species = 1:numSpecies
                cols = cellfun(@(U) U(:, species), varargin, 'UniformOutput', false);
                obj.solvers{species}.setStateHistory(cols{:});
            end
        end

        function out = bdf1Step(obj, t, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            out = stepColumns(obj, t, forcing, NeuCoeffFunc, DirCoeffFunc, bc, "bdf1Step");
        end

        function out = bdf2Step(obj, t, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            out = stepColumns(obj, t, forcing, NeuCoeffFunc, DirCoeffFunc, bc, "bdf2Step");
        end

        function out = bdf3Step(obj, t, forcing, NeuCoeffFunc, DirCoeffFunc, bc)
            out = stepColumns(obj, t, forcing, NeuCoeffFunc, DirCoeffFunc, bc, "bdf3Step");
        end

        function out = returnsDistributedState(~)
            out = false;
        end

        function out = getOutputRange(obj)
            if isempty(obj.solvers)
                n = obj.Domain.getNumIntBdryNodes();
                out = [1, n];
            else
                out = obj.solvers{1}.getOutputRange();
            end
        end

        function out = getOutputNodes(obj)
            if isempty(obj.solvers)
                out = obj.Domain.getIntBdryNodes();
            else
                out = obj.solvers{1}.getOutputNodes();
            end
        end
    end
end

function ensureSolvers(obj, numSpecies)
if numSpecies <= 0
    error('kp:solvers:NoSpecies', 'MultiSpeciesPUDiffusionSolver requires at least one species.');
end
if isempty(obj.solvers)
    obj.solvers = cell(1, numSpecies);
    for species = 1:numSpecies
        solver = kp.solvers.PUDiffusionSolver();
        solver.init(obj.Domain, obj.xi, obj.dt, obj.nu, obj.num_omp_threads);
        obj.solvers{species} = solver;
    end
    obj.num_species = numSpecies;
    return;
end
if obj.num_species ~= numSpecies
    error('kp:solvers:SpeciesMismatch', 'MultiSpeciesPUDiffusionSolver was initialized for a different number of species.');
end
end

function out = stepColumns(obj, t, forcing, NeuCoeffFunc, DirCoeffFunc, bc, stepName)
if isempty(obj.solvers)
    error('kp:solvers:MissingInitialState', 'MultiSpeciesPUDiffusionSolver requires setInitialState() before stepping.');
end
for species = 1:obj.num_species
    forcingSpecies = @(nu, time, X) pickColumn(forcing(nu, time, X), species);
    bcSpecies = @(NeuCoeffs, DirCoeffs, nr, time, Xb) pickColumn(bc(NeuCoeffs, DirCoeffs, nr, time, Xb), species);
    col = obj.solvers{species}.(stepName)(t, forcingSpecies, NeuCoeffFunc, DirCoeffFunc, bcSpecies);
    if species == 1
        out = zeros(numel(col), obj.num_species);
    end
    out(:, species) = col;
end
end

function col = pickColumn(M, species)
col = M(:, species);
end
