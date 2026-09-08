classdef DualNodeDomainGenerator < handle
    %DUALNODEDOMAINGENERATOR Paired velocity/pressure domain generator.

    properties (SetAccess = private)
        velocityGenerator kp.nodes.DomainNodeGenerator = kp.nodes.DomainNodeGenerator()
        pressureGenerator kp.nodes.DomainNodeGenerator = kp.nodes.DomainNodeGenerator()
    end

    methods
        function generateSmoothDomainNodes(obj, geometry, velocityRadius, pressureRadius, varargin)
            obj.velocityGenerator.buildDomainDescriptorFromGeometry(geometry, velocityRadius, varargin{:});
            obj.pressureGenerator.buildDomainDescriptorFromGeometry(geometry, pressureRadius, varargin{:});
        end

        function result = generateSmoothDomainNodesAutoPressure( ...
                obj, geometry, velocityRadius, maxPressureFraction, minPressureNodes, varargin)
            if nargin < 5 || isempty(minPressureNodes)
                minPressureNodes = 0;
            end
            if ~(maxPressureFraction > 0 && maxPressureFraction < 1)
                error('kp:nodes:BadPressureFraction', ...
                    'maxPressureFraction must lie strictly between 0 and 1.');
            end

            obj.velocityGenerator.buildDomainDescriptorFromGeometry(geometry, velocityRadius, varargin{:});
            velocityDescriptor = obj.velocityGenerator.getDomainDescriptor();
            velocityNodes = velocityDescriptor.getNumIntBdryNodes();
            pressureTarget = max(1, floor(maxPressureFraction * velocityNodes));
            if pressureTarget < minPressureNodes
                error('kp:nodes:PressureFractionTooSmall', ...
                    ['Requested pressure-node fraction is incompatible with the stencil-size requirement. ' ...
                     'Increase the fraction, decrease the required pressure stencil size, or refine the velocity cloud.']);
            end

            lowRadius = max(velocityRadius, 1.0e-8);
            highRadius = 1.25 * lowRadius;
            pressureCount = inf;
            expansionIters = 0;
            while pressureCount > pressureTarget && expansionIters < 30
                obj.pressureGenerator.buildDomainDescriptorFromGeometry(geometry, highRadius, varargin{:});
                pressureCount = obj.pressureGenerator.getDomainDescriptor().getNumIntBdryNodes();
                if pressureCount > pressureTarget
                    lowRadius = highRadius;
                    highRadius = 1.2 * highRadius;
                end
                expansionIters = expansionIters + 1;
            end

            if pressureCount > pressureTarget
                error('kp:nodes:PressureSelectionFailed', ...
                    'Unable to select a pressure node cloud satisfying the requested density ratio.');
            end

            bestRadius = highRadius;
            bestCount = pressureCount;
            for iter = 1:20
                trialRadius = 0.5 * (lowRadius + highRadius);
                obj.pressureGenerator.buildDomainDescriptorFromGeometry(geometry, trialRadius, varargin{:});
                trialCount = obj.pressureGenerator.getDomainDescriptor().getNumIntBdryNodes();
                if trialCount <= pressureTarget
                    bestRadius = trialRadius;
                    bestCount = trialCount;
                    highRadius = trialRadius;
                else
                    lowRadius = trialRadius;
                end
            end

            obj.pressureGenerator.buildDomainDescriptorFromGeometry(geometry, bestRadius, varargin{:});
            while obj.pressureGenerator.getDomainDescriptor().getNumIntBdryNodes() > pressureTarget
                bestRadius = 1.02 * bestRadius;
                obj.pressureGenerator.buildDomainDescriptorFromGeometry(geometry, bestRadius, varargin{:});
            end
            bestCount = obj.pressureGenerator.getDomainDescriptor().getNumIntBdryNodes();

            if bestCount < minPressureNodes
                error('kp:nodes:PressureTooSmall', ...
                    'Automatic pressure coarsening produced too few pressure nodes for the requested stencil.');
            end

            result = struct( ...
                'pressure_radius', bestRadius, ...
                'num_velocity_nodes', velocityNodes, ...
                'num_pressure_nodes', bestCount);
        end

        function dual = createDualNodeDomainDescriptor(obj, stripPressureGhosts)
            if nargin < 2 || isempty(stripPressureGhosts)
                stripPressureGhosts = true;
            end

            dual = kp.domain.DualNodeDomainDescriptor();
            dual.setVelocityDomain(obj.velocityGenerator.getDomainDescriptor());

            pressureDomain = obj.pressureGenerator.getDomainDescriptor();
            if stripPressureGhosts
                pressureGhosts = zeros(0, pressureDomain.getDim());
                pressureDomain.setNodes( ...
                    pressureDomain.getInteriorNodes(), ...
                    pressureDomain.getBdryNodes(), ...
                    pressureGhosts);
                pressureDomain.setNormals(pressureDomain.getNrmls());
                pressureDomain.setSepRad(pressureDomain.getSepRad());
                pressureDomain.setOuterLevelSet(pressureDomain.getOuterLevelSet());
                pressureDomain.setBoundaryLevelSets(pressureDomain.getBoundaryLevelSets());
                pressureDomain.buildStructs();
            end
            dual.setPressureDomain(pressureDomain);
        end

        function out = getVelocityGenerator(obj)
            out = obj.velocityGenerator;
        end

        function out = getPressureGenerator(obj)
            out = obj.pressureGenerator;
        end
    end
end
