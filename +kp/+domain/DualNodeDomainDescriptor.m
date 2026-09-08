classdef DualNodeDomainDescriptor < handle
    %DUALNODEDOMAINDESCRIPTOR KernelPack-style paired velocity/pressure domains.

    properties (Access = private)
        velocityDomain_ kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
        pressureDomain_ kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
        hasVelocity_ (1,1) logical = false
        hasPressure_ (1,1) logical = false
    end

    methods
        function setVelocityDomain(obj, domain)
            obj.velocityDomain_ = domain;
            obj.hasVelocity_ = true;
        end

        function setPressureDomain(obj, domain)
            obj.pressureDomain_ = domain;
            obj.hasPressure_ = true;
        end

        function buildStructs(obj)
            if obj.hasVelocity_
                obj.velocityDomain_.buildStructs();
            end
            if obj.hasPressure_
                obj.pressureDomain_.buildStructs();
            end
        end

        function out = hasVelocityDomain(obj)
            out = obj.hasVelocity_;
        end

        function out = hasPressureDomain(obj)
            out = obj.hasPressure_;
        end

        function out = getVelocityDomain(obj)
            if ~obj.hasVelocity_
                error('kp:domain:MissingVelocityDomain', ...
                    'DualNodeDomainDescriptor does not contain a velocity domain.');
            end
            out = obj.velocityDomain_;
        end

        function out = getPressureDomain(obj)
            if ~obj.hasPressure_
                error('kp:domain:MissingPressureDomain', ...
                    'DualNodeDomainDescriptor does not contain a pressure domain.');
            end
            out = obj.pressureDomain_;
        end

        function out = getDim(obj)
            if obj.hasVelocity_
                out = obj.velocityDomain_.getDim();
                return;
            end
            if obj.hasPressure_
                out = obj.pressureDomain_.getDim();
                return;
            end
            out = 0;
        end
    end
end
