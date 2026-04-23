function [x_min, x_max] = boundingBoxExtents(geometry, preferUniform)
%BOUNDINGBOXEXTENTS Axis-aligned bounds from a geometry bounding-box cloud.

    if nargin < 2
        preferUniform = true;
    end

    box = [];
    if preferUniform && ismethod(geometry, 'getUniformBoundingBox')
        box = geometry.getUniformBoundingBox();
        if isempty(box) && ismethod(geometry, 'computeUniformBoundingBox')
            geometry.computeUniformBoundingBox();
            box = geometry.getUniformBoundingBox();
        end
    end

    if isempty(box) && ismethod(geometry, 'getBoundingBox')
        box = geometry.getBoundingBox();
        if isempty(box) && ismethod(geometry, 'computeBoundingBox')
            geometry.computeBoundingBox();
            box = geometry.getBoundingBox();
        end
    end

    if isempty(box)
        if preferUniform && ismethod(geometry, 'getUniformSampleSites')
            box = geometry.getUniformSampleSites();
        elseif ismethod(geometry, 'getSampleSites')
            box = geometry.getSampleSites();
        elseif ismethod(geometry, 'getBdryNodes')
            box = geometry.getBdryNodes();
        else
            error('kp:nodes:NoBoundingData', 'Geometry object does not expose bounding-box or boundary samples.');
        end
    end

    x_min = min(box, [], 1);
    x_max = max(box, [], 1);
end
