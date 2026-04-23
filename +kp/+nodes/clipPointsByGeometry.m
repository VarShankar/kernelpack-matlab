function [Xkeep, keepMask, phi] = clipPointsByGeometry(X, geometry, varargin)
%CLIPPOINTSBYGEOMETRY Keep points inside or outside a geometry level set.

    parser = inputParser();
    parser.addRequired('X', @(x) validateattributes(x, {'numeric'}, {'2d', 'real', 'finite'}));
    parser.addRequired('geometry');
    parser.addParameter('Keep', 'inside', @(x) any(strcmpi(x, {'inside', 'outside'})));
    parser.addParameter('Tolerance', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite'}));
    parser.addParameter('AutoBuildLevelSet', true, @(x) islogical(x) || isnumeric(x));
    parser.parse(X, geometry, varargin{:});
    opts = parser.Results;

    levelSet = ensureGeometryLevelSet(geometry, logical(opts.AutoBuildLevelSet));
    phi = levelSet.Evaluate(X);

    switch lower(opts.Keep)
        case 'inside'
            keepMask = phi <= opts.Tolerance;
        case 'outside'
            keepMask = phi >= -opts.Tolerance;
        otherwise
            error('kp:nodes:BadKeepMode', 'Keep must be ''inside'' or ''outside''.');
    end

    Xkeep = X(keepMask, :);
end

function levelSet = ensureGeometryLevelSet(geometry, autoBuild)
    if ~ismethod(geometry, 'getLevelSet')
        error('kp:nodes:MissingLevelSet', 'Geometry object must provide getLevelSet().');
    end

    levelSet = geometry.getLevelSet();
    needsBuild = isempty(levelSet) || ~isa(levelSet, 'kp.geometry.RBFLevelSet') || levelSet.n == 0;
    if needsBuild
        if ~autoBuild
            error('kp:nodes:MissingBuiltLevelSet', 'Geometry level set is not built.');
        end
        if ismethod(geometry, 'buildLevelSetFromGeometricModel')
            geometry.buildLevelSetFromGeometricModel([]);
        elseif ismethod(geometry, 'buildLevelSet')
            geometry.buildLevelSet();
        else
            error('kp:nodes:CannotBuildLevelSet', 'Geometry object does not expose a recognized level-set builder.');
        end
        levelSet = geometry.getLevelSet();
    end
end
