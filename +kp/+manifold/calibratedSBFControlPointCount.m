function [count, info] = calibratedSBFControlPointCount(geom, xi, varargin)
%CALIBRATEDSBFCONTROLPOINTCOUNT Rate-based reduced SBF control count.
%   The theory count is the first M whose material fill distance H_M satisfies
%   H_M^q <= safety*h^xi.  The returned count applies a calibration factor
%   and a small floor so low-order runs do not use under-resolved geometry.

parser = inputParser();
parser.addParameter('Degree', 7, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('NormalOrder', NaN, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('BalanceSafety', 0.1, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('ControlPointScale', 1 / 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('MinControlPointCount', 48, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('MaxControlPointCount', Inf, @(x) isnumeric(x) && isscalar(x));
parser.parse(varargin{:});
opts = parser.Results;

normalOrder = opts.NormalOrder;
if ~isfinite(normalOrder)
    normalOrder = opts.Degree + 1;
end

n = size(geom.X, 1);
targetFillDistance = max(opts.BalanceSafety, eps)^(1 / normalOrder) * ...
    max(geom.h, eps)^(xi / normalOrder);
[~, theoryInfo] = kp.manifold.selectSBFControlPoints(geom, ...
    'FillDistanceTarget', targetFillDistance, ...
    'MaxCount', opts.MaxControlPointCount);

theoryCount = theoryInfo.controlPointCount;
minCount = min(n, max(1, round(opts.MinControlPointCount)));
scaledCount = ceil(opts.ControlPointScale * theoryCount);
count = min(n, max(minCount, scaledCount));

[~, selectedInfo] = kp.manifold.selectSBFControlPoints(geom, 'Count', count);
info = selectedInfo;
info.theoryControlPointCount = theoryCount;
info.controlPointScale = opts.ControlPointScale;
info.minControlPointCount = minCount;
info.normalOrder = normalOrder;
info.balanceSafety = opts.BalanceSafety;
info.targetFillDistance = targetFillDistance;
info.predictedRMSAngleDegrees = rad2deg(info.fillDistance ^ normalOrder);
end
