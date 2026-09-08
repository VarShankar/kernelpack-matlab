function inds = resampleClosedCurveByArcLength(curve, targetCount)
%RESAMPLECLOSEDCURVEBYARCLENGTH Choose nearly uniform samples by arc length.

    n = size(curve, 1);
    if n == 0
        inds = zeros(0, 1);
        return;
    end
    if n == 1 || targetCount <= 1
        inds = uint32(1);
        return;
    end

    shifted = curve([2:end, 1], :);
    segLens = sqrt(sum((shifted - curve).^2, 2));
    total = sum(segLens);
    if total <= eps
        inds = uint32(round(linspace(1, n, min(n, targetCount))).');
        inds = unique(inds, 'stable');
        return;
    end

    cumLen = [0; cumsum(segLens)];
    targets = (0:targetCount-1).' * (total / targetCount);
    inds = zeros(targetCount, 1, 'uint32');
    for k = 1:targetCount
        t = targets(k);
        idx = find(cumLen <= t, 1, 'last');
        if isempty(idx)
            idx = 1;
        end
        idx = min(idx, n);
        nextIdx = mod(idx, n) + 1;
        if idx < n
            if abs(cumLen(idx + 1) - t) < abs(cumLen(idx) - t)
                idx = idx + 1;
            end
        else
            if abs(total - t) < abs(cumLen(idx) - t)
                idx = nextIdx;
            end
        end
        inds(k) = uint32(idx);
    end
    inds = unique(inds, 'stable');
end
