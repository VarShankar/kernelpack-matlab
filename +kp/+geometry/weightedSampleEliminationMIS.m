function keep = weightedSampleEliminationMIS(X, h, varargin)
%WEIGHTEDSAMPLEELIMINATIONMIS KernelPack-style weighted sample elimination.

parser = inputParser();
parser.addParameter('UseParallel', true, @(x) islogical(x) || isnumeric(x));
parser.parse(varargin{:});
useParallel = logical(parser.Results.UseParallel);

n = size(X, 1);
if n == 0
    keep = false(0, 1);
    return;
end

searcher = buildSearcher(X);
weights = zeros(n, 1);
conflictNeighbors = cell(n, 1);
useParfor = localShouldUseParfor(n, useParallel);

if useParfor
    parfor i = 1:n
        [weightIdx, weightDist] = rangeQuery(searcher, X, X(i, :), 2 * h);
        maskWeight = weightIdx ~= i;
        c = min(weightDist(maskWeight), 2 * h);
        weights(i) = sum((1 - c ./ (2 * h)).^8);
        [conflictIdx, ~] = rangeQuery(searcher, X, X(i, :), h);
        conflictNeighbors{i} = conflictIdx(conflictIdx ~= i);
    end
else
    for i = 1:n
        [weightIdx, weightDist] = rangeQuery(searcher, X, X(i, :), 2 * h);
        maskWeight = weightIdx ~= i;
        c = min(weightDist(maskWeight), 2 * h);
        weights(i) = sum((1 - c ./ (2 * h)).^8);
        [conflictIdx, ~] = rangeQuery(searcher, X, X(i, :), h);
        conflictNeighbors{i} = conflictIdx(conflictIdx ~= i);
    end
end

active = true(n, 1);
winners = zeros(n, 1);
nWinners = 0;

while true
    roundWinner = false(n, 1);
    for i = 1:n
        if ~active(i)
            continue;
        end
        dominates = true;
        nbrs = conflictNeighbors{i};
        for jj = 1:numel(nbrs)
            j = nbrs(jj);
            if ~active(j)
                continue;
            end
            if higherPriority(weights(j), j, weights(i), i)
                dominates = false;
                break;
            end
        end
        roundWinner(i) = dominates;
    end

    idx = find(roundWinner);
    if isempty(idx)
        break;
    end

    nw = numel(idx);
    winners(nWinners + (1:nw)) = idx;
    nWinners = nWinners + nw;
    active(idx) = false;

    for ii = 1:nw
        nbrs = conflictNeighbors{idx(ii)};
        active(nbrs) = false;
    end
end

winners = winners(1:nWinners);
if isempty(winners)
    keep = false(n, 1);
    return;
end

[~, order] = sortrows([-weights(winners), winners], [1 2]);
winners = winners(order);
keep = false(n, 1);
keep(winners) = true;

function tf = higherPriority(lhsWeight, lhsIndex, rhsWeight, rhsIndex)
if lhsWeight > rhsWeight
    tf = true;
elseif lhsWeight < rhsWeight
    tf = false;
else
    tf = lhsIndex < rhsIndex;
end

function tf = localShouldUseParfor(nPts, useParallel)
tf = false;
if ~useParallel
    return;
end
if nPts < 64
    return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
    return;
end
tf = true;

function searcher = buildSearcher(X)
searcher = [];
if exist('KDTreeSearcher', 'class') == 8 && exist('rangesearch', 'file') == 2
    searcher = KDTreeSearcher(X);
end

function [idx, dist] = rangeQuery(searcher, X, xq, radius)
if ~isempty(searcher)
    [idxCell, distCell] = rangesearch(searcher, xq, radius);
    idx = idxCell{1}(:);
    dist = distCell{1}(:);
else
    d = sqrt(sum((X - xq).^2, 2));
    idx = find(d <= radius);
    dist = d(idx);
end
