function [patch_ids_per_query, alphaMat] = puQueryPatchWeights(patchData, Xq)
%PUQUERYPATCHWEIGHTS Query active PU patches and normalized patch weights.

if isempty(Xq)
    patch_ids_per_query = cell(0, 1);
    alphaMat = sparse(0, size(patchData.centers, 1));
    return;
end

radius = patchData.radius;
centers = patchData.centers;
patch_ids_per_query = queryPatchIds(patchData, Xq, radius);

I = cell(size(Xq, 1), 1);
J = cell(size(Xq, 1), 1);
S = cell(size(Xq, 1), 1);
for q = 1:numel(patch_ids_per_query)
    patch_ids = patch_ids_per_query{q};
    if isempty(patch_ids)
        continue;
    end
    center_dist = sqrt(sum((centers(patch_ids, :) - Xq(q, :)).^2, 2));
    alpha = kp.solvers.puPatchWeight(center_dist ./ radius);
    alpha_sum = sum(alpha);
    if alpha_sum <= 1.0e-14
        alpha = ones(size(alpha));
        alpha_sum = sum(alpha);
    end
    alpha = alpha / alpha_sum;
    I{q} = q * ones(numel(patch_ids), 1);
    J{q} = patch_ids(:);
    S{q} = alpha(:);
end

alphaMat = sparse(vertcat(I{:}), vertcat(J{:}), vertcat(S{:}), size(Xq, 1), size(centers, 1));
end

function patch_ids_per_query = queryPatchIds(patchData, Xq, radius)
tree = patchData.center_tree;
centers = patchData.centers;
if tree.HasSearcher
    patch_ids_per_query = rangesearch(tree.Searcher, Xq, radius);
else
    D = kp.geometry.distanceMatrix(Xq, centers);
    patch_ids_per_query = cell(size(Xq, 1), 1);
    for q = 1:size(Xq, 1)
        patch_ids_per_query{q} = find(D(q, :) < radius);
    end
end

for q = 1:numel(patch_ids_per_query)
    if isempty(patch_ids_per_query{q})
        if isempty(centers)
            continue;
        end
        d = sqrt(sum((centers - Xq(q, :)).^2, 2));
        [~, nearest_patch] = min(d);
        patch_ids_per_query{q} = nearest_patch;
    else
        patch_ids_per_query{q} = patch_ids_per_query{q}(:).';
    end
end
end
