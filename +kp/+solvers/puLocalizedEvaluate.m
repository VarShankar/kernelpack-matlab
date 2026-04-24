function values = puLocalizedEvaluate(domain, patchData, ~, coeffs, Xq)
%PULOCALIZEDEVALUATE Partition-of-unity evaluation on arbitrary query points.

if isempty(Xq)
    values = zeros(0, size(coeffs, 2));
    return;
end

radius = patchData.radius;
centers = patchData.centers;
node_ids = patchData.node_ids;
nc = size(coeffs, 2);
stencils = patchData.stencils;
sp = patchData.stencil_props;
patch_ids_per_query = queryPatchIds(patchData, Xq, radius);

values = zeros(size(Xq, 1), nc);
weight_sum = zeros(size(Xq, 1), 1);
for q = 1:size(Xq, 1)
    patch_ids = patch_ids_per_query{q};
    if isempty(patch_ids)
        continue;
    end
    center_dist = sqrt(sum((centers(patch_ids, :) - Xq(q, :)).^2, 2));
    alpha = puPatchWeight(center_dist ./ radius);
    alpha_sum = sum(alpha);
    if alpha_sum <= 1.0e-14
        alpha = ones(size(alpha));
        alpha_sum = sum(alpha);
    end
    alpha = alpha / alpha_sum;
    for k = 1:numel(patch_ids)
        p = patch_ids(k);
        ids = node_ids{p};
        stencil = stencils{p};
        vloc = stencil.EvalStencil(sp, Xq(q, :), coeffs(ids, :), false);
        values(q, :) = values(q, :) + alpha(k) * vloc;
    end
    weight_sum(q) = 1.0;
end

missing = weight_sum <= 1.0e-14;
if any(missing)
    [idx, ~] = domain.queryKnn("all", Xq(missing, :), 1);
values(missing, :) = coeffs(idx(:, 1), :);
    weight_sum(missing) = 1.0;
end

values = values ./ weight_sum;
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

function w = puPatchWeight(r)
w = zeros(size(r));
mask = r < 1;
t = 1 - r(mask);
rm = r(mask);
w(mask) = t.^8 .* (32 * rm.^3 + 25 * rm.^2 + 8 * rm + 1);
end
