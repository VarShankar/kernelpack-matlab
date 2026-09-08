function values = puLocalizedEvaluate(domain, patchData, ~, coeffs, Xq)
%PULOCALIZEDEVALUATE Partition-of-unity evaluation on arbitrary query points.

if isempty(Xq)
    values = zeros(0, size(coeffs, 2));
    return;
end

node_ids = patchData.node_ids;
nc = size(coeffs, 2);
stencils = patchData.stencils;
sp = patchData.stencil_props;
[~, alphaMat] = kp.solvers.puQueryPatchWeights(patchData, Xq);

values = zeros(size(Xq, 1), nc);
weight_sum = zeros(size(Xq, 1), 1);
for p = 1:numel(node_ids)
    qIdx = find(alphaMat(:, p));
    if isempty(qIdx)
        continue;
    end
    alpha = full(alphaMat(qIdx, p));
    ids = node_ids{p};
    stencil = stencils{p};
    vloc = stencil.EvalStencil(sp, Xq(qIdx, :), coeffs(ids, :), false);
    values(qIdx, :) = values(qIdx, :) + alpha .* vloc;
    weight_sum(qIdx) = weight_sum(qIdx) + alpha;
end

missing = weight_sum <= 1.0e-14;
if any(missing)
    [idx, ~] = domain.queryKnn("all", Xq(missing, :), 1);
values(missing, :) = coeffs(idx(:, 1), :);
    weight_sum(missing) = 1.0;
end

values = values ./ weight_sum;
end
