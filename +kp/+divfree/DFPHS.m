function [dfrbf, lap_dfrbf, hess_rbf, full_rbf] = DFPHS(dim, m_in)
%DFPHS Divergence-free PHS kernel blocks in 2D or 3D.
%   [dfrbf, lap_dfrbf, hess_rbf, full_rbf] = kp.divfree.DFPHS(dim, m)
%   returns cell arrays of numeric function handles matching the older
%   DFPHS prototype API.

validateattributes(dim, {'numeric'}, {'scalar', 'integer', '>=', 2, '<=', 3});
validateattributes(m_in, {'numeric'}, {'scalar', 'integer', 'positive'});
if mod(m_in, 2) == 0
    error('kp:divfree:EvenDegree', 'Divergence-free PHS currently expects an odd degree.');
end

m = m_in;
dfrbf = cell(dim, dim);
lap_dfrbf = cell(dim, dim);
hess_rbf = cell(dim, dim);
full_rbf = cell(dim, dim);

for i = 1:dim
    for j = 1:dim
        dfrbf{i, j} = makeKernelEntry(dim, m, i, j, "divfree");
        lap_dfrbf{i, j} = makeKernelEntry(dim, m, i, j, "lap_divfree");
        hess_rbf{i, j} = makeKernelEntry(dim, m, i, j, "hessian");
        full_rbf{i, j} = makeKernelEntry(dim, m, i, j, "full");
    end
end
end

function fh = makeKernelEntry(dim, m, i, j, mode)
switch dim
    case 2
        if i == j
            fh = @(r, a, b) evaluateEntry(dim, m, i, j, mode, r, makeDiffTensor(dim, r, i, a - b));
        else
            fh = @(r, ax, bx, ay, by) evaluateEntry(dim, m, i, j, mode, r, ...
                makeDiffTensor(dim, r, 1, ax - bx, 2, ay - by));
        end
    case 3
        if i == j
            fh = @(r, a, b) evaluateEntry(dim, m, i, j, mode, r, makeDiffTensor(dim, r, i, a - b));
        elseif any([i, j] == 1) && any([i, j] == 2)
            fh = @(r, ax, bx, ay, by) evaluateEntry(dim, m, i, j, mode, r, ...
                makeDiffTensor(dim, r, 1, ax - bx, 2, ay - by));
        elseif any([i, j] == 1) && any([i, j] == 3)
            fh = @(r, ax, bx, az, bz) evaluateEntry(dim, m, i, j, mode, r, ...
                makeDiffTensor(dim, r, 1, ax - bx, 3, az - bz));
        else
            fh = @(r, ay, by, az, bz) evaluateEntry(dim, m, i, j, mode, r, ...
                makeDiffTensor(dim, r, 2, ay - by, 3, az - bz));
        end
    otherwise
        error('kp:divfree:BadDimension', 'Only 2D and 3D are supported.');
end
end

function out = evaluateEntry(dim, m, i, j, mode, r, diffParts)
re = r + eps;
diff = diffParts;

lap_scalar = m * (m + dim - 2) * re.^(m - 2);
full_kernel = zeros(size(r));
if i == j
    full_kernel = -lap_scalar;
end

if i == j
    hess_term = m * (m - 2) * re.^(m - 4) .* diff(:, :, i) .* diff(:, :, j) + ...
        m * re.^(m - 2);
else
    hess_term = m * (m - 2) * re.^(m - 4) .* diff(:, :, i) .* diff(:, :, j);
end

divfree_kernel = full_kernel + hess_term;
lap_divfree_kernel = zeros(size(r));
if i == j
    lap_divfree_kernel = m * (m - 2) * re.^(m - 4) .* ...
        (2 - (m + dim - 3) * (m + dim - 4));
end
lap_divfree_kernel = lap_divfree_kernel + ...
    m * (m - 2) * (m - 4) * (m + dim - 2) * re.^(m - 6) .* diff(:, :, i) .* diff(:, :, j);

switch mode
    case "divfree"
        out = divfree_kernel;
    case "lap_divfree"
        out = lap_divfree_kernel;
    case "hessian"
        out = hess_term;
    case "full"
        out = full_kernel;
    otherwise
        error('kp:divfree:BadMode', 'Unknown DFPHS kernel mode.');
end
end

function diff = makeDiffTensor(dim, r, varargin)
diff = zeros([size(r), dim]);
for k = 1:2:numel(varargin)
    coord = varargin{k};
    values = varargin{k + 1};
    diff(:, :, coord) = values;
end
end
