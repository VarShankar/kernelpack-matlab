classdef MeshfreeGeometricMultilevelSolver < handle
    %MESHFREEGEOMETRICMULTILEVELSOLVER MGM integration shell.
    %
    % Keeps the KernelPack-facing configuration and input plumbing for the
    % meshfree geometric multilevel method of Wright, Jones, and Shankar,
    % SIAM J. Sci. Comput. 45 (2023), A312-A337. The executable MGM
    % hierarchy/V-cycle algorithm is intentionally not included here.

    properties
        LevelsData struct = struct([])
        DomainVolume (1,1) double = 0
        CoarseningFactor (1,1) double = 4
        Nmin (1,1) double = 250
        PreSmooth (1,1) double = 1
        PostSmooth (1,1) double = 1
        MaxIters (1,1) double = 100
        HasConstNullspace (1,1) logical = false
        StencilSize (1,1) double = 3
        PolynomialDegree (1,1) double = 0
        TransferOperator (1,:) char = 'rbf'
        Verbose (1,1) logical = false
        CoarseningSeed (1,1) double = 17
        UseParallel (1,1) logical = true
    end

    methods
        function obj = MeshfreeGeometricMultilevelSolver(Lh, X, domainVolume, varargin)
            if nargin == 0
                return;
            end

            parser = inputParser();
            parser.addRequired('Lh', @(x) isnumeric(x) && ismatrix(x));
            parser.addRequired('X', @(x) isnumeric(x) && ismatrix(x));
            parser.addRequired('domainVolume', @(x) isempty(x) || (isscalar(x) && isnumeric(x) && isfinite(x) && x > 0));
            parser.addParameter('HasConstNullspace', obj.HasConstNullspace, @(x) islogical(x) || isnumeric(x));
            parser.addParameter('CoarseningFactor', obj.CoarseningFactor, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite', '>', 1}));
            parser.addParameter('Nmin', obj.Nmin, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 2}));
            parser.addParameter('PreSmooth', obj.PreSmooth, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 0}));
            parser.addParameter('PostSmooth', obj.PostSmooth, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 0}));
            parser.addParameter('MaxIters', obj.MaxIters, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 1}));
            parser.addParameter('StencilSize', obj.StencilSize, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 1}));
            parser.addParameter('PolynomialDegree', obj.PolynomialDegree, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 0}));
            parser.addParameter('TransferOperator', obj.TransferOperator, @(x) any(strcmpi(string(x), ["rbf", "gmls"])));
            parser.addParameter('Verbose', obj.Verbose, @(x) islogical(x) || isnumeric(x));
            parser.addParameter('CoarseningSeed', obj.CoarseningSeed, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite'}));
            parser.addParameter('UseParallel', obj.UseParallel, @(x) islogical(x) || isnumeric(x));
            parser.parse(Lh, X, domainVolume, varargin{:});

            if size(Lh, 1) ~= size(Lh, 2) || size(Lh, 1) ~= size(X, 1)
                error('kp:solvers:MGMSizeMismatch', 'Lh must be square with one row per point in X.');
            end

            obj.HasConstNullspace = logical(parser.Results.HasConstNullspace);
            obj.CoarseningFactor = parser.Results.CoarseningFactor;
            obj.Nmin = parser.Results.Nmin;
            obj.PreSmooth = parser.Results.PreSmooth;
            obj.PostSmooth = parser.Results.PostSmooth;
            obj.MaxIters = parser.Results.MaxIters;
            obj.StencilSize = parser.Results.StencilSize;
            obj.PolynomialDegree = parser.Results.PolynomialDegree;
            obj.TransferOperator = lower(char(parser.Results.TransferOperator));
            obj.Verbose = logical(parser.Results.Verbose);
            obj.CoarseningSeed = parser.Results.CoarseningSeed;
            obj.UseParallel = logical(parser.Results.UseParallel);

            if isempty(domainVolume)
                obj.DomainVolume = estimateDomainVolume(X);
            else
                obj.DomainVolume = domainVolume;
            end
            obj.build(Lh, X);
        end

        function build(obj, Lh, X)
            validateTransferStencil(obj.StencilSize, obj.PolynomialDegree, size(X, 2));
            obj.LevelsData = struct( ...
                'nodes', double(X), ...
                'stencilSize', obj.StencilSize, ...
                'polyDegree', obj.PolynomialDegree, ...
                'Lh', sparse(Lh), ...
                'DLh', [], ...
                'I', [], ...
                'R', [], ...
                'Mhf', [], ...
                'Nhf', [], ...
                'Mhb', [], ...
                'Nhb', [], ...
                'w', ones(size(X, 1), 1));
        end

        function varargout = solve(obj, varargin) %#ok<INUSD>
            algorithmUnavailable();
            varargout = cell(1, nargout);
        end

        function varargout = precondition(obj, varargin) %#ok<INUSD>
            algorithmUnavailable();
            varargout = cell(1, nargout);
        end

        function y = applyConstrainedOperator(obj, x)
            if isempty(obj.LevelsData)
                error('kp:solvers:MGMNotInitialized', 'MGM solver has not been initialized with a matrix and point cloud.');
            end
            u = x(1:end-1);
            lambda = x(end);
            level = obj.LevelsData(1);
            y = [level.Lh * u + level.w * lambda; level.w.' * u];
        end
    end
end

function validateTransferStencil(stencilSize, polyDegree, dim)
polyDim = size(kp.poly.total_degree_indices(dim, polyDegree), 1);
if stencilSize <= polyDim
    error('kp:solvers:MGMStencilTooSmall', ...
        'MGM interpolation stencil size %d must be larger than polynomial dimension %d.', stencilSize, polyDim);
end
end

function volume = estimateDomainVolume(X)
if size(X, 1) < 2
    volume = 1;
    return;
end
if exist('KDTreeSearcher', 'class') == 8 && exist('knnsearch', 'file') == 2
    tree = KDTreeSearcher(X);
    [~, d] = knnsearch(tree, X, 'K', 2);
else
    d = zeros(size(X, 1), 2);
    for i = 1:size(X, 1)
        distances = sqrt(sum((X - X(i, :)).^2, 2));
        distances = sort(distances, 'ascend');
        d(i, :) = distances(1:2).';
    end
end
spacing = mean(d(:, 2));
volume = max(size(X, 1) * spacing^size(X, 2), eps);
end

function algorithmUnavailable()
error('kp:solvers:MGMAlgorithmUnavailable', ...
    ['MGM configuration plumbing is present, but the executable meshfree ' ...
     'geometric multilevel hierarchy and V-cycle algorithm are not included.']);
end
