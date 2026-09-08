function T = convergenceResultsTable(results)
%CONVERGENCERESULTSTABLE Flatten convergence results into publication rows.
%   T has one row for every (xi, N) pair. Numeric diagnostics stored as
%   vectors in the result structure are copied to columns, while cell and
%   structure-valued implementation diagnostics remain in the MAT file.

arguments
    results (1,:) struct
end

requiredFields = {'xi', 'N', 'h', 'relerr'};
for k = 1:numel(requiredFields)
    assert(isfield(results, requiredFields{k}), ...
        'kp:manifold:MissingConvergenceField', ...
        'Convergence results must contain the field "%s".', requiredFields{k});
end

rowCounts = arrayfun(@(result) numel(result.N), results);
numRows = sum(rowCounts);
xi = zeros(numRows, 1);
exactStartup = false(numRows, 1);
N = zeros(numRows, 1);
sqrtN = zeros(numRows, 1);
fitRate = nan(numRows, 1);

rowStart = 1;
for k = 1:numel(results)
    rows = rowStart:(rowStart + rowCounts(k) - 1);
    xi(rows) = results(k).xi;
    if isfield(results, 'exactStartup')
        exactStartup(rows) = logical(results(k).exactStartup);
    end
    N(rows) = results(k).N(:);
    sqrtN(rows) = sqrt(results(k).N(:));

    valid = isfinite(results(k).N(:)) & results(k).N(:) > 0 & ...
        isfinite(results(k).relerr(:)) & results(k).relerr(:) > 0;
    if nnz(valid) >= 2
        fit = polyfit(log(sqrt(results(k).N(valid))), ...
            log(results(k).relerr(valid)), 1);
        fitRate(rows) = -fit(1);
    end
    rowStart = rows(end) + 1;
end

T = table(xi, exactStartup, N, sqrtN, fitRate);

excludedFields = {'xi', 'exactStartup', 'N', 'rate', ...
    'rearrangementTimes', 'updateStats', 'solveStats', 'solution'};
resultFields = fieldnames(results);
for k = 1:numel(resultFields)
    field = resultFields{k};
    if ismember(field, excludedFields)
        continue;
    end

    values = nan(numRows, 1);
    rowStart = 1;
    isTabularField = true;
    for j = 1:numel(results)
        rows = rowStart:(rowStart + rowCounts(j) - 1);
        value = results(j).(field);
        if ~(isnumeric(value) || islogical(value)) || ...
                ~(isscalar(value) || numel(value) == rowCounts(j))
            isTabularField = false;
            break;
        end
        if isscalar(value)
            values(rows) = double(value);
        else
            values(rows) = double(value(:));
        end
        rowStart = rows(end) + 1;
    end
    if isTabularField
        T.(field) = values;
    end
end

if isfield(results, 'rate')
    localRate = nan(numRows, 1);
    rowStart = 1;
    for k = 1:numel(results)
        rows = rowStart:(rowStart + rowCounts(k) - 1);
        localRate(rows) = results(k).rate(:);
        rowStart = rows(end) + 1;
    end
    T.localRate = localRate;
end
end
