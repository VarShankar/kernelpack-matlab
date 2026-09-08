function a = total_degree_indices(d, k)
%TOTAL_DEGREE_INDICES Total-degree multi-index set up to order k in d dimensions.
%
% Ported from the polynomial-support routines used alongside Akil Narayan's
% Jacobi evaluation code.

n = (1:nchoosek(d + k, d));
n = n - 1;
n = n(:);
n = n + 1;
n_max = max(n);

N = 0;
while nchoosek(N + d, d) < n_max
    N = N + 1;
end

a = zeros(nchoosek(N + d, d), d);
row_id = 1;
a(row_id, :) = zeros(1, d);
row_id = 2;

if N == 0
    a = zeros(size(n, 1), d);
else
    for q = 1:N
        current_row = zeros(1, d);
        current_row(1) = q;
        a(row_id, :) = current_row;
        row_id = row_id + 1;

        onesman_home = 1;
        onesman_location = 1;
        finished = false;

        while ~finished
            [a, row_id, current_row, onesman_home, onesman_location, finished] = ...
                onesmanPilgrimage(a, row_id, current_row, onesman_home, onesman_location, q, d);
        end
    end

    a = a(n, :);
end
end

function [a, row_id, current_row, onesman_home, onesman_location, finished] = ...
    onesmanPilgrimage(a, row_id, current_row, onesman_home, onesman_location, q, d)

while onesman_location < d
    onesman_location = onesman_location + 1;
    current_row(onesman_location - 1) = current_row(onesman_location - 1) - 1;
    current_row(onesman_location) = current_row(onesman_location) + 1;
    a(row_id, :) = current_row;
    row_id = row_id + 1;
end

if onesman_home + 1 == d
    while current_row(onesman_home) > 0
        current_row(end) = current_row(end) + 1;
        current_row(end - 1) = current_row(end - 1) - 1;
        a(row_id, :) = current_row;
        row_id = row_id + 1;
    end
end

if current_row(end) == q
    finished = true;
    return;
end

columns = find(current_row, 2, 'last');
current_row(columns(1)) = current_row(columns(1)) - 1;
current_row(columns(1) + 1) = current_row(end) + 1;
current_row(end) = 0;
a(row_id, :) = current_row;
row_id = row_id + 1;

onesman_home = columns(1) + 1;
onesman_location = columns(1) + 1;
finished = false;
end
