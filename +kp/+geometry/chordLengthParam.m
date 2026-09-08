function t = chordLengthParam(X, isClosed)
%CHORDLENGTHPARAM Chord-length parameter values on [0,1].

if nargin < 2
    isClosed = false;
end

n = size(X, 1);
if n == 1
    t = 0;
    return;
end

d = sqrt(sum(diff(X, 1, 1).^2, 2));
if isClosed
    d = [d; norm(X(1, :) - X(end, :), 2)];
    s = [0; cumsum(d(1:end-1))];
    total = sum(d);
else
    s = [0; cumsum(d)];
    total = s(end);
end

if total <= eps
    t = linspace(0, 1, n).';
else
    t = s / total;
end

