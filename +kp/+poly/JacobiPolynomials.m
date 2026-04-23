classdef JacobiPolynomials
    %JACOBIPOLYNOMIALS Jacobi polynomial recurrence and evaluation helpers.
    %
    % This module is based on the Jacobi polynomial recurrence and
    % orthonormal evaluation pathway used in Akil Narayan's MATLAB routines.

    methods (Static)
        function [a, b] = recurrence(N, alph, bet)
            %RECURRENCE Return the first N Jacobi three-term recurrence coefficients.
            N = max(N(:));
            n = (1:N).' - 1;

            a = (bet^2 - alph^2) * ones(size(n));
            b = ones(size(n));

            flags0 = (n == 0);
            if any(flags0)
                a(flags0) = (bet - alph) / (alph + bet + 2);
                b(flags0) = exp((alph + bet + 1) * log(2) + gammaln(alph + 1) + ...
                    gammaln(bet + 1) - gammaln(alph + bet + 2));
            end

            flags1 = (n == 1);
            if any(flags1)
                a(flags1) = a(flags1) ./ ((2 + alph + bet) * (4 + alph + bet));
                b(flags1) = 4 * (1 + alph) * (1 + bet) ./ ...
                    ((2 + alph + bet)^2 * (3 + alph + bet));
            end

            flags = ~(flags0 | flags1);
            if any(flags)
                a(flags) = a(flags) ./ ...
                    ((2 * n(flags) + alph + bet) .* (2 * n(flags) + alph + bet + 2));
                b(flags) = 4 * n(flags) .* (n(flags) + alph) .* (n(flags) + bet) .* ...
                    (n(flags) + alph + bet);
                b(flags) = b(flags) ./ ...
                    ((2 * n(flags) + alph + bet).^2 .* ...
                    (2 * n(flags) + alph + bet + 1) .* ...
                    (2 * n(flags) + alph + bet - 1));
            end
        end

        function p = evaluate(a, b, x, N, d)
            %EVALUATE Evaluate orthonormal polynomials and derivatives from recurrence data.
            if nargin < 5
                d = 0;
            end

            assert(d >= 0, 'Derivative order d must be nonnegative.');
            assert(N >= 0, 'Polynomial degree N must be nonnegative.');
            assert(N <= length(a), 'N must not exceed the recurrence length.');
            assert(N <= length(b), 'N must not exceed the recurrence length.');

            nx = numel(x);
            p = zeros(nx, N + 1);
            xf = x(:);

            p(:, 1) = ones(nx, 1) / sqrt(b(1));
            if N > 0
                p(:, 2) = ((xf - a(1)) .* p(:, 1)) / sqrt(b(2));
            end

            for q = 2:N
                p(:, q + 1) = (xf - a(q)) .* p(:, q) - sqrt(b(q)) * p(:, q - 1);
                p(:, q + 1) = p(:, q + 1) / sqrt(b(q + 1));
            end

            if d == 0
                return;
            end

            for qd = 1:d
                pd = zeros(size(p));
                for q = qd:N
                    if q == qd
                        pd(:, q + 1) = exp(gammaln(qd + 1) - 0.5 * sum(log(b(1:q + 1))));
                    else
                        pd(:, q + 1) = (xf - a(q)) .* pd(:, q) - sqrt(b(q)) * pd(:, q - 1) + qd * p(:, q);
                        pd(:, q + 1) = pd(:, q + 1) / sqrt(b(q + 1));
                    end
                end
                p = pd;
            end
        end

        function p = tensorEvaluate(x, alpha, recurrenceHandle, d)
            %TENSOREVALUATE Evaluate tensor-product orthonormal polynomials.
            [m, dim] = size(x);
            if nargin < 4
                d = zeros(1, dim);
            end

            assert(size(alpha, 2) == dim, 'x and alpha must have the same number of columns.');
            assert(size(d, 2) == dim, 'x and d must have the same number of columns.');

            n = size(alpha, 1);
            pcount = size(d, 1);

            [a, b] = recurrenceHandle(max(alpha(:)) + 1);
            p = ones(m, n, pcount) / sqrt(b(1)^dim);
            activeAlpha = alpha > 0;

            for qd = 1:pcount
                for qdim = 1:dim
                    temp = kp.poly.JacobiPolynomials.evaluate(a, b, x(:, qdim), ...
                        max(alpha(:, qdim)), d(qd, qdim));
                    if d(qd, qdim) > 0
                        p(:, :, qd) = p(:, :, qd) .* temp(:, alpha(:, qdim) + 1) * sqrt(b(1));
                    else
                        mask = activeAlpha(:, qdim);
                        p(:, mask, qd) = p(:, mask, qd) .* ...
                            temp(:, alpha(mask, qdim) + 1) * sqrt(b(1));
                    end
                end
            end
        end
    end
end
