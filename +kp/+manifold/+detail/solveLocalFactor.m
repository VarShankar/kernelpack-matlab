function X = solveLocalFactor(factor, B)
%SOLVELOCALFACTOR Solve using a factor from factorLocalAugmentedSystem.

switch string(factor.kind)
    case "lu"
        X = factor.U \ (factor.L \ B(factor.P, :));

    case "qr"
        Y = zeros(size(factor.R, 2), size(B, 2));
        if factor.rank > 0
            Y(1:factor.rank, :) = factor.R(1:factor.rank, 1:factor.rank) \ ...
                (factor.Q(:, 1:factor.rank).' * B);
        end

        if isvector(factor.E)
            X = zeros(size(Y));
            X(factor.E, :) = Y;
        else
            X = factor.E * Y;
        end

    otherwise
        error('kp:manifold:BadLocalFactor', ...
            'Unknown local factorization kind "%s".', factor.kind);
end
end
