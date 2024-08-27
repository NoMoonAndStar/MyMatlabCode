function J = jordan_matrix(alpha, n)

    arguments , alpha(1, 1), n(1, 1){mustBePositiveInteger};
    end

    v = ones(1, n - 1);
    J = -alpha * eye(n) + diag(v, 1);
end
