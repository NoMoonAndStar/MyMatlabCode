function A = nonsingularmatrix(n)

    arguments , n(1, 1){mustBePositiveInteger};
    end

    while true
        A = randi([0 1], n);

        if rank(A) == n
            break
        end

    end

end
