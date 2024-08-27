function C = kronfuns(A, B, fun)
    [ma, na] = size(A);
    C = [];

    for i = 1:ma
        c = fun(A(i, 1), B)

        for j = 2:na
            c = [c, fun(A(i, j), B)];
        end

        C = [C; c];
    end

end
