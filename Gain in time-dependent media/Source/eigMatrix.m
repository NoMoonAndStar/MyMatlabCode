function M = eigMatrix(Nt, g, omega, Omega, eps, alpha)
    mu0 = 4 * pi * 1e-7;
    Nnum = length(Nt);
    MEE = -diag(Nt * g);
    MHH = -diag(Nt * g);
    MEH = -mu0 * diag(omega + Nt * Omega);
    MHE = -eps * diag(omega + Nt * Omega);

    for i = 1:Nnum

        for j = 1:Nnum

            if i - j == 1
                MHE(i, j) = -eps * (alpha * (omega + Nt(i) * Omega) + Omega * alpha);
            end

            if i - j == -1
                MHE(i, j) = -eps * (alpha * (omega + Nt(i) * Omega) - Omega * alpha);
            end

        end

    end

    M = [MEE, MEH; MHE, MHH];
end
