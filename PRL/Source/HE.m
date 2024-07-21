function M = HE(a, w, n, g ,Omega, eps0, mu0)
    [MEE, MEH, MHH, MHE]=deal(zeros(2*n+1,2*n+1));
    for i=-n:n
        MEE(i+n+1,i+n+1)=-i*g;
        MHH(i+n+1,i+n+1)=-i*g;
        MEH(i+n+1,i+n+1)=-mu0*(w+i*Omega);
    end
    MHE=diag(-eps0*(w+(-n:n)*Omega));
    for i=-n:n
        for j=-n:n
            if i-j==1
                MHE(i+n+1,j+n+1)=-a*eps0*(w+(i)*Omega);
            end
            if i-j==-1
                MHE(i+n+1,j+n+1)=-a*eps0*(w+(i)*Omega);
            end
        end
    end
    M=[MEE, MEH; MHE, MHH];
end