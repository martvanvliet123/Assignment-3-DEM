
function Unew=solveP(f)
    tol=1E-9;
    Unew=0;
    Uold=Unew+2*tol;
    while max(abs(Uold-Unew))>tol
        Uold=Unew;
        Unew=f(Uold);
    end
end