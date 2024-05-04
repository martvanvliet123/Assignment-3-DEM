
function Uout=bellman(U,xi)
global rho
    Uout=xi+rho*log(1+exp(U));
end