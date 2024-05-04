
function nFirms=simdata(gCheck,tCheck,mu,sigma,delta)
global rho
    beta=mu+sigma*randn(2,gCheck);
    betaH=max(beta);
    betaL=min(beta);

    pH=logit(solveP(@(U)(bellman(U,betaH))));
    pL=logit(solveP(@(U)(bellman(U,betaL+delta))));
    
    unif = rand(tCheck,gCheck);
    stayL = unif<pL;
    stayH = unif<pH;

    nFirms = 2*ones(1,gCheck);
    for t=1:tCheck
        S = nFirms(t,:);
        S(S==2) = 1+stayL(t,S==2);
        S(S==1) = stayH(t,S==1);
        nFirms=[nFirms;S];
    end
end