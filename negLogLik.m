
function nll=negLogLik(data,parq,parx,pard)
    gCheck=size(data,2);
    x=cumsum([parx(1);exp(parx(2:end))]);
    q=[exp(parq);1];
        q=q/sum(q);
        q=q.*q';
        q=q(:);
    nHet=size(x,1);
    v=repmat(x,1,nHet);
    xH=max(x,x');
        xH=xH(:)';
    xL=min(x,x');
        xL=xL(:)';
    delta=exp(pard);    
    pL=logit(solveP(@(U)(bellman(U,xL+delta))));
    pH=logit(solveP(@(U)(bellman(U,xH))));
    p22=pL;
    p21=max(pH-pL,0);
    p20=1-p22-p21;
    p11=pH;
    p10=1-pH;
    nll=0;    
    for n=1:gCheck
        ly=data(1:end-1,n);
        nw=data(2:end,n);
        lik=(ly==2).*((nw==2).*p22+(nw==1).*p21+(nw==0).*p20)+...
            (ly==1).*((nw==1).*p11+(nw==0).*p10)+...
            (ly==0).*ones(1,nHet^2);
        nll=nll-log(prod(lik)*q);
    end    
end