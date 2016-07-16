function [pPL, resus] = delVcalc(lambda, e, N)
global Isp g0 Mpay W Udes

%N = round(N);
if size(lambda,2) == 1,
    pPL = (lambda/(1+lambda));
    num = (pPL^(1/N))*(1-e) + e;
    U = log(1/num);
    resus = U*Isp*g0*N;
end
if size(lambda,2)>1,
    pPL = (lambda./(1+lambda));
    num = (pPL.^(1./N)).*(1-e) + e;
    U = log(1./num);
    resus = U.*Isp*g0*N;
end


end
