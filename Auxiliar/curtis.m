clc
clear all

e = 0.15%0.4968;
pPL = 0.05;
Isp = 350;
g = 0.00981;
N = 3;
Mpay = 1000;

lambda = (pPL^(1/N))/(1 - (pPL^(1/N)))
%lambda = 0.4219; %0.5833

pPL = (lambda/(1+lambda))^N;
num = (pPL^(1/N))*(1-e) + e;
U = log(1/num);
delV = U*Isp*g*N

Nmax = round(N);
for i =0:1:Nmax - 1
    pPL = (lambda/(1+lambda))^Nmax;
    mE(i+1) = (1 - pPL^(1/Nmax))*e*Mpay/(pPL^((Nmax - i)/Nmax));
    mp(i+1) = (1 - pPL^(1/Nmax))*(1-e)*Mpay/(pPL^((Nmax - i)/Nmax));
end
Mvec = sum(mE) + sum(mp) + Mpay;

for n = 1:Nmax,
    if n == 1,
        m(n) = Mvec - mp(n);
        v(n) = Isp*g*log(Mvec/m(n));
    else
        m(n-1) = m(n-1) - mE(n-1);
        m(n) = m(n-1) - mp(n);
        v(n) = Isp*g*log(m(n-1)/m(n));
    end
end
v
sum(v)

% m01 = Mvec - mp(1);
% v2(1) = Isp*g*log(Mvec/m01);
% 
% m01 = m01 - mE(1);
% m02 = m01 - mp(2);
% v2(2) = Isp*g*log(m01/m02);
% 
% m02 = m02 - mE(2);
% m03 = m02 - mp(3) ;
% v2(3) = Isp*g*log(m02/m03);
% 
% sum(v2)

%vfinal = Isp*g*log(Mvec/m03)