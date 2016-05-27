function val = fitn(x)
global cut Mpay W1 W2 g0 Isp Udes DEBUG type test

A = 0.6:-0.1:0.1;
B = 15:2:25;
W2 = 200000;

if strcmp(type,'bitString'),
    
    out = divVec(x,cut);
    N = de2re(out{1},test.Nmin,test.Nmax);
    lambda = de2re(out{2},test.Lmin,test.Lmax);
    e = de2re(out{3},test.Emin,test.Emax);
    
    if DEBUG,
        for i =1:size(out,2),
            fprintf('%g ', out{i});
        end
        
    end
    
elseif strcmp(type,'doubleVector'),
    N = round(x(1));
    lambda = x(2);
    e = x(3);
end

%=======================
%          U
%=======================
U = log((1+lambda)/(e + lambda));
delV = abs(Udes-g0*Isp*U)/Udes;


%=======================
%         Mvec
%=======================
pPL = lambda/(1+ lambda);
mf = Mpay*((lambda+e)/lambda);

num = (pPL^(1/N))*(1-e) + e;
n = (1/num)^N;

Mvec = mf*n;

%Nmax = round(N);
%for i =0:1:Nmax
%mE(i+1) = (1 - pPL^(1/Nmax))*e*Mpay/(pPL^((Nmax - i)/Nmax));
%mp(i+1) = (1 - pPL^(1/Nmax))*(1-e)*Mpay/(pPL^((Nmax - i)/Nmax));
%end
%Mvec = sum(mE)+sum(mp) + Mpay;
%Mvecfit = Mvec/W2;
%Mvec = Mpay*(lambda^(1/N) + 1)/(lambda^(1/N));

%=======================
%        Cost
%=======================
%Cost =


%=======================
%         Val
%=======================
Mvecfit =0;
val = delV; %+ Mvecfit;

if DEBUG
    fprintf('\nLamb  \t estr \t DelvaV \t Mvec \t\t Cost')
    fprintf('\n%f  \t %f \t %f \t %f \t %f \n\n',lambda,e,N,delV,Mvec)
end
%
% if DEBUG
%     g = sprintf('%d ', x);
%     fprintf('Answer: %s\n', g)
%     fprintf('N: %f \n',N)
%     fprintf('lambda: %f \n',lambda)
%     fprintf('obj(1): %f \n',obj(1))
%     fprintf('obj(2): %f \n',obj(2))
%     fprintf('val: %f \n',val)
%     fprintf('\n\n')
% end

